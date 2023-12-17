#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process filtool {
    label 'filtool'
    container "${params.fold_singularity_image}"
    //scratch "${params.tmp_dir}"

    input:
    val(input_data) 
    val rfi_filter
    val threads
    val telescope

    output:
    path("*.fil")

    script:
    def outputFile = "${input_data[3].trim()}_${input_data[5].trim()}_${input_data[4].trim()}"
    """
    filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${input_data[2]} -s ${input_data[3]}
    """
}

process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "${params.fold_singularity_image}"

    input:
    path(fil_file)

    output:
    path "nearest_power_of_two.txt"

    script:
    """
    #!/bin/bash
   
    output=\$(readfile ${fil_file})
    echo "\$output"

    value=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')

    log2=\$(echo "l(\$value)/l(2)" | bc -l)
    decimal_part=\$(echo "\$log2" | awk -F"." '{print "0."\$2}')
    rounded_log2=\$(echo "\$log2" | awk -F"." '{if (\$2 >= 35) print \$1+1; else print \$1}')

    nearest_power_of_2=\$((2**\$rounded_log2))
    echo \$nearest_power_of_2 > nearest_power_of_two.txt
    """
}
process peasoup {
    label 'peasoup'
    container "${params.search_singularity_image}"
    publishDir "SEARCH/${target_name}/${utc}/", pattern: "**/*.xml", mode: 'copy'

    input:
    path(fil_file)
    path(dm_file) 
    val(fft_size)
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    val(target_name)
    val(utc)

    output:
    path("**/*.xml")

    script:
    """
    #!/bin/bash
   
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus}

    """
}

process create_phase_predictor {
    label 'create_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    path(xml_file)
    val(period_start)
    val(period_end)
    val(target_name)

    output:
    path "predictor_candidate_*"

    script:
    """
    python ${params.phase_predictor} -i ${xml_file} -ps ${period_start} -pe ${period_end} -t ${target_name}
    """
}

process fold_apsuse {
    label 'fold_apsuse' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/00_APSUSE_FOLDS/", pattern: "*.ar", mode: 'copy'

    input:
    path(fil_file)
    path(phase_predictor)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)
    val(target_name)

    output:
    path "*.ar"

    script:
    """
    output_filename=\$(basename ${fil_file} .fil)
    dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename} ${fil_file}
    """
}

process fold_ptuse {
    label 'fold_ptuse' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/01_PTUSE_FOLDS/", pattern: "*.ar", mode: 'copy'

    input:
    tuple val(target_name), val(utc), val(psrfits_file_path)
    path(ephemeris_file)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    output:
    path "*.ar"

    script:
    """
    # Find at least one file that matches the wildcard pattern
    found_file=\$(ls -v ${psrfits_file_path} 2> /dev/null | head -n 1)

    if [[ -z "\$found_file" ]]; then
        echo "Error: No matching PTUSE files found."
        exit 1
    fi

    # Construct the output file name
    output_filename="${target_name}_${utc}"

    # Run the dspsr command
    dspsr -scloffs -k ${telescope} -t${threads} -E ${ephemeris_file} -L ${subint_length} -b${bins} -A -O \${output_filename} \$(ls -v ${psrfits_file_path})
    """
}


process clfd_apsuse {
    label 'clfd' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/00_APSUSE_FOLDS/", pattern: "*.clfd", mode: 'copy'

    input:
    path(fold_archive)
    val(target_name)
    val(qmask)
    val(qspike)
    val(processes)

    output:
    path "*.clfd"

    script:
    """
    clfd --fmt psrchive ${fold_archive} --features std ptp lfamp --qmask ${qmask} --despike --qspike ${qspike} --processes ${processes} -o \${PWD}
    """
}

//Nextflow requires to have separate processes if you need to re-use them. Very strange!
process clfd_ptuse {
    label 'clfd' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/01_PTUSE_FOLDS/", pattern: "*.clfd", mode: 'copy'

    input:
    path(fold_archive)
    val(target_name)
    val(qmask)
    val(qspike)
    val(processes)

    output:
    path "*.clfd"

    script:
    """
    clfd --fmt psrchive ${fold_archive} --features std ptp lfamp --qmask ${qmask} --despike --qspike ${qspike} --processes ${processes} -o \${PWD}
    """
}

process pdmp_apsuse {
    label 'pdmp' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/00_APSUSE_FOLDS/", pattern: "*.png", mode: 'copy'

    input:
    path(fold_archive)
    val(target_name)
    val(nchan)
    val(nsubint)
    val(nbin)

    output:
    path "*.png"

    script:
    """
    output_filename=\$(basename "${fold_archive}" | sed 's/\\.[^.]*\$//')
    pdmp -mc ${nchan} -ms ${nsubint} -mb ${nbin} -g \${output_filename}.png/PNG ${fold_archive}
    """

}

process pdmp_ptuse {
    label 'pdmp' 
    container "${params.fold_singularity_image}"
    publishDir "TIMING/${target_name}/01_PTUSE_FOLDS/", pattern: "*.png", mode: 'copy'

    input:
    path(fold_archive)
    val(target_name)
    val(nchan)
    val(nsubint)
    val(nbin)

    output:
    path "*.png"

    script:
    """
    output_filename=\$(basename "${fold_archive}" | sed 's/\\.[^.]*\$//')
    pdmp -mc ${nchan} -ms ${nsubint} -mb ${nbin} -g \${output_filename}.png/PNG ${fold_archive}
    """

}

process query_db {
    label 'query_db'
    container "${params.sql_query_image}"
    input:
    val(target_name)
    val(beam_name)
    val(utc_start)
    val(utc_end)

    output:
    path "grouped_files.txt"

    script:
    def beamArg = beam_name ? "-b ${beam_name}" : ""
    def utc_arg = utc_start ? "-us ${utc_start} -ue ${utc_end}" : ""
    """
    python ${params.db_query} -t ${target_name} ${beamArg} ${utc_arg}
    """

}

workflow {

    grouped_query_output = query_db(params.target_name, params.beam_name, params.utc_start, params.utc_end)
    processed_data_channel = grouped_query_output.splitText()
                                  .toList()
                                  .flatMap { it.toList()[1..-1] } // Skip the first line (header)
                                  .map { line ->
                                      def parts = line.split('\t')

                                      def pointing_id = parts[0]
                                      def beam_id = parts[1]
                                      def filterbank_files = parts[2]
                                      def target = parts[3]
                                      def beam_num = parts[4].replace(" ", "")
                                      // Replace space with underscore in utc_start
                                      def utc_start = parts[5].replace(" ", "-")


                                      return tuple(pointing_id, beam_id, filterbank_files, target, beam_num, utc_start)
         
                                  }
    utc_current = processed_data_channel.map { tuple -> tuple[5].trim() }
    
    filtool_output = filtool(processed_data_channel, params.rfi_filter, params.threads, params.telescope)
    nearest_two_output = nearest_power_of_two_calculator(filtool_output)
    fft_size_value = nearest_two_output.map{ file -> file.text.trim() }

    peasoup_output = peasoup(filtool_output, params.dm_file, fft_size_value, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.target_name, utc_current)
    phase_predictor_output = create_phase_predictor(peasoup_output, params.period_start, params.period_end, params.target_name)
    if (params.APSUSE_FOLDS == 1) {

        apsuse_folds = fold_apsuse(filtool_output, phase_predictor_output, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins, params.target_name)
        if (params.use_clfd == 1) {
            clfd_output = clfd_apsuse(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            pdmp_output = pdmp_apsuse(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
        else {
            pdmp_output = pdmp_apsuse(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
        }

    }
    if (params.PTUSE_FOLDS == 1) {

        processed_data_channel.branch { tuple ->
        def year = tuple[5].split("-")[0] as int
        before2023: year < 2023
        after2023: year >= 2023
        }.set { yearBranch }

        // Before 2023 PTUSE data is kept in a different directory.
        ptuse_data_before2023 = yearBranch.before2023.map { tuple ->
        def parts = tuple[5].split(/-|:/)
        def dateWithHour = parts[0..3].take(4).join("-")
        def psrfits_file = "${params.PTUSE1}/${dateWithHour}*/${tuple[3]}/**/*.sf"

        return [tuple[3].trim(), tuple[5].trim(), psrfits_file]
    }
        // After 2023 case
        yearBranch.after2023.map { tuple ->
        def parts = tuple[5].split(/-|:/) // Splitting by both hyphens and colons
        def dateWithHour = parts[0..3].take(4).join("-") // Now takes only the first 4 elements
        def psrfits_file = "${params.PTUSE2}/${dateWithHour}*/${tuple[3]}/**/*.sf"

        // Create a new tuple that includes the psrfits_file
            return [tuple[3].trim(), tuple[5].trim(), psrfits_file]
        }.set { ptuse_data_after2023 }

        // Now combine both channels

        all_ptuse_data = ptuse_data_before2023.mix(ptuse_data_after2023)

        ptuse_folds = fold_ptuse(all_ptuse_data, params.ephemeris_file, params.dspsr_ptuse_threads, params.telescope, params.dspsr_ptuse_subint_length, params.dspsr_ptuse_bins)
        if (params.use_clfd == 1) {
            clfd_ptuse_output = clfd_ptuse(ptuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            pdmp_output = pdmp_ptuse(clfd_ptuse_output, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
        else {
            pdmp_output = pdmp_ptuse(ptuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
    }


}

