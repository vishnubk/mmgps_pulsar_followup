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

    output:
    path("**/*.xml")

    script:
    """
    #!/bin/bash
   
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus}

    """
}

process create_phase_predictor {
    label 'create_phase_predictor' // Replace with appropriate label
    container "${params.fold_singularity_image}"
    //container "appropriateContainer" // Replace with the appropriate container if needed

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

    input:
    path(fil_file)
    path(phase_predictor)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    // output:
    // path "predictor_candidate_*"

    script:
    """
    output_filename=\$(basename ${fil_file} .fil)
    dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename} ${fil_file}
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

    // Create channels for each line in the grouped output
    //grouped_query_output.splitText().view()
    filtool_output = filtool(processed_data_channel, params.rfi_filter, params.threads, params.telescope)
    nearest_two_output = nearest_power_of_two_calculator(filtool_output)
    fft_size_value = nearest_two_output.map{ file -> file.text.trim() }

    peasoup_output = peasoup(filtool_output, params.dm_file, fft_size_value, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus)
    phase_predictor_output = create_phase_predictor(peasoup_output, params.period_start, params.period_end, params.target_name)
    apsuse_folds = fold_apsuse(filtool_output, phase_predictor_output, params.dspsr_threads, params.telescope, params.dspsr_subint_length, params.dspsr_apsuse_bins)
}

