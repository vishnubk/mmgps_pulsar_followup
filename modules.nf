process filtool {
    label 'filtool'
    container "${params.pulsarx_singularity_image}"

    input:
    val(input_data) 
    val(metadata)
    val rfi_filter
    val threads
    val telescope

    output:
    path("*.fil")

    script:
    def outputFile = "${metadata[3].trim()}_${metadata[5].trim()}_${metadata[4].trim()}"
    """
    # Get the first file from the input_data string
    first_file=\$(echo ${input_data} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ \${file_extension} == "sf" ]]; then
        filtool -psrfits --scloffs -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${input_data} -s ${metadata[3]}
    else 
        filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${input_data} -s ${metadata[3]}
    fi

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
    val(input_file) //filterbank file used for grouping later.

    output:

    tuple (val(input_file), path("predictor_candidate_*"), optional: true)


    script:
    """
    python ${params.phase_predictor} -i ${xml_file} -ps ${period_start} -pe ${period_end} -t ${target_name}
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

process dspsr_fold_phase_predictor {
    label 'fold_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    tuple val(fil_file), path(phase_predictor)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)
    val(target_name)

    output:
    path "*.ar"

    script:
    """
    # Extract the basename of the filterbank file
    output_filename=\$(basename ${fil_file} .fil)

    # Extract the basename of the phase_predictor file
    phase_predictor_basename=\$(basename ${phase_predictor})

    # Extract the candidate number from the phase_predictor basename
    candidate_id=\$(echo \${phase_predictor_basename} | sed 's/predictor_candidate_//')

    # Append the candidate_id to the output filename
    output_filename_with_candidate_id=\${output_filename}_candidate_id_\${candidate_id}

    # Run dspsr with the modified output filename
    dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename_with_candidate_id} ${fil_file}
    """
}


process dspsr_fold_ephemeris {
    label 'fold_ephemeris' 
    container "${params.fold_singularity_image}"

    input:
    val(target_name) 
    val(utc)
    val(file_path)
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
    found_file=\$(ls -v ${file_path} 2> /dev/null | head -n 1)

    if [[ -z "\$found_file" ]]; then
        echo "Error: No matching files found."
        exit 1
    fi

    # Extracting basename and filename without extension from ephemeris_file
    ephemeris_basename=\$(basename ${ephemeris_file})
    ephemeris_name_no_ext=\${ephemeris_basename%.*}

    # Construct the output file name with ephemeris file suffix
    output_filename="${target_name}_${utc}_\${ephemeris_name_no_ext}"

    # Check the file extension of the found file
    file_extension="\${found_file##*.}"
    echo \${file_extension}

    # Run the dspsr command
    if [[ \${file_extension} == "sf" ]]; then
        dspsr -scloffs -k ${telescope} -t${threads} -E ${ephemeris_file} -L ${subint_length} -b${bins} -A -O \${output_filename} \$(ls -v ${file_path})
    else
        dspsr -k ${telescope} -t${threads} -E ${ephemeris_file} -L ${subint_length} -b${bins} -A -O \${output_filename} \$(ls -v ${file_path})
    fi
    """
}


process clfd {
    label 'clfd' 
    container "${params.fold_singularity_image}"

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

process pdmp {
    label 'pdmp' 
    container "${params.fold_singularity_image}"

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
    output_filename=\$(basename ${fold_archive} | sed 's/\\.[^.]*\$//')
    pdmp -mc ${nchan} -ms ${nsubint} -mb ${nbin} -g \${output_filename}.png/PNG ${fold_archive}
    """

}

process digifil {
    label 'digifil'
    container "${params.digifil_singularity_image}"

    input:
    val(input_data)
    val(metadata)
    val(nbits)
    val(nthreads)
    val(time_decimate_factor)

    output:
    path "*.fil"

    script:
    def outputFile = "${metadata[3].trim()}_${metadata[5].trim()}_${metadata[4].trim()}.fil"
    """
    # Get the first file from the input_data string
    first_file=\$(echo ${input_data} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"

    # Determine the -scloffs flag based on file extension
    scloffs_flag=""
    if [[ "\${file_extension}" == "sf" ]]; then
        scloffs_flag="-scloffs"
    fi

    # Conditional execution of digifil based on time_decimate_factor
    if [ ${time_decimate_factor} -eq 1 ]; then
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -d 1 -o ${outputFile} \$(ls -v ${input_data})
    else
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -t ${time_decimate_factor} -d 1 -o ${outputFile} \$(ls -v ${input_data})
    fi
    """
}



