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

process peasoup {
    label 'peasoup'
    container "${params.search_singularity_image}"
    // This will only publish the XML files
    publishDir "RESULTS/${POINTING}/${UTC_OBS}/${BAND}/${BEAM}/03_SEARCH/", pattern: "**/*.xml", mode: 'copy'



    input:
    tuple path(fil_file), val(POINTING), val(BAND), val(UTC_OBS), val(BEAM)
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
    tuple path("**/*.xml"), path(fil_file), val(POINTING), val(BAND), val(UTC_OBS), val(BEAM)



    script:
    """
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus} 

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
    filtool(processed_data_channel, params.rfi_filter, params.threads, params.telescope)
}

