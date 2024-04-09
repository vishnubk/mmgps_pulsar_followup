#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool_apsuse } from './modules'
include { filtool_ptuse } from './modules'

include { digifil as digifil_apsuse } from './modules'
include { digifil_ptuse } from './modules'



include { create_phase_predictor as create_phase_predictor_apsuse } from './modules'
include { create_phase_predictor as create_phase_predictor_ptuse } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'
include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_ptuse } from './modules'


//include { dspsr_fold_phase_predictor as apsuse_fold_phase_predictor } from './modules'
include { dspsr_fold_phase_predictor_parallel as apsuse_fold_phase_predictor_parallel } from './modules'
include { dspsr_fold_phase_predictor_serial as apsuse_fold_phase_predictor_serial } from './modules'

include { dspsr_fold_phase_predictor as ptuse_fold_phase_predictor } from './modules'
include { dspsr_fold_ephemeris_apsuse as apsuse_fold_ephemeris } from './modules'
include { dspsr_fold_ephemeris_ptuse_updated as ptuse_fold_ephemeris } from './modules'

include { clfd as clfd_apsuse_predictor } from './modules'
include { clfd as clfd_apsuse_eph } from './modules'
include { clfd as clfd_ptuse_predictor } from './modules'
include { clfd as clfd_ptuse_eph } from './modules'

include { pdmp as pdmp_apsuse_predictor } from './modules'
include { pdmp as pdmp_apsuse_eph } from './modules'
include { pdmp as pdmp_ptuse_predictor } from './modules'
include { pdmp as pdmp_ptuse_eph } from './modules'



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

process create_ptuse_paths {
    label 'query_db'
    container "${params.sql_query_image}"
    input:
    path(filename)

    output:
    path "ptuse_data_paths.txt"

    script:
    
    """
    python ${params.ptuse_path_creater} -f ${filename}
    """

}

process peasoup_apsuse {
    label 'peasoup'
    container "${params.search_singularity_image}"
    publishDir "SEARCH/${params.target_name}/APSUSE/${utc}/${beam_name}/", pattern: "**/*.xml", mode: 'copy'

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size)
    path(dm_file) 
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    

    output:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size), path("**/*.xml")

    script:
    """
    #!/bin/bash
   
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus}

    """
}

process peasoup_ptuse {
    label 'peasoup'
    container "${params.search_singularity_image}"
    publishDir "SEARCH/${params.target_name}/PTUSE/${utc}/${beam_name}/", pattern: "**/*.xml", mode: 'copy'
    

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size)
    path(dm_file) 
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

process FIND_PTUSE_DATA {
    // Inputs
    input:
    val target_name

    // Outputs
    output:
    path "ptuse_data.csv"

    script:
    """
    #!/bin/bash
    echo "PTUSE_PATH,target,beam,utc" > ptuse_data.csv
    find /beegfs/DATA/{MeerTIME,PTUSE}/SCI-20200703-MK-{01,02}/search -name "${target_name}" | while read line; do
        ptuse_utc=\$(echo \$line | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{2}:[0-9]{2}:[0-9]{2}')
        echo "\${line}/**/*.sf,${target_name},ptuse,\${ptuse_utc}" >> ptuse_data.csv
    done
    """
}

workflow {

    if (params.APSUSE_SEARCH == 1 || params.APSUSE_EPH_FOLD == 1 ) {

    grouped_query_output = query_db(params.target_name, params.beam_name, params.utc_start, params.utc_end)
    
    // Selecting filterbank_files, target, beam_num, utc_start
    filterbank_channel_with_metadata = grouped_query_output.splitText()
                                  .toList()
                                  .flatMap { it.toList()[1..-1] } // Skip the first line (header)
                                  .map { line ->
                                      def parts = line.split('\t')

                                      def filterbank_files = parts[2]
                                      def target = parts[3]
                                      def beam_num = parts[4].replace(" ", "")
                                      // Replace space with underscore in utc_start
                                      def utc_start = parts[5].replace(" ", "-")
                                      target = target.trim()
                                      beam_num = beam_num.trim()
                                      utc_start = utc_start.trim()
                                      return tuple(filterbank_files, target, beam_num, utc_start)
         
                                  }
        
    }
    if (params.APSUSE_SEARCH == 1) {

        if (params.use_filtool_apsuse == 1){
        
            processed_filterbank = filtool_apsuse(filterbank_channel_with_metadata, params.filtool_rfi_filter, params.filtool_threads, params.telescope)
        }
        else
        {
            processed_filterbank = digifil_apsuse(filterbank_channel_with_metadata, params.nbits, params.digifil_threads, params.digifil_decimate_time_apsuse)

        }
        //Add the filterbank file from filtool/digifil along with metadata to a single channel.
        updated_filterbank_channel = processed_filterbank.map { metadata, filepath ->
        def (raw_filterbank, target, beam_name, utc_start) = metadata
        return tuple(filepath, target, beam_name, utc_start)
                                                              }
       
        nearest_two_output = nearest_power_of_two_calculator_apsuse(updated_filterbank_channel)
        peasoup_output = peasoup_apsuse(nearest_two_output, params.dm_file, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus)
        phase_predictor_output = create_phase_predictor_apsuse(peasoup_output, params.period_start, params.period_end)

        if (params.parallel_fold == 1) {
            restructured_phase_predictor_output = phase_predictor_output.flatMap { fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictors ->
                // Check if phase_predictors is a list (or tuple)
                if (phase_predictors instanceof List || phase_predictors.getClass().isArray()) {
                    phase_predictors.collect { phase_predictor ->
                        return [fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictor]
                    }
                } 
                else {
                    // If it's a single element, return it as it is
                    return [[fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictors]]
                }
            }
            
            apsuse_folds = apsuse_fold_phase_predictor_parallel(restructured_phase_predictor_output, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)
        }
        //Folding candidates in serial. Use this when you have only a few spin period cands.
        else{

            apsuse_folds = apsuse_fold_phase_predictor_serial(phase_predictor_output, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)
        }

        if (params.use_clfd == 1) {
            clfd_output = clfd_apsuse_predictor(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            pdmp_output = pdmp_apsuse_predictor(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
        else {
            pdmp_output = pdmp_apsuse_predictor(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
        }

    
    }
       // User asked for APSUSE_EPH_FOLD.
    if (params.APSUSE_EPH_FOLD == 1) {
        all_par_files_apsuse = Channel.fromPath("${params.ephemeris_files_dir}/*.par")
        // Combine par file and filterbank file channel using a cartesian product. Each par file will apply on all filterbank files.
        combined_channel_apuse_eph_fold = filterbank_channel_with_metadata.combine(all_par_files_apsuse)  

        apsuse_folds = apsuse_fold_ephemeris(combined_channel_apuse_eph_fold, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)

        if (params.use_clfd == 1) {
            clfd_output = clfd_apsuse_eph(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            pdmp_output = pdmp_apsuse_eph(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
        else {
            pdmp_output = pdmp_apsuse_eph(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
        }

    }

    if (params.PTUSE_EPH_FOLD == 1) {

        ptuse_path_creater = FIND_PTUSE_DATA(params.target_name)
        ptuse_data = ptuse_path_creater.splitText()
                                  .toList()
                                  .flatMap { it.toList()[1..-1] } // Skip the first line (header)
                                  .map { line ->
                                      def parts = line.split(',')
                                      def input_data = parts[0]
                                      def target = parts[1].trim()
                                      def beam_name = parts[2].trim()
                                      def ptuse_utc_start = parts[3].trim()
                                      
                                      return tuple(input_data, target, beam_name, ptuse_utc_start)
         
                                  }
        
        all_par_files_ptuse = Channel.fromPath("${params.ephemeris_files_dir}/*.par")
        
        // Combine par file and filterbank file channel using a cartesian product. Each par file will apply on all filterbank files.
        combined_channel_ptuse_eph_fold = ptuse_data.combine(all_par_files_ptuse)
        ptuse_folds_eph = ptuse_fold_ephemeris(combined_channel_ptuse_eph_fold, params.dspsr_ptuse_threads, params.telescope, params.dspsr_ptuse_subint_length, params.dspsr_ptuse_bins)
        if (params.use_clfd == 1) {
            clfd_output_ptuse_eph = clfd_ptuse_eph(ptuse_folds_eph, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            pdmp_output = pdmp_ptuse_eph(clfd_output_ptuse_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
        }
        else {
            pdmp_output = pdmp_ptuse_eph(ptuse_folds_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
        }


         }
        
}