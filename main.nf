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

include { dspsr_fold_ephemeris as apsuse_fold_ephemeris } from './modules'
//include { dspsr_fold_ephemeris as ptuse_fold_ephemeris } from './modules'
include { dspsr_fold_ephemeris_ptuse as ptuse_fold_ephemeris } from './modules'

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

workflow {

    grouped_query_output = Channel.fromPath("${params.filterbank_list}")
    // Selecting filterbank_files, target, beam_num, utc_start
    filterbank_channel_with_metadata = grouped_query_output.splitText()
                                  .toList()
                                  .flatMap { it.toList()[1..-1] } // Skip the first line (header)
                                  .map { line ->
                                      def parts = line.split('\t')

                                    //   def pointing_id = parts[0]
                                    //   def beam_id = parts[1]
                                      def filterbank_files = parts[2]
                                      def target = parts[3]
                                      def beam_num = parts[4].replace(" ", "")
                                      // Replace space with underscore in utc_start
                                      def utc_start = parts[5].replace(" ", "-")


                                      return tuple(filterbank_files, target, beam_num, utc_start)
         
                                  }

    if (params.APSUSE_SEARCH == 1 || params.APSUSE_EPH_FOLD == 1) {

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
        // Trim the string values
        target = target.trim()
        beam_name = beam_name.trim()
        utc_start = utc_start.trim()
        return tuple(filepath, target, beam_name, utc_start)
                                                              }
       
        if (params.APSUSE_SEARCH == 1) {
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
                    } else {
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

       
        }}
    
        
        // User asked for APSUSE_EPH_FOLD.
        if (params.APSUSE_EPH_FOLD == 1) {
            all_par_files_apsuse = Channel.fromPath("${params.ephemeris_files_dir}/*.par")
            
            // Combine par file and filterbank file channel using a cartesian product. Each par file will apply on all filterbank files.
            combined_channel_apuse_eph_fold = updated_filterbank_channel.combine(all_par_files_apsuse)  

            apsuse_folds = apsuse_fold_ephemeris(combined_channel_apuse_eph_fold, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)

            if (params.use_clfd == 1) {
                clfd_output = clfd_apsuse_eph(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_apsuse_eph(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_apsuse_eph(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
            }

        }
    

    if (params.PTUSE_SEARCH == 1 || params.PTUSE_EPH_FOLD == 1) {

        filterbank_channel_with_metadata.branch { tuple ->
        def year = tuple[3].split("-")[0] as int
        before2023: year < 2023
        after2023: year >= 2023
        }.set { yearBranch }

        // Before 2023 PTUSE data is kept in a different directory.
        ptuse_data_before2023 = yearBranch.before2023.map { tuple ->
        def parts = tuple[3].split(/-|:/)
        def dateWithHour = parts[0..3].take(4).join("-")
        def psrfits_file = "${params.PTUSE1}/${dateWithHour}*/${params.target_name}/**/*.sf"

        return [psrfits_file, tuple[1].trim(), "ptuse", dateWithHour]
    }
        // After 2023 case
        yearBranch.after2023.map { tuple ->
        def parts = tuple[3].split(/-|:/) // Splitting by both hyphens and colons
        def dateWithHour = parts[0..3].take(4).join("-") // Now takes only the first 4 elements
        def psrfits_file = "${params.PTUSE2}/${dateWithHour}*/${params.target_name}/**/*.sf"

        // Create a new tuple that includes the psrfits_file
            return [psrfits_file, tuple[1].trim(), "ptuse", dateWithHour]
        }.set { ptuse_data_after2023 }

        // Now combine both channels
        
        all_ptuse_data = ptuse_data_before2023.mix(ptuse_data_after2023)
        // PTUSE EPH FOLD CASE
        if (params.PTUSE_EPH_FOLD == 1)

        {
            all_par_files_ptuse = Channel.fromPath("${params.ephemeris_files_dir}/*.par")
            // Combine par file and filterbank file channel using a cartesian product. Each par file will apply on all filterbank files.
            combined_channel_ptuse_eph_fold = all_ptuse_data.combine(all_par_files_ptuse) 
            //combined_channel_ptuse_eph_fold.view() 
            ptuse_folds_eph = ptuse_fold_ephemeris(combined_channel_ptuse_eph_fold, params.dspsr_ptuse_threads, params.telescope, params.dspsr_ptuse_subint_length, params.dspsr_ptuse_bins)
            if (params.use_clfd == 1) {
                clfd_output_ptuse_eph = clfd_ptuse_eph(ptuse_folds_eph, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_ptuse_eph(clfd_output_ptuse_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_ptuse_eph(ptuse_folds_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
            }


        }

        if (params.PTUSE_SEARCH == 1) {

        if (params.use_filtool_ptuse == 1){
            merged_filterbank_ptuse = filtool_ptuse(all_ptuse_data, params.filtool_rfi_filter, params.filtool_threads, params.telescope)
            //merged_filterbank_ptuse = digifil_ptuse(ptuse_data, processed_data_channel, params.nbits, params.digifil_threads, params.digifil_decimate_time_ptuse)
        }
        else {
            merged_filterbank_ptuse = digifil_ptuse(all_ptuse_data, params.nbits, params.digifil_threads, params.digifil_decimate_time_ptuse)

        }
    
    }
    merged_filterbank_ptuse.view()




        }
        
        


 
    



     

}

