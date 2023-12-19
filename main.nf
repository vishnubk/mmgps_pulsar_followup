#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool_apsuse } from './modules'
include { filtool as filtool_ptuse } from './modules'

include { digifil as digifil_apsuse } from './modules'
include { digifil as digifil_ptuse } from './modules'

include { peasoup as peasoup_apsuse } from './modules'
include { peasoup as peasoup_ptuse } from './modules'

include { create_phase_predictor as create_phase_predictor_apsuse } from './modules'
include { create_phase_predictor as create_phase_predictor_ptuse } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'
include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_ptuse } from './modules'


include { dspsr_fold_phase_predictor as apsuse_fold_phase_predictor } from './modules'
include { dspsr_fold_phase_predictor as ptuse_fold_phase_predictor } from './modules'
include { dspsr_fold_ephemeris as apsuse_fold_ephemeris } from './modules'
include { dspsr_fold_ephemeris as ptuse_fold_ephemeris } from './modules'

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
    raw_data = processed_data_channel.map { tuple -> tuple[2].trim() }
    if (params.APSUSE_SEARCH == 1 || params.APSUSE_EPH_FOLD == 1) {

        if (params.use_filtool_apsuse == 1){
        
            merged_filterbank = filtool_apsuse(raw_data, processed_data_channel, params.filtool_rfi_filter, params.filtool_threads, params.telescope)
        }
        else
        {
            merged_filterbank = digifil_apsuse(raw_data, processed_data_channel, params.nbits, params.digifil_threads, params.digifil_decimate_time_apsuse)

        }
        utc_current = processed_data_channel.map { tuple -> tuple[5].trim() }

        if (params.APSUSE_SEARCH == 1) {
            nearest_two_output = nearest_power_of_two_calculator_apsuse(merged_filterbank)
            fft_size_value = nearest_two_output.map{ file -> file.text.trim() }
            peasoup_output = peasoup_apsuse(merged_filterbank, params.dm_file, fft_size_value, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.target_name, utc_current)
            phase_predictor_output = create_phase_predictor_apsuse(peasoup_output, params.period_start, params.period_end, params.target_name)
            collected_phase_predictors = phase_predictor_output.collect()
            apsuse_folds = apsuse_fold_phase_predictor(merged_filterbank, collected_phase_predictors, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins, params.target_name)

            if (params.use_clfd == 1) {
                clfd_output = clfd_apsuse_predictor(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_apsuse_predictor(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_apsuse_predictor(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
            }

        }
        // User asked for APSUSE_EPH_FOLD.
        else {
        
            apsuse_folds = apsuse_fold_ephemeris(params.target_name,utc_current, merged_filterbank, params.ephemeris_file, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)
            if (params.use_clfd == 1) {
                clfd_output = clfd_apsuse_eph(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_apsuse_eph(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_apsuse_eph(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
            }

        }
    }

    if (params.PTUSE_SEARCH == 1 || params.PTUSE_EPH_FOLD == 1) {

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

        return [dateWithHour, psrfits_file]
    }
        // After 2023 case
        yearBranch.after2023.map { tuple ->
        def parts = tuple[5].split(/-|:/) // Splitting by both hyphens and colons
        def dateWithHour = parts[0..3].take(4).join("-") // Now takes only the first 4 elements
        def psrfits_file = "${params.PTUSE2}/${dateWithHour}*/${tuple[3]}/**/*.sf"

        // Create a new tuple that includes the psrfits_file
            return [dateWithHour, psrfits_file]
        }.set { ptuse_data_after2023 }

        // Now combine both channels

        all_ptuse_data = ptuse_data_before2023.mix(ptuse_data_after2023)
        ptuse_utc = all_ptuse_data.map { tuple -> tuple[0].trim() }
        ptuse_data = all_ptuse_data.map { tuple -> tuple[1].trim() }

        if (params.use_filtool_ptuse == 1){
            // Temporarily running filtool on PTUSE data is disabled until the bug is resolved. So we run digifil and no cleaning in both cases!
            //merged_filterbank_ptuse = filtool_ptuse(ptuse_data, processed_data_channel, params.filtool_rfi_filter, params.filtool_threads, params.telescope)
            merged_filterbank_ptuse = digifil_ptuse(ptuse_data, processed_data_channel, params.nbits, params.digifil_threads, params.digifil_decimate_time_ptuse)
        }
        else {
            merged_filterbank_ptuse = digifil_ptuse(ptuse_data, processed_data_channel, params.nbits, params.digifil_threads, params.digifil_decimate_time_ptuse)

        }
        if (params.PTUSE_SEARCH == 1){
    
        
            nearest_two_output_ptuse = nearest_power_of_two_calculator_ptuse(merged_filterbank_ptuse)
            fft_size_value_ptuse = nearest_two_output_ptuse.map{ file -> file.text.trim() }
            peasoup_output_ptuse = peasoup_ptuse(merged_filterbank_ptuse, params.dm_file, fft_size_value_ptuse, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.target_name, ptuse_utc)
            phase_predictor_output_ptuse = create_phase_predictor_ptuse(peasoup_output_ptuse, params.period_start, params.period_end, params.target_name)
            collected_phase_predictors_ptuse = phase_predictor_output_ptuse.collect()

            ptuse_folds_pred = ptuse_fold_phase_predictor(merged_filterbank_ptuse, collected_phase_predictors_ptuse, params.dspsr_ptuse_threads, params.telescope, params.dspsr_ptuse_subint_length, params.dspsr_ptuse_bins, params.target_name)

            if (params.use_clfd == 1) {
                clfd_output_ptuse_pred = clfd_ptuse_predictor(ptuse_folds_pred, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_ptuse_predictor(clfd_output_ptuse_pred, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_ptuse_predictor(ptuse_folds_pred, params.target_name, params.nchan, params.nsubint, params.nbins)
            }

          

        }
        // PTUSE EPH FOLD CASE
        else 

        {
        
            ptuse_folds_eph = ptuse_fold_ephemeris(params.target_name, ptuse_utc, merged_filterbank_ptuse, params.ephemeris_file, params.dspsr_ptuse_threads, params.telescope, params.dspsr_ptuse_subint_length, params.dspsr_ptuse_bins)
            if (params.use_clfd == 1) {
                clfd_output_ptuse_eph = clfd_ptuse_eph(ptuse_folds_eph, params.target_name, params.qmask, params.qspike, params.clfd_processes)
                pdmp_output = pdmp_ptuse_eph(clfd_output_ptuse_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
            }
            else {
                pdmp_output = pdmp_ptuse_eph(ptuse_folds_eph, params.target_name, params.nchan, params.nsubint, params.nbins)
            }


        }
    



    }

}

