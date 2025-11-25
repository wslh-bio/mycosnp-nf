include { PRE_MYCOSNP_WF } from '../workflows/pre_mycosnp'
include { MYCOSNP } from '../workflows/mycosnp'

include { WSLH_SUMMARY } from '../modules/local/wslh_summary'

workflow BOTH {

    //
    // WORKFLOW: Run pre-mycosnp pipeline
    //
    PRE_MYCOSNP_WF ()

    //
    // WORKFLOW: Run main nf-core/mycosnp analysis pipeline
    //
    MYCOSNP ()

    //
    // MODULE: Generate WSLH specific summary report
    
    WSLH_SUMMARY (
        MYCOSNP.out.qc_stats,
        params.runname,
        MYCOSNP.out.fks1_combined,
        PRE_MYCOSNP_WF.out.wslh_results
    )
}