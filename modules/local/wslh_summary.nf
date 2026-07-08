process WSLH_SUMMARY {
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path qc_stats
    val runname
    path fks1_combined
    path clade_designation

    output:
    path("*_qc_report.csv"), emit: qc_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mycosnp_results_summary.py \\
    -qc ${qc_stats} \\
    -r ${runname} \\
    -f ${fks1_combined} \\
    -c ${clade_designation} \\
    -wv ${workflow.version} \\
    """
}