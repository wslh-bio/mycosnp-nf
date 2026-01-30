process REPORT_REJECTED_SAMPLES {
    tag "rejected_samples"
    label 'process_single'

    publishDir "${params.outdir}/rejected_samples", mode: params.publish_dir_mode

    input:
    path csv_file
    val workflow

    output:
    path "*empty_samples.csv"

    script:
    """
    cp ${csv_file} ${workflow}_empty_samples.csv
    """
}