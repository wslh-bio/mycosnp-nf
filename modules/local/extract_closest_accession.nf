process EXTRACT_CLOSEST_ACCESSION {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    tuple val(meta), path(gambit)
    path mapping_reference_file 

    output:
    tuple val(meta), env(closest_accession), env(s3_uri), emit: accession_and_uri

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract closest accession from gambit file
    closest = \$(cat "${gambit}" | grep -v 'closest.description' | cut -f 7 -d ',')
    closest_accession = \$(echo "\${closest}" | cut -f 1 -d ' ' | tr -d '[] \\t\\n\\r')
    
    # Look up the S3 URI from the mapping file
    if [[ -n "\${closest_accession}" ]]; then
        s3_uri=\$(awk -v acc="\${closest_accession}" '\$1 == acc {print \$2}' ${mapping_reference_file})
        
        if [[ -z "\${s3_uri}" ]]; then
            echo "Warning: No S3 URI found for accession \${closest_accession}" >&2
            s3_uri="NONE"
        fi
    else
        closest_accession="NONE"
        s3_uri="NONE"
    fi
    
    echo "Found accession: \${closest_accession}" >&2
    echo "Mapped to S3 URI: \${s3_uri}" >&2
    """

    stub:
    """
    # Testing
    closest_accession="GCF_000001234"
    s3_uri="s3://test-bucket/genomes/GCF_000001234.00_ABCh12.P12_genomic.fna"
    
    echo "STUB MODE: Mock accession: \${closest_accession}" >&2
    echo "STUB MODE: Mock S3 URI: \${s3_uri}" >&2
    """
}
