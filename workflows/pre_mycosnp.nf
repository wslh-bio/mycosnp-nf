/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ] // params.snpeffdb
if (params.skip_samples_file) { // check for skip_samples_file
    checkPathParamList.add(params.skip_samples_file)
}
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
sra_list = []
sra_ids = [:]
ch_input = null;
if (params.input) 
{ 
    ch_input = file(params.input) 
}
if(params.add_sra_file)
{
    sra_file = file(params.add_sra_file, checkIfExists: true)
    allLines  = sra_file.readLines()
    for( line : allLines ) 
    {
        row = line.split(',')
        if(row.size() > 1)
        {
            println "Add SRA ${row[1]} => ${row[0]}"
            sra_list.add(row[1])
            sra_ids[row[1]] = row[0]
        } else
        {
            if(row[0] != "")
            {
                println " ${row[0]} => ${row[0]}"
                sra_list.add(row[0])
                sra_ids[row[0]] = row[0]
            }
        }
    }
}

vcf_file_list = []
vcfidx_file_list = []
if(params.add_vcf_file)
{
    vcf_file = file(params.add_vcf_file, checkIfExists: true)
    allLines  = vcf_file.readLines()
    for( line : allLines ) 
    {
        if(line != "")
        {
            println " Add VCF => $line"
            t_vcf = file(line)
            t_idx = file(line + ".tbi")
            vcf_file_list.add(t_vcf)
            vcfidx_file_list.add(t_idx)
        }
    }
}

if(! ( params.input || params.add_sra_file || params.add_vcf_file ) ) { exit 1, 'Input samplesheet, sra file, or vcf file not specified!' }


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { SRA_FASTQ_SRATOOLS       } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { SEQKIT_PAIR              } from '../modules/nf-core/modules/seqkit/pair/main'
include { FAQCS                    } from '../modules/nf-core/modules/faqcs/main'
include { GAMBIT_QUERY             } from '../modules/local/gambit'
include { SUBTYPE                  } from '../modules/local/subtype'
include { EXTRACT_CLOSEST_ACCESSION} from '../modules/local/extract_closest_accession'
include { PRE_MYCOSNP_INDV_SUMMARY } from '../modules/local/pre_mycosnp_indv_summary'
include { PRE_MYCOSNP_COMB_SUMMARY } from '../modules/local/pre_mycosnp_comb_summary'
include { REJECTED_SAMPLES         } from '../modules/local/rejected_samples'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/modules/fastqc/main'
include { SHOVILL as SHOVILL          } from '../modules/nf-core/modules/shovill/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []


workflow PRE_MYCOSNP_WF {

    ch_versions = Channel.empty()

    // Create empty channels for reference files required byt the main workflow
    fas_file = Channel.empty()
    fai_file = Channel.empty()
    bai_file = Channel.empty()
    dict_file = Channel.empty()
    meta_val = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_all_reads = Channel.empty()
    ch_sra_reads = Channel.empty()
    ch_sra_list  = Channel.empty()
    if(params.add_sra_file)
    {   
        ch_sra_list = Channel.fromList(sra_list)
                             .map{valid -> [ ['id':sra_ids[valid],single_end:false], valid ]}
        SRA_FASTQ_SRATOOLS(ch_sra_list)
        ch_all_reads = ch_all_reads.mix(SRA_FASTQ_SRATOOLS.out.reads)
    }
    
    if(params.input)
    {
        INPUT_CHECK (
            ch_input
        )

        INPUT_CHECK.out.reads
            .branch{ meta, file -> 
                single_end: meta.single_end
                paired_end: !meta.single_end
                }
            .set{ ch_filtered }

        ch_filtered.paired_end
            .map{ meta, file ->
                [meta, file, file[0].countFastq(), file[1].countFastq()]}
            .branch{ meta, file, count1, count2 ->
                pass: count1 > 0 && count2 > 0
                fail: count1 == 0 || count2 == 0 || count1 == 0 && count2 == 0
            }
            .set{ ch_paired_end }

        ch_paired_end.pass
            .map { meta, file, count1, count2 -> 
                [meta, file]
                }
            .set{ ch_filtered }

        ch_paired_end.fail
            .map { meta, file, count1, count2 ->
                [meta.id]
                }
            .set{ ch_paired_end_fail }

        ch_paired_end_fail
            .flatten()
            .set{ ch_failed }

        ch_failed
            .ifEmpty ('NO_EMPTY_SAMPLES')
            .collectFile(
                name: 'empty_samples.csv',
                newLine: true
            )
            .set{ ch_rejected_file }

        REJECTED_SAMPLES (
                ch_rejected_file,
                "Pre_mycosnp"
            )

        ch_all_reads = ch_all_reads.mix(ch_filtered)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    //
    // MODULE: Run Pre-FastQC 
    //
    if (params.workflow != "BOTH") {
        FASTQC_RAW (
            ch_all_reads
        )
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)
    }
    //
    // MODULE: Run seqkit to remove unpaired reads
    //
    SEQKIT_PAIR(
        ch_all_reads
    )
    ch_versions = ch_versions.mix(SEQKIT_PAIR.out.versions)

    //
    // MODULE: Run FAQCs - no downsampling option because a reference cannot be supplied before knowing the species
    //
    FAQCS(
        SEQKIT_PAIR.out.reads
    )
    ch_versions = ch_versions.mix(FAQCS.out.versions)

    //
    // MODULE: Run Shovill
    //
    SHOVILL (
        FAQCS.out.reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions)

    //
    // MODULE: Run Gambit
    //
    GAMBIT_QUERY(
        SHOVILL.out.contigs,
        params.gambit_db,
        params.gambit_h5_dir
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions)

    //
    // MODULE: Subtype
    //

    // Join the GAMBIT output and the spades assembly into a single channel   
    SHOVILL.out.contigs.map{ meta, contigs -> [meta, contigs] }.set{ ch_contigs }
    GAMBIT_QUERY.out.taxa.map{ meta, gambit -> [meta, gambit] }.join(ch_contigs).set{ ch_gambit_assembly }

    // Define path to subtyper files
    SUBTYPE(
        ch_gambit_assembly,
        params.subtype_db
    )
    ch_versions = ch_versions.mix(SUBTYPE.out.versions)

    //
    // MODULE: Create line summary for each sample
    //

    // Combine trimmed reads and the QC reference into single channel
    FAQCS.out.txt.map{ meta, txt -> [meta, txt] }.set{ ch_faqcs_txt }
    GAMBIT_QUERY.out.taxa.map{ meta, gambit -> [meta, gambit] }.set{ ch_gambit }
    SUBTYPE.out.subtype.map{ meta, subtype -> [meta, subtype] }.set{ ch_subtype }
    SHOVILL.out.contigs.map{ meta, contigs -> [meta, contigs] }.join(ch_faqcs_txt).join(ch_gambit).join(ch_subtype).set{ ch_line_summary_input }

    EXTRACT_CLOSEST_ACCESSION(
        ch_gambit,
        params.closest_accession_map_file
    )

    EXTRACT_CLOSEST_ACCESSION
        .out
        .acc_and_uri
        .map{ meta, s3_uri, closest_accession -> 
            [meta.id, s3_uri] }
        .set{ ch_s3_uri }

    PRE_MYCOSNP_INDV_SUMMARY(
        ch_line_summary_input,
        ch_s3_uri
    )

    //
    // MODULE: Combine line summaries into single output
    //
    PRE_MYCOSNP_COMB_SUMMARY(
        PRE_MYCOSNP_INDV_SUMMARY.out.result.map{ meta, result -> [result] }.collect()
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (params.workflow != "BOTH") {
        workflow_summary    = WorkflowMycosnp.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
        pre_mycosnp_summary = PRE_MYCOSNP_COMB_SUMMARY.out.wslh_results
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
