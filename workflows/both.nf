include { PRE_MYCOSNP_WF } from '../workflows/pre_mycosnp'
include { MYCOSNP } from '../workflows/mycosnp'

include { WSLH_SUMMARY } from '../modules/local/wslh_summary'
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// params.snpeffdb = WorkflowMain.getGenomeAttribute(params, 'snpeffdb')
params.snpeffconfig = WorkflowMain.getGenomeAttribute(params, 'snpeffconfig')


// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ] // params.snpeffdb
if (params.skip_samples_file) { // check for skip_samples_file
    checkPathParamList.add(params.skip_samples_file)
}
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// if (params.snpeffdb == null) { exit 1, 'Input path to snpeffdb not specified!' }
// if (params.snpeffconfig == null) { exit 1, 'Input snpeff config file not specified' }


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
include { SRA_FASTQ_SRATOOLS } from '../subworkflows/local/sra_fastq_sratools'
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { BWA_PREPROCESS     } from '../subworkflows/local/bwa-pre-process'
include { BWA_REFERENCE      } from '../subworkflows/local/bwa-reference'
include { GATK_VARIANTS      } from '../subworkflows/local/gatk-variants'
include { CREATE_PHYLOGENY   } from '../subworkflows/local/phylogeny'
include { SNPEFF_BUILD       } from '../subworkflows/local/snpeff_build'
include { SNPEFF             } from '../subworkflows/local/snpeff'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/modules/fastqc/main'
include { QC_REPORTSHEET              } from '../modules/local/qc_reportsheet.nf'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { GATK4_HAPLOTYPECALLER       } from '../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GATK4_COMBINEGVCFS          } from '../modules/nf-core/modules/gatk4/combinegvcfs/main'
include { SEQKIT_REPLACE              } from '../modules/nf-core/modules/seqkit/replace/main'
include { SNPDISTS                    } from '../modules/nf-core/modules/snpdists/main'
include { GATK4_LOCALCOMBINEGVCFS     } from '../modules/local/gatk4_localcombinegvcfs.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []


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
        PRE_MYCOSNP_WF.out.pre_mycosnp_summary
    )
}