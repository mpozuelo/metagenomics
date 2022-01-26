#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/metagenomics
========================================================================================
 mpozuelo/metagenomics
 #### Homepage / Documentation
 https://github.com/mpozuelo/metagenomics
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/metagenomics --input '*.txt' -profile docker

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.

    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }



// Validate inputs
samplesheet_ch = file(params.input, checkIfExists: true)
nucleotide_DB = params.ntDB
protein_DB = params.ppDB


if (!params.outdir) {
  params.outdir = params.run
}




// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Name'] = custom_runName ?: workflow.runName
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-metagenomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/metagenomics Workflow Summary'
    section_href: 'https://github.com/mpozuelo/metagenomics'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version &> v_ngi_rnaseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
Get the ONLY the samplesheet for getting information per sample for the Trimming
Also get information about the cycles for the demultiplexing
*/


/*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing
 It contains the following columns:
0- Sample_ID (mandatory for Illumina)
1- Sample_Name (mandatory for Illumina)
2- index (mandatory for Illumina)
3- index2 (mandatory for Illumina)
4- SampleProject
5- SampleSource
6- Organism
7- Method
8- Library
9- Platform
10- Run
11- Date
12- User
13- Coverage
*/


Channel
  .from( s )
  .splitCsv(header:false, sep:',')
  .map { it = ["${it[1]}", "${it[2]}", "${it[3]}",
  ]}
  .set { ch_samplesheet }


/*
Trimmomatic
*/

process kneaddata_index {
  container "biobakery/kneaddata"
  tag "$sample"
  label "process_low"

  output:
  path "human" into ch_index

  script:
  """
  mkdir human
  kneaddata_database --download human_genome bowtie2 human
  """

}


/* Bowtie2 */

process kneaddata_trim_rmHost {
  tag "$sample"
  label "process_medium"

  input:
  set val(sample), path(fastq1), path(fastq2) from ch_samplesheet
  path index from ch_index

  output:
  path "rmHost/*_paired_*.fastq" into ch_remove_host
  path "rmHost"


  script:
  """
  mkdir rmHost
  kneaddata --input1 ${fastq1}--input2 ${fastq2} --reference-db ${index} --output rmHost --run-fastqc-start --run-fastqc-end
  """


}


/*
 * STEP MERGE P1 AND P2


process concatenate {
  tag "$sample"
  label "process_medium"

  input:
  set val(sample), path(fastq) from ch_remove_host

  output:
  path "*.gz" into

  script:
  output = sample".fq"
  """
  cat ${fastq[0]} > ${output}
  cat ${fastq[1]} >> ${output}

  pigz -p 6 ${output}
  """
}


/*
 * STEP HUMAnN3.0


process humann3 {

  container 'biobakery/humann:3.0.0'
  tag "$sample"
  label 'process_high'
  publishDir "results/humann3", mode: 'copy'

  input:
  path fastq from ch_fastq_merged

  output:
  path "*.tsv" into ch_gene_out
  path "*_temp/*bugs_list.tsv" into ch_tax_out
  path "*_temp/*.log" into ch_logs

  script:
  """
  humann --input ${fastq} --output . --threads 10 --nucleotide-database ${nucleotide_DB} --protein-database ${protein_DB}
  """

}


workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/metagenomics] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/metagenomics] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/metagenomics]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/metagenomics]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/metagenomics v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
