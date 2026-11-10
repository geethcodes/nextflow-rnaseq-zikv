nextflow.enable.dsl=2

params.samplesheet = "metadata/samplesheet.csv"

process DOWNLOAD_SRA {

    tag "$sample"
    publishDir "${projectDir}/data/raw/SRA", mode: 'copy'
    input:
    tuple val(sample), val(accession)

    output:
    tuple val(sample), path("${accession}"), emit: sra_files

    script:
    """
    prefetch --output-directory . ${accession}
    """
}

process SRA_TO_FASTQ {
   tag "$sample"
   publishDir "${projectDir}/data/raw/FASTQ" , mode: 'copy'
   
   input:
   tuple val(sample), path(sra_dir)
  
   output:
    tuple val(sample),
          path("${sra_dir.getName()}_1.fastq"),
          path("${sra_dir.getName()}_2.fastq"),
          emit: fastq_files

    script:
    """
    mkdir -p tmp

    fasterq-dump \
        --split-files \
        --threads 4 \
        --temp tmp \
        --outdir . \
        ${sra_dir}
    """
}

 process FASTQC {
   tag "$sample"
   publishDir "${projectDir}/results/fastqc", mode: 'copy'

   input:
   tuple val(sample), path(read1), path(read2)

   output:
tuple val(sample),
      path("*_fastqc.html"),
      path("*_fastqc.zip"),
      emit: fastqc_reports

    script:
    """
    fastqc \
        --threads 2 \
        ${read1} ${read2}
    """
}
   

process FASTP {

    tag "$sample"

    publishDir "${projectDir}/results/fastp", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_R1.trim.fastq.gz"),
          path("${sample}_R2.trim.fastq.gz"),
          path("${sample}.fastp.html"),
          path("${sample}.fastp.json"),
          emit: trimmed_reads

    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \
        --out1 ${sample}_R1.trim.fastq.gz \
        --out2 ${sample}_R2.trim.fastq.gz \
        --detect_adapter_for_pe \
        --thread 4 \
        --html ${sample}.fastp.html \
        --json ${sample}.fastp.json
    """
}


process POST_TRIM_FASTQC {

    tag "$sample"

    publishDir "${projectDir}/results/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("*_fastqc.html"),
          path("*_fastqc.zip"),
          emit: trimmed_fastqc_reports

    script:
    """
    fastqc \
        --threads 2 \
        ${read1} ${read2}
    """
}



process SALMON_INDEX {

    tag "Mouse transcriptome index"

    module 'gnu/12.2.0', 'salmon/1.10.1'

    publishDir "${projectDir}/results/salmon/index", mode: 'copy'

    input:
    path transcriptome

    output:
    path "salmon_index", emit: salmon_index

    script:
    """
    salmon index \
        -t ${transcriptome} \
        -i salmon_index \
        --threads 8
    """
}


process SALMON_QUANT {

    tag "$sample"

    module 'gnu/12.2.0', 'salmon/1.10.1'

    publishDir "${projectDir}/results/salmon/quant", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2), path(index)

    output:
    tuple val(sample),
          path("${sample}"),
          emit: quant_results

    script:
    """
    salmon quant \
        -i ${index} \
        -l A \
        -1 ${read1} \
        -2 ${read2} \
        --validateMappings \
        --seqBias \
        --gcBias \
        -p 8 \
        -o ${sample}
    """
}

process TXIMPORT {

    tag "tximport"

beforeScript 'export R_LIBS_USER=/N/slate/gkarutur/Rlibs'

    publishDir "${projectDir}/results/tximport", mode: 'copy'

    input:
    path quant_dirs
    path tx2gene
    path tximport_script

    output:
    path "txi.rds", emit: txi

    script:
    """
    Rscript ${tximport_script} \
        ${tx2gene} \
        txi.rds \
        ${quant_dirs.join(' ')}
    """
}


process DESEQ2 {

     tag "Differential expression: ZIKV vs PBS"

    module 'gnu/14.2.0', 'r/4.5.1'

    beforeScript 'export R_LIBS_USER=/N/slate/gkarutur/Rlibs'

    publishDir "${projectDir}/results/deseq2", mode: 'copy'

    input:
    path txi_file
    path samplesheet
    path deseq2_script

    
    output:
    path "DESeq2_results.csv", emit: deseq2_results
    path "volcano_plot.png", emit: volcano_plot

    script:
    """
    Rscript ${deseq2_script} \
        ${txi_file} \
        ${samplesheet} \
        DESeq2_results.csv
    """
}


  
workflow {

    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, row.accession)
        }
     transcriptome_ch = Channel.fromPath(
        "${projectDir}/data/reference/Mus_musculus.GRCm39.cdna.all.fa.gz",
        checkIfExists: true
    )
           
    DOWNLOAD_SRA(samples_ch)
    SRA_TO_FASTQ(DOWNLOAD_SRA.out.sra_files)

    FASTQC(SRA_TO_FASTQ.out.fastq_files) 
    FASTP(SRA_TO_FASTQ.out.fastq_files)


    trimmed_fastq_ch = FASTP.out.trimmed_reads.map {
    sample, read1, read2, html, json ->
    tuple(sample, read1, read2)
}

POST_TRIM_FASTQC(trimmed_fastq_ch)

 SALMON_INDEX(transcriptome_ch)

// Connect trimmed reads with the Salmon index
salmon_input_ch = trimmed_fastq_ch
    .combine(SALMON_INDEX.out.salmon_index)

// Run Salmon quantification
SALMON_QUANT(salmon_input_ch)


// Collect all Salmon quantification directories
quant_dirs_ch = SALMON_QUANT.out.quant_results
    .map { sample, quant_dir -> quant_dir }
    .collect()

// Reference mapping and R script
tx2gene_ch = Channel.fromPath(
    "${projectDir}/data/reference/tx2gene.tsv",
    checkIfExists: true
)

tximport_script_ch = Channel.fromPath(
    "${projectDir}/scripts/tximport.R",
    checkIfExists: true
)

// Run transcript-to-gene import
TXIMPORT(quant_dirs_ch, tx2gene_ch, tximport_script_ch)

// Samplesheet and DESeq2 script
samplesheet_ch = Channel.fromPath(
    params.samplesheet,
    checkIfExists: true
)

deseq2_script_ch = Channel.fromPath(
    "${projectDir}/scripts/deseq2.R",
    checkIfExists: true
)

// Run differential expression analysis
DESEQ2(
    TXIMPORT.out.txi,
    samplesheet_ch,
    deseq2_script_ch
)




}
