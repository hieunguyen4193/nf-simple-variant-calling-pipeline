// ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
// N E X T F L O W pipeline for ctDNA project - GS 
// =====  ==== ===== ===== ===== ===== ===== ===== ===== =====
params.CPU_BWA = ""
// ----- ----- ----- FOLDERS STRUCTURE  ----- ----- ----- 
// input parameters.
params.BWA_REPO= "" 

params.INPUT_DIR= ""
params.OUTPUT_DIR= "" 

params.BEDFILE=""
bedfile=file(params.BEDFILE)

params.HG38DIR=""
hg38=file(params.HG38DIR)

params.RESOURCES=""
params.dbsnp="${params.RESOURCES}/genomes/dbSNP_INDELS/dbsnp_146.hg38.vcf"
params.ADAPTER="${params.RESOURCES}/adapters/TruSeq3-PE-2.fa"
ADAPTER=file(params.ADAPTER)
dbsnp=file(params.dbsnp)
// 
params.AF_THR="0.00001"
params.INPUT_PAIRS= "$params.INPUT_DIR/*_R{1,2}*" 
// 
// params.WORK=""
// work=file(params.WORK)

params.dir_cache=""
dir_cache=file(params.dir_cache)

params.fasta=""
fasta=file(params.fasta)

params.clinvar=""
clinvar=file(params.clinvar)
params.clinvar_index=""
clinvar_index=file(params.clinvar_index)

params.COSMIC=""
COSMIC=file(params.COSMIC)
params.COSMIC_INDEX=""
COSMIC_INDEX=file(params.COSMIC_INDEX)


// ----- ----- ----- CHANNEL ----- ----- -----
Channel
    .fromFilePairs( params.INPUT_PAIRS )
    .ifEmpty { error "Cannot find any reads matching: ${params.INPUT_PAIRS}"  }
    // .view()
    .into {fastqc_ch; trim_ch}

//  ----- ----- ----- FASTQC ----- ----- -----
process fastqc {
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/0_fastqc_results", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5

    input:
        tuple sample_id, file(READS) from fastqc_ch 
    output:
        tuple sample_id, '*_fastqc.{zip,html}' into fastqc_results
    
    script:
        """
        fastqc --quiet --outdir . ${READS[0]} ${READS[1]}
        """
}

//  ----- ----- ----- TRIM  ----- ----- -----
process trim {
    cache "deep"; 
    tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/1_trimmed_fastq", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5
    // cpus    params.CPU_trim
    // memory  params.RAM_trim
    input:
        tuple sample_id, file(READS) from trim_ch

    output:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R3_paired.fastq.gz" into bwa_ch, fastqc_postTRIM

    script:
        """
        trimmomatic PE -phred33 -threads 12 \
        ${READS} \
        ${sample_id}_R1_paired.fastq.gz ${sample_id}_R1_unpaired.fastq.gz \
        ${sample_id}_R3_paired.fastq.gz ${sample_id}_R3_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTER}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 CROP:75 MINLEN:36
        """
}

//  ----- ----- ----- Alignment ----- ----- -----
process align{
    cache "deep"; tag "$sample_id"
    // publishDir "$params.OUTPUT_DIR/alignment_bwamem", mode: 'copy'
    // errorStrategy 'retry'
    maxForks 5

    input:
        tuple sample_id, "${sample_id}_R1_paired.fastq.gz", "${sample_id}_R3_paired.fastq.gz" from bwa_ch
    output:
        tuple sample_id, "alignBWA_${sample_id}.sam" into sam2bam_ch
    script:
        """
        bwa mem -t ${params.CPU_BWA} -M -R \
            '@RG\\tID:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\tPM:MINISEQ\\tSM:${sample_id}' \\
            ${params.BWA_REPO}\
            ${sample_id}_R1_paired.fastq.gz ${sample_id}_R3_paired.fastq.gz > alignBWA_${sample_id}.sam
        """
}

//  ----- ----- ----- SAM to BAM, mapped and sorted reads ----- ----- -----
process SortSAM{
    cache "deep"; tag "$sample_id"
    errorStrategy 'retry'
    maxRetries 10
    publishDir "$params.OUTPUT_DIR/2_SortedSAM", mode: "copy"
    // maxForks 5
    input:
        tuple sample_id, "alignBWA_${sample_id}.sam" from sam2bam_ch

    output:
        tuple sample_id, "${sample_id}_sorted_mapped.bam" into markdup_ch
        
    script:
        """
        samtools view -b -F 4 alignBWA_${sample_id}.sam | samtools sort -@ 10 | samtools view -@ 10 -bS -o ${sample_id}_sorted_mapped.bam; 
        samtools index -@ 10 ${sample_id}_sorted_mapped.bam
        """
}

//  ----- ----- ----- Picard: Sort Bam and Mark Duplicates ----- ----- -----
process picard_MarkDuplicates{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/3_MarkDup", mode: "copy"
    errorStrategy 'retry'
    maxRetries 1
    // maxForks 5   
    input:
        tuple sample_id, "${sample_id}_sorted_mapped.bam" from markdup_ch
    output:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.dedup.bai" into vardict, mutect2
        file "${sample_id}.metrics.txt" into picard_metrics_markdup
    script:
        """
        picard SortSam --INPUT ${sample_id}_sorted_mapped.bam --SORT_ORDER coordinate --OUTPUT ${sample_id}.sorted.bam;
        picard MarkDuplicates --INPUT ${sample_id}.sorted.bam --OUTPUT ${sample_id}.dedup.bam --METRICS_FILE ${sample_id}.metrics.txt --BARCODE_TAG RX;
        picard BuildBamIndex --INPUT ${sample_id}.dedup.bam
        """
}

process vardict{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/9_vardict_output", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 3 
    input:
        tuple sample_id, "${sample_id}.dedup.bam", "${sample_id}.dedup.bai" from vardict
    output:
        tuple sample_id, "vardict_single_${sample_id}.vcf" into vardict_output
    script:
        """
        samtools index ${sample_id}.dedup.bam;
        vardict-java -t -G ${hg38} -q 20 -f ${params.AF_THR} -N $sample_id -b ${sample_id}.dedup.bam \
        -c 1 -S 2 -E 3 -g 4 ${bedfile} \
        | teststrandbias.R | var2vcf_valid.pl -N ${sample_id} -E -f ${params.AF_THR} -r 3 > vardict_single_${sample_id}.vcf
        """
}

process annotate_vep_vardict{
    cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUT_DIR/10_annotated_vcf_vardict", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 5
    input:
        tuple sample_id, "vardict_single_${sample_id}.vcf" from vardict_output
        file dir_cache 
        file fasta 
        file clinvar 
        file COSMIC
        file clinvar_index
        file COSMIC_INDEX
    output:
        tuple sample_id, "vardict_single_${sample_id}_annotated.vcf" into final_output_vardict
    script:
    """
    vep --stats_file ${sample_id}.vep_summary.html --cache\
    --dir_cache $dir_cache --refseq --force_overwrite --format vcf \
    --fasta $fasta --everything --pick   \
    --custom $clinvar,ClinVar,vcf,exact,0,ALLELEID,CLNSIG,GENEINFO,CLNREVSTAT,CLNDN,CLNHGVS \
    --custom $COSMIC,COSMIC,vcf,exact,0,GENE,STRAND,LEGACY_ID,SNP,CDS,HGVSC,HGVSG,CNT \
    -i vardict_single_${sample_id}.vcf --check_existing --symbol --vcf \
    -o vardict_single_${sample_id}_annotated.vcf --offline;
    """
}
