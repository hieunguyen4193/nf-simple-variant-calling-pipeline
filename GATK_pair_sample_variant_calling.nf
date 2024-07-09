// ======================================================================================
// PAIR-SAMPLE VARIANT CALLING WITH MUTECT2 AND VARDICT PIPELINE
// Author: Trong-Hieu Nguyen
// hieunguyen@genesolutions.vn  
// ======================================================================================
println """\
        ===================================================================
                    S O M A T I C - V A R I A N T - C A L L I N G  
                                ----- ----- -----
                   W I T H - G A T K - M U T E C T 2 - & - V A R D I C T 
                             
         ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** 
        
         Author: TRONG - HIEU NGUYEN
         Email: hieunguyen@genesolutions.vn
         
        ===================================================================
         """
         .stripIndent()

params.METADATA                 = "" 
// metadata file containing paths to paired BAM files.
/* Run the following BASH script to create metadata.csv 
echo -e sampleID "\t" tumor_file "\t" normal_file >> sample_metadata.csv && \\ 
for file in $files;do sampleid=${file%F_*} && \\ 
sampleid=${sampleid#*-} && echo $sampleid && fileF=$(ls *${sampleid}F*.bam) && \\
fileW=$(ls *${sampleid}W*.bam) && echo -e $sampleid "\t" $(pwd)/${fileF} "\t" $(pwd)/${fileW} "\n" \\
>> sample_metadata.csv;done
*/

params.OUTPUTDIR                = ""
params.BWA_REPO                 = ""
BWA_REPO                        = file(params.BWA_REPO)
params.GERMLINE_RESOURCE        = ""
GERMLINE_RESOURCE               = file(params.GERMLINE_RESOURCE)
params.TARGET_GENES             = ""
params.GENE_PANEL               = ""
params.AF_THR                   = 0.00001;

Channel.fromPath(params.METADATA)
        .splitCsv(header:true)
        .map{ row -> tuple(row.sampleID, file(row.tumor_file), file(row.normal_file)) }
        .view()
        .into { Pair_samples_GATK_ch; Pair_samples_GATK_PON_ch; Pair_samples_VARDICT_ch}


process Vardict{
    // cache "deep"; tag "$sample_id"
    publishDir "$params.OUTPUTDIR/Vardict_SomaticCalling", mode: 'copy'
    maxForks 6
    input: 
        set sample_id, file(FFPE_file), file(NORMAL_file) from Pair_samples_VARDICT_ch
        file(BWA_REPO)
    output: 
        tuple sample_id, "Vardict_${sample_id}.vcf" into somatic_var_ch
    shell: 
    '''
    FFPE_file_name=!{FFPE_file};
    FFPE_file_name=${FFPE_file_name%.bam*};
    NORMAL_file_name=!{NORMAL_file};
    NORMAL_file_name=${NORMAL_file_name%.bam*};

    vardict-java -G ./HG38DIR/hg38_selected.fa -f !{params.AF_THR} -N ${FFPE_file_name} -b "!{FFPE_file}|!{NORMAL_file}" -c 1 -S 2 -E 3 -g 4 \
    !{params.TARGET_GENES}/!{params.GENE_PANEL}.targets.hg38.bed | testsomatic.R | var2vcf_paired.pl -N "$FFPE_file_name|$NORMAL_file_name" -f !{params.AF_THR} > Vardict_!{sample_id}.vcf
    '''
}

process Variant_calling_for_NORMAL_samples {
    cache "deep"; tag "$sample_id"
    publishDir "${params.OUTPUTDIR}/Normal_samples_vcf"

    input:
        set sample_id, file(FFPE_file), file(NORMAL_file) from Pair_samples_GATK_PON_ch
        file(BWA_REPO)
    output: 
        set sample_id, file("${sample_id}_normal.vcf.gz"), file("${sample_id}_normal.vcf.gz.tbi") into PON_ch
    shell: 
    '''
    gatk Mutect2 -R ./HG38DIR/hg38_selected.fa -I !{NORMAL_file} -max-mnp-distance 0 -O !{sample_id}_normal.vcf.gz;
    gatk IndexFeatureFile -I !{sample_id}_normal.vcf.gz
 
    '''
}

process CreatePanelOfNormals {
    cache "deep"; 
    publishDir "${params.OUTPUTDIR}/PanelOfNormals", mode: 'copy'

    input: 
        file "*" from PON_ch.collect()
        file(BWA_REPO)
    output: 
        set file("PON.vcf.gz"), file("PON.vcf.gz.tbi") into pon_mutect2_ch
    shell: 
    '''
    normal_vcfs=$(ls *_normal.vcf.gz);
    for file in ${normal_vcfs};do echo -e ${file%_normal.vcf.gz*}"\t"${file} >> normal_vcfs_to_import_PON.txt;done
    
    gatk GenomicsDBImport -R ./HG38DIR/hg38_selected.fa -L !{params.TARGET_GENES}/!{params.GENE_PANEL}.targets.hg38.interval_list \
       --genomicsdb-workspace-path pon_db \
       --sample-name-map normal_vcfs_to_import_PON.txt;

    gatk CreateSomaticPanelOfNormals -R ./HG38DIR/hg38_selected.fa -V gendb://pon_db -O PON.vcf.gz --disable-tool-default-read-filters --min-sample-count 1
    gatk IndexFeatureFile -I PON.vcf.gz

    '''
}

process SomaticCalling_GATK {
    cache "deep"; tag "$sample_id"
    publishDir "${params.OUTPUTDIR}/MUTECT2_raw_variants", mode: 'copy'

    input:
        set sample_id, file(FFPE_file), file(NORMAL_file) from Pair_samples_GATK_ch
        set file("PON.vcf.gz"), file("PON.vcf.gz.tbi") from pon_mutect2_ch
        file(BWA_REPO)
        file GERMLINE_RESOURCE
    output: 
        set sample_id, file("Mutect2_raw_${sample_id}.vcf.gz"), file("Mutect2_raw_${sample_id}.vcf.gz.tbi"), file("Mutect2_raw_${sample_id}.vcf.gz.stats") into filter_ch
    shell: 
    '''
    FFPE_file_name=!{FFPE_file};
    FFPE_file_name=${FFPE_file_name%.bam*};
    NORMAL_file_name=!{NORMAL_file};
    NORMAL_file_name=${NORMAL_file_name%.bam*};

    gatk Mutect2 \
     -R ./HG38DIR/hg38_selected.fa \
     -I !{FFPE_file} \
     -I !{NORMAL_file} \
     -normal $NORMAL_file_name \
     --germline-resource !{GERMLINE_RESOURCE}/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
     --panel-of-normals PON.vcf.gz \
     --af-of-alleles-not-in-resource 0.0000025 \
     -O Mutect2_raw_!{sample_id}.vcf.gz;
    gatk IndexFeatureFile -I Mutect2_raw_!{sample_id}.vcf.gz

    '''
}

process Filter_variants{
    cache "deep"; tag "$sample_id"
    publishDir "${params.OUTPUTDIR}/MUTECT2_filtered_variants", mode: 'copy'

    input:
        file(BWA_REPO) 
        set sample_id, file("Mutect2_raw_${sample_id}.vcf.gz"), file("Mutect2_raw_${sample_id}.vcf.gz.tbi"), file("Mutect2_raw_${sample_id}.vcf.gz.stats") from filter_ch


    output: 
        set sample_id, file("Mutect2_filtered_${sample_id}.vcf.gz") into output_ch

    script: 
    """
    gatk FilterMutectCalls \
   -R ./HG38DIR//hg38_selected.fa \
   -V Mutect2_raw_${sample_id}.vcf.gz \
   -O Mutect2_filtered_${sample_id}.vcf.gz
    """
}

// end of pipeline ======================================================================
