METADATA=/media/tronghieu/HNSD01/WORKING_DATA/pvacseq/pvacseq_samples/align/sample_metadata.csv;
OUTPUTDIR=/media/tronghieu/HNSD01/WORKING_DATA/pvacseq/pvacseq_samples/test_output;
BWA_REPO=/media/tronghieu/GSHD_HN01/OFFICIAL/annotation_resources/HG38DIR;
GERMLINE_RESOURCE=/media/tronghieu/HNSD01/WORKING_DATA/pvacseq/germline_resources_GATK;
TARGET_GENES=/media/tronghieu/HNSD01/work/target_genes;
GENE_PANEL=onco08;
nextflow=/media/tronghieu/HNSD01/work/nextflow;

$nextflow run GATK_pair_sample_variant_calling.nf -w ./work \
-c GATK_pair_sample_variant_calling.config \
-with-report test_report.html \
--METADATA $METADATA \
--OUTPUTDIR $OUTPUTDIR \
--BWA_REPO $BWA_REPO \
--GERMLINE_RESOURCE $GERMLINE_RESOURCE \
--TARGET_GENES $TARGET_GENES \
--GENE_PANEL $GENE_PANEL \
-resume