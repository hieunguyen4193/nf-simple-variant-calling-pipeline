# Run N E X T F L O W pipeline for ctDNA project.
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
# Define input parameters.
# ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
# 
# 
# 
RUN=H004;
maindir=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD;
srcdir=$maindir/src_$RUN;
# 
# 
# 
for group in liver colon;do \
CPU_BWA=4;
re=-resume;

BWA_REPO=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/annotation_resources/HG38DIR/hg38_selected;
# --------
INPUT_DIR=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/ctDNA/fastq/ffpe/$group;
OUTPUT_DIR=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/ECD/output/$RUN/$group;

BEDFILE_colon=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/annotation_resources/target_genes/new_version/oncoTTYSH.targets.hg38.bed;
BEDFILE_liver=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/annotation_resources/target_genes/new_version/onco08.targets.hg38.bed;
if [ "$group" = "colon" ]; then BEDFILE=$BEDFILE_colon
else BEDFILE=$BEDFILE_liver
fi;
# --------
echo -e "Using bedfile for" $group "\n" >> $srcdir/log.txt
echo $BEDFILE >> $srcdir/log.txt;
HG38DIR=$BWA_REPO.fa;

RESOURCES=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/annotation_resources;

dir_cache=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/dir_cache_VEP99;

fasta=$dir_cache/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz;
clinvar=$dir_cache/clinvar/clinvar_20201026.vcf.gz;
clinvar_index=$dir_cache/clinvar/clinvar_20201026.vcf.gz.tbi
COSMIC=$dir_cache/COSMIC/CosmicCodingMuts.normal.vcf.gz;
cosmic_index=$dir_cache/COSMIC/CosmicCodingMuts.normal.vcf.gz.tbi;

nextflow=/media/genesolutions/DELL02_HD01/DataResearch/TrongHieu-2021/data/nextflow_exe/nextflow;

chmod -R a+rwx $RESOURCES;
chmod -R a+rwx $dir_cache;
chmod -R a+rwx $OUTPUT_DIR;

$nextflow $srcdir/ECD_pipeline_noUMIsamples_ffpe.nf -c $srcdir/ECD_pipeline_noUMIsamples_ffpe.config \
-w $srcdir/work \
--BWA_REPO $BWA_REPO \
--CPU_BWA $CPU_BWA \
--INPUT_DIR $INPUT_DIR \
--OUTPUT_DIR $OUTPUT_DIR \
--BEDFILE $BEDFILE \
--HG38DIR $HG38DIR \
--RESOURCES $RESOURCES \
--dir_cache $dir_cache \
--fasta $fasta \
--clinvar_index $clinvar_index \
--clinvar $clinvar \
--COSMIC $COSMIC \
--COSMIC_INDEX $cosmic_index \
$re;done