while read SM LB PU
do
# use "_" to sperate SM,LB and PU tag, the newfastq name is same with gropu name.
group_prefix="${SM}_${LB}_$PU"
fastq_prefix=$group_prefix
bam_input="$bam_input -i sorted_$group_prefix.bam"

# check
if [ ! -e "sorted_$group_prefix.bam.bai" ];then
echo "rerun[sorted_$group_prefix.bam]"
( $BWA mem -Y -M -R "@RG\tID:$group_prefix\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$fastq_prefix$i$fastq_suffix_1 $fastq_folder/$fastq_prefix$i$fastq_suffix_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort -r $fasta -o sorted_$group_prefix.bam -t 32 --sam2bam -i -

else
echo "algo bwa-sort [$group_prefix] done!"
fi

done < ../$CFG	# 1028: next using absolute path for cfg file

# ******************************************
# 2. Metrics on the multiple sorted BAM files
# ******************************************

# ******************************************
# 3. Remove Duplicate Reads on the multiple sorted BAM files
# ******************************************

if [ ! -e "score.txt.idx" ]; then
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt $bam_input --algo LocusCollector --fun score_info score.txt
else
echo "algo Dedup.sore.done"
fi

if [ ! -e "deduped.bam.bai" ]; then
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt $bam_input --algo Dedup --score_info score.txt --metrics dedup_metrics.txt deduped.bam
else
echo "algo Dedup done"
fi

# ******************************************
# 4. Indel realigner removed by gatk4
# ******************************************


# ******************************************
# 5. Base recalibration
# ******************************************

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table

# tnhp2/tnsp do not nedd output --algo ReadWriter recaled.bam , but gatk4 m2 need
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post --algo ReadWriter recaled.bam && echo "DONE:BQSR for M2"
# ******************************************
#  QC bam
# ******************************************
if [ ! -e "bamqc.information.tsv" ]; then
$basepath/NCbamInfo -b $interval_file -r $fasta -t -p bamqc deduped.bam
else
echo "bamqc done"
fi

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

