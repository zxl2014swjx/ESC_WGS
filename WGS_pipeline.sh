#creat WGS-sample.txt 
perl  /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/latitude/data_path_ln.pl -f  WGS/ -l WGS-sample.txt -o output/
# get bam
for i in *;do dir=/mnt/X500/farmers/songmm/project/Lung_Cheng/output/ echo "python /mnt/NL200/chenzhx/pipeline/single-WES/bin/wes.py getbam -i $dir/$i/input -o $dir/$i/output -c nc_g -n 5 -s $i -l 150 > $dir/$i/$i.log 2>&1 " > $dir/$i/$i.sh;done
for i in *-*; do echo "qsub $i.sh" >> WGS_getbam.sh ; done
sh WGS_getbam.sh
# creat snv & indel needed config file
for i in *-*; do echo "/mnt/ssd/local/repo/strelka/src/python/bin/configureStrelkaSomaticWorkflow.py --config /usr/local/bin/configureStrelkaSomaticWorkflow.py.ini --normalBam=$PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.bam --tumorBam=$PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.bam --referenceFasta=/refgenomes/hs37d5/hs37d5.fa --runDir $PWD/$i/output/" >> creat_snv_indel_config.sh; done
sh creat_snv_indel_config.sh
# call snv & indel
for i in *-*; do echo "python /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/output/runWorkflow.py -m sge -q baseline.q -j 30 -g 40 > /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/${i}_call_snv_indel.log 2>&1" > /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/${i}_call_snv_indel.sh; done
for i in *-*; do echo "sh $i/${i}_call_snv_indel.sh" >> call_snv_indel.sh; done
sh call_snv_indel.sh
# fix vcf
for i in *-* ; do echo "python /mnt/NL200/chenzhx/pipeline/AKso/DevOps/chenzhx/tools/fixStrelkaVcf.py -i $PWD/$i/output/results/variants/somatic.indels.vcf.gz -o $PWD/$i/output/results/variants/somatic_standard.indels.vcf.gz && python /mnt/NL200/chenzhx/pipeline/AKso/DevOps/chenzhx/tools/fixStrelkaVcf.py -i $PWD/$i/output/results/variants/somatic.snvs.vcf.gz -o $PWD/$i/output/results/variants/somatic_standard.snvs.vcf.gz" >> fix_snv_indel.sh ; done
sh fix_snv_indel.sh
# anno vcf
for i in *-* ; do echo "perl /repos/chenzhx/single-vcf/bin/uniAnno.pl /repos/chenzhx/single-vcf/lib/NoahCare/config/nc_g/anno/config_nc.pl -q -n 5 -b 500 -i $PWD/$i/output/results/variants/somatic_standard.indels.vcf.gz -o $PWD/$i/output/results/variants/somatic_standard.indels.tsv && perl /repos/chenzhx/single-vcf/bin/uniAnno.pl /repos/chenzhx/single-vcf/lib/NoahCare/config/nc_g/anno/config_nc.pl -q -n 5 -b 500 -i $PWD/$i/output/results/variants/somatic_standard.snvs.vcf.gz -o $PWD/$i/output/results/variants/somatic_standard.snvs.tsv" >> anno_snv_indel.sh; done
sh anno_snv_indel.sh
# create sv config file
for i in *-*; do rm_str=""; if [ -f "$i/output/sv/runWorkflow.py" ]; then rm_str="rm $i/output/sv/runWorkflow.py* && ";fi;echo "$rm_str /repos/chenzhx/manta/bin/configManta.py --config /repos/chenzhx/manta/bin/configManta.py.ini --normalBam=$PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.bam --tumorBam=$PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.bam --referenceFasta=/refgenomes/hs37d5/hs37d5.fa --runDir $PWD/$i/output/sv" >> creat_sv_config.sh;done
sh creat_sv_config.sh
# call sv
for i in *-*; do echo "python /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/output/sv/runWorkflow.py -m sge -q baseline.q -j 30 -g 40 > /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/${i}_call_sv.log 2>&1" > /mnt/X500/farmers/songmm/project/Lung_Cheng/output/$i/${i}_call_sv.sh; done
for i in *-*; do echo "sh $i/${i}_call_sv.sh" >> call_sv.sh; done
sh call_sv.sh
# merge snvs & indels of all sample 
perl /mnt/X500/farmers/songmm/pipeline/extract_snv_indel.pl -f /mnt/X500/farmers/songmm/project/Lung_Cheng/output/ -o /mnt/X500/farmers/songmm/project/Lung_Cheng/output/all
# filter snv & indel
cd all && perl /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/filter_SNV.pl -l /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/filter.list -s all/snv_all.txt  -ex /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/ncbi_anno_rel104_db_hg19_EXN_Primary.tab && perl /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/filter_SNV.pl -l /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/filter.list -s ./indel_all.txt  -ex /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/ncbi_anno_rel104_db_hg19_EXN_Primary.tab
# cnv
for i in *-*; do echo "/mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/Longitude/CNV_sh/100wesFacets.sh $PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.bam $PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.bam  $PWD/$i/output/100wesFacets $i" > $PWD/$i/${i}_100wesFacets.sh; done
ls *-*/*_100wesFacets.sh|awk '{print "qsub "$0}' > 100wesFacets.sh
sh 100wesFacets.sh
# msi
for i in *-*; do echo "/repos/yangxx/Akso/BNC/bin/msisensor msi -d /repos/yangxx/Akso/BNC/program/NoahCare/db/analyze/msi/hs37d5.microsatellites.list -n $PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.bam -t $PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.bam -e  /mnt/NL200/suyao/tools/BNC/BNC/program/NoahCare/config/nc_e/regs/flk50_chip.bed -c 15 -o $PWD/$i/output/$i.msi" > $PWD/$i/msi_$i.sh ; done
ls *-*/msi_*.sh|awk '{print "qsub "$0}' > wgs_msi.sh
# extract sample snv gene
for i in *-*; do echo "awk '{if(\$1~/$i/) print \$0}' all/snv_all.txt.filter |cut -d $'\t' -f22 > all/venn/$i.snv" >> extract_svn_gene.sh ; done
# stat mutation gene in cosmic top20
perl /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/merge_SNV/cosmic_topgene/count_Mutation_gene_Freq.pl -l cosmic_top20.list -asample_mutationGene.list -g group.txt -o top20_vs_WGS
# 
awk  -F "\t" '{print $1"\t"$22"\t"$32"\t"$27"\t"$6"\t"$8"\t"$7"\t"$9}' snv_all.txt.filter > snv_forSig
# pin pu
Rscript /mnt/X500/farmers/suyao/NL200/analysis2/module/oncoprint/v2.0/mut_oncoprint.R -m snv_forheatmap --groupfile group2.txt -a LC --height 10 -o mutation_heatmap.pdf --cut 2 --colname T
# signature
export PATH="/usr/local/MATLAB/R2017b/bin/:$PATH" && perl /mnt/X500/farmers/zhaogf/signature/select_mutsig_number.pl -i snv_forheatmap -o sig -t lung
perl /mnt/X500/farmers/zhaogf/signature/mutsigframeSigNum.pl -i sig -n lung -t 2
# CNV 
for i in *-*; do echo -e "samtools mpileup -f /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa -q 20 -Q 20 -B -C 50 $PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.bam > $PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.mpileup\nsamtools mpileup -f /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa -q 20 -Q 20 -B -C 50 $PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.bam > $PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.mpileup" > $PWD/$i/mpileup_$i.sh; done
# creat control_freec config
for i in *-*; do echo "perl /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/latitude/control-freec/control_freec-config-500k.pl -c $PWD/$i/output/cancer/4_recal_bam/${i}_case_sort_markdup_recald.mpileup -n $PWD/$i/output/normal/4_recal_bam/${i}_normal_sort_markdup_recald.mpileup -o $PWD/$i/output/CNV -s XY -f $PWD/$i/output/CNV/$i.control_freec.config" >> control_freec_config.sh; done
#
for i in *-*; do echo "/mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/src/freec -conf $PWD/$i/output/CNV/$i.control_freec.config > $PWD/$i/output/CNV/$i.control_freec.log 2>&1" > $PWD/$i/$i.freec.sh; done
#
for i in *-*; do echo "cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/assess_significance.R | R --slave --args $PWD/$i/output/CNV/${i}_case_sort_markdup_recald.mpileup_CNVs $PWD/$i/output/CNV/${i}_case_sort_markdup_recald.mpileup_ratio.txt" > $PWD/$i/assess_significance.$i.sh; done
#
for i in *-*; do echo "cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/makeGraph.R | R --slave --args 2 $PWD/$i/output/CNV/${i}_case_sort_markdup_recald.mpileup_ratio.txt $PWD/$i/output/CNV/${i}_case_sort_markdup_recald.mpileup_BAF.txt" >$PWD/$i/makeGraph.$i.sh; done
#
for i in *-*; do echo "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/control_freec_exon.R $PWD/$i/output/CNV/${i}_case_sort_markdup_recald.mpileup_CNVs.p.value.txt /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/ncbi_anno_rel104_db_hg19_EXN_Primary.tab $PWD/$i/output/CNV/$i.CNV.anno.xls" > $PWD/$i/control_freec_exon.$i.sh; done
# anno CNV
perl control_freeC_anno_500k.pl -f /mnt/X500/farmers/songmm/project/Lung_Cheng/output/ -cyto /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/merge_CNV/cytoBand_hg19_grch37.modify.txt -exon /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/merge_CNV/ncbi_anno_rel104_db.bed
# check CNV
perl  ../../../pipeline/check_control_freeC_anno.pl -f /mnt/X500/farmers/songmm/project/Lung_Cheng/output -k lung
