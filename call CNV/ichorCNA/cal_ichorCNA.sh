#!/bin/bash

if [ $# -lt 5 ] ; then
	echo -e "\n\tsh $0 <tumor_bam> <normal_bam> <tumor_id> <normal_id> <out_dir> \n"
	exit
fi

tumor_bam=$1
normal_bam=$2
tumor_id=$3
normal_id=$4
out_dir=$5

if [ ! -d $out_dir ] ; then
	mkdir -p $out_dir
fi

########step_1
if [ ! -f $out_dir/normal.wig ] ; then
	/mnt/X500/farmers/linhx/bin/cna/ichorCNA_bamtools/ichorCNA_bamtools readcount 999999999 30 \
	/mnt/X500/farmers/linhx/bin/ref/region_bed/combine.black.merge.complement.nochr.bed \
	$normal_bam \
	$out_dir/normal.wig \
	1000000 
fi

########step_2
if [ ! -f $out_dir/tumor.wig ] ; then
	/mnt/X500/farmers/linhx/bin/cna/ichorCNA_bamtools/ichorCNA_bamtools readcount 999999999 30 \
	/mnt/X500/farmers/linhx/bin/ref/region_bed/combine.black.merge.complement.nochr.bed \
	$tumor_bam \
	$out_dir/tumor.wig \
	1000000
fi

########step_3
/mnt/X500/farmers/linhx/bin/miniconda3/envs/lohhla/bin/Rscript /mnt/X500/farmers/linhx/bin/miniconda3/envs/py37_withR/bin/runIchorCNA.R \
--id ${normal_id}_${tumor_id} \
--WIG $out_dir/tumor.wig \
--NORMWIG $out_dir/normal.wig \
--ploidy "c(2,3)" \
--normal "c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)" \
--maxCN 5  \
--gcWig /mnt/X500/farmers/linhx/bin/miniconda3/envs/py37_withR/share/r-ichorcna-0.1.0.20180710-0/extdata/gc_hg19_1000kb.wig \
--mapWig /mnt/X500/farmers/linhx/bin/miniconda3/envs/py37_withR/share/r-ichorcna-0.1.0.20180710-0/extdata/map_hg19_1000kb.wig  \
--centromere /mnt/X500/farmers/linhx/bin/miniconda3/envs/py37_withR/share/r-ichorcna-0.1.0.20180710-0/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
--exons.bed /mnt/X500/farmers/linhx/bin/ref/region_bed/combine.black.merge.complement.nochr.bed \
--rmCentromereFlankLength 1e5 \
--includeHOMD False \
--estimateNormal True \
--estimatePloidy True \
--estimateScPrevalence True \
--scStates "c(1,3)" \
--txnE 0.9999999 \
--txnStrength 1e7 \
--outDir $out_dir \
--libdir /mnt/X500/farmers/linhx/bin/cna/ichorCNA/R \
--chrNormalize "c(1:22)" \
--chrTrain "c(1:22)" \
--chrs "c(1:22)"


