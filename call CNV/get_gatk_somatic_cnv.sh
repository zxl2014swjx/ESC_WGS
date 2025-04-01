#!/bin/bash

if [ $# -lt 5 ] ; then
	echo -e "\n\tsh $0 <tumor_bam> <normal_bam> <tumor_name> <PON.file|no_pon> <out_dir>\n"
	exit
fi

#parameters
tumor_bam=$1
normal_bam=$2
tumor_name=$3
pon_file=$4
out_dir=$5

if [ ! -d $out_dir ] ; then
	mkdir -p $out_dir
fi
cd $out_dir

###config
interval_list="/mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/cnv_result/gatk_somatic_cnv/hs37d5.fa_autosomal.interval_list"
gc_file="/mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/cnv_result/gatk_somatic_cnv/hs37d5.fa_autosomal.interval_list.GC.tsv"
ref_dict="/mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/cnv_result/gatk_somatic_cnv/hg19_ref/hs37d5.dict"
ref_fa="/mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/cnv_result/gatk_somatic_cnv/hg19_ref/hs37d5.fa"

#step_1 get interval read_counts
if [ ! -f $out_dir/normal_bam.counts.hdf5 ] ; then 
	gatk CollectReadCounts \
	-I $normal_bam \
	-L $interval_list \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O $out_dir/normal_bam.counts.hdf5
fi

if [ ! -f $out_dir/tumor_bam.counts.hdf5 ] ; then
	gatk CollectReadCounts \
	-I $tumor_bam \
	-L $interval_list \
	--interval-merging-rule OVERLAPPING_ONLY \
	-O $out_dir/tumor_bam.counts.hdf5
fi
echo -e "CollectReadCounts done ..."


#step_2 get pon
if [ ! -f $pon_file ] ; then 
	pon_file="$out_dir/cnv_pon.hdf5"
	gatk CreateReadCountPanelOfNormals \
	-I $out_dir/normal_bam.counts.hdf5 \
	--minimum-interval-median-percentile 5.0 \
	--annotated-intervals $gc_file \
	-O $pon_file
	echo -e "PON_file:$pon_file"
else 
	echo -e "PON_file:$pon_file"
fi

#step_3 DenoiseReadCounts
if [ ! -f $out_dir/tumor_bam.counts.hdf5.denoisedCR.tsv ] ; then
	gatk DenoiseReadCounts \
	-I $out_dir/tumor_bam.counts.hdf5 \
	--count-panel-of-normals $pon_file \
	--annotated-intervals $gc_file \
	--standardized-copy-ratios $out_dir/tumor_bam.counts.hdf5.standardizedCR.tsv  \
	--denoised-copy-ratios $out_dir/tumor_bam.counts.hdf5.denoisedCR.tsv 
fi

echo -e "DenoiseReadCounts done ..."

if [ ! -f $out_dir/$tumor_name.denoisedLimit4.png ] ; then
#step_4 PlotDenoisedCopyRatios
	gatk PlotDenoisedCopyRatios \
	--standardized-copy-ratios $out_dir/tumor_bam.counts.hdf5.standardizedCR.tsv \
	--denoised-copy-ratios $out_dir/tumor_bam.counts.hdf5.denoisedCR.tsv \
	--sequence-dictionary $ref_dict \
	--minimum-contig-length 1000000 \
	--output $out_dir \
	--output-prefix $tumor_name
fi
echo -e "PlotDenoisedCopyRatios done ..."

#step_5 CollectAllelicCounts
if [ ! -f normal_allelicCounts.tsv ] ; then
	gatk CollectAllelicCounts \
	-L $interval_list \
	-I $normal_bam \
	-R $ref_fa \
	-O normal_allelicCounts.tsv
fi
if [ ! -f tumor_allelicCounts.tsv ] ; then
	gatk CollectAllelicCounts \
	-L $interval_list \
	-I $tumor_bam \
	-R $ref_fa \
	-O tumor_allelicCounts.tsv
fi
echo -e "CollectAllelicCounts done ..."

#step_6 ModelSegments
if [ ! -f $out_dir/$tumor_name.cr.seg ] ; then
	gatk  ModelSegments \
	--denoised-copy-ratios $out_dir/tumor_bam.counts.hdf5.denoisedCR.tsv \
        --allelic-counts tumor_allelicCounts.tsv \
        --normal-allelic-counts normal_allelicCounts.tsv \
	--output $out_dir \
	--output-prefix $tumor_name
fi

echo -e "ModelSegments done ..."

#step_7 CallCopyRatioSegments
if [ ! -f $out_dir/$tumor_name.cr.seg.called.seg ] ; then
	gatk CallCopyRatioSegments \
	--input $out_dir/$tumor_name.cr.seg \
	--output $out_dir/$tumor_name.cr.seg.called.seg
fi
echo  -e "CallCopyRatioSegments done ..."

#step_8 PlotModeledSegments
if [ ! -f $out_dir/$tumor_name.modeled.png ] ; then
	gatk PlotModeledSegments \
	--denoised-copy-ratios $out_dir/tumor_bam.counts.hdf5.denoisedCR.tsv  \
	--allelic-counts tumor_allelicCounts.tsv \
	--segments $out_dir/$tumor_name.modelFinal.seg \
	--sequence-dictionary $ref_dict  \
	--minimum-contig-length 1000000 \
	--output $out_dir  \
	--output-prefix $tumor_name
fi
echo  -e "PlotModeledSegments done ..."

