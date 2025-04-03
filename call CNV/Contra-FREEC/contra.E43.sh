perl /mnt/X500/farmers/zhuxl/bin/getcnv/control_freec-config-1M.pl -c /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup -n /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43ND_normal_sort_markdup_recald.mpileup -o /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43 -s XY -f /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD.control_freec.config 
/mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/src/freec -conf /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD.control_freec.config > /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD.control_freec.log 2>&1 
cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/assess_significance.R | R --slave --args /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup_CNVs /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup_ratio.txt
cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/makeGraph.R | R --slave --args 2 /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup_ratio.txt /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup_BAF.txt
Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/control_freec_exon.R /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD_case_sort_markdup_recald.mpileup_CNVs.p.value.txt /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/ncbi_anno_rel104_db_hg19_EXN_Primary.tab /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/15_cancer_analysis/contra/E43/E43TD.CNV.anno.xls
