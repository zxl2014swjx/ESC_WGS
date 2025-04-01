The pipeline and analysis codes of low depths whole genome sequence (WGS) for 28 esophageal cancer (ESC) samples.

## The copy number variations (CNV) were detected utilizing GATK following steps below：
step_1_get_tumor_readcont

step_2_getPON

step_3_get_DenoiseReadCounts

step_4_PlotDenoisedCopyRatios

step_5_get_allelic_work

step_6_ModelSegments

step_7_CallCopyRatioSegments

step_8_PlotModeledSegments

## The structure variations (SV) were detected utilizing Meerkat following steps below：
step1_pre_process

step2_meerkat

step3_mechanism

step4_somatic_filting

step 5 germline_calling

step_6 annotation

step_7 transform2vcf

step_8 Calculate allele frequency

# The GeneFuse was employed to call gene fuse events.
# The Strelka was utilized to call somatic mutations.
# The ichorCNA was used to calculated the ploidy of tumor samples.
# The Facets was used to call fragment variations.







