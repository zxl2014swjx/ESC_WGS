#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"i=s",'s=s');

my $usage =<< "USAGE";
Program: $0 -i InputDir -s SampleGenderFile
            -i input directory (output directory of BNC)
            -s sample gender, 'XY' for man and 'XX' for woman, one sample per line,separated by tab.
               example:
                   sample1  XY
                   sample2  XX
            -h help
USAGE

die $usage if(!$opts{'i'} || $opts{'h'});

my $input_dir = $opts{'i'};
my $gender_file = $opts{'s'};
my $ref_seq = '/mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa';
my $ref_arm = '/mnt/X500/farmers/suyao/NL200/tools/VarScan2/scripts/hg19.ref.arm';
my $chip_bed = '/mnt/NL200/suyao/tools/BNC/BNC/program/NoahCare/config/nc_e/regs/flk50_chip.bed';
my $exclude_bed = '/mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/wgEncodeDacMapabilityConsensusExcludable.txt';
my $exon_file = '/mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/ncbi_anno_rel104_db_hg19_EXN_Primary.tab';
my $exon_intron_file = '/mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/ncbi_anno_rel104_db_hg19_EXN_IVS_primary.tab';

my %gender;
open(G, $gender_file);
while(<G>){
    chomp;
    my ($sample,$sex) = split /\t/;
    $gender{$sample} = $sex;
#    print "$sample\t$sex\n";
}
close(G);

opendir I, $input_dir;
while(my $sample = readdir I){
    next if($sample =~ /^\.+$/);
    my $sample_dir = $input_dir . '/' . $sample;
    my $cancer_dir = $sample_dir . '/output/cancer';
    my $normal_dir = $sample_dir . '/output/normal';
    next unless(-d $cancer_dir);
    next unless(-d $normal_dir);
 #   print "$sample\n";
    if(!exists($gender{$sample})){
        print "No gender information for $sample!\n";
        next;
    }
    my $sample_gender = $gender{$sample};
    
    my $cancer_recal_dir = $cancer_dir . '/4_recal_bam';
    my $normal_recal_dir = $normal_dir . '/4_recal_bam';
#Modified by baijing at 20180427
#    my $cancer_recal_bam = $cancer_recal_dir . '/' . $sample . '_cancer_sort_markdup_realign_recal.bam';
#    my $normal_recal_bam = $normal_recal_dir . '/' . $sample . '_normal_sort_markdup_realign_recal.bam';
    my $cancer_recal_bam = $cancer_recal_dir . '/' . $sample . '_case_sort_markdup_recald.bam';
    my $normal_recal_bam = $normal_recal_dir . '/' . $sample . '_normal_sort_markdup_recald.bam';
    my $cancer_recal_mpileup = $cancer_recal_dir . '/' . $sample . '_case_sort_markdup_recald.mpileup';
    my $normal_recal_mpileup = $normal_recal_dir . '/' . $sample . '_normal_sort_markdup_recald.mpileup';
    
    if(!-e $cancer_recal_bam){
        print "$cancer_recal_dir not exists!\n";
        next;
    }
    if(!-e $normal_recal_bam){
        print "$normal_recal_dir not exists!\n";
        next;
    }
    #my $wes_dir = $sample_dir . '/WES';
#    my $snv_dir = $sample_dir . '/SNV';
    my $cnv_dir = $sample_dir . '/CNV';
#    my $sv_dir = $sample_dir . '/SV';
    #mkdir $wes_dir;
#   mkdir $snv_dir unless(-d $snv_dir);
    mkdir $cnv_dir unless(-d $cnv_dir);
#   mkdir $sv_dir unless(-d $sv_dir);
#    my $snv_varscan2_dir = $snv_dir . '/VarScan2';
#    mkdir $snv_varscan2_dir;
#   my $cnv_varscan2_dir = $cnv_dir . '/VarScan2';
#   mkdir $cnv_varscan2_dir;
    #my $cnv_Platypus_dir = $cnv_dir . '/Platypus';
    #mkdir $cnv_Platypus_dir;
    
#    my $varscan2_somatic = $snv_varscan2_dir . '/' . $sample . '.somatic';
    #my $varscan2_copynumber = $cnv_varscan2_dir . '/cnv';
    my $snv_sh = $sample_dir . '/' . $sample . '.500kcnv.sh';
    
    open(SNVSH, ">$snv_sh");
    print SNVSH "#!/bin/bash\n";
=pod
    print SNVSH "samtools mpileup -f $ref_seq -q 20 -Q 20 -B -C 50 $cancer_recal_bam > $cancer_recal_mpileup\n";
    print SNVSH "samtools mpileup -f $ref_seq -q 20 -Q 20 -B -C 50 $normal_recal_bam > $normal_recal_mpileup\n";
#-C --adjust-MQ INT 用于降低比对质量的系数，如果reads中含有过多的错配。不能设置为零。BWA推荐值为50。
    print SNVSH "################## SNV ###########################\n\n\n";
    
    print SNVSH "################ SNV VarScan2 ###############\n\n";
    my $varscan2_somatic_log = $varscan2_somatic . '.log';
    print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar somatic $normal_recal_mpileup $cancer_recal_mpileup $varscan2_somatic --output-vcf 1 --min-coverage-normal 10 --min-var-freq 0.01 --tumor-purity 0.5 > $varscan2_somatic_log 2>&1 \n";
#--output-vcf 1 生成结果vcf
    my $varscan2_somatic_indel = $varscan2_somatic . '.indel.vcf';
    my $varscan2_somatic_snp = $varscan2_somatic . '.snp.vcf';
    my $varscan2_somatic_indel_log = $varscan2_somatic . '.indel.log';
    my $varscan2_somatic_snp_log = $varscan2_somatic . '.snp.log';
    print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar processSomatic --min-tumor-freq 0.07 --max-normal-freq 0.02 --p-value 0.05 $varscan2_somatic_indel > $varscan2_somatic_indel_log 2>&1 \n";
    print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar processSomatic --min-tumor-freq 0.07 --max-normal-freq 0.02 --p-value 0.05 $varscan2_somatic_snp > $varscan2_somatic_snp_log 2>&1 \n";
#processSomatic          Isolate Germline/LOH/Somatic calls from output
    my $varscan2_somatic_snp_hc = $varscan2_somatic . '.snp.Somatic.hc.vcf';
    my $varscan2_somatic_snp_hc_filter = $varscan2_somatic . '.snp.Somatic.hc.filter.vcf';
    my $varscan2_somatic_snp_hc_filter_log = $varscan2_somatic . '.snp.Somatic.hc.filter.log';
    print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar somaticFilter $varscan2_somatic_snp_hc --indel-file $varscan2_somatic_indel --output-file $varscan2_somatic_snp_hc_filter > $varscan2_somatic_snp_hc_filter_log 2>&1 \n";
#somaticFilter           Filter somatic variants for clusters/indels
    my $varscan2_somatic_snp_hc_filter_readcount = $varscan2_somatic . '.snp.Somatic.hc.filter.readcount';
    print SNVSH "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/somatic_snp_bed.R $varscan2_somatic_snp_hc_filter\n";
    my $varscan2_somatic_snp_hc_filter_pos = $varscan2_somatic . '.snp.Somatic.hc.filter.vcf.pos';
    print SNVSH "bam-readcount -f $ref_seq --site-list $varscan2_somatic_snp_hc_filter_pos $cancer_recal_bam > $varscan2_somatic_snp_hc_filter_readcount 2>/dev/null\n";
    my $varscan2_somatic_snp_hc_filter_fpfilter = $varscan2_somatic . '.snp.Somatic.hc.filter.fpfilter.vcf';
    my $varscan2_somatic_snp_hc_filter_fpfilter_log = $varscan2_somatic . '.snp.Somatic.hc.filter.fpfilter.log';
    print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar fpfilter $varscan2_somatic_snp_hc_filter $varscan2_somatic_snp_hc_filter_readcount --output-file $varscan2_somatic_snp_hc_filter_fpfilter > $varscan2_somatic_snp_hc_filter_fpfilter_log 2>&1 \n";
#fpfilter                Apply the false-positive filter
    print SNVSH "### SNV Screeen ###\n";
    print SNVSH "python /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/screen_somatic_snv_vcf_20180611.py -i $sample_dir\n";
    my $merge_vcf = $varscan2_somatic . '.final.vcf';
    my $merge_final_file = $varscan2_somatic . '.filt.final';
    #my $merge_file_vcf = $varscan2_somatic . '.filt.final.vcf';
    my $final_file_vcf = $snv_dir . '/' . $sample . '.final.snv.vcf';
    print SNVSH "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/somatic_screen.R $merge_vcf $chip_bed $exclude_bed\n";
    print SNVSH "cat <(grep '^#' $merge_vcf) $merge_final_file > $final_file_vcf\n";
    my $annobed = $final_file_vcf;
    my $annolog = $final_file_vcf;
    $annobed =~ s/\.vcf/.anno.bed/;
    $annolog =~ s/\.vcf/.anno.bed.log/;
    print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/bin/perl /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/bin/anno/NCanno.pl /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/config/nc_e/anno/config_nc.pl $final_file_vcf -t vcf -f -q -n 5 -b 500 -o $annobed > $annolog\n";
     
    print SNVSH "################## CNV ##############################\n\n\n";
    
    #print SNVSH "############## CNV VarScan2 #######################\n\n";
    #my $cancer_recal_bam_stat = $cancer_recal_bam;
    #my $cancer_recal_bam_read_length_stat = $cancer_recal_bam;
    #my $normal_recal_bam_stat = $normal_recal_bam;
    #my $normal_recal_bam_read_length_stat = $normal_recal_bam;
    #$cancer_recal_bam_stat =~ s/bam$/stat/;
    #$normal_recal_bam_stat =~ s/bam$/stat/;
    #$cancer_recal_bam_read_length_stat =~ s/bam$/read.length/;
    #$normal_recal_bam_read_length_stat =~ s/bam$/read.length/;
    #print SNVSH "samtools flagstat $cancer_recal_bam > $cancer_recal_bam_stat\n";
    #print SNVSH "samtools flagstat $normal_recal_bam > $normal_recal_bam_stat\n";
    #print SNVSH "samtools view $cancer_recal_bam | head -n1000 | awk -F'\t' '$6!~/H/' | cut -f10 | awk '{print length($0)}' > $cancer_recal_bam_read_length_stat\n";
    #print SNVSH "samtools view $normal_recal_bam | head -n1000 | awk -F'\t' '$6!~/H/' | cut -f10 | awk '{print length($0)}' > $normal_recal_bam_read_length_stat\n";
    #my $cancer_unique_reads = &uni_reads($cancer_recal_bam_stat);
    #my $normal_unique_reads = &uni_reads($normal_recal_bam_stat);
    #my $cancer_read_length = &read_len($cancer_recal_bam_read_length_stat);
    #my $normal_read_length = &read_len($normal_recal_bam_read_length_stat);
    #my $data_ratio = $normal_unique_reads * $normal_read_length / ($cancer_unique_reads * $cancer_read_length);
    #my $varscan2_copynumber_log = $cnv_varscan2_dir . '/cnv.log';
    #print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar copynumber $normal_recal_mpileup $cancer_recal_mpileup $varscan2_copynumber --min-coverage=8 --min-segment-size=50 --data-ratio $data_ratio > $varscan2_copynumber_log 2>&1 \n";
    #my $varscan2_copynumber_output = $varscan2_copynumber . '.copynumber';
    #my $varscan2_copyCaller_output = $varscan2_copynumber . '.called';
    #my $varscan2_copyCaller_log = $varscan2_copynumber . '.called.log';
    #print SNVSH "java -jar /mnt/X500/farmers/suyao/NL200/tools/VarScan2/VarScan.v2.3.9.jar copyCaller $varscan2_copynumber_output --output-file $varscan2_copyCaller_output > $varscan2_copyCaller_log 2>&1 \n";
    #print SNVSH "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/run_DNAcopy.R $varscan2_copyCaller_output\n";
    #my $varscan2_copyCaller_DNAcopy = $varscan2_copynumber . '.called.DNAcopy';
    #my $varscan2_copyCaller_Merge = $varscan2_copynumber . '.called.MergeSeg';
    #my $varscan2_copyCaller_Merge_log = $varscan2_copynumber . '.called.MergeSeg.log';
    #print SNVSH "perl /mnt/X500/farmers/suyao/NL200/tools/VarScan2/scripts/varscan-copynumber-mergeSegments.pl $varscan2_copyCaller_DNAcopy  --ref-arm-sizes $ref_arm --output-basename $varscan2_copyCaller_Merge >$varscan2_copyCaller_Merge_log 2>&1 \n";
    #my $varscan2_copyCaller_Merge_out = $varscan2_copynumber . '.called.MergeSeg.events.tsv';
=cut    
    print SNVSH "################# CNV Control-FREEC #####################################\n\n";
    my $cnv_Control_FREEC_dir = $cnv_dir . '/Control_FREEC_500k';
    mkdir $cnv_Control_FREEC_dir;
    my $control_freec_config = $cnv_Control_FREEC_dir . '/' . $sample . '.control_freec.config';
    my $control_freec_log = $cnv_Control_FREEC_dir . '/' . $sample . '.control_freec.log';
    my $control_freec_CNV = $cnv_Control_FREEC_dir . '/' . $sample . '_case_sort_markdup_recald.mpileup_CNVs';#GP10-1B_case_sort_markdup_recald.mpileup_CNVs
    my $control_freec_ratio = $cnv_Control_FREEC_dir . '/' . $sample . '_case_sort_markdup_recald.mpileup_ratio.txt';#GP10-1B_case_sort_markdup_recald.mpileup_ratio.txt
    my $control_freec_BAF = $cnv_Control_FREEC_dir . '/' . $sample . '_case_sort_markdup_recald.mpileup_BAF.txt';#GP10-1B_case_sort_markdup_recald.mpileup_BAF.txt
    my $control_freec_CNV_pvalue = $cnv_Control_FREEC_dir . '/' . $sample . '_case_sort_markdup_recald.mpileup_CNVs.p.value.txt';#GP10-1B_case_sort_markdup_recald.mpileup_CNVs.p.value.txt
    my $control_freec_CNV_anno = $cnv_dir . '/' . $sample . '.CNV.anno.xls';
    print SNVSH "perl /mnt/X500/farmers/zhaogf/Project/PanCancer_Cheng/Gastric_WGS/latitude/control-freec/control_freec-config-500k.pl -c $cancer_recal_mpileup -n $normal_recal_mpileup -o $cnv_Control_FREEC_dir -s $sample_gender -f $control_freec_config\n";
    print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/src/freec -conf $control_freec_config > $control_freec_log 2>&1 \n";
    print SNVSH "cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/assess_significance.R | R --slave --args $control_freec_CNV $control_freec_ratio\n";
    print SNVSH "cat /mnt/X500/farmers/suyao/NL200/tools/Control-FREEC/FREEC-10.6/scripts/makeGraph.R | R --slave --args 2 $control_freec_ratio $control_freec_BAF\n";
    print SNVSH "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/control_freec_exon.R $control_freec_CNV_pvalue $exon_file $control_freec_CNV_anno\n";
    
=pod    
    #print SNVSH "################ CNV Platypus ########################\n\n";
    
    #my $cnv_Platypus_cancer_output = $cnv_dir . '/Platypus/' . $sample . '.Platypus.cancer.snp';
    #my $cnv_Platypus_cancer_output_log = $cnv_dir . '/Platypus/' . $sample . '.Platypus.cancer.snp.log';
    #my $cnv_Platypus_normal_output = $cnv_dir . '/Platypus/' . $sample . '.Platypus.normal.snp';
    #my $cnv_Platypus_normal_output_log = $cnv_dir . '/Platypus/' . $sample . '.Platypus.normal.snp.log';
    
    #print SNVSH "python /mnt/X500/farmers/suyao/NL200/tools/Platypus/Platypus_0.8.1/Platypus.py callVariants --refFile=$ref_seq --bamFiles=$cancer_recal_bam --genIndels=0 --output=$cnv_Platypus_cancer_output > $cnv_Platypus_cancer_output_log 2>&1 \n";
    #print SNVSH "python /mnt/X500/farmers/suyao/NL200/tools/Platypus/Platypus_0.8.1/Platypus.py callVariants --refFile=$ref_seq --bamFiles=$normal_recal_bam --genIndels=0 --output=$cnv_Platypus_normal_output > $cnv_Platypus_normal_output_log 2>&1 \n";
       
    print SNVSH "################## SV ##############################\n\n\n";
    my $delly_dir = $sv_dir . '/delly';
    mkdir $delly_dir;
    #my @svTypes = qw/DEL DUP INV BND INS/;
    my @svTypes = qw/BND/;
    my $hg19_excl = "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/excludeTemplates/human.hg19.excl.tsv";
    my $delly_sample_tsv = "/mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/samples.tsv";
    foreach my $svType (@svTypes){
        my $svType_dir = $delly_dir . '/' . $svType;
        mkdir $svType_dir;
        my $bcffile = $svType_dir . '/' . $sample . '.' . $svType . '.bcf';
        my $bcffile_log = $svType_dir . '/' . $sample . '.' . $svType . '.bcf.log';
        print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly call -t $svType -x $hg19_excl -o $bcffile -g $ref_seq $cancer_recal_bam $normal_recal_bam > $bcffile_log 2>&1 \n";
        my $prebcffile = $svType_dir . '/' . $sample . '.' . $svType . '.pre.bcf';
        my $prebcffile_log = $svType_dir . '/' . $sample . '.' . $svType . '.pre.bcf.log';
        print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly filter -t $svType -f somatic -o $prebcffile -s $delly_sample_tsv $bcffile > $prebcffile_log 2>&1 \n";
        my $genobcffile = $svType_dir . '/' . $sample . '.' . $svType . '.geno.bcf';
        my $genobcffile_log = $svType_dir . '/' . $sample . '.' . $svType . '.geno.bcf.log';
        print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly call -t $svType -g $ref_seq -v $prebcffile -o $genobcffile -x $hg19_excl $cancer_recal_bam $normal_recal_bam > $genobcffile_log 2>&1 \n";
        my $somaticbcffile = $svType_dir . '/' . $sample . '.' . $svType . '.somatic.bcf';
        my $somaticbcffile_log = $svType_dir . '/' . $sample . '.' . $svType . '.somatic.bcf.log';
        print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly filter -t $svType -f somatic -o $somaticbcffile -s $delly_sample_tsv $genobcffile > $somaticbcffile_log 2>&1 \n";
        my $somaticvcffile = $svType_dir . '/' . $sample . '.' . $svType . '.somatic.vcf';
        my $somaticvcffile_log = $svType_dir . '/' . $sample . '.' . $svType . '.somatic.vcf.log';
        print SNVSH "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/bcftools/bcftools view $somaticbcffile > $somaticvcffile 2> $somaticvcffile_log \n";
    }
    
    my $bnd_vcf = $delly_dir . '/BND/' . $sample . '.BND.somatic.vcf';
    my $bnd_tab = $delly_dir . '/BNC/' . $sample . '.BND.somatic.tab';
    my $bnd_anno = $sv_dir . '/' . $sample . '.SV.anno.xls';
    print SNVSH "perl /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/delly2_BND_parse.pl -i $bnd_vcf -o $bnd_tab\n";
    print SNVSH "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/delly2_BND_exon.R $bnd_tab $exon_intron_file $bnd_anno\n";

    close(SNVSH); 
=cut
}
closedir I;

sub uni_reads(){
    my $f = shift;
    my ($total, $mapped, $dup);
    open(F, $f);
    while(<F>){
        chomp;
        my @tmp = split /\s/;
        if(/in total/){
            $total = $tmp[0];
        }elsif(/mapped\s\(/){
            $mapped = $tmp[0];
        }elsif(/duplicates/){
            $dup = $tmp[0];
        }
    }
    close(F);
    my $uni_reads_count = $mapped * (1 - $dup / $total);
    return $uni_reads_count;
}

sub read_len(){
    my $f = shift;
    my $read_length;
    open(F, $f);
    while(<F>){
        chomp;
        $read_length = $_;
        last;
    }
    close(F);
    return $read_length;
}
