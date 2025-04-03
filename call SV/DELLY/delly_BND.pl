use File::Basename;

my @dir = glob("$ARGV[0]/*");
foreach $d(@dir){
	next if(!-d $d);
	my $s = basename $d;
	open OUT,">P$s.SV.sh";
	my $case_bam = (glob "$d/case/5_recal_bam/*.bam")[0];
	my $control_bam = (glob "$d/normal/5_recal_bam/*.bam")[0];
	mkdir -p "$d/SV/delly/BND/" if(!-d "$d/SV/delly/BND/");
	print OUT "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly call -t BND -x /mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/excludeTemplates/human.hg19.excl.tsv -o $d/SV/delly/BND/$s.BND.bcf -g /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa $case_bam $control_bam  > $d/SV/delly/BND/$s.BND.bcf.log 2>&1\n";
	print OUT "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly filter -t BND -f somatic -o $d/SV/delly/BND/$s.BND.pre.bcf -s /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/samples.tsv $d/SV/delly/BND/$s.BND.bcf\n";
	print OUT "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly call -t BND -g /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/tgp_phase2_flat/hs37d5.fa -v $d/SV/delly/BND/$s.BND.pre.bcf -o $d/SV/delly/BND/$s.BND.geno.bcf -x /mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/excludeTemplates/human.hg19.excl.tsv $case_bam $control_bam > $d/SV/delly/BND/$s.BND.geno.bcf.log 2>&1\n";
	print OUT "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/delly filter -t BND -f somatic -o $d/SV/delly/BND/$s.BND.somatic.bcf -s /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/samples.tsv $d/SV/delly/BND/$s.BND.geno.bcf > $d/SV/delly/BND/$s.BND.somatic.bcf.log 2>&1\n";
	print OUT "/mnt/X500/farmers/suyao/NL200/tools/Delly2/delly/src/bcftools/bcftools view $d/SV/delly/BND/$s.BND.somatic.bcf > $d/SV/delly/BND/$s.BND.somatic.vcf 2> $d/SV/delly/BND/$s.BND.somatic.vcf.log\n";
	print OUT "perl /mnt/X500/farmers/baijing/Program/Share/wes/V1.0/delly2_BND_parse.pl -i $d/SV/delly/BND/$s.BND.somatic.vcf -o $d/SV/delly/BND/$s.BND.somatic.tab\n";
	print OUT "Rscript /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/delly2_BND_exon.R $d/SV/delly/BND/$s.BND.somatic.tab /mnt/X500/farmers/suyao/NL200/analysis2/wes/V1.0/ncbi_anno_rel104_db_hg19_EXN_IVS_primary.tab $d/SV/delly/BND/$s.SV.anno.xls\n";
	close OUT;
	
}

