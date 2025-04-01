mkdir -p /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND || exit 1
ln -sf /mnt/X500/farmers/lijun/prj/baijin_ES_cancer_Cheng/01.MP/E43TD/recaled.bam /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD || exit 1
ln -sf /mnt/X500/farmers/lijun/prj/baijin_ES_cancer_Cheng/01.MP/E43ND/recaled.bam /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND || exit 1
ln -sf /mnt/X500/farmers/lijun/prj/baijin_ES_cancer_Cheng/01.MP/E43TD/recaled.bam.bai /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD || exit 1
ln -sf /mnt/X500/farmers/lijun/prj/baijin_ES_cancer_Cheng/01.MP/E43ND/recaled.bam.bai /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND || exit 1
export LD_LIBRARY_PATH="/mnt/X500/farmers/linhx/bin/software/Meerkat/Meerkat/src/mybamtools/lib/:/mnt/X500/farmers/linhx/bin/lib/zlib-1.2.8/local/lib"
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/pre_process.pl -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.bam -l 0 -I /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/hg19_bwa_idx -A /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19 -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -W /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -k 1500 -u 0 -s 20 -q 15" > /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_pre.sh
export LD_LIBRARY_PATH="/mnt/X500/farmers/linhx/bin/software/Meerkat/Meerkat/src/mybamtools/lib/:/mnt/X500/farmers/linhx/bin/lib/zlib-1.2.8/local/lib"
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/pre_process.pl -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam -l 0 -I /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/hg19_bwa_idx -A /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19 -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -W /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -k 1500 -u 0 -s 20 -q 15" >> /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_pre.sh
cat /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_pre.sh |xargs -i -P 2 sh -c "{}" || exit 1
mv /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.blacklist.gz /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.blacklist.gz.backup || exit 1
ln -fs /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.blacklist.gz /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/ || exit 1
export LD_LIBRARY_PATH="/mnt/X500/farmers/linhx/bin/software/Meerkat/Meerkat/src/mybamtools/lib/:/mnt/X500/farmers/linhx/bin/lib/zlib-1.2.8/local/lib"
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts//meerkat.pl -l 0 -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.bam -Q 10 -s 20 -m 0 -F /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/hg19_fasta -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -W /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -B /mnt/X500/farmers/linhx/bin/software/Meerkat//blast-legacy-2.2.26-1/bin -d 5 -p 3 -o 1" > /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_meerkat.sh
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts//meerkat.pl -l 0 -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam -Q 10 -s 20 -m 0 -F /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/hg19_fasta -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -W /mnt/X500/farmers/linhx/bin/software/Meerkat/opt -B /mnt/X500/farmers/linhx/bin/software/Meerkat//blast-legacy-2.2.26-1/bin -d 5 -p 3 -o 1" >> /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_meerkat.sh
cat /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_meerkat.sh |xargs -i -P 2 sh -c "{}" || exit 1
export LD_LIBRARY_PATH="/mnt/X500/farmers/linhx/bin/software/Meerkat/Meerkat/src/mybamtools/lib/:/mnt/X500/farmers/linhx/bin/lib/zlib-1.2.8/local/lib"
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/mechanism.pl -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.bam -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19" > /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_mechanism.sh
echo "perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/mechanism.pl -b /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19" >> /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_mechanism.sh
cat /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/multi_mechanism.sh |xargs -i -P 2 sh -c "{}" || exit 1
grep -v MT /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.variants | grep -v GL000 |grep -v NC_007605 |grep -v hs37d5 > /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.variants || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.variants -F /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.discord -l 1000 -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_b.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -I /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.isinfo -n 1 -D 5 -Q 10 -B /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_b.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_c.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -u 1 -Q 10 -B /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_c.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_d.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -f 1 -Q 10 -B /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43ND/recaled.bam || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_d.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_e.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -e 1 -D 5 -Q 10 -B /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.bam -I /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.isinfo || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_e.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_f.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -z 1 || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/somatic_sv.pl -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt  -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter_f.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter.final.variants -R /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/rmsk-hg19.txt -d 40 -t 20 || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts//fusions.pl -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD//autoXY.recaled.somatic.filter.final.variants -G /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/refGene_hg19_sorted.txt || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/meerkat2vcf.pl -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter.final.variants  -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter.final.variants.vcf -H /mnt/X500/farmers/linhx/bin/software/Meerkat//head.vcf -F /mnt/X500/farmers/linhx/bin/software/Meerkat/db/hg19/hg19_fasta || exit 1
perl /mnt/X500/farmers/linhx/bin/software/Meerkat//Meerkat/scripts/discon.pl -d 5 -Q 10 -i /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter.final.variants -o /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/autoXY.recaled.somatic.filter.final.variants.rp -B /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.bam  -C /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.clusters  -I /mnt/X500/farmers/luozw/factory/20190821_15WGS_sv_cnv/E43TD_E43ND/E43TD/recaled.isinfo  -S /mnt/X500/farmers/linhx/bin/software/Meerkat/opt || exit 1
