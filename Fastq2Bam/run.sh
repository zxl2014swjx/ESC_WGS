cat lib.txt|xargs -I {} echo 'name=`echo "{}" |cut -d " " -f1 `;echo  "{}" >cfg.$name'|sh &
for i in  E10ND E10TD E15ND E15TD E18ND E18TD E1ND E1TD E23ND E23TD E24ND E24TD E26ND E26TD E2ND E2TD E31ND E31TD E33ND E33TD E35ND E35TD E47ND E47TD E51ND E51TD;do echo /repos/lijun/pip/mp/sentieon/MP4Sentieon cfg.$i\ $i /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/13_normal_fastq >>map.smaple.sh ;done
perl /mnt/X500/farmers/wuaw/software/bin/qsub-sge.pl --maxjob 40 --jobprefix qsubjob --resource vf=4G:p=2  --getmem map.smaple.sh &
wait