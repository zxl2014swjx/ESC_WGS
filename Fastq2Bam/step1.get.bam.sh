less   sample.list |while read i; do readlink -f 01.MP/$i[TN]D/recaled.bam|xargs ; done|awk '{print $2"\t"$1}'|less  > bam.list
