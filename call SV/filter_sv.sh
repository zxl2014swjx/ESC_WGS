for i in E1 E10 E15 E18 E2 E23 E24 E26 E31 E33 E35 E47 E51 E11 E21 E25 E43 E44 E52 E53 EL07 EL08 EL12 EL16 EL18 EL19 EL26 EL29;do echo "/mnt/X500/farmers/luozw/bin/anaconda/envs/python27/bin/python /products/repos/prod/Akso/OncoWESuper/bin/variantCalling/filtersv.py -c /products/repos/prod/Akso/OncoWESuper/program/NoahCare/config/nc_we4/sv/sv.cfg -g /products/repos/prod/Akso/OncoWESuper/program/NoahCare/config/nc_we4/anno/tr_nc.lst /mnt/X500/farmers/zhuxl/Project/ESC_Cheng_WGS/result/casesv/$i\ND_$i\TD.callsv.anno.csv" >>filter.sh ;done &
