# median_filte

Программа для получения и применения медианных коррекций в площадке
Пример median2df gsc_scos_psc_65557.dat gsc_scos_psc_65557.dat x1=11 y1=12 dra1=25 ddec1=26 x2=11 y2=12 ra2=1 dec2=2 --nx=1000 --ny=1000 --dx=3 --dy=3 --cs=5 -v > result.txt
#!/bin/bash
PID=196610
while [ $PID -lt 197431 ]
do
echo "PID=$PID"
median2df /mnt/work/gsc_scosmos/N/gsc_BRI_psc/gsc_BRI_psc_$PID.dat /mnt/work/gsc_scosmos/N/gsc_BRI/pairs_gsc23_$PID.dat x1=11 y1=12 dra1=43 ddec1=44 x2=11 y2=12 ra2=1 dec2=2 --nx=1000 --ny=1000 --dx=7 --dy=7 --cs=12 > /mnt/work/gsc_scosmos/N/GScor/gsc_cor_$PID.dat
#median2df --help
PID=$[$PID+1]
done

Программа для получения медианных коррекций в площадке, принимает ra,dec, x,y, dra.ddec возращает на 6 + 2 столбца cdra, cddec
cut -f 1,2,11,12,25,26 gsc_scos_psc_65557.dat |median2d --nx=1000 --ny=1000 --dx=7 --dy=7 --cs=20
plot '< cut -f 1,2,11,12,25,26 gsc_scos_psc_65557.dat |median2d --nx=1000 --ny=1000 --dx=7 --dy=7 --cs=12' u 1:8 w d

