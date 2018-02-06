#!/bin/bash

rm -rf errorbond.info
rm -rf errorbond.xyz

for i in `seq 1 24`
do
    echo " " >> errorbond.info
    idn=$(awk -v n=$i '{if(NR==n) {print $1;exit;}}' error)
    awk -v n=$i '{if(NR==n) {print $0;exit;}}' error >> errorbond.info 
    x=$(awk -v k=$idn '{if($1==k) {print $3;exit;}}' trueglasszif4_3.xyz)
    y=$(awk -v k=$idn '{if($1==k) {print $4;exit;}}' trueglasszif4_3.xyz)
    z=$(awk -v k=$idn '{if($1==k) {print $5;exit;}}' trueglasszif4_3.xyz)
    #echo "idn=$idn x=$x y=$y z=$z"
    awk -v x=$x -v y=$y -v z=$z '{
    if(x==$3 && y==$4 && z==$5){print $1,"C_"$2,$3,$4,$5;}
    else if((x-$3)*(x-$3)+(y-$4)*(y-$4)+(z-$5)*(z-$5)<10){print $0;}
    else {;}
    }' trueglasszif4_3.xyz >> errorbond.xyz
    awk -v n=$idn '{if($2==n || $3==n) print "row:"NR, $0;}' glassbondlist >> errorbond.info

done
