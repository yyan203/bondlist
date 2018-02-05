#! /bin/bash

sed "1,2d" zif4.xyz > zif4_2.xyz

INPUTFILE="zif4_2.xyz" # delete first two lines please
OUT="out"
OUTFILE="out2"
# awk 'BEGIN{
# Zn=C=N=H=1;
# id=1;
# }{
#  name=$1
#  if(name=="Zn"){print id, $1Zn;Zn++;id++;}
#  else if(name=="C"){}
#  else if(name=="N"){}
#  else {}
# }' $INPUTFILE > idtoatom

awk 'BEGIN {
Zn=C=N=H=1;
id=1;
while (getline < "'"$INPUTFILE"'")
    {
        split($0, ft, " " );
        name=ft[1];
        if(name=="Zn"){idtoname[id]=ft[1]Zn;Zn++;id++;}
        else if(name=="C"){idtoname[id]=ft[1]C;C++;id++;}
        else if(name=="N"){idtoname[id]=ft[1]N;N++;id++;}
        else {idtoname[id]=ft[1]H;H++;id++;}
    }
close ("'"$INPUTFILE"'")
 ##print idtoname[270];

b1=b2=b3=b4=1;
 while (getline < "'"$OUT"'")
     {
         if($1==1)    {print "$bond:Zn_N"b1,"$atom:"idtoname[$2],"$atom:"idtoname[$3]> "'"$OUTFILE"'"; b1++;}
         else if($1==2){print "$bond:H_C"b2,"$atom:"idtoname[$2],"$atom:"idtoname[$3]> "'"$OUTFILE"'"; b2++;}
         else if($1==3){print "$bond:C_N"b3,"$atom:"idtoname[$2],"$atom:"idtoname[$3]> "'"$OUTFILE"'"; b3++;}
         else          {print "$bond:C_C"b4,"$atom:"idtoname[$2],"$atom:"idtoname[$3]> "'"$OUTFILE"'"; b4++;} }
}'
