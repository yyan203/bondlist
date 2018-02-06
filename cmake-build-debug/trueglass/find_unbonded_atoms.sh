#!/bin/bash

#sed "1,2d" trueglasszif4.xyz > trueglasszif4_2.xyz

INPUTFILE="trueglasszif4_2.xyz" # with index and without title
INPUTFILE2="glassbondlist" # with index
#OUT="trueglassout"
#OUTFILE="trueglassout2"
nmols=$(wc -l < $INPUTFILE)
echo $nmols "atoms"

awk -v l=$nmols 'BEGIN {
for (i=1;i<=l;i++){
bondnum[i]=0;
}
line=0;
while (getline < "'"$INPUTFILE"'")
    {
        line++;
        split($0, ft, " " );
        name=ft[1];
        if(name=="Zn")    {bondnum[line]=4;atomtype[line]="Zn";}
        else if(name=="C"){bondnum[line]=3;atomtype[line]="C";}
        else if(name=="N"){bondnum[line]=3;atomtype[line]="N";}
        else              {bondnum[line]=1;atomtype[line]="H";}
    }
close ("'"$INPUTFILE"'");

while (getline < "'"$INPUTFILE2"'")
    {
        split($0, ft, "  ");
        atom1=ft[2];atom2=ft[3];
        bondnum[atom1]--;bondnum[atom2]--;
    }
print "D";
close ("'"$INPUTFILE2"'");

for(i=1;i<=l;i++){
    if(bondnum[i]<0) print i,atomtype[i],">bonds",-bondnum[i];
    if(bondnum[i]>0) print i,atomtype[i],"<bonds",bondnum[i];
}
}'
