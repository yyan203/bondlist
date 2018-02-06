#! /bin/bash

#sed "1,2d" trueglasszif4.xyz > trueglasszif4_2.xyz

INPUTFILE="trueglasszif4_3.xyz" # with index and without title
INPUTFILE2="trueglassout_1.2" # with index
#OUT="trueglassout"
#OUTFILE="trueglassout2"
nmols=$(wc -l < $INPUTFILE)
echo $nmols "atoms"

awk -v l=$nmols 'BEGIN {
for (i=1;i<=l;i++){
bondnum[i]=0;
}

while (getline < "'"$INPUTFILE"'")
    {
        split($0, ft, " " );
        name=ft[2];
        if(name=="Zn")    {bondnum[ft[1]]=4;atomtype[ft[1]]="Zn";}
        else if(name=="C"){bondnum[ft[1]]=3;atomtype[ft[1]]="C";}
        else if(name=="N"){bondnum[ft[1]]=3;atomtype[ft[1]]="N";}
        else              {bondnum[ft[1]]=1;atomtype[ft[1]]="H";}
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
