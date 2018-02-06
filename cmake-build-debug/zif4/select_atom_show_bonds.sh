#! /bin/bash
if [  ! $# -eq 2 ]
then
    echo "example: bondfile  atomindex"
else

awk -v idn=$2 ' {
if($2==idn) print $2,$3;
if($3==idn) print $3,$2;
}' $1
fi
