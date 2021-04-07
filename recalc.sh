#!/bin/bash

basePath=/Volumes/MyPass2019/MD_simulation/Polaron/200130/CollectionsXYZ
outPath=SASArecalc
for i in XSmetaD1a_1.txt XSmetaD2_1.txt XSmetaD2_2.txt XSmetaD2_3.txt XSmetaD3_1.txt XSmetaD3_2.txt XSmetaD3_3.txt
do
    echo $i
    ./a.out $basePath/$i 10000 1.8 $outPath/$i.SASArecalc #> SASArecalc/$i.SASArecalc
done

echo "00_eq.txt" 
./a.out $basePath/00_eq.txt 2001 1.8 $outPath/00_eq.txt.SASArecalc #> SASArecalc/00_eq.txt.SASArecalc
