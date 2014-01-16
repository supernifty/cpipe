#!/bin/bash
echo "Renaming to flagship names ..."
for i in *.bed; do mv $i `echo $i | sed 's/^.*_//'`; done

echo "Renaming chromosomes ..."
for i in *.bed; do sed -i.bak 's/^/chr/' $i; done

echo "Renaming chrMT to chrM ..."
for i in *.bed; do sed -i.bak 's/^chrMT/chrM/' $i; done

rm *.bak

echo "Flattening ..."
mkdir -p flattened
for i in *.bed; do
    sortBed -i $i | bedtools merge -nms -i - | sed 's/;.*$//' > flattened/$i
done
rm -rf unflattened
mkdir -p unflattened
mv *.bed unflattened
mv flattened/*.bed .
rmdir flattened

echo "Sorting ..."
mkdir sorted
for i in *.bed; do 
    ../../tools/IGVTools/2.3.6/igvtools.lowmem sort $i sorted/$i
done

mkdir -p unsorted
mv *.bed unsorted
mv sorted/*.bed .
rmdir sorted

