#!/bin/bash
echo "Renaming to flagship names ..."
for i in *.bed; do mv $i `echo $i | sed 's/^.*_//'`; done

echo "Renaming chromosomes ..."
for i in *.bed; do sed -i.bak 's/^/chr/' $i; done

echo "Renaming chrMT to chrM ..."
for i in *.bed; do sed -i.bak 's/^chrMT/chrM/' $i; done

echo "Sorting ..."

mkdir sorted
for i in *.bed; do 
    ../../tools/IGVTools/2.3.6/igvtools.lowmem sort $i sorted/$i
done

mkdir -p unsorted
mv *.bed unsorted
mv sorted/*.bed .
rmdir sorted

