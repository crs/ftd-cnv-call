#!/bin/bash

for i in `ls -1 processing_allsnps`;
do
	echo -n "Intersection $i "
	awk 'FNR==NR{f[$0];next}($1 in f)' snplist/SNPs_Common.txt processing_allsnps/$i > processing_intersect/$i
	echo "done."
done
