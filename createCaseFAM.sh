#!/bin/bash

input=$1
output=$2

#printf ${input}


cat ${input} | awk '{print $1}' |sort -u | awk '{printf "%s %s 0 0 0 2\n",$1,$1}' > ${output}
sed -i '1d' ${output}
