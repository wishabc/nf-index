#!/bin/bash

#Path to relevant Masterlist Files
masterlist=$1
repeats=$2
outfile=$3

#############################
#Parse and Map Repeat Masker#
#############################

cut -f1-3 ${masterlist} \
> tmp.masterlist.bed3


zcat ${repeats} \
    | bedmap --echo --echo-map --fraction-map .5 --echo-overlap-size --echo-map-size --ec --skip-unmapped tmp.masterlist.bed3 - \
    > repeated_mapped.bed

########################
#Choose Best Annotation#
########################
biggest=0
col=0
fraction=0

cat repeated_mapped.bed \
    | awk -F'|' -v OFS='\t' -v f=$fraction -v b=$biggest -v c=$col '{
        line=$3
        split(line,a,";")
        mapped=$2
        split(mapped,m,";")
        size=$4
        split(size,s,";")
        
        if (length(a) == 1) {
            c=1;
        }
        else {
            for(i=1;i<=NF;i++) {
                if (a[i] > b) {
                        b=a[i];
                        c=i;
                        f=a[i]/s[i];
                }
                else if (a[i] == b) {
                    if(a[i]/s[i] > f) {
                        b=a[i];
                        c=i;
                    }
                } 
            }      
        }
        print $1, m[c];
        b=0;      
    }' \
    |  awk -v OFS='\t' \
        '{print $1,$2,$3,$9,$10,$8}' \
    | sort-bed - > dhs_annotated_all-repeats.bed


###################################
#Map Annotation BACK to Masterlist#
###################################
bedmap --echo-map --fraction-both 1 tmp.masterlist.bed3 dhs_annotated_all-repeats.bed \
    | cut -f4-6 \
    | awk -F'\t' -v OFS='\t' '{if($1 == "") print "", "", ""; else print}' \
    > repeats.txt

echo -e "repeat_class\trepeat_family\trepeat_name" > ${outfile}
paste repeats.txt >> ${outfile}

echo "Finished Repeat Annotations"
