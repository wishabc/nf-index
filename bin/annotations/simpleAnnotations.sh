#!/bin/bash

#Path to relevant Masterlist Files
masterlist=$1
encode3=$2
gwas_catalog=$3
repeats=$4

#Bedops command to print whether or not the masterlist file overlaps an encode3 DHS and by how much
bedmap --echo --echo-map --bp-ovr 1 --indicator --bases-uniq-f ${masterlist} ${encode3} \
| awk -F'|' '{print $(NF-1)"\t"$NF}'> is_encode3.txt


#Double check number of rows are the same
#Can do later

echo "Finished encode3 annotation"

####################
#Parse gwas_catalog#
####################
#Got rid of some weird chr5 X 12, chr1 x 3 in the first column
cut -f7,12,13,21,28,31 ${gwas_catalog} \
    | awk -F'\t' '{print "chr"$2"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6"\t"$1}' \
    | awk -F' ' '{if($2 != "x") print}' \
    | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/' \
    | sort-bed - \
> catalog_parsed.tsv

################################
#Map Masterlist to gwas_catalog#
################################

bedmap --count ${masterlist} catalog_parsed.tsv > gwas_catalog_count.txt
bedmap --echo --echo-map ${masterlist} catalog_parsed.tsv > gwas_catalog_mapped.bed


#Double check number of rows are the same
#Can do later
#Also maybe have a mapped file to show the overlapping trait(s)

echo "Finished gwas_catalog annotations"

#############################
#Parse and Map Repeat Masker#
#############################

cut -f1-3 ${masterlist} \
> tmp.masterlist.bed3


zcat ${repeats} \
| tail -n +2 \
| cut -f6-8,10-13 \
| sort-bed - \
| grep -v LTR? | grep -v DNA? | grep -v RC? | grep -v SINE? \
| bedmap --echo --echo-map --fraction-map .5 --echo-overlap-size --echo-map-size --ec --skip-unmapped tmp.masterlist.bed3 - \
> repeated_mapped.bed

########################
#Choose Best Annotation#
########################
biggest=0
col=0
fraction=0

awk -F'|' -v f=$fraction -v b=$biggest -v c=$col '{
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
        print $1"\t"m[c];
        b=0;      
}'  repeated_mapped.bed  > overlap-answer.txt


awk -v OFS='\t' \
    '{print $1,$2,$3,$9,$10,$8}'\
     overlap-answer.txt \
    | sort-bed - > dhs_annotated_all-repeats.bed


###################################
#Map Annotation BACK to Masterlist#
###################################
bedmap --echo-map --fraction-both 1 tmp.masterlist.bed3 dhs_annotated_all-repeats.bed \
| cut -f4-6 \
| awk -F'\t' '{if($1 == "") print """\t""""\t"""; else print}' \
> repeats.txt

paste is_encode3.txt gwas_catalog_count.txt repeats.txt > simpleAnnotations.txt
echo -e "is_encode3\tencode3_ovr-fraction\tnum_gwasCatalog_variants\trepeat_class\trepeat_family\trepeat_name" > simpleAnnotations_header.txt

echo "Finished Repeat Annotations"

