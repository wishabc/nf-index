#!/bin/bash

#Path to relevant Masterlist Files
masterlist=$1
encode3=$2
gencode=$3
gwas_catalog=$4
repeats=$5

#Bedops command to print whether or not the masterlist file overlaps an encode3 DHS and by how much
bedmap --echo --echo-map --bp-ovr 1 --indicator --bases-uniq-f ${masterlist} ${encode3} \
| awk -F'|' '{print $(NF-1)"\t"$NF}'> is_encode3.txt


#Double check number of rows are the same
#Can do later

echo "Finished encode3 annotation"

###############
#Parse Gencode#
###############

zcat ${gencode} \
| awk -F'\t' '{if($4 != $5) print}' \
| awk -F'\t' '{
        if($3 == "transcript") {
                if($7 == "+") {
                        print $1"\t"$4"\t"$4+1"\t"$9;
                }
                else if($7 == "-") {
                         print $1"\t"$5-1"\t"$5"\t"$9;
                }
        }
}' - \
| grep -v chrM | grep -v Selenocysteine | grep -v codon \
| sort-bed - \
| awk -F';' '{print $1"\t"$3}' \
| awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$5}' \
| sed 's/gene_name//g' \
| sed 's/\"//g' \
| sed 's/ //g' \
> tss.bed

###########################
#Closest-Features to Genes#
##########################

closest-features --closest --no-ref --dist ${masterlist} tss.bed \
| awk -F'\t' '{print $4}' \
| awk -F'|' '{print $2"\t"$1}' \
> dist_gene.txt


echo "Finished Distance to TSS"


####################
#Parse gwas_catalog#
####################
touch tmp2.catalog_parsed.tsv
rm tmp2.catalog_parsed.tsv
touch tmp2.catalog_parsed.tsv

#Get columns I want from the gwas catalog
#Got rid of some weird chr5 X 12, chr1 x 3 in the first column
cut -f7,12,13,21,28,31 ${gwas_catalog} \
| awk -F'\t' '{print "chr"$2"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6"\t"$1}' \
| awk -F' ' '{if($2 != "x") print}' \
> tmp.catalog_parsed.tsv

#Only extract chr1-22,X,Y
for i in {1..22} X Y
do
        awk -v chrom=$i '{if($1 =="chr"chrom) print}' tmp.catalog_parsed.tsv >> tmp2.catalog_parsed.tsv
done

cat tmp2.catalog_parsed.tsv \
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


awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$10"\t"$8}' overlap-answer.txt \
| sort-bed - \
> dhs_annotated_all-repeats.bed


###################################
#Map Annotation BACK to Masterlist#
###################################
bedmap --echo-map --fraction-both 1 tmp.masterlist.bed3 dhs_annotated_all-repeats.bed \
| cut -f4-6 \
| awk -F'\t' '{if($1 == "") print """\t""""\t"""; else print}' \
> repeats.txt

echo "Finished Repeat Annotations"

