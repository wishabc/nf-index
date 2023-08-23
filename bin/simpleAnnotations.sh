#!/bin/bash

#Path to relevant Masterlist Files
masterlist=$1
encode3=$2
gencode=$3
gwas_catalog=$4

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

