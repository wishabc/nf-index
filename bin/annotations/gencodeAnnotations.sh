#!/bin/bash

#Path to relevant Masterlist Files
masterlist=$1
gencode=$2
chromSize=$3

##################################
#Parse Gencode File (utr updated)#
##################################
chromInfo="chrom_sizes.bed"
awk -v OFS='\t' \
    '{ print \$1,0,\$2 }' \
    ${chromSize} \
    | sort-bed - > ${chromInfo}

#Remove row if start = end
zcat ${gencode} \
    | awk -F'\t' -v OFS='\t' '{
        if($4 != $5) {
            print $1,$4,$5,$3,$7
        }
    }' \
    | sort-bed -  \
    > gencode.filtered.gtf

#Expand the transcription region to say promoter. +/- 1KB of TSS
awk -F'\t' -v OFS='\t' '{
        if($4 == "transcript") {
            if ($5 == "+") {
                print $1,$2,$2+1000,"promoter";
            } else {
                print $1,$3-1000,$3,"promoter";
            }
        } else  {
            print $1,$2,$3,$4;
        }
    }' gencode.filtered.gtf \
    | grep -v chrM \
    | grep -v Selenocysteine \
    | grep -v codon \
    | sort-bed - \
> gencode.bed

#Need to find the INTRONS. Difference between gene and (CDS + PROMOTER + UTR) 
awk '{if($4 == "gene") print}' gencode.bed > gene.bed
awk '{if($4 == "exon") print}' gencode.bed > exon.bed
awk '{if($4 == "CDS") print}' gencode.bed > cds.bed
awk '{if($4 == "promoter") print}' gencode.bed > promoter.bed
awk '{if($4 == "five_prime_utr" || $4 == "three_prime_utr") print}' gencode.bed > utr.bed

bedops --ec -m utr.bed exon.bed promoter.bed cds.bed \
    | bedops --ec -d gene.bed - \
    | awk -v OFS='\t' '{print $1,$2,$3,"intron"}' tmp.intron.bed > intron.bed

#Need to find the Intergenic region. Difference between Genome and gene-body + promoter region
bedops --ec -d \
    ${chromInfo} \
    gene.bed \
    promoter.bed \
    awk -v OFS='\t' \
        '{print $1,$2,$3,"intergenic"}' > intergenic.bed

#Clean Up
rm tmp.intron.bed
rm tmp.intergenic.bed

#Unite promoter, exon, intron, and intergenic regions in one bed file
#Map united bed file and map to DHS_Index.bed
bedops --ec -u \
    promoter.bed \
    exon.bed \
    intron.bed \
    intergenic.bed \
    | sort-bed - \
    | bedmap --ec --echo \
        --echo-map \
        --skip-unmapped \
        --echo-overlap-size \
        <(cut -f1-3 ${masterlist}) - > gencode_mapped.bed

#Filter Initial Gencode file based on protein coding and non-protein coding regions
zcat ${gencode} \
    | grep protein_coding \
    | awk -v OFS='\t' '{
        if($3 == "transcript") {
            if ($7 == "+") {
                print $1,$4,$4+1000,"PC";
            } else {
                print $1,$5-1000,$5,"PC";
            }
        }
        else {
            print $1,$4,$5,"PC"
        }
    }' \
    | awk -F'\t' '{if($2 != $3) print}' - \
    | sort-bed - \
    > PC.bed

zcat ${gencode} \
    | grep -v protein_coding \
    | awk '{
            if($3 == "transcript") {
                    if ($7 == "+") {
                        print $1"\t"$4"\t"$4+1000"\t""NPC";
                    }
                    else {
                        print $1"\t"$5-1000"\t"$5"\t""NPC";
            }
            }
            else {

                    print $1"\t"$4"\t"$5"\t""NPC"
            }
    }' \
    | awk -F'\t' '{if($2 != $3) print}' - \
    | sort-bed - \
    > NPC.bed

bedops -u PC.bed NPC.bed > PC-NPC-gencode.bed


####################################
#Choose the best Gencode Annotation#
####################################

biggest=0
col=0

#Promoter > Exon > Intron > Intergenic
sed 's/intergenic/1/g' gencode_mapped.bed \
| sed 's/intron/2/g' \
| sed 's/exon/3/g' \
| sed 's/promoter/4/g' \
> choose_best_annotation.bed

awk -F'|' -v b=$biggest -v c=$col '{
        line=$3
        split(line,a,";")

        mapped=$2
        split(mapped,m,";")
        
    if (length(a) == 1) {
        print $1"\t"$2
    }
    else {
        for(i=1;i<=NF;i++) {
            if (a[i] > b) {
                b=a[i];
                c=i;
            }
            else if (a[i] == b) {
                old=m[c];
                split(old,o,"\t");
                new=m[i];
                split(new,n,"\t");
                if (o[4] < n[4]) {
                    b=a[i];
                    c=i;
                }
            }
        }
    print $1"\t"m[c];
    b=0;
    }

}' choose_best_annotation.bed > best_annotation.bed

awk -v OFS='\t' \
    '{print $1,$2,$3,$NF}' best_annotation.bed \
    | sort-bed - \
    | awk -v OFS='\t' '{
        if($4 == 1) {
            print $1,$2,$3,"intergenic";
        } else if($4 == 2) {
            print $1,$2,$3,"intron";
        } else if($4 == 3) {
            print $1,$2,$3,"exon";
        } else if($4 == 4) {
            print $1,$2,$3,"promoter"
        }' - \
    > dhs_annotated.bed

#####################################
#Update Exon/Promoter/NPC/PC Regions#
#####################################

#Map Exon regions to DHS Annotations
awk '{if($4 == "exon") print}' dhs_annotated.bed > dhs-exon.bed
bedops -u utr.bed cds.bed > utr-cds-gencode.bed
bedmap --echo --echo-map --echo-overlap-size --echo-map-size --ec dhs-exon.bed  utr-cds-gencode.bed \
> exon_mapped.bed

#Choose the element with the largest overlap or the largest fraction of overlap
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
        else if(length(a) > 1) {
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
        else {
            c=1;
        }
        print $1"\t"m[c];
        b=0;      
}' exon_mapped.bed > best_exon_mapped.bed


awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' best_exon_mapped.bed \
| sort-bed - \
| awk -F'\t' '{if($5 == "") print $1"\t"$2"\t"$3"\t"$4"\t""NPC"; else print}' \
> dhs_annotated_exon.bed

#Map Promoter
awk '{if($4 == "promoter") print}' dhs_annotated.bed > dhs-promoter.bed

bedmap --echo --echo-map --echo-overlap-size --echo-map-size --skip-unmapped --ec dhs-promoter.bed PC-NPC-gencode.bed \
> promoter_mapped.bed


#Pick the element with the largest overlap or the largest fraction of overlap

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
}' promoter_mapped.bed > best_promoter_mapped.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' best_promoter_mapped.bed \
| sort-bed - \
> dhs_annotated_promoter.bed

#Map PC/NPC regions to Intronic DHSs
awk -F'\t' '{if($4 == "intron") print}' dhs_annotated.bed > dhs-intron.bed

bedmap --echo --echo-map --echo-overlap-size --echo-map-size --skip-unmapped --ec dhs-intron.bed PC-NPC-gencode.bed \
> intron_mapped.bed


biggest=0
col=0
fraction=0

#Choose Best Annotation
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
}' intron_mapped.bed > best_intron_mapped.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' best_intron_mapped.bed \
| sort-bed - \
> dhs_annotated_intron.bed


#Paste the Annotations together

#gene_body
cut -f4 dhs_annotated.bed > gene_body.txt

#exon_subfamily
bedmap --echo-map --fraction-both 1 ${masterlist} dhs_annotated_exon.bed \
    | awk -F'\t' '{if($5 == "NPC") print $1"\t"$2"\t"$3"\t"$4"\t"""; else print}' \
    | cut -f5 \
    | sed 's/;.*//' \
    > exon_subfamily.txt

#PC/NP
cat dhs_annotated_exon.bed dhs_annotated_intron.bed dhs_annotated_promoter.bed \
    | sort-bed - \
    | bedmap --echo-map --fraction-both 1 ${masterlist} - \
    | awk -F'\t' '{if($5 == "CDS" || $5 == "three_prime_utr" || $5 == "five_prime_utr") print $1"\t"$2"\t"$3"\t"$4"\t""PC"; else print}' \
    | awk -F'\t' '{if($5 == "") print $1"\t"$2"\t"$3"\t"$4"\t""NPC"; else print}' \
    | cut -f5 \
    > is_coding.txt


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
    | awk -F'|' -v OFS='\t' '{print $2,$1}' \
    > dist_gene.txt


echo "Finished Distance to TSS"

paste dist_gene.txt gene_body.txt exon_subfamily.txt is_coding.txt > gencode_annotations.txt
echo -e "dist_tss\tgene\tgene_body\texon_subgroup\tis_coding" > gencodeAnnotations_header.txt

echo "Finished Gencode Annotation"
