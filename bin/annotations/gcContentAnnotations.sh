masterlist=$1
genome_fasta=$2
mappable_file=$3
outfile=$4

echo -e 'n_gc\tpercent_gc\tn_mappable' > ${outfile}

faidx -i nucleotide -b ${masterlist} ${genome_fasta} \
    | awk -v OFS="\t" \
        'NR>1 { 
            total=$4+$5+$6+$7+$8;
            cg=$6+$7;
            print $1, $2-1, $3, cg, cg/total; }' \
    | bedmap \
        --delim "\t" --echo \
        --bases-uniq - ${mappable_file} \
    | cut -f4- \
    >> ${outfile}