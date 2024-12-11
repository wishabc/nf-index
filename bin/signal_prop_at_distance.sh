#!/bin/bash
PEAKS_FILE="$1"
CUTCOUNTS_FILE="$2"
DISTANCE="$3"
SUM_CUTCOUNTS="$4"
SAMPLE_ID="$5"
CHROM_SIZES="$6"
DFOUT="$7"


if ! [[ "$DISTANCE" =~ ^[0-9]+$ ]]; then
    echo "Error: Distance must be a non-negative integer."
    exit 1
fi

echo -e "#length\tsignal\tdistance\tsum_cutcounts\tsample_id" > $DFOUT

zgrep -v '^#' "$PEAKS_FILE" \
    | awk -v OFS="\t" '{print $1, $7, $7 + 1}' \
    | bedtools slop -i stdin -g $CHROM_SIZES -b $DISTANCE \
    | bedops --merge - \
    | bedmap --echo --sum - \
        <(zgrep -v '^#' "$CUTCOUNTS_FILE" | awk -v OFS="\t" '{print $1, $2, $3, 0, $4}') \
    | awk -v OFS='\t' \
        -v total_cutcounts=$SUM_CUTCOUNTS \
        -v sample_id=$SAMPLE_ID \
        -v distance=$DISTANCE \
        '{length_sum += $3 - $2; signal_sum += $4}
         END {print length_sum, signal_sum, distance, total_cutcounts, sample_id}'
