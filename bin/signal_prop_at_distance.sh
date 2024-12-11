#!/bin/bash

SAMPLE_ID="$1"
DISTANCE="$2"

PEAKS_FILE="$3"
CUTCOUNTS_FILE="$4"

CHROM_SIZES="$5"
DFOUT="$6"


if ! [[ "$DISTANCE" =~ ^[0-9]+$ ]]; then
    echo "Error: Distance must be a non-negative integer."
    exit 1
fi

echo -e "length\tsignal\tdistance\tsample_id" > $DFOUT

zgrep -v '^#' "$PEAKS_FILE" \
    | awk -v OFS="\t" '{print $1, $7, $7 + 1}' \
    | bedtools slop -i stdin -g $CHROM_SIZES -b $DISTANCE \
    | bedops --merge - \
    | bedmap --delim '\t' --echo --sum - \
        <(zgrep -v '^#' "$CUTCOUNTS_FILE" | awk -v OFS="\t" '{print $1, $2, $3, 0, $4}') \
    | awk -v OFS='\t' \
        -v total_cutcounts=$SUM_CUTCOUNTS \
        -v sample_id=$SAMPLE_ID \
        -v distance=$DISTANCE \
        '{length_sum += $3 - $2; signal_sum += $4}
         END {print length_sum, signal_sum, distance, sample_id}' >> $DFOUT
