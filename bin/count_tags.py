#!/usr/bin/env python3
"""
Counts sequencing reads within a set of regions
"""

import sys
import pysam

from genome_tools import bed

bam = pysam.AlignmentFile(sys.argv[1], 'r')

for interval in bed.bed3_iterator(sys.stdin):
	print(bam.count(interval.chrom, interval.start, interval.end, read_callback='all'))

