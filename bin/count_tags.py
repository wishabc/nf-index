#!/usr/bin/env python3
"""
Counts sequencing reads within a set of regions
"""

import sys
import pysam

bam = pysam.AlignmentFile(sys.argv[1], 'r')
def bed3_iterator(filehandle):
	"""
	Generator that parses BED3 format from a string iterator

	Returns:
		genomic_interval
	"""
	for line in filehandle:
		
		fields = line.strip().split()
		
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])

		yield (chrom, start, end)


for chrom, start, end in bed3_iterator(sys.stdin):
	print(bam.count(chrom, start, end, read_callback='all'))

