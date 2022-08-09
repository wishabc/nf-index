#!/usr/bin/env python3
"""
Counts sequencing reads within a set of regions
"""

import sys
import pysam


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

with pysam.AlignmentFile(sys.argv[1], 'rc') as bam:
	for index, (chrom, start, end) in enumerate(bed3_iterator(sys.stdin)):
		try:
			print(bam.count(chrom, start, end, read_callback='all'))
		except ValueError:
			if index == 0:
				pysam.index(sys.argv[1])
				print(bam.count(chrom, start, end, read_callback='all'))
			else:
				raise
