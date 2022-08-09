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
	for fields in filehandle:
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])

		yield (chrom, start, end)


with pysam.AlignmentFile(sys.argv[1], 'rc') as bam,\
 pysam.TabixFile(sys.argv[2]) as bed:
	for index, (chrom, start, end) in enumerate(bed3_iterator(bed)):
		try:
			print(bam.count(chrom, start, end, read_callback='all'))
		except Exception as e:
			print(f'Problems with {chrom}:{start}-{end}')
			raise e

