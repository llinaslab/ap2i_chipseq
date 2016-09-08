#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description: Creating a bed file with regions before and after the first and last genes
of every chromosome

Version = 1.0
Author  = "Philipp Ross"
"""

# import modules
from __future__ import division
import os
import sys
import argparse
import subprocess

def zero_or_greater(num):
	if num < 0:
		num = 1
	return num

def main():

	parser = argparse.ArgumentParser(description = "Create a bed file with regions before and after the first and last genes in every chromosome")
	parser.add_argument("-f", "--gff", dest = "gff", help = "GFF3 genes feature file")
	parser.add_argument("-g", "--genome", dest = "genome", help = "Chromosome lengths file")
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	lengths = dict()
	with open(args.genome) as genome:
		for line in genome:
			line =line.rstrip().split("\t")
			lengths[line[0]] = line[1]

	pos_edges = dict()
	neg_edges = dict()
	with open(args.gff) as gff:
		for line in gff:
			if line.startswith("#"):
				continue
			columns = line.rstrip().split("\t")
			chrom = columns[0]
			start = int(columns[3])
			stop = int(columns[4])
			strand = columns[6]
			if strand == "+":
				if chrom in pos_edges.keys():
					if pos_edges[chrom]["first"] > start:
						pos_edges[chrom]["first"] = zero_or_greater(start - 2000)
					elif pos_edges[chrom]["last"] < stop:
						pos_edges[chrom]["last"] = stop
				else:
					pos_edges[chrom] = {"first" : start, "last" : stop}
			if strand == "-":
				if chrom in neg_edges.keys():
					if neg_edges[chrom]["first"] > start:
						neg_edges[chrom]["first"] = start
					elif neg_edges[chrom]["last"] < stop:
						neg_edges[chrom]["last"] = stop + 2000
				else:
					neg_edges[chrom] = {"first" : start, "last" : stop}
			
	# I'm subtracting and adding 2000 basepairs to each strand on one particular edge so as not to exclude possible regulatory
	# sequences
	for key in pos_edges.keys():
		if key != "PfDd2_MT" and key != "PfDd2_API":
			print >> sys.stdout, key + "\t" + str(0) + "\t" + str(pos_edges[key]["first"]) + "\t" + "5_prime_edge" + "\t" + "." + "\t" + "+" + "\n" + key + "\t" + str(pos_edges[key]["last"]) + "\t" + lengths[key] + "\t" + "3_prime_edge" + "\t" + "." + "\t" + "+"
	for key in neg_edges.keys():
		if key != "PfDd2_MT" and key != "PfDd2_API":
			print >> sys.stdout, key + "\t" + str(0) + "\t" + str(neg_edges[key]["first"]) + "\t" + "5_prime_edge" + "\t" + "." + "\t" + "-" + "\n" + key + "\t" + str(neg_edges[key]["last"]) + "\t" + lengths[key] + "\t" + "3_prime_edge" + "\t" + "." + "\t" + "-"


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt! Ciao!\n")
		sys.exit(0)
