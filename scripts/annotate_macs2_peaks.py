#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description: Annotate a MACS2 narrowPeaks file

Version = 0.2
Author  = Philipp Ross
"""

# import modules
import os
import re
import sys
import math
import argparse
import random
import subprocess
import pandas as pd
import pybedtools as pybed

# define main function
def main():

	parser = argparse.ArgumentParser(description = "Script to annotate ChIP-seq peaks")
	parser.add_argument("-f", "--features", dest = "features", required = True, help = "Features GFF3 file")
	parser.add_argument("-g", "--genome", dest = "genome", required = True, help = "BED genome file or chromosome lengths")
	parser.add_argument("-p", "--peaks", dest = "peaks", required = True, help = "MACS2 formatted narrowPeak file")
	parser.add_argument("-m", "--pbm-motifs", dest = "pbm_motifs", help = "GFF3 or BED formatted motif locations file")
	parser.add_argument("-d", "--dreme-motifs", dest = "dreme_motifs", help = "GFF3 or BED formatted motif locations file")
	parser.add_argument("-a", "--random-motifs", dest = "random_motifs", help = "GFF3 or BED formatted motif locations file")
	parser.add_argument("-i", "--random-intervals", dest = "random_intervals", help = "BED formatted interval file")
	parser.add_argument("-r", "--motif-region", dest = "motif_region", default = 100, type = int, help = "Region left and right of peak summit to look for motif matches")
	parser.add_argument("-o", "--output-dir", dest = "out_dir", default = "macs2", help = "Output directory")
	parser.add_argument("-n", "--name", dest = "name", default = "name", help = "File name prefix")
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	# create output directory
	if not os.path.exists(args.out_dir):
		os.makedirs(args.out_dir)

	# create bedtool objects
	features = pybed.BedTool(args.features)
	peak_summits = pybed.BedTool(args.peaks)

	# create dictionary of every peak name
	peak_dict = dict()
	for peak in peak_summits:
		peak_dict[peak.name] = {\
				"seqid"     : peak.chrom,\
				"source"    : "annotation",\
				"type"      : "macs2_peak",\
				"start"     : peak.start + 1,\
				"stop"      : peak.stop,\
				"score"     : peak.score,\
				"strand"    : ".",\
				"phase"     : ".",\
				"attrs"     : {\
						"name"         : peak.name,\
						"length"       : peak.length,\
						"ugene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : "", "dist" : "-1"},\
						"dgene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : "", "dist" : "-1"},\
						"ogene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : ""},\
						"pbm_motifs"   : [],\
						"dreme_motifs" : [],\
						"pvalue"       : peak[7],\
						"qvalue"       : peak[8],\
						"source"       : int(peak[9])\
						}
				}

	# find overlapping and closest genes and create new interval files
	peak_summits = peak_summits.sort()
	feature_overlap = peak_summits.intersect(features, wb=True).saveas("tmpa.bed")
	upstream_features = peak_summits.closest(features, D="ref", id=True, io=True, g=args.genome).saveas("tmpb.bed")
	downstream_features = peak_summits.closest(features, D="ref", iu=True, io=True, g=args.genome).saveas("tmpc.bed")

	# annotate peaks that overlap ORFS
	for interval in feature_overlap:
		if interval.name in peak_dict.keys():
			if interval[18] != ".":
				attributes = interval[18].split(";")
				peak_dict[interval.name]["attrs"]["ogene"]["strand"] = interval[16]
				if attributes[0].split("=")[0] == "ID":
					peak_dict[interval.name]["attrs"]["ogene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						peak_dict[interval.name]["attrs"]["ogene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							peak_dict[interval.name]["attrs"]["ogene"]["name"] = attributes[2].split("=")[1]
	
	# annotate peaks that are upstream of negative stranded genes
	for interval in upstream_features:
		if interval.name in peak_dict.keys():
			if len(interval.fields) != 20:
				print >> sys.stderr, "Upstream features"
				print >> sys.stderr, interval
			if interval[18] != ".":
				attributes = interval[18].split(";")
				peak_dict[interval.name]["attrs"]["ugene"]["strand"] = interval[16]
				peak_dict[interval.name]["attrs"]["ugene"]["dist"] = int(interval[19]) * -1
				if attributes[0].split("=")[0] == "ID":
					peak_dict[interval.name]["attrs"]["ugene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						peak_dict[interval.name]["attrs"]["ugene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							peak_dict[interval.name]["attrs"]["ugene"]["name"] = attributes[2].split("=")[1]

	# annotate peaks that are upstream of positive stranded genes
	for interval in downstream_features:
		if interval.name in peak_dict.keys():
			if len(interval.fields) != 20:
				print >> sys.stderr, "Downstream features"
				print >> sys.stderr, interval
			if interval[18] != ".":
				attributes = interval[18].split(";")
				peak_dict[interval.name]["attrs"]["dgene"]["strand"] = interval[16]
				peak_dict[interval.name]["attrs"]["dgene"]["dist"] = int(interval[19]) * -1
				if attributes[0].split("=")[0] == "ID":
					peak_dict[interval.name]["attrs"]["dgene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						peak_dict[interval.name]["attrs"]["dgene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							peak_dict[interval.name]["attrs"]["dgene"]["name"] = attributes[2].split("=")[1]

	if(args.random_intervals):
		random_replicates = 100
		# create random summits for random intervals file
		random_intervals = dict()
		with open(args.random_intervals) as ri:
			for line in ri:
				line = line.split("\t")
				start = int(line[1])
				stop = int(line[2])
				name = line[3]
				length = stop - start
				random_sources = []
				for i in range(random_replicates):
					random_sources.append(random.randrange(math.floor(length * (1.0/3.0)), math.floor(length * (2.0/3.0)))) 
				if name not in random_intervals.keys():
					random_intervals[name] = random_sources
				else:
					print >> sys.stderr, "This shouldn't be happening...Check random interval motif searching"
					sys.exit(1)

		# check if motif hits for random intervals fall within motif region
		random_occurences = dict()
		with open(args.random_motifs) as rm:
			next(rm)
			for line in rm:
				line = line.split("\t")
				motif_name = line[8].split(";")[0].split("=")[1]
				random_occurences[motif_name] = [0] * 100

		for i in range(random_replicates):
			with open(args.random_motifs) as rm:
				next(rm)
				for line in rm:
					line = line.split("\t")
					interval_name = line[0]
					motif_start = int(line[3])
					motif_name = line[8].split(";")[0].split("=")[1]
					if interval_name in random_intervals.keys():
						left_bound = random_intervals[interval_name][i] - args.motif_region
						if left_bound < 1:
							left_bound = 1
						right_bound = random_intervals[interval_name][i] + args.motif_region
						if motif_start > left_bound and motif_start < right_bound:
							if motif_name in random_occurences.keys():
								random_occurences[motif_name][i] += 1
	
		# write random occurences data frame to file
		pd.DataFrame(random_occurences).to_csv(\
				path_or_buf = args.out_dir + "/random.occurences", \
				sep="\t", \
				header=True, \
				index_label="replicates")

	# annotate pbm motifs that fall within motif region
	if args.pbm_motifs:
		pbm_motif_occ = dict()
		with open(args.pbm_motifs) as pmf:
			next(pmf)
			for line in pmf:
				line = line.rstrip().split("\t")
				peak_name = line[0]
				motif_start = int(line[3])
				motif_name = line[8].split(";")[0].split("=")[1]
				if peak_name in peak_dict.keys():
					peak_source = peak_dict[peak_name]["attrs"]["source"]
					peak_length = peak_dict[peak_name]["attrs"]["length"]
					left_bound = peak_source - args.motif_region
					if left_bound < 1:
						left_bound = 1
					right_bound = peak_source + args.motif_region
					if right_bound > peak_length:
						right_bound = peak_length
					if motif_start > left_bound and motif_start < right_bound:
						if motif_name in pbm_motif_occ.keys():
							pbm_motif_occ[motif_name] += 1
						else:
							pbm_motif_occ[motif_name] = 1
						peak_dict[peak_name]["attrs"]["pbm_motifs"].append('"' + motif_name + '"')
	
	# write peak occurences data frame to file
		pd.DataFrame(pbm_motif_occ, index=[1]).to_csv( \
				path_or_buf = args.out_dir + "/peak.occurences", \
				sep="\t", \
				header=True, \
				index=False)

	# annotate enriched motifs that fall within motif region
	if args.dreme_motifs:
		with open(args.dreme_motifs) as dmf:
			next(dmf)
			for line in dmf:
				line = line.rstrip().split("\t")
				peak_name = line[0]
				motif_start = int(line[3])
				motif_name = line[8].split(";")[0].split("=")[1]
				if peak_name in peak_dict.keys():
					peak_source = peak_dict[peak_name]["attrs"]["source"]
					peak_length = peak_dict[peak_name]["attrs"]["length"]
					left_bound = peak_source - args.motif_region
					if left_bound < 1:
						left_bound = 1
					right_bound = peak_source + args.motif_region
					if right_bound > peak_length:
						right_bound = peak_length
					if motif_start > left_bound and motif_start < right_bound:
						peak_dict[peak_name]["attrs"]["dreme_motifs"].append('"' + motif_name + '"')

	# open output files
	csv_out = open(args.out_dir + "/" + args.name + "_annotation.csv", "w")
	gff_out = open(args.out_dir + "/" + args.name + "_annotation.gff", "w")

	# print header
	print >> csv_out, "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27}".format(\
			"seqid",\
			"source",\
			"type",\
			"start",\
			"stop",\
			"score",\
			"strand",\
			"phase",\
			"name",\
			"length",\
			"overlap_id",\
			"overlap_ortho_id",\
			"overlap_name",\
			"overlap_strand",\
			"left_hand_id",\
			"left_hand_ortho",\
			"left_hand_name",\
			"left_hand_dist",\
			"left_hand_strand",\
			"right_hand_id",\
			"right_hand_ortho",\
			"right_hand_name",\
			"right_hand_dist",\
			"right_hand_strand",\
			"pbm_motifs",\
			"dreme_motifs",\
			"p-value",\
			"q-value")

	# print annotations to output files
	for key in peak_dict.keys():
		seqid   = peak_dict[key]["seqid"]
		source  = peak_dict[key]["source"]
		type    = peak_dict[key]["type"]
		start   = str(peak_dict[key]["start"])
		stop    = str(peak_dict[key]["stop"])
		score   = str(peak_dict[key]["score"])
		strand  = peak_dict[key]["strand"]
		phase   = peak_dict[key]["phase"]
		attrs   = "Name=" + peak_dict[key]["attrs"]["name"] + ";" +\
			"Length=" + str(peak_dict[key]["attrs"]["length"]) + ";" +\
			"OverlapID=" + peak_dict[key]["attrs"]["ogene"]["id"] + ";" +\
			"OverlapOrtho=" + peak_dict[key]["attrs"]["ogene"]["ortho"] + ";" +\
			"OverlapName=" + peak_dict[key]["attrs"]["ogene"]["name"] + ";" +\
			"OverlapStrand=" + peak_dict[key]["attrs"]["ogene"]["strand"] + ";" +\
			"LeftHandID=" + peak_dict[key]["attrs"]["ugene"]["id"] + ";" +\
			"LeftHandOrtho=" + peak_dict[key]["attrs"]["ugene"]["ortho"] + ";" +\
			"LeftHandName=" + peak_dict[key]["attrs"]["ugene"]["name"] + ";" +\
			"LeftHandStrand=" + peak_dict[key]["attrs"]["ugene"]["strand"] + ";" +\
			"LeftHandDist=" + str(peak_dict[key]["attrs"]["ugene"]["dist"]) + ";" +\
			"RightHandID=" + peak_dict[key]["attrs"]["dgene"]["id"] + ";" +\
			"RightHandOrtho=" + peak_dict[key]["attrs"]["dgene"]["ortho"] + ";" +\
			"RightHandName=" + peak_dict[key]["attrs"]["dgene"]["name"] + ";" +\
			"RightHandStrand=" + peak_dict[key]["attrs"]["dgene"]["strand"] + ";" +\
			"RightHandDist=" + str(peak_dict[key]["attrs"]["dgene"]["dist"]) + ";" +\
			"PBMMotifs=" + ",".join(set(peak_dict[key]["attrs"]["pbm_motifs"])) + ";" +\
			"DremeMotifs=" + ",".join(set(peak_dict[key]["attrs"]["dreme_motifs"])) + ";" +\
			"Pvalue=" + str(peak_dict[key]["attrs"]["pvalue"]) + ";" +\
			"Qvalue=" + str(peak_dict[key]["attrs"]["qvalue"])
		print >> gff_out, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(\
			seqid,\
			source,\
			type,\
			start,\
			stop,\
			score,\
			strand,\
			phase,\
			attrs)
		print >> csv_out, "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27}".format(\
			seqid,\
			source,\
			type,\
			start,\
			stop,\
			score,\
			strand,\
			phase,\
			peak_dict[key]["attrs"]["name"],\
			str(peak_dict[key]["attrs"]["length"]),\
			peak_dict[key]["attrs"]["ogene"]["id"],\
			peak_dict[key]["attrs"]["ogene"]["ortho"],\
			peak_dict[key]["attrs"]["ogene"]["name"],\
			peak_dict[key]["attrs"]["ogene"]["strand"],\
			peak_dict[key]["attrs"]["ugene"]["id"],\
			peak_dict[key]["attrs"]["ugene"]["ortho"],\
			peak_dict[key]["attrs"]["ugene"]["name"],\
			str(peak_dict[key]["attrs"]["ugene"]["dist"]),\
			peak_dict[key]["attrs"]["ugene"]["strand"],\
			peak_dict[key]["attrs"]["dgene"]["id"],\
			peak_dict[key]["attrs"]["dgene"]["ortho"],\
			peak_dict[key]["attrs"]["dgene"]["name"],\
			str(peak_dict[key]["attrs"]["dgene"]["dist"]),\
			peak_dict[key]["attrs"]["dgene"]["strand"],\
			" ".join(set(peak_dict[key]["attrs"]["pbm_motifs"])),\
			" ".join(set(peak_dict[key]["attrs"]["dreme_motifs"])),\
			str(peak_dict[key]["attrs"]["pvalue"]),\
			str(peak_dict[key]["attrs"]["qvalue"]))

	csv_out.close()
	gff_out.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("\nUser interrupt! Ciao!\n")
		sys.exit(0)
