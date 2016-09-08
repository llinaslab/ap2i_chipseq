#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Description: Annotate a multi-intersect file created from ChIP-seq peak overlaps

Version = 0.1
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
import pybedtools as pybed

# define main function
def main():

	parser = argparse.ArgumentParser(description = "Script to annotate ChIP-seq peaks")
	parser.add_argument("-m", "--multi-inter", dest = "minter", required = True, help = "Multi-intersect BED file")
	parser.add_argument("-f", "--features", dest = "features", required = True, help = "Features GFF3 file")
	parser.add_argument("-g", "--genome", dest = "genome", required = True, help = "BED genome file or chromosome lengths")
	parser.add_argument("-n", "--name", dest = "name", default = "mintersect", help = "File name prefix")
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	# create bedtool objects
	minter = pybed.BedTool(args.minter).sort()
	features = pybed.BedTool(args.features)

	# create dictionary of every peak name
	minter_dict = dict()
	for interval in minter:
		num = interval.score.split("_")[0]
		reps = interval.score.split("_")[1]
		minter_dict[interval.name] = {\
				"seqid"     : interval.chrom,\
				"source"    : "annotation",\
				"type"      : "intersect_peaks",\
				"start"     : interval.start + 1,\
				"stop"      : interval.stop,\
				"score"     : ".",\
				"strand"    : ".",\
				"phase"     : ".",\
				"attrs"     : {\
						"name"         : interval.name,\
						"length"       : interval.length,\
						"ugene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : "", "dist" : "-1"},\
						"dgene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : "", "dist" : "-1"},\
						"ogene"        : {"id" : ".", "ortho" : ".", "name" : ".", "strand" : ""},\
						"reps"         : {"num" : num, "reps" : reps }
						}
				}

	# find overlapping and closest genes and create new interval files
	feature_overlap = minter.intersect(features, wb=True).saveas("tmpa.bed")
	upstream_features = minter.closest(features, D="ref", id=True, io=True, g=args.genome).saveas("tmpb.bed")
	downstream_features = minter.closest(features, D="ref", iu=True, io=True, g=args.genome).saveas("tmpc.bed")

	# annotate peaks that overlap ORFS
	for interval in feature_overlap:
		if interval.name in minter_dict.keys():
			if interval[13] != ".":
				attributes = interval[13].split(";")
				minter_dict[interval.name]["attrs"]["ogene"]["strand"] = interval[11]
				if attributes[0].split("=")[0] == "ID":
					minter_dict[interval.name]["attrs"]["ogene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						minter_dict[interval.name]["attrs"]["ogene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							minter_dict[interval.name]["attrs"]["ogene"]["name"] = attributes[2].split("=")[1]
	
	# annotate peaks that are upstream of negative stranded genes
	for interval in upstream_features:
		if interval.name in minter_dict.keys():
			if interval[13] != ".":
				attributes = interval[13].split(";")
				minter_dict[interval.name]["attrs"]["ugene"]["strand"] = interval[11]
				minter_dict[interval.name]["attrs"]["ugene"]["dist"] = int(interval[14]) * -1
				if attributes[0].split("=")[0] == "ID":
					minter_dict[interval.name]["attrs"]["ugene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						minter_dict[interval.name]["attrs"]["ugene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							minter_dict[interval.name]["attrs"]["ugene"]["name"] = attributes[2].split("=")[1]

	# annotate peaks that are upstream of positive stranded genes
	for interval in downstream_features:
		if interval.name in minter_dict.keys():
			if interval[13] != ".":
				attributes = interval[13].split(";")
				minter_dict[interval.name]["attrs"]["dgene"]["strand"] = interval[11]
				minter_dict[interval.name]["attrs"]["dgene"]["dist"] = int(interval[14]) * -1
				if attributes[0].split("=")[0] == "ID":
					minter_dict[interval.name]["attrs"]["dgene"]["id"] = attributes[0].split("=")[1]
				if len(attributes) > 1:
					if attributes[1].split("=")[0] == "ratt_ortholog":
						minter_dict[interval.name]["attrs"]["dgene"]["ortho"] = attributes[1].split("=")[1]
					if len(attributes) > 2:
						if attributes[2].split("=")[0] == "Name":
							minter_dict[interval.name]["attrs"]["dgene"]["name"] = attributes[2].split("=")[1]

	# open output files
	csv_out = open(args.name + "_annotation.csv", "w")
	gff_out = open(args.name + "_annotation.gff", "w")

	# print header
	print >> csv_out, "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25}".format(\
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
			"num_replicates",\
			"replicates")

	# print annotations to output files
	for key in minter_dict.keys():
		seqid   = minter_dict[key]["seqid"]
		source  = minter_dict[key]["source"]
		type    = minter_dict[key]["type"]
		start   = str(minter_dict[key]["start"])
		stop    = str(minter_dict[key]["stop"])
		score   = str(minter_dict[key]["score"])
		strand  = minter_dict[key]["strand"]
		phase   = minter_dict[key]["phase"]
		attrs   = "Name=" + minter_dict[key]["attrs"]["name"] + ";" +\
			"Length=" + str(minter_dict[key]["attrs"]["length"]) + ";" +\
			"OverlapID=" + minter_dict[key]["attrs"]["ogene"]["id"] + ";" +\
			"OverlapOrtho=" + minter_dict[key]["attrs"]["ogene"]["ortho"] + ";" +\
			"OverlapName=" + minter_dict[key]["attrs"]["ogene"]["name"] + ";" +\
			"OverlapStrand=" + minter_dict[key]["attrs"]["ogene"]["strand"] + ";" +\
			"LeftHandID=" + minter_dict[key]["attrs"]["ugene"]["id"] + ";" +\
			"LeftHandOrtho=" + minter_dict[key]["attrs"]["ugene"]["ortho"] + ";" +\
			"LeftHandName=" + minter_dict[key]["attrs"]["ugene"]["name"] + ";" +\
			"LeftHandStrand=" + minter_dict[key]["attrs"]["ugene"]["strand"] + ";" +\
			"LeftHandDist=" + str(minter_dict[key]["attrs"]["ugene"]["dist"]) + ";" +\
			"RightHandID=" + minter_dict[key]["attrs"]["dgene"]["id"] + ";" +\
			"RightHandOrtho=" + minter_dict[key]["attrs"]["dgene"]["ortho"] + ";" +\
			"RightHandName=" + minter_dict[key]["attrs"]["dgene"]["name"] + ";" +\
			"RightHandStrand=" + minter_dict[key]["attrs"]["dgene"]["strand"] + ";" +\
			"RightHandDist=" + str(minter_dict[key]["attrs"]["dgene"]["dist"]) + ";" +\
			"numReplicates=" + minter_dict[key]["attrs"]["reps"]["num"] + ";" +\
			"Replicates=" + minter_dict[key]["attrs"]["reps"]["reps"]
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
		print >> csv_out, "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25}".format(\
			seqid,\
			source,\
			type,\
			start,\
			stop,\
			score,\
			strand,\
			phase,\
			minter_dict[key]["attrs"]["name"],\
			str(minter_dict[key]["attrs"]["length"]),\
			minter_dict[key]["attrs"]["ogene"]["id"],\
			minter_dict[key]["attrs"]["ogene"]["ortho"],\
			minter_dict[key]["attrs"]["ogene"]["name"],\
			minter_dict[key]["attrs"]["ogene"]["strand"],\
			minter_dict[key]["attrs"]["ugene"]["id"],\
			minter_dict[key]["attrs"]["ugene"]["ortho"],\
			minter_dict[key]["attrs"]["ugene"]["name"],\
			str(minter_dict[key]["attrs"]["ugene"]["dist"]),\
			minter_dict[key]["attrs"]["ugene"]["strand"],\
			minter_dict[key]["attrs"]["dgene"]["id"],\
			minter_dict[key]["attrs"]["dgene"]["ortho"],\
			minter_dict[key]["attrs"]["dgene"]["name"],\
			str(minter_dict[key]["attrs"]["dgene"]["dist"]),\
			minter_dict[key]["attrs"]["dgene"]["strand"],\
			minter_dict[key]["attrs"]["reps"]["num"],\
			minter_dict[key]["attrs"]["reps"]["reps"])

	csv_out.close()
	gff_out.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("\nUser interrupt! Ciao!\n")
		sys.exit(0)
