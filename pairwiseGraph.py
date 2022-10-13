#! /usr/bin/env python3

# Script Programming Project
# GARCIA Damien, M2BB

# Please consider taking a look at my GitHub repository for in depth explanation of the code :
# https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph


### Library imports ###

import os
import sys
import time
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from random import randint, choice
from netgraph import Graph
from matplotlib import pyplot


### Functions ###

# Help function that is printed if script returns an error in execution.
def help(ask4help) :
	if os.path.exists("asciiArt") and ask4help == False :
		print(easter_egg())
		print('[pairwiseGraph] [Error] : pairwiseGraph.py execution error : Fasta file not found.\n')
	
	print('Syntax is : ./pairwiseGraph [optional arguments] ')
	print('\tjobName   : set name for the directory containing all results')
	print('\tsequences : path to fasta file containing sequences.')
	print('\tsaveAlign : save all alignment files in a dedicated repertory "align" (-yes or -no)')
	print('\tthreshold : set a threshold score for edges to visualize on graph (float)')
	print('\nFor more informations, please take a look at my GitHub repository :\nhttps://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph')
	exit(1)


# Easter egg function printing an ascii art on error message if you tried to execute the script without providing a fasta file
def easter_egg() :
	asciiFile = choice(os.listdir("asciiArt"))
	asciiString = ""
	with open(f"asciiArt/{asciiFile}") as f :
		for line in f :
			asciiString += line
	return asciiString


# Returns current time which is used for file naming purpose
def get_time() :
	return time.strftime('%y%m%d_%H%M%S')


# Extracts parameters from a plain text file (scriptParam.txt) and returns a dictionary
def read_param_file() :
	dicParam = {}
	with open("scriptParam.txt") as paramFile :
		for line in paramFile :
			line = line.rstrip('\n')
			if not line.startswith('#') :
				words = line.split(';')
				words = words[0].split('=')
				if words[0] == "file" :
					dicParam["file"] = words[1]
					continue
				elif words[0] == "job_name" and not (words[1] == 'None' or words[1] == '') :
					dicParam["job_name"] = words[1]
					continue
				elif words[0] == "job_name" :
					if words[1] in ['None',''] :
						dicParam['job_name'] = get_time()
						continue
					dicParam['job_name'] = words[1]
					continue
				elif words[0] == "number_of_sequences" :
					dicParam["number_of_sequences"] = int(words[1])
					continue
				elif words[0] == "minimum_length" :
					dicParam["minimum_length"] = int(words[1])
					continue
				elif words[0] == "maximum_length" :
					dicParam["maximum_length"] = int(words[1])
					continue
				elif words[0] == "save_alignment_files" :
					if words[1] == 'Yes' :
						dicParam["save_alignment_files"] = True
						continue
					dicParam["save_alignment_files"] = False
					continue
				elif words[0] == "threshold" :
					dicParam["threshold"] = float(words[1])
					continue
	return dicParam


# Generates a fasta file containing nucleic acid sequences
def fasta_gen(jobName, nbSeq, minLen, maxLen) :
	tablebase = 'ATCG'
	sequences = []

	# Generation of fasta sequences
	for i in range(nbSeq) :
		seq = ''
		for j in range(randint(minLen, maxLen)) :
			seq += tablebase[randint(0,3)]
		sequences.append(seq)

	# Writes generated sequences in fasta file
	with open(f"results/{jobName}/sequences.fasta", 'w+') as file :
		for i in range(len(sequences)) :
			file.write(f">random{i}\n{sequences[i]}\n")


# Calculates the percentage of GC in nucleic acid sequence
def GCcontent(seq) :
	seqLen = len(seq)
	count = 0
	for i in range(seqLen) :
		if seq[i] in ['G','C'] :
			count += 1
	return (count/seqLen)*100


# Reads a fasta file and creates a dictionary containing ID, sequence, and GC percentage of every entry
def read_fasta(file) :
	dicFasta = {}
	with open(file) as f :
		seq = ''
		for line in f :
			line = line.rstrip('\n')
			if line.startswith('>') :
				seq = ''
				id = line[1:]
				continue
			else :
				seq += line
			dicFasta[id] = {"sequence" : seq, "GCcontent" : GCcontent(seq)}
	return dicFasta


# Prints alignment results in a plain text format file. Every alignment is systematicaly named and stored in 'align' directory
def beautiful_print(jobName, result, seq1, seq2) :
	maxSize = 80
	p = len(result[0]) // maxSize
	r = len(result[0]) % maxSize
	if r != 0 :
		p += 1
	
	os.makedirs(f"results/{jobName}/align", exist_ok=True)
	with open(f"results/{dicParam['job_name']}/align/{seq1}_{seq2}.txt", "w+") as f :
		f.write(f"{seq1} - {seq2} result alignment :\n\n")
		for i in range(p) :
			f.write(f"{result[0][i*maxSize:(i+1)*maxSize]}\n")
			f.write(f"{result[1][i*maxSize:(i+1)*maxSize]}\n")
			f.write(f"{result[2][i*maxSize:(i+1)*maxSize]}\n\n")
		f.write(result[3])


# Generates RGBA color format considering weight and threshold parameters.
def colors_RGBA(weight, threshold) :
	R = min((weight*(1/threshold)-(2*threshold)), 1)
	G = 0
	B = 1 - R
	A = R
	return [R,G,B,A]


# Generates a csv format file with compact result of every alignement 
def write_csv(jobName, seq1, seq2, GCcontent1, GCcontent2, score) :
	# Opens file in append mode
	with open(f"results/{jobName}/compact_results.csv", "a") as f :
		if os.path.getsize(f"results/{jobName}/compact_results.csv") == 0:
			f.write("seq1;seq2;GCcontent1;GCcontent2;score\n")
		f.write(f"{seq1};{seq2};{GCcontent1};{GCcontent2};{score}\n")


# Generates a graph visualization of pairwise sequence alignments.
def graph_gen(dicFasta, threshold, saveAlign, jobName) :
	lstID = list(dicFasta) # List of nodes name
	
	# Creates edges weighted by percentage of alignments between sequences
	seqLen = []
	edges2draw = []
	for i in range(len(lstID)) :
		for j in range(len(lstID)) :
			if i < j :
				alignments = pairwise2.align.globalxx(
					dicFasta[lstID[i]]['sequence'],
					dicFasta[lstID[j]]['sequence'],
					one_alignment_only=True)
				edges2draw.append((lstID[i], lstID[j], alignments[0][2]))

				# Writes alignments in plain text format files, in directory 'align/'
				if saveAlign :
					for alignment in alignments :
						result = (pairwise2.format_alignment(*alignment)).split('\n')
						beautiful_print(jobName, result, lstID[i], lstID[j])
				
				# List containing length of longest sequence between the two sequences used for pairwise alignment
				seqLen.append(max(len(dicFasta[lstID[i]]['sequence']), len(dicFasta[lstID[j]]['sequence'])))

				# Writes a compact result file in csv format
				score = alignments[0][2]/seqLen[i]
				write_csv(jobName, lstID[i], lstID[j], dicFasta[lstID[i]]['GCcontent'], dicFasta[lstID[j]]['GCcontent'], score)

	# Nodes parameters
	nodesSize = {}
	for i in range(len(lstID)) :
		nodesSize[lstID[i]] = (0.5+(dicFasta[lstID[i]]['GCcontent'] / 100))*4

	# Edges parameters
	weights = []
	edgesColor = {}
	edgesAlpha = {}
	for i in range(len(edges2draw)) :
		weights.append((edges2draw[i][2]/seqLen[i])) # Weight between 0 an 1
		edgesColor[edges2draw[i][0], edges2draw[i][1]] = colors_RGBA(weights[i], threshold)
		if edgesColor[edges2draw[i][0],edges2draw[i][1]][3] < 0 :		# Case where weight < threshold
			edgesColor[edges2draw[i][0],edges2draw[i][1]] = [0,0,0,0]
			edgesAlpha[edges2draw[i][0],edges2draw[i][1]] = 0 
		else :															# Case where weight >= threshold
			edgesAlpha[edges2draw[i][0],edges2draw[i][1]] = 0.8

	# Graph generation
	Graph(
		edges2draw,
		edge_color=edgesColor,
		edge_alpha=edgesAlpha,
		edge_layout='curved',
		node_layout='circular',
		node_color='mediumpurple',
		node_size=nodesSize,
		node_edge_width=0.5,
		node_labels=True, 
		node_label_fontdict=dict(size=9),
		node_label_offset=0.14)
	pyplot.tight_layout()
	pyplot.savefig(f"results/{dicParam['job_name']}/graph.png", format='png')
	#pyplot.show()


### Script execution ###

if __name__ == '__main__' :

	if len(sys.argv) == 2 and sys.argv[1].lower() == 'help' :
		ask4help = True
		help(ask4help)
	ask4help = False

	if not os.path.exists('results/') :
		os.makedirs('results/')

	# Extraction of parameters from "script_parameters.txt"
	print("[pairwiseGraph] Reading script parameters")
	dicParam = {}
	dicParam = read_param_file()

	# Extraction of command line parameters : Parameters from command line have a higher priority which means that they will replace the ones set in scriptParam.txt.
	if len(sys.argv) >= 2 :
		for i in range(len(sys.argv)-1) :
			if sys.argv[i].lower() == "jobname" :
				dicParam['job_name'] = sys.argv[i+1]
				continue
			if sys.argv[i] == "sequences" :
				dicParam['file'] = sys.argv[i+1]
				continue
			if sys.argv[i] == "saveAlign" and (sys.argv[i+1].upper() == "-YES" or sys.argv[i+1].upper() == "-Y") :
				dicParam['save_alignment_files'] = True
				continue
			if sys.argv[i] == "threshold" :
				dicParam['threshold'] = float(sys.argv[i+1])
				continue
			if sys.argv[i] == "generation" :
				dicParam['number_of_sequences'] = int(sys.argv[i+1])
				dicParam['minimum_length'] = int(sys.argv[i+2])
				dicParam['maximum_length'] = int(sys.argv[i+3])
				continue
	
	# Creation of repertory to job results
	#print(dicParam['job_name'])
	os.makedirs(f"results/{dicParam['job_name']}/", exist_ok=True)

	# If a fasta file is given in script parameter file
	if dicParam["file"] not in ['None',''] and os.path.exists(dicParam["file"]) :
		print("[pairwiseGraph] Extraction of data from fasta file")
		dicFasta = read_fasta(dicParam["file"])
	# If no fasta file is provided
	elif dicParam["file"] in ['None',''] :
		print("[pairwiseGraph] Generating fasta sequences")
		fasta_gen(dicParam['job_name'], dicParam['number_of_sequences'], dicParam['minimum_length'], dicParam['maximum_length'])
		dicFasta = read_fasta(f"results/{dicParam['job_name']}/sequences.fasta")
	else :
		help(ask4help)

	# Creation and visualization of graph
	print("[pairwiseGraph] Generating graph file")
	graph_gen(dicFasta, dicParam['threshold'], dicParam['save_alignment_files'], dicParam['job_name'])