#! /usr/bin/env python3

# Script Programming Project
# GARCIA Damien, M2BB

# Please consider reading the README file for explanation of the code.


### Library imports ###

import os
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from random import randint
from netgraph import Graph
from matplotlib import pyplot


### Functions ###

def help() :
	print('''
Error : pairwiseGraph.py execution error : Fasta file not found.
Syntax is : ./pairwiseGraph [optional arguments] 
\tsequences : path to fasta file containing sequences.
\tgraphName : set file name for the graph representation (example.png)
\tsaveAlign : save all alignment files in a dedicated repertory "align" (-yes or -no)
\ttreshold  : set a treshold score for edges to visualize on graph (float)

Go to : https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph for more informations.
	''')
	exit(1)


# Extracts parameters from a plain text file and returns a dictionary
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
				elif words[0] == "number_of_sequences" :
					dicParam["number_of_sequences"] = int(words[1])
					continue
				elif words[0] == "minimum_length" :
					dicParam["minimum_length"] = int(words[1])
					continue
				elif words[0] == "maximum_length" :
					dicParam["maximum_length"] = int(words[1])
					continue
				elif words[0] == "filename" :
					dicParam["filename"] = words[1]
					continue
				elif words[0] == "save_alignment_files" :
					if words[1] == 'Yes' :
						dicParam["save_alignment_files"] = True
						continue
					dicParam["save_alignment_files"] = False
					continue
				elif words[0] == "graph_name" :
					dicParam["graph_name"] = words[1]
					continue
				elif words[0] == "treshold" :
					dicParam["treshold"] = float(words[1])
					continue
	return dicParam


# Generates a fasta file containing nucleic acid sequences
def fasta_gen(filename, nbSeq, minLen, maxLen) :
	tablebase = 'ATCG'
	sequences = []

	# Generation of fasta sequences
	for i in range(nbSeq) :
		seq = ''
		for j in range(randint(minLen, maxLen)) :
			seq += tablebase[randint(0,3)]
		sequences.append(seq)

	# Writing of generated sequences in fasta file
	with open(filename, 'w+') as file :
		for i in range(len(sequences)) :
			file.write(f">random{i}\n{sequences[i]}\n")


# Calculates the percentage of GC in nucleic acid sequence
def GCcontent(seq) :
	seqLen = len(seq)
	count = 0
	for i in range(seqLen) :
		if seq[i] in ['G','C'] :
			count += 1
	return ((count/seqLen)*100)


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
def beautiful_print(result, seq1, seq2) :
	maxSize = 80
	p = len(result[0]) // maxSize
	r = len(result[0]) % maxSize
	if r != 0 :
		p += 1
	if not os.path.exists("align") :
		os.makedirs("align")
	with open(f"align/{seq1}_{seq2}.align", "w+") as f :
		f.write(f"{seq1} - {seq2} result alignment :\n\n")
		for i in range(p) :
			f.write(f"{result[0][i*maxSize:(i+1)*maxSize]}\n")
			f.write(f"{result[1][i*maxSize:(i+1)*maxSize]}\n")
			f.write(f"{result[2][i*maxSize:(i+1)*maxSize]}\n\n")
		f.write(result[3])


# Generates RGBA color format considering weight and treshold parameters.
def colors_RGBA(weight, treshold) :
	R = min((weight*(1/treshold)-(2*treshold)), 1)
	G = 0
	B = 1 - R
	A = R
	return [R,G,B,A]


# Generates a graph visualization of pairwise sequence alignments.
def graph_gen(dicFasta, treshold, saveAlign, graphName) :
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
						beautiful_print(result, lstID[i], lstID[j])
				
				# List containing length of longest sequence between the two sequences used for pairwise alignment
				seqLen.append(max(len(dicFasta[lstID[i]]['sequence']), len(dicFasta[lstID[j]]['sequence'])))

	# Nodes parameters
	nodesSize = {}
	for i in range(len(lstID)) :
		nodesSize[lstID[i]] = (0.5+(dicFasta[lstID[i]]['GCcontent'] / 100))*4

	# Edges parameters
	weights = []
	edgesColor = {}
	edgesAlpha = {}
	for i in range(len(edges2draw)) :
		weights.append((edges2draw[i][2]/seqLen[i])) #Poids entre 0 et 1
		edgesColor[edges2draw[i][0], edges2draw[i][1]] = colors_RGBA(weights[i], treshold)
		if edgesColor[edges2draw[i][0],edges2draw[i][1]][3] < 0 :		# Case where weight < treshold
			edgesColor[edges2draw[i][0],edges2draw[i][1]] = [0,0,0,0]
			edgesAlpha[edges2draw[i][0],edges2draw[i][1]] = 0 
		else :															# Case where weight >= treshold
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
	pyplot.savefig(graphName, format='png')
	pyplot.show()


### Script execution ###

if __name__ == '__main__' :
	# Extraction of parameters from "script_parameters.txt"
	dicParam = {}
	dicParam['save_alignment_files'] = False
	dicParam = read_param_file()

	# Extraction of command line parameters : Parameters from command line have a higher priority which means that they will replace the ones set in scriptParam.txt.
	if len(sys.argv) >= 2 :
		for i in range(len(sys.argv)-1) :
			if sys.argv[i] == "sequences" :
				dicParam['file'] = sys.argv[i+1]
				continue
			if sys.argv[i] == "saveAlign" and (sys.argv[i+1].upper() == "-YES" or sys.argv[i+1].upper() == "-Y") :
				dicParam['save_alignment_files'] = True
				continue
			if sys.argv[i] == "graphName" :
				dicParam['graph_name'] = sys.argv[i+1]
				continue
			if sys.argv[i] == "treshold" :
				dicParam['treshold'] = float(sys.argv[i+1])
				continue

	# # If a fasta file is given in command line parameter
	# if len(sys.argv) == 2 and os.path.exists(sys.argv[1]) :
	# 	dicFasta = read_fasta(sys.argv[1])

	# If a fasta file is given in script parameter file
	if not dicParam["file"] == "None" and os.path.exists(dicParam["file"]):
		dicFasta = read_fasta(dicParam["file"])
	# If no fasta file is provided
	elif dicParam["file"] == "None" :
		fasta_gen(dicParam['filename'], dicParam['number_of_sequences'], dicParam['minimum_length'], dicParam['maximum_length'])
		dicFasta = read_fasta(dicParam['filename'])
	else :
		help()

	# Creation and visualization of graph
	graph_gen(dicFasta, dicParam['treshold'], dicParam['save_alignment_files'], dicParam['graph_name'])
