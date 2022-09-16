#! /usr/bin/env python3

# Script Programming Project
# GARCIA Damien, M2BB

# Please consider reading the README file for detailed explanation of the code.

### Library imports ###

import os 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from random import randint
from netgraph import Graph
import matplotlib.pyplot as plt


### Functions ###

def help() :
	terminalSize = os.get_terminal_size()
	print(terminalSize)


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
				elif words[0] == "number_of_sequences" :
					dicParam["number_of_sequences"] = int(words[1])
				elif words[0] == "minimum_length" :
					dicParam["minimum_length"] = int(words[1])
				elif words[0] == "maximum_length" :
					dicParam["maximum_length"] = int(words[1])
				elif words[0] == "filename" :
					dicParam["filename"] = words[1]
				elif words[0] == "treshold" :
					dicParam["treshold"] = float(words[1])
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


# Creates a 
def colors_RGBA(weight, treshold) :
	R = min((weight*(1/treshold)-(2*treshold)),1)
	G = 0
	B = 1 - R
	A = R 
	return [R,G,B,A]



#
def graph_gen(dicFasta, treshold) :
	lstID = list(dicFasta) # List of nodes name
	
	# Creates edges weighted by percentage of alignments between sequences
	seqLen = []
	edges2draw = []
	for i in range(len(lstID)) :
		for j in range(len(lstID)) :
			if i < j :
				alignments = pairwise2.align.globalxx(dicFasta[lstID[i]]['sequence'],
													  dicFasta[lstID[j]]['sequence'],
													  one_alignment_only=True)
				edges2draw.append((lstID[i], lstID[j], alignments[0][2]))

				# Writes alignments in plain text format files, in directory 'align/'
				for alignment in alignments :
					result = (pairwise2.format_alignment(*alignment)).split('\n')
					beautiful_print(result, lstID[i], lstID[j])
				
				# List containing length of longest sequence between the two sequences used for pairwise alignment
				seqLen.append(max(len(dicFasta[lstID[i]]['sequence']), len(dicFasta[lstID[j]]['sequence'])))

	# Nodes parameters
	nodesSize = {}
	for i in range(len(lstID)) :
		nodesSize[lstID[i]] = pow((0.5+(dicFasta[lstID[i]]['GCcontent'] / 100))*4,1)

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
			edgesAlpha[edges2draw[i][0],edges2draw[i][1]] = 0.7


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
		node_label_fontdict=dict(size=8),
		node_label_offset=0.14
	)
	plt.show()


### Script execution ###

if __name__ == '__main__' :
	# Extraction of parameters from "script_parameters.txt"
	dicParam = read_param_file()

	# If no fasta file is given in parameter
	if dicParam["file"] == "None" :
		fasta_gen(dicParam['filename'], dicParam['number_of_sequences'], dicParam['minimum_length'], dicParam['maximum_length'])
		dicFasta = read_fasta(dicParam['filename'])

	# If a fasta file is given in parameter
	elif not dicParam["file"] == "None" and os.path.exists(dicParam["file"]) :
		dicFasta = read_fasta(dicParam["file"])

	# Error
	else :
		print(help())
		exit(1)

	# Creation and visualization of graph
	graph_gen(dicFasta, dicParam['treshold'])