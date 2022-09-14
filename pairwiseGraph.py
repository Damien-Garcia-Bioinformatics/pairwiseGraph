#! /usr/bin/env python3

import sys
import os 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from random import randint
import networkx as nx
import matplotlib.pyplot as plt

### Functions ###

def help() :
	terminalSize = os.get_terminal_size()
	print(terminalSize)

def read_param_file() :
	dicParam = {}
	with open("script_parameters.txt") as paramFile :
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
			file.write(f">seq{i}\n{sequences[i]}\n")

def GCcontent(seq) :
	seqLen = len(seq)
	count = 0
	for i in range(seqLen) :
		if seq[i] in ['G','C'] :
			count += 1
	return ((count/seqLen)*100)

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

# Generates a graph from pairwise alignment of fasta sequences
def graph_gen(dicFasta, treshold) :
	nodes2draw = []
	edges2draw = []

	G = nx.Graph()

	# for id in dicFasta.items() :
	# 	nodes2draw.append((id[0], id[1]["GCcontent"]))

	for id in dicFasta.items() :
		G.add_nodes_from(dicFasta.keys())

	# for i in range(len(nodes2draw)) :
	# G.add_nodes_from([(id, {'sequence': seq}) for (id, seq) in dicFasta.items()])

	lstID = list(dicFasta)
	for i in range(len(lstID)) :
		for j in range(len(lstID)) :
			if not i >= j :
				minScore = max(len(dicFasta[lstID[i]]['sequence']), len(dicFasta[lstID[j]]['sequence']))*treshold
				alignments = pairwise2.align.globalxx(dicFasta[lstID[i]]['sequence'],
													  dicFasta[lstID[j]]['sequence'],
													  one_alignment_only=True)
				if alignments[0][2] > minScore :
					edges2draw.append((lstID[i], lstID[j], alignments[0][2]))
					for alignment in alignments :
						result = (pairwise2.format_alignment(*alignment)).split('\n')
						beautiful_print(result, lstID[i], lstID[j])

	# lstID = list(dicFasta.keys())
	# for i in range(len(lstID)) :
	# 	for j in range(len(lstID)) :
	# 		if not i >= j :
	# 			treshold = int((max(len(dicFasta[lstID[i]]), len(dicFasta[lstID[j]]))) * 0.61)
	# 			alignments = pairwise2.align.globalxx(dicFasta[lstID[i]], dicFasta[lstID[j]], one_alignment_only=True)
	# 			if alignments[0][2] > treshold :
	# 				G.add_edge(lstID[i], lstID[j], weight=alignments[0][2])
	# 				for alignment in alignments :
	# 					result = (pairwise2.format_alignment(*alignment)).split('\n')
	# 					beautiful_print(result, lstID[i], lstID[j])

	pos = nx.spring_layout(G)
	nx.draw(G, pos, edgelist=edges2draw, with_labels=True) # nodelist=nodes2draw,
	nx.write_graphml(G, 'pairwise.graphml', encoding='utf-8', prettyprint=True, named_key_ids=False)
	plt.show()

	# nx.draw(G, nodelist=nodes2draw, with_labels=True) #edgelist=edges2draw
	# nx.write_graphml(G, 'graph.graphml', encoding='utf-8', prettyprint=True, infer_numeric_types=False, named_key_ids=False)
	# plt.show()


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