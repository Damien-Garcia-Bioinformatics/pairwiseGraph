# pairwiseGraph
Script programming project in 2nd year of Master's Degree.

## What is pairwiseGraph?
This script aims to generate a graph to easily visualize pairwise alignment distances between sequences. Every node represent a sequence and every edge is the alignment score. The node size varies with percentage of 'G' and 'C' bases throughout the sequence. The visible edges on the graph are the ones that score higher than the choosen treshold. The color goes from blue to red to represent the actual alignment score. If the edge is blue, the alignment score is near the treshold score. The edge is red if the alignment score is near 100%.

## Results format


## How to use the script?


## Requirements


## Roadmap
[x] Fasta reader function
[x] Construction of dictionary structure containing fasta sequence data
[x] GCcontent calculation function
[x] Graph construction considering GCcontent and pairwise distance values
[x] Function to save alignments files
[x] File containing script parameters
[ ] Caption and graph legend
