# pairwiseGraph
Script programming project in 2nd year of Master's Degree.

### What is pairwiseGraph?
This script aims to generate a graph to easily visualize pairwise alignment distances between sequences.

## Results format
### Graph result example
![Graph result image](https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph/blob/main/generated_graph.png)

Every node represent a sequence and every edge is the alignment score.

The node size varies with percentage of 'G' and 'C' bases throughout the sequence. The minimum and maximum size of a node is represented by minGC and maxGC.

The visible edges on the graph are the ones that score higher than the choosen treshold. The color goes from blue to red to represent the actual alignment score. If the edge is blue, the alignment score is near the treshold score. The edge is red if the alignment score is near 100% which can be seen with the edge exactAlign1 - exactAlign2.

## How to use the script?

### Installation
```bash
git clone https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph.git
pip install -r requirements
python3 pairwiseGraph.py [file]
```

### Code explanation


## Roadmap
- [x] Fasta reader function
- [x] Construction of dictionary structure containing fasta sequence data
- [x] GCcontent calculation function
- [x] Graph construction considering GCcontent and pairwise distance values
- [x] Function to save alignments files
- [x] File containing script parameters
- [ ] Caption and graph legend

## Requirements
Running the script requires matplotlib, pyqt5, netgraph, biopython and sub-requirements from previously cited libraries.

You can use ``` pip install matplotlib pyqt5 netgraph bio ``` or ```pip install -r requirements.txt``` to install everything needed to run the script.

### Full list of requirements
```bash
bio==1.4.0
biopython==1.79
biothings-client==0.2.6
certifi==2022.9.14
charset-normalizer==2.1.1
contourpy==1.0.5
cycler==0.11.0
fonttools==4.37.2
grandalf==0.7
idna==3.4
kiwisolver==1.4.4
matplotlib==3.6.0
mygene==3.2.2
netgraph==4.9.6
numpy==1.23.3
packaging==21.3
Pillow==9.2.0
pyparsing==3.0.9
PyQt5==5.15.7
PyQt5-Qt5==5.15.2
PyQt5-sip==12.11.0
python-dateutil==2.8.2
rectangle-packer==2.0.1
requests==2.28.1
scipy==1.9.1
six==1.16.0
tqdm==4.64.1
urllib3==1.26.12
```

## Author
Damien GARCIA, M2BB
