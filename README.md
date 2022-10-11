# pairwiseGraph
Script programming project in 2nd year of Master's Degree.


### What is pairwiseGraph?
This script aims to generate a graph to easily visualize pairwise alignment distances between sequences.
Other files are created :
- A compact csv file containing names of sequences aligned, the corresponding align and GC contents of sequences.
- A directory containing the formated alignment.


## Results format
### Graph result example
![Graph result image](https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph/blob/main/figureReadme.png)

- Every node represent a sequence and every edge is the alignment score.

- The node size varies with percentage of 'G' and 'C' bases throughout the sequence. The minimum and maximum size of a node is represented by minGC and maxGC.

- The visible edges on the graph are the ones that score higher than the choosen treshold. The color goes from blue to red to represent the actual alignment score. If the edge is blue, the alignment score is near the treshold score. The edge is red if the alignment score is near 100% which can be seen with the edge exactAlign1 - exactAlign2.


## Installation
- Cloning the repository
```bash
git clone https://github.com/Damien-Garcia-Bioinformatics/pairwiseGraph.git
```

- Creation of virtual environment to run the script
```bash
virtualenv venvPairwiseGraph -p python3 &&
source venvPairwiseGraph/bin/activate &&
pip install -r requirements.txt
```
- Running the script
```plaintext
python3 pairwiseGraph.py [optional arguments]
  sequences : path to fasta file containing sequences.
  graphName : set file name for the graph representation (example.png)
  saveAlign : save all alignment files in a dedicated repertory "align" (-yes or -no)
  threshold  : set a treshold score for edges to visualize on graph (float)
```


## Running the script
You can use ```python3 pairwiseGraph.py help``` to print execution help sheet.

To execute the script, multiple ways of providing a fasta file are possible.
- Provide a file name in 'scriptParam.txt'. After cloning the repository, parameters are set to execute the script with the example sequences.
This option was created for non-programmers users to be able to easily use this script and tweak the options without having to look through the code.

- Provide a file name through command line parameter : This method will override every other parameters. This is the prefered option if you want to use this script inside of a pipeline.

- If no file is provided to the script, random sequences are generated. You can change the parameters of the random sequence generator in the 'scriptParam.txt' file.


## Roadmap
### Done
- [x] Fasta reader function
- [x] Construction of dictionary structure containing fasta sequence data
- [x] GCcontent calculation function
- [x] Graph construction considering GCcontent and pairwise distance values
- [x] Function to save alignments files
- [x] File containing script parameters
- [x] Add parameters and options to tweak script execution (having the option to change the name of generated graph, etc)
- [x] Help function
- [x] An easter egg because programmation should be fun ;)
### Work in progress
- [ ] Caption and graph legend


## Easter Egg

The directory "asciiArt" is just here as an easter egg and is used if the user doesn't provide a fasta file for the script to run. It is of course not a necessary directory for to script to be fully working. You can delete this directory (asciiArt) at any time... :'(


## Full list of requirements
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
