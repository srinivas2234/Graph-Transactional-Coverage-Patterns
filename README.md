Graph Transactional Coverage Patterns:

This file provides details about dataset and steps to execute graph transactional coverage patterns.
Dataset Format:
SMILES Format: C1CCCC1C
Graph Transaction Format: A Graph transaction is represented as follows:

t # 0 --> represents graph transaction id 0

v 0 a --> represents vertex id 0 and vertex label a

v 1 b --> represents vertex id 1 and vertex label b

v 2 b --> represents vertex id 2 and vertex label b
v 3 c --> represents vertex id 3 and vertex label c
e 0 1 x --> represents edge between vertices 0 and 1 with label x
e 1 2 y --> represents edge between vertices 1 and 2 with label y
e 2 3 x --> represents edge between vertices 2 and 3 with label x
e 3 0 z --> represents edge between vertices 3 and 0 with label z

If your dataset is in SMILES Format which represent chemical compounds, Please convert it into graph transactions format shown above using the following steps:
Create a new conda environment using:
$ conda create -n cigin
Installing RDKit:
$ conda install -c rdkit rdkit==2019.03.1
Installing other dependencies:
$ conda install -c pytorch pytorch
$ pip install dgl (Please check here for installing for different cuda builds)
$ pip install numpy
$ pip install pandas
$conda activate cigin
$python molecular_graph_Final.py

This will read Smilesdata.txt from /Dataset/ and generates Smilesdata_to_GT.txt in same directory which contain graph transactions corresponding to smiles data.

Once you have graph transactional dataset:
Step-1: Compute FID based flat transactions using "python3 -m gspan_mining -s support -p True ./Dataset/filename.txt". /*This will produce flat transactions with name "filename_support_Flat_tra.txt" in ./Dataset/ directory
Step-2: Compute graph transactional coverage patterns using "python SetCoverProblem_gSpan.py minTC minTPC maxOR filename_support_Flat_tra 'True' ". /*This will produce graph transactional coverage patterns with name "filename_support_Flat_tra.txt" in ./Dataset/ directory.
If the fifth argument is 'True' it will write all GTCP to "filename_support_Flat_tra_GTCPs.txt"

The file filename_support_Flat_tra_Results.txt provides information about execution time, No.of candidate patterns and number of GTCPs.

Note:In Extracted GTCPs each line is GTCP with transactional pattern coverage and overlap ratio.
For example: The GTCP [['355', '516'], 1.0, 0.22] represents transactional pattern coverage = 1.0 and overlap ratio=0.22"
