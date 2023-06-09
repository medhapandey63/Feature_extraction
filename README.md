# Feature_extraction
>This program will calculate different features based on protein sequence and associated mutations 

This is a python based program which can be used to calculate sequence based features from protein sequences and also be used to derive properties based on mutations.

You can download the librarirs using "requirements.txt" file with following command: <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; pip install -r requirements.txt
   
This program uses several other methods for obtaining features such as: <br />
	**1.** Psi-BLAST for Position Spcifica Scoring matrices (PSSM) (https://www.ncbi.nlm.nih.gov/books/NBK279690/) <br />
	**2.** AACon for obtaining conservation scores (https://www.compbio.dundee.ac.uk/aacon/docs/library.html) <br />
	**3.** NetsurfP for predicted secondary structures (https://services.healthtech.dtu.dk/services/NetSurfP-2.0/) <br />
   	For example file, I have provided the relevant features to be calculated in the data directory. User can see the detailed description of the above mentioned programs in the provide links.
   	
This program calculates following properties: <br />
	**1.** Sequence based: <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Molecular weight, secondary structures, physicochemical properties, compositions of amino acids, encoded feature <br />
	**2.** Properties for the mutations: <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Physicochemical properties <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;PSSMs <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Conservation scores <br />
		
**User requirements:** <br />
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**1.** An input file containing the space separated information on fasta files name and comma separated mutations (input_file.txt) <br />
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**2.** Fasta file for the respective protein in "example" directory (protein.fasta)<br />
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**3.** Required feature files for the protein of interest in "data" directory (protein.features, protein.pssm, protein_netsurfp.csv) <br />
	
		
**How to use the code:** <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**1.** Extract the zipped file and keep the directories in the exact paths
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**2.** To run the program, use the following command from the terminal: <br />
	$ python sample_code.py --input_file input_file.txt
