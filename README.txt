This is a python based program which can be used to calculate sequence based features from protein sequences and also be used to derive properties based on mutations.

You can download the librarirs using "requirements.txt" file with following command:
   pip install -r requirements.txt
   
This program uses several other methods for obtaining features such as:
	1. Psi-BLAST for Position Spcifica Scoring matrices (PSSM) (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
	2. AACon for obtaining conservation scores (https://www.compbio.dundee.ac.uk/aacon/docs/library.html)
	3. NetsurfP for predicted secondary structures (https://services.healthtech.dtu.dk/services/NetSurfP-2.0/)
   	For example file, I have provided the relevant features to be calculated in the data directory. User can see the detailed description of the above mentioned programs in the provide links.
   	
This program calculates following properties:
	1. Sequence based:
		Molecular weight, secondary structures, physicochemical properties, compositions of amino acids, encoded feature
	2. Properties for the mutations:
		Physicochemical properties
		PSSMs
		Conservation scores
		
User requires:
	1. An input file containing the space separated information on fasta files name and comma separated mutations (input_file.txt)
	2. Fasta file for the respective protein in "example" directory (protein.fasta)
	3. Required feature files for the protein of interest in "data" directory (protein.features, protein.pssm, protein_netsurfp.csv)
	
		
HOw to use the code:
1. Extract the zipped file and keep the directories in the exact paths
2. To run the program, use the following command from the terminal:
	python sample_code.py --input_file input_file.txt



