# SCAProtein
Search, concatenate and alignment of protein sequences

## **Prerequired softwares**
Make sure python and biopython are installed.

To install biopython:

		pip install biopython

and if you do not have sudo authority, try:

		pip install --user biopython 

## **Download and Run**

After download and decompress, run the following command to make sure the multiple sequence alignment tool can run:

		chmod +x muscle3.8.31_i86linux64


Configuration (set the parameters at config.txt):
		
		BatchSize 5000     #Batch size when fetch
		
		RetMax 10000       #How many records in all will be fetched from NCBI
		
		N_Threads 4        #Number of threads
		
		WindowSize 15      #When kick out species, the window size to check gap
		
		N_Species_Gap 1    #If there is only N_Species_Gap species have value, other are all gap, then these N_Species_Gap species will be removed

To run the program:

		python main.py Options subtrate_term_file subunit_term_file subtrate_cleave_file subunit_cleave_file

		Options:

							fetch          Fetch from NCBI by give list
							
							cleave         Cleave pre-seq
							
							convert        Convert plain txt to fasta
							
							cat            Concatenate sequences

							align          Run multiple sequence alignment
							
							kick           Kick out the badly aligned sequences
							
							realign        Re-run multiple sequence alignment after kicking
							
							all            Run the whole pipeline

							clear          Clear all old files
							
Sample usage:

	python main.py all subtrate.list subunit.list subtrate_cleave.list subunit_cleave.list