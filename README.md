# Genome_Assembly_Pipeline
Python wrapper automating the execution of genome assembly for HCMV patients 2 days and 6 days post infection (dpi).

Tools included in this pipeline are:

bowtie2 for indexing HCMV and mapping reads to the index

SPAdes to produce a single assembly from the transcriptomes

Blast+ to check contig(s) against reference Betaherpesvirinae database

Biopython to parse FASTQ and FASTA files 


#Running the code
After downloading the repository, the full pipeline can be run by running the 'Main_Script.py' python script. Test data is available to quickly check the effectiveness of the pipeline. Otherwise, full sample data can be used to produce results.
Running 'Main_Script.py' will prompt user input for either 'full' or 'test' data to be run. Please specify which you would like to run by typing 'full' or 'test'. The prompt IS case sensitive, so respond exactly as the words appear.

The script will generate a directory 'PipelineProject_Alec_Loftus' where all files created during the pipeline will be generated. If this directory already exists in the path, 'Main_Script.py' will fail to execute. Please either change directories to one where it can be generated or delete the preexisting one using 'rm -r PipelineProject_Alec_Loftus/'. 


#Betaherpesvirinae database
For the Blast conducted, a reference database of Betaherpesvirinae has been preconstructed for time and replicability purposes. The code used to generate it can be found in a block comment in step 5 of the Main_Script.py. There were 15,065 fasta sequences used for the database at the time of creation. If you want to run the pipeline gathering the most updated list of Betaherpesvirinae subfamily sequences, you can remove the block comments and run this step.

