import os

#STEP 1#######################################################################################################################################

#create a directory for all files to be added to
directory = "PipelineProject_Alec_Loftus"
#the current working directory for the user must be grabbed in order to add the new directory for generated files to go into
parent_dir = os.path.abspath(os.getcwd())

#create the path to move into of the parent directory and new directory
path = parent_dir + "/" + directory
os.mkdir(path)

#visual indicator of the directory that has been created
print("Directory '%s' created" % directory)
#move into this directory so all outfiles from each step can be organized within the folder
os.chdir(directory)
#Visual indicator of the full path for all written files in this script, for convenience
print('All outfiles will be written within the path ' + os.path.abspath(os.getcwd()))

#user will be asked to input which dataset they are running to determine which steps the script will take
#needs user input because full data is for generating results, test data is for testing efficiency and troubleshooting
#assigns the user input to a variable that then determines which path of the conditional we take
dataset_type = input("Which dataset are you running?(full/test): ")

#if running the full dataset after typing 'full'
if dataset_type == "full":
    #all of the SRA files to be used, grabbed using 'copy link' on the web pages for the files
    list_of_links = ["https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"]

    #generate list for all file names to be tracked
    #these file names are the last 10 characters of the web link (the ID)
    filenames = []
    for link in list_of_links:
        #add each ID to the list of file names
        filenames.append(link[-10:])
        #linux wget command to download the contents of each link
        command = 'wget ' + link
        os.system(command)
    #generate the paired end files using the IDs grabbed in the last loop
    for file in filenames:
        #lets the user know which file is being dumped while it's happening so they're not left in limbo
        print('Now fastq dumping file ' + file)
        #runs paired end fastq-dump for each file
        paired_command = "fastq-dump -I --split-files " + file
        os.system(paired_command)
    #since we're running the full dataset, the path is the current working directory, which is where these fastq-dumped files are located
    dataset_path = os.path.abspath(os.getcwd()) + "/"

#if dataset is 'test', then the files are the test files located in the 'test_data' directory within the main repository
elif dataset_type == "test":
    #assigns the prefixes + _test to be the filenames
    filenames = ['SRR5660030_test','SRR5660033_test','SRR5660044_test','SRR5660045_test']
    #the path to the test data is going back to the parent directory and then moving into test_data from there
    dataset_path = '../test_data/'

#Just a checkpoint to make sure you are in the directory for the dataset you want to be
print(dataset_path)
#STEP 2#############################################################################################################################################################

#logging file is needed for the outputs
#set to the .INFO level for simplicity
#formatted to provide the time of the log, the level that it was logged at, and the message attached
#written to PipelineProject.log
#filemode is 'w' because it is written only once (never called in other scripts where it would need to be appended to instead
import logging
logging.basicConfig(
    level = logging.INFO,
    format = '{asctime} {levelname:<8} {message}',
    style = '{',
    filename = 'PipelineProject.log',
    filemode = 'w')

#function wc created to check how many lines (therefore reads) are in each fastq-dump file before and after bowtie2 filtering/mapping
def wc(filename):
    #check_output runs the linux command and allows it to be assigned to the variable it is returned to
    return int(check_output(['wc','-l', filename]).split()[0])
#subprocess is where we get check_output from
import subprocess
from subprocess import check_output

#create empty list for the number of reads before bowtie2 can be appended
reads_before = []
for file in filenames:
    #full path and command grabbing of the fastq file (both full and test data files are properly grabbed here)
    file_command = dataset_path + file + '_1.fastq'
    #assign read length to the wc function above
    read_length = wc(file_command)
    #this is kind of an unnecessary step but ensures that the format is integer
    before_length = int(read_length)
    #adds the number of reads to the list in order in which the files are processed
    #linux wc -l command grabs number of lines, but there are 4 lines per read in fastq files, so divide by 4 to get number of reads
    reads_before.append(int(before_length/4))
print(reads_before)


#pull reference genome from repository file and create HCMV index using bowtie2-build command
#HCMV reference genome was grabbed using accession number 'NC_006273.2' in the NCBI nucleotide database and put in the repository already
bowtie2_build_command = 'bowtie2-build ../HCMV_reference_genome.fasta HCMV'
os.system(bowtie2_build_command)

for file in filenames:
    #runs a bowtie2 for each set of paired-end reads
    #only writes to the outfile paired-end reads that mapped to the reference genome
    bowtie2_command = 'bowtie2 --quiet -x HCMV -1 ' + dataset_path + file + '_1.fastq -2 ' + dataset_path + file + '_2.fastq -s ' + file + 'map.sam --al-conc-gz ' + file + '_mapped_%.fq.gz'
    os.system(bowtie2_command)

#same process for reading in and grabbing the number of reads after bowtie2 as above for before

reads_after = []
for file in filenames:
    #the pathing no longer matters, as both datasets after bowtie2 are written to the PipelineProject_Alec_Loftus directory
    #grabs mapped reads from the _1 paired-end file (we care about read PAIRS, so we dont need to double it or grab both mapped files per Donor and dpi
    file_command = file +  '_mapped_1.fq.gz'
    read_length = wc(file_command)
    after_length = int(read_length)
    #linux command grabs number of lines, but fastq has 4 lines per read, so divide by 4
    reads_after.append(int(after_length/4))

#Since we know the order of files, we can log directly to each Donor dpi the before and after bowtie2 reads from the lists easily
logging.info("Donor 1 (2dpi) had " + str(reads_before[0]) + " read pairs before Bowtie2 filtering and " + str(reads_after[0]) + " read pairs after.\n")

logging.info("Donor 1 (6dpi) had " + str(reads_before[1]) + " read pairs before Bowtie2 filtering and " + str(reads_after[1]) + " read pairs after.\n")

logging.info("Donor 3 (2dpi) had " + str(reads_before[2]) + " read pairs before Bowtie2 filtering and " + str(reads_after[2]) + " read pairs after.\n")

logging.info("Donor 3 (6dpi) had " + str(reads_before[3]) + " read pairs before Bowtie2 filtering and " + str(reads_after[3]) + " read pairs after.\n")



#STEP 3#############################################################################################################################################################

#spades to assemble all 4 transcriptomes together

#os.system(SPAdes_assembly_script.py)
spades_command = 'spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 ' + filenames[0] + '_mapped_1.fq.gz --pe-2 1 ' + filenames[0] + '_mapped_2.fq.gz --pe-1 2 ' + filenames[1] + '_mapped_1.fq.gz --pe-2 2 ' + filenames[1] + '_mapped_2.fq.gz --pe-1 3 ' + filenames[2] + '_mapped_1.fq.gz --pe-2 3 ' + filenames[2] + '_mapped_2.fq.gz --pe-1 4 ' + filenames[3] + '_mapped_1.fq.gz --pe-2 4 ' + filenames[3] + '_mapped_2.fq.gz -o SPAdes_assembly/'
os.system(spades_command)

logging.info('spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 ' + filenames[0] + '_mapped_1.fq.gz --pe-2 1 ' + filenames[0] + '_mapped_2.fq.gz --pe-1 2 ' + filenames[1] + '_mapped_1.fq.gz  --pe-2 2 ' + filenames[1] + '_mapped_2.fq.gz --pe-1 3 ' + filenames[2] + '_mapped_1.fq.gz --pe-2 3 ' + filenames[2] + '_mapped_2.fq.gz --pe-1 4 ' + filenames[3] + '_mapped_1.fq.gz --pe-2 4 ' + filenames[3] + '_mapped_2.fq.gz -o SPAdes_assembly/')

#STEP 4#############################################################################################################################################################

#os.system(total_number_of_contigs.py)
contigs_file = 'SPAdes_assembly/contigs.fasta'
import Bio
from Bio import SeqIO
large_contigs = []

#record the largest contig id for use in Blast+ later
largest_contig_id = ""
#record the largest contig for use in Blast+ later
largest_contig_sequence = ""
#counter variable will count each contig >1000bp in length
counter = 0
#the total number of base pairs in all contigs >1000bp will be added to this variable
bps_in_assembly = 0

#SeqIO as part of the Bio package in Biopython enables parsing of fasta files to grab each record in the file
for seq_record in SeqIO.parse("SPAdes_assembly/contigs.fasta", "fasta"):
    #only consider sequences (.seq) that are larger than 1000bp
    if len(seq_record.seq) > 1000:
        #counter moves up for each record this conditional is true for
        counter+=1
        #total bps adds the length of each sequence (the number of bases) for each record
        bps_in_assembly+= len(seq_record.seq)
        #of the records >1000bp, make note of the size of the current largest record and replace it if any subsequent sequence is longer
        if len(seq_record.seq) > len(largest_contig_sequence):
            largest_contig_sequence = str(seq_record.seq)
            largest_contig_id = str(seq_record.id)

#turn integer of number of large contigs into a string for writing to the log file
counter = str(counter)
logging.info("There are %s contigs > 1000 bp in the assembly." %counter)

#turn integer of base pairs in contigs >1000bp into string for the log file
bps_in_assembly = str(bps_in_assembly)
logging.info("There are %s bp in the assembly." %bps_in_assembly)

#STEP 5#############################################################################################################################################################

#largest contig is used in Blast+ against the reference sequences in the repository

#take the longest contig we found earlier and write it to a fasta file that can be used by the blast command as the query
outfile = open("longest_contig.fasta","w")
outfile.write(">" + largest_contig_id + "\n" + largest_contig_sequence)
outfile.close()


'''
The Betaherpesvirinae reference sequences and subsequent database were created in advance using the code in this block text
At the time of gathering, there were 15065 sequences in the Betaherpesvirinae subfamily
For exact replication, use this number of sequences. To produce more accurate results with the most updated subfamily, inclusive of new sequences added, use all NCBI ids

#Entrez from Biopython allows access to the ID numbers to fetch their fasta sequences from NCBI
from Bio import Entrez
Entrez.email = "aloftus1@luc.edu"
handle = Entrez.esearch(db='nucleotide',term='betaherpesvirinae[organism]', retmax ='20000')
record = Entrez.read(handle)
string_of_ids = ""
id_list = record["IdList"]

for item in id_list:
    string_of_ids+= item + ','

string_of_ids = string_of_ids[:-1]
print('Now fetching reference sequences for Betaherpesvirinae.')

handle = Entrez.efetch(db='nucleotide',id=string_of_ids,rettype='fasta', retmode='text')

outfile = open("Betaherpesvirinae_refseqs.fasta","w")
outfile.write(handle.read())
outfile.close()


refseq_file = "Betaherpesvirinae_refseqs.fasta"
makeblast_command = 'makeblastdb -in ' + refseq_file + ' -out ../Betaherpesvirinae/Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl'
os.system(makeblast_command)
'''
#assemble the blast command to take the longest contig as the query file
input_file = 'longest_contig.fasta'
#name the output file whatever you want in .csv format
output_file = 'HCMV_longest_contig_blast.csv'

#blast command takes the longest contig and the pregenerated Betaherpesvirinae database and creates a csv outfile
#-max_hsps 1 flag is used to remove duplicate titles and ensures only the best hit for each HSP is recorded
#out format (outfmt) 6 is used to tab delimit the items with the subsequent further specifications indicating what should be written to the .csv file
blast_command = 'blastn -query ' + input_file + ' -db ../Betaherpesvirinae/Betaherpesvirinae -max_hsps 1 -out ' + output_file + ' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_command)

#function takes the .csv outfile from blast and grabs the top 10 hits using check_output again
def blast_results(output_file):
    return check_output(['head',output_file])
#assigns the return from blast_results function to blast_output (but as a byte)
blast_output = blast_results(output_file)
#decodes the byte file into a proper string that can be logged
blast_output = blast_output.decode('utf-8')
#logs the headers for the .csv file followed by the top 10 blast hits for the longest contig to the PipelineProject.log file
logging.info('\nsacc\t\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore evalue\tstitle\t\n'+blast_output)

#A little message to thank you for running the pipeline and let you know it's done :)
print('Blast+ has been completed, thank you for running the pipeline!')
