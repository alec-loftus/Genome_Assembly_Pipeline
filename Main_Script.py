import os
#create a directory for all files to be added to
directory = "PipelineProject_Alec_Loftus"
parent_dir = os.path.abspath(os.getcwd())

path = parent_dir + "/" + directory
os.mkdir(path)
print("Directory '%s' created" % directory)
os.chdir(directory)
print('All outfiles will be written within the path ' + os.path.abspath(os.getcwd()))

dataset_type = input("Which dataset are you running?(full/test): ")

if dataset_type == "full":
    #all of the SRA files to be used
    list_of_links = ["https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044","https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"]

    #do i need to do this step here or can i include the wget files
    #in my repository

    #generate list for all file names to be tracked
    filenames = []
    for link in list_of_links:
        filenames.append(link[-10:])
        command = 'wget ' + link
        os.system(command)
    #generate the paired end files
    for file in filenames:
        print('Now fastq dumping file ' + file)
        paired_command = "fastq-dump -I --split-files " + file
        os.system(paired_command)
    dataset_path = os.path.abspath(os.getcwd())
    #test fastq-dump for how many reads to take
    #or cut however many lines (multiple of 4 left)
elif dataset_type == "test":
    filenames = ['SRR5660030_test','SRR5660033_test','SRR5660044_test','SRR5660045_test']
    dataset_path = '../test_data/'

print(dataset_path)
#STEP 2#############################################################################################################################################################

import logging
logging.basicConfig(
    level = logging.INFO,
    format = '{asctime} {levelname:<8} {message}',
    style = '{',
    filename = 'PipelineProject.log',
    filemode = 'w')

def wc(filename):
    return int(check_output(['wc','-l', filename]).split()[0])

import subprocess
from subprocess import check_output

reads_before = []
for file in filenames:
    file_command = dataset_path + file + '_1.fastq'
    read_length = wc(file_command)
    before_length = int(read_length)
    reads_before.append(int(before_length/4))
print(reads_before)


#pull reference genome from repository file and create HCMV index using bowtie2-build command
bowtie2_build_command = 'bowtie2-build ../HCMV_reference_genome.fasta HCMV'
os.system(bowtie2_build_command)

for file in filenames:
    bowtie2_command = 'bowtie2 --quiet -x HCMV -1 ' + dataset_path + file + '_1.fastq -2 ' + dataset_path + file + '_2.fastq -s ' + file + 'map.sam --al-conc-gz ' + file + '_mapped_%.fq.gz'
    os.system(bowtie2_command)


reads_after = []
for file in filenames:
    file_command = file +  '_mapped_1.fq.gz'
    read_length = wc(file_command)
    after_length = int(read_length)
    reads_after.append(int(after_length/4))


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
            largest_contig_sequence = seq_record.seq
            largest_contig_id = seq_record.id

#turn integer of number of large contigs into a string for writing to the log file
counter = str(counter)
logging.info("There are %s contigs > 1000 bp in the assembly." %counter)

#turn integer of base pairs in contigs >1000bp into string for the log file
bps_in_assembly = str(bps_in_assembly)
logging.info("There are %s bp in the assembly." %bps_in_assembly)

#STEP 5#############################################################################################################################################################

#largest contig is used in Blast+ against the reference sequences in the repository

makeblast_command = 'makeblastdb -in ' + refseq_file + ' -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl'
os.system(makeblast_command)

#take the longest contig we found earlier and write it to a fasta file that can be used by the blast command as the query
outfile.open("longest_contig.fasta","w")
outfile.write(largest_contig_id,"\n",largest_contig_sequence)
outfile.close()

#Entrez from Biopython allows accessing NCBI accession numbers and 
from Bio import Entrez
Entrez.email = "aloftus1@luc.edu"
handle = Entrez.efetch(db='nucleotide',id='',rettype='fasta',
#assemble the blast command to take the longest contig as the query file
input_file = 'longest_contig.fasta'
#name the output file whatever you want in .csv format
output_file = 'HCMV_longest_contig_blast.csv'

blast_command = 'blastn -query ' + input_file + ' -db Betaherpesvirinae -out ' + output_file + ' -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_command)

