import os
import sys
#create a directory for all files to be added to
directory = "PipelineProject_Alec_Loftus"
parent_dir = os.path.abspath(os.getcwd())

path = parent_dir + "/" + directory
os.mkdir(path)
print("Directory '%s' created" % directory)
path = parent_dir + "/" + directory
os.mkdir(path)
print("Directory '%s' created" % directory)

os.chdir(directory)
print('All outfiles will be written within the path ' + os.path.abspath(os.getcwd()))

dataset_type = ""
while dataset_type != "full" or dataset_type != "test":
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
    filenames = ['SRR5660030','SRR5660033','SRR5660044','SRR5660045']
    dataset_path = '../test_data/'
#STEP 2#############################################################################################################################################################

import logging
logging.basicConfig(
    level = logging.INFO,
    format = '{asctime} {levelname:<8} {message}',
    style = '{',
    filename = 'PipelineProject.log',
    filemode = 'w')


reads_before = []
for file in filenames:
    file_length_command = 'wc -l ' + directory_path + file + '_1.fastq'
    read_length = os.system(file_length_command)
    reads_before.append(int(read_length)/4)





#pull reference genome from repository file
HCMV_genome = "HCMV_reference_genome.fasta"
#create HCMV index using bowtie2-build command
bowtie2_build_command = 'bowtie2-build ' + HCMV_genome + ' HCMV_index'
os.system(bowtie2_build_command)

for file in filenames:
    bowtie2_command = 'bowtie2 --quiet -x HCMV_index -1 ' + dataset_path + file + '_1.fastq -2 ' + dataset_path + file + '_2.fastq -s ' + file + 'map.sam --al-conc-gz ' + file + '_mapped_%.fq.gz'
    os.system(bowtie2_command)


reads_after = []
for file in filenames:
    file_length_command = 'wc -l ' + file +  '_mapped_1.fq.gz'
    read_length = os.system(file_length_command)
    reads_after.append(int(read_length)/4)


logging.info("Donor 1 (2dpi) had " + reads_before[1] + " read pairs before Bowtie2 filtering and " + reads_after[1] + " read pairs after.\n")

logging.info("Donor 1 (6dpi) had " + reads_before[2] + " read pairs before Bowtie2 filtering and " + reads_after[2] + " read pairs after.\n")

logging.info("Donor 3 (2dpi) had " + reads_before[3] + " read pairs before Bowtie2 filtering and " + reads_after[3] + " read pairs after.\n")

logging.info("Donor 3 (6dpi) had " + reads_before[4] + " read pairs before Bowtie2 filtering and " + reads_after[4] + " read pairs after.\n")



#STEP 3#############################################################################################################################################################

#spades to assemble all 4 transcriptomes together

#os.system(SPAdes_assembly_script.py)



#STEP 4#############################################################################################################################################################

#os.system(total_number_of_contigs.py)


#STEP 5#############################################################################################################################################################

