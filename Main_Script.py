import os

#create a directory for all files to be added to
directory = "PipelineProject_Alec_Loftus"
parent_dir = os.path.abspath(os.getcwd())

path = parent_dir + "/" + directory
os.mkdir(path)
print("Directory '%s' created" % directory)

os.chdir(directory)
print(os.path.abspath(os.getcwd()))

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

#test fastq-dump for how many reads to take
#or cut however many lines (multiple of 4 left)

#STEP 2#############################################################################################################################################################

#pull reference genome from repository file?
#HCMV_reference_genome = ../Genome_Assembly_Pipeline

#creating the HCMV index
#os.system(bowtie2-build HCMV_reference_genome HCMV_index)


#run bowtie2 for Donor1_2dpi
#is there a way to loop this????
#os.system(bowtie2 --quiet -x HCMV_index -1 Donor1/SRR5660030_2dpi_1.fastq -2 Donor1/SRR5660030_2dpi_2.fastq -s Donor1_2dpi_map.sam --al-conc-gz Donor1/SRR5660030_2dpi_mapped_%.fq.gz)


#STEP 3#############################################################################################################################################################

#spades to assemble all 4 transcriptomes together

#os.system(SPAdes_assembly_script.py)



#STEP 4#############################################################################################################################################################

#os.system(total_number_of_contigs.py)


#STEP 5#############################################################################################################################################################

