import logging
logging.basicConfig(
    level = logging.INFO,
    format = "{asctime} {levelname:<8} {message}",
    style = "{",
    filename= "PipelineProject.log",
    filemode = "w")

#before files are the raw SRR56600XX.fastq(fasta?)
#find a character or series of characters that is only found once in each read
#use grep -c to find it



#after files are the _mapped_ files
#again, find a character or series of characters that is only found once in each read
#use grep -c to find it



logging.info("Donor 1 (2dpi) had %s" % before_reads " read pairs before Bowtie2 filtering and %s" % after_reads " read pairs after.\n")

logging.info("Donor 1 (6dpi) had %s" % before_reads " read pairs before Bowtie2 filtering and %s" % after_reads " read pairs after.\n")





logging.info("Donor 3 (2dpi) had %s" % before_reads " read pairs before Bowtie2 filtering and %s" % after_reads " read pairs after.\n")

logging.info("Donor 3 (6dpi) had %s" % before_reads " read pairs before Bowtie2 filtering and %s" % after_reads " read pairs after.\n")
