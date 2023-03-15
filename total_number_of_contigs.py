import logging
logging.basicConfig(
    level = logging.INFO,
    format = "{asctime} {levelname:<8} {message}",
    style = "{",
    filename= "PipelineProject.log",
    filemode = "a")

#input file is the SPAdes_assembly file to grab the total number of contigs
#contigs.fasta might not use the "^>" grep check
os.system(grep -c "^>" contigs.fasta | wc -l)


import Bio
from Bio import SeqIO

large_contigs = []
for seq_record in SeqIO.parse("SPAdes_assembly/contigs.fasta", "fasta"):
    if len(seq_record.seq) > 1000:
        large_contigs.append(seq_record.seq)

num_large_contigs = len(large_contigs)

os.system(grep -c "" contigs.fasta | wc -l)

logging.info("There are %s contigs > 1000 bp in the assembly." %num_large_contigs)


bps_in_assembly = 0
for contig in large_contigs:
    bps_in_assembly+= len(contig)

logging.info("There are %s bp in the assembly." %bps_in_assembly)
