import Bio
from Bio import SeqIO

longest_contig = ""
for seq_record in SeqIO.parse("SPAdes_assembly/contigs.fasta", "fasta"):
    if len(seq_record.seq) > longest_contig:
        longest_contig_id = seq_record.id
        longest_contig = seq_record.seq

outfile = open("longest_contig.fasta","w")
outfile.write(longest_contig_id, "\n", longest_contig)
outfile.close()

import os

file_name = "Genome_Assembly_Pipeline/Betaherpesvirinae_refseqs.fasta"
db_name = "Betaherpesvirinae"
makeblast_command='makeblastdb -in '+file_name+' -out '+db_name+' -title '+db_name+' -dbtype nucl' 

os.system(makeblast_command)

input_file = 'longest_contig.fasta'
output_file = 'Betaherpesvirinae_blast_results.csv'
blast_command = 'blastn -query '+input_file+' -db Betaherpesvirinae -out '+output_file' + -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"

os.system(blast_command)

