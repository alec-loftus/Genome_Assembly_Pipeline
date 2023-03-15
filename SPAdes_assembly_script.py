import logging
logging.basicConfig(
    level = logging.INFO,
    format = "{asctime} {levelname:<8} {message}",
    style = "{",
    filename= "PipelineProject.log",
    filemode = "a")

os.system(spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 SRR5660030__mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz --pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz --pe-1 3 SRR5660044_mapped_1.fq.gz --pe-2 3 SRR5660044_mapped_2.fq.gz --pe-1 4 SRR5660045_mapped_1.fq.gz --pe-2 4 SRR5660045_mapped_2.fq.gz -o SPAdes_assembly/)
logging.info(spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 SRR5660030__mapped_1.fq.gz --pe-2 1 SRR5660030_mapped_2.fq.gz --pe-1 2 SRR5660033_mapped_1.fq.gz --pe-2 2 SRR5660033_mapped_2.fq.gz --pe-1 3 SRR5660044_mapped_1.fq.gz --pe-2 3 SRR5660044_mapped_2.fq.gz --pe-1 4 SRR5660045_mapped_1.fq.gz --pe-2 4 SRR5660045_mapped_2.fq.gz -o SPAdes_assembly/)
