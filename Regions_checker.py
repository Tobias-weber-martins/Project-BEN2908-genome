from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import sys
import os

cline= NcbiblastnCommandline(cmd='blastn', query='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/All_regions.fasta' , subject="78-Pyelo.fa",
                             outfmt = 5, out='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/tst78')

os.system(str(cline))

print (cline)

with open ('tst78','r') as result_read:
    for region, record in enumerate(NCBIXML.parse(result_read)):
        print ('REGION ' + str(region + 1) + '\n')
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print (hsp)
                input()
