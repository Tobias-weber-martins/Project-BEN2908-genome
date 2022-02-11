from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import csv
import sys
import os

list_regions = [(261763, 291329), (315134, 333163), (700643, 710124), (907553, 915900), (1162719, 1213653), (1314505, 1367765), (1475467,1484180), (1528787, 1536017), (1993833, 2061837), (2074427, 2107866), (2111619, 2152727), (2171638, 2186768), (2528973, 2570476), (2704116, 2715523), (2832786, 2890922), (3049735, 3085267), (3235318, 3252268), (3318817, 3326487), (3339011, 3344062), (3354829, 3357997), (3563844, 3572274), (3723195, 3732610), (3768643, 3772952), (3816306, 3829417), (4002999, 4056152), (4224069, 4244202), (4331263, 4362554), (4388561, 4395751), (4527238, 4559762), (4656646, 4700398), (4843537, 4898144), (4916549, 4936724), (4983971, 5029411)]


for file in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/Strains"):
    file_split = file.split('.')
    if file_split[1] == 'fa':
        '''#Descomentar para blastar e produzir os arquivos .xml
        cline= NcbiblastnCommandline(cmd='blastn', query='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/Regions/All_regions.fasta',
                                   subject='Strains/' + file, reward=1, penalty = -2, gapopen = 5, gapextend = 2, outfmt = 5,
                                   out='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/xml_' +  file_split[0])
        os.system(str(cline))
        print (cline)
        '''
        print (file_split)
        with open ('/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/xml_' +  file_split[0],'r') as result_read:
            for region, record in enumerate(NCBIXML.parse(result_read)):
                print ('REGION ' + str(region + 1) + '\n')
                #input()
                for alignment in record.alignments:
                    for hsp in alignment.hsps:   
                        #print (hsp)
                        if (hsp.align_length >= 450 and hsp.expect == 0) or hsp.align_length >= 1000:
                            region_st = list_regions[region][0]
                            #region_end = list_regions[region][1]
                            query_st = hsp.query_start
                            query_end = hsp.query_end
                            Qbeg = query_st - 1 + region_st
                            Qend = query_end - 1 + region_st
                            Sbeg = hsp.sbjct_start
                            Send = hsp.sbjct_end
                            scor = hsp.score
                            expct = hsp.expect
                            print ('tamanho do alinhamento: ' + str(hsp.align_length))
                            print ('score: ' + str(scor))
                            print ('expect: ' + str(expct))
                            print ('Query: ' + str(Qbeg), str(Qend))
                            print ('Subject: ' + str(Sbeg), str(Send) + '\n')
                            for record in SeqIO.parse('Strains/' + file_split[0] + '_gbk.gbk', "genbank"):
                                features = record.features
                            for feat in features:
                                seq_position_start = feat.location.start + 1
                                seq_position_end = feat.location.end
                                #print ('ALbeg =' + str(Sbeg), 'ALend =' + str(Send))
                                #print ('gbk_strt = ' + str(seq_position_start), 'gbk_end = ' +  str(seq_position_end))
                                #input()
                                ''' ARRUMAR, os condicionais ainda não estão printando exclusivamente as regiões 
                                if Sbeg < Send:                                    
                                    if (Sbeg <= seq_position_start) and (Send >= seq_position_end):
                                        print ('Caso Sbeg < Send')
                                        print ('ALi =' + str(Sbeg), 'ALf =' + str(Send))
                                        print ('Começo do gene = ' + str(seq_position_start), 'Fim do gene = ' +  str(seq_position_end))
                                        print (feat.qualifiers)#['product'])
                                        print ('\n\n')
                                        input()
                                else:
                                    if (Sbeg <= seq_position_end) and (Send >= seq_position_start):
                                        print ('Caso Sbeg > Sen[else]')
                                        print ('ALi =' + str(Send), 'ALf =' + str(Sbeg))
                                        print ('Começo do gene = ' + str(seq_position_start), 'Fim do gene = ' +  str(seq_position_end))
                                        print (feat.qualifiers)#['product'])
                                        print ('\n\n')
                                        input()
                                '''
