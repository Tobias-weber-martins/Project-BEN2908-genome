from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import csv
import sys
import os

list_regions = [(261763, 291329), (315134, 333163), (700643, 710124), (907553, 915900), (1162719, 1213653), (1314505, 1367765), (1475467,1484180), (1528787, 1536017), (1993833, 2061837), (2074427, 2107866), (2111619, 2152727), (2171638, 2186768), (2528973, 2570476), (2704116, 2715523), (2832786, 2890922), (3049735, 3085267), (3235318, 3252268), (3318817, 3326487), (3339011, 3344062), (3354829, 3357997), (3563844, 3572274), (3723195, 3732610), (3768643, 3772952), (3816306, 3829417), (4002999, 4056152), (4224069, 4244202), (4331263, 4362554), (4388561, 4395751), (4527238, 4559762), (4656646, 4700398), (4843537, 4898144), (4916549, 4936724), (4983971, 5029411)]

for file in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/Strains"):
        list_files = os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/Strains")
        num_files = len(list_files)
        file_split = file.split('.')
        '''
        if file_split[1] == 'fa':
            #Descomentar para blastar e produzir os arquivos .xml
            cline= NcbiblastnCommandline(cmd='blastn', query='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/Regions/' + rgn,
                                       subject='Strains/' + file, reward=1, penalty = -2, gapopen = 5, gapextend = 2, outfmt = 5,
                                       out='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/xml_' +  file_split[0] + '_' + rgn)
            os.system(str(cline))
            print (cline)
           '''
        for num_rgn, item in enumerate(list_regions):
            linha_csv_Region = ['Region ' + str(num_rgn+1)]
            linha_csv_CDS = ['Region ' + str(num_rgn+1)]
            print (linha_csv_Region)
            for xml_rgn_str in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/xml_Regions"):
                Albeg = (xml_rgn_str.split('_')[3].split('..')[0]) #inicio da regiao
                Alend = (xml_rgn_str.split('_')[3].split('..')[1].split('.')[0]) # final da regiao
                if ((int(Albeg),int(Alend))) == item:
                    with open ('/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Python_regions_check/xml_Regions/' +  xml_rgn_str,'r') as result_read:
                         for record in NCBIXML.parse(result_read):
                             cont =1
                             for alignment in record.alignments:
                                 for hsp in alignment.hsps:
                                     for record in SeqIO.parse('Strains/' + file_split[0] + '_gbk.gbk', "genbank"):
                                        features = record.features
                                     if (hsp.align_length >= 1000 and hsp.expect == 0) or hsp.align_length >= 1500:
                                        try:
                                            region_st = list_regions[cont-1][0]
                                        except IndexError:
                                             print('More than 33 alignments')
                                        #region_end = list_regions[region][1]
                                        query_st = hsp.query_start
                                        query_end = hsp.query_end
                                        Qbeg = query_st - 1 + region_st
                                        Qend = query_end - 1 + region_st
                                        Sbeg = hsp.sbjct_start
                                        Send = hsp.sbjct_end
                                        al_length = hsp.align_length
                                        scor = hsp.score
                                        expct = hsp.expect
                                        if cont != 1:
                                            linha_csv_Region[-1] = linha_csv_Region[-1] + '\n' + 'AL' + str(cont) + '(' + str(al_length) + ' bp)' + '(BEN2908:' + str(Qbeg) + ', BEN2908:' +  str(Qend) + ')' + '(' + xml_rgn_str.split('_')[1] + ':' +  str(Sbeg)+ ', ' + xml_rgn_str.split('_')[1] + ':' + str(Send) + ')'
                                        else:
                                            linha_csv_Region.append('AL' + str(cont) + '(' + str(al_length) + ' bp)' + '(BEN2908:' + str(Qbeg) + ', BEN2908:' +  str(Qend) + ')' + '(' + xml_rgn_str.split('_')[1] + ':' +  str(Sbeg)+ ', ' + xml_rgn_str.split('_')[1] + ':' + str(Send) + ')')
                                        print (xml_rgn_str.split('_')[1])
                                        cont_CDS = 1
                                        list_cont_CDS=[]
                                        for feat in features:
                                                CDS_position_start = feat.location.start + 1
                                                CDS_position_end = feat.location.end
                                                if (CDS_position_start and CDS_position_end in range(Sbeg, Send)) or (CDS_position_start and CDS_position_end in range(Send, Sbeg)):
                                                    gene_length = abs(CDS_position_end - CDS_position_start + 1)
                                                    list_cont_CDS.append(cont_CDS)
                                                    try:
                                                        if cont != 1:
                                                           if cont_CDS != list_cont_CDS[0]: 
                                                                linha_csv_CDS[-1] = linha_csv_CDS[-1] + ', ' + '(' + str(gene_length) + ' bp)' + feat.qualifiers['product'][0] + '[' + str(CDS_position_start) + '..' + str(CDS_position_end)  + ']\n'
                                                           else:
                                                                linha_csv_CDS[-1] = linha_csv_CDS[-1] + '\n' +  'AL ' + str(cont) + '_' + xml_rgn_str.split('_')[1] + ':' + '(' + str(gene_length) + ' bp)' + feat.qualifiers['product'][0] + '[' + str(CDS_position_start) + '..' + str(CDS_position_end) + ']\n'                                           
                                                        else:
                                                           if cont_CDS != list_cont_CDS[0]:
                                                                linha_csv_CDS[-1] = linha_csv_CDS[-1] + ', ' + '(' + str(gene_length) + ' bp)' + feat.qualifiers['product'][0] + '[' + str(CDS_position_start) + '..' + str(CDS_position_end) + ']\n'
                                                           else:
                                                                linha_csv_CDS.append('AL ' + str(cont) + '_' + xml_rgn_str.split('_')[1] + ':' + '(' + str(gene_length) + ' bp)' + feat.qualifiers['product'][0] + '[' + str(CDS_position_start) + '..' + str(CDS_position_end) + ']\n')                                                 
                                                    except KeyError:
                                                        print ('KeyError passed')
                                                cont_CDS= cont_CDS + 1
                                        cont = cont +1
            with open ('Demo_Regions.csv', 'a', encoding='UTF8', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    #print (linha_csv_Region)
                    #print (linha_csv_CDS)
                    #input('prescrita')
                    writer.writerow(linha_csv_Region)
                    writer.writerow(linha_csv_CDS)
                    print('escreveu')























