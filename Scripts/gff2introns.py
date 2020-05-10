#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/bin/python
# generate intron and exon annotations from gff file

import getopt, re, os, sys, logging, time, datetime, copy
from collections import defaultdict



# gff file
options, args = getopt.getopt(sys.argv[1:], 'g:',['GFF='])
for opt, arg in options:
    if opt in ('--GFF'):
        gff = arg

        
        
# generate Exon.gff file
cmd = "less %s | awk '{if($3==\"exon\"){print $0}}' > Exon_%s" % (gff, gff)
os.system(cmd)




# load files
def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

    
    
    
# gff to introns
def gff2introns(gffFiles):
    exon_s, exon_e = defaultdict(lambda: []), defaultdict(lambda: [])
    parse_line_start = 0
    for i in range(0, len(gffFiles)):
        if gffFiles[i].startswith('#') or gffFiles[i].split('\t')[2] != 'exon':
            continue
        parse_line_start = i
        break
    for i in range(parse_line_start, len(gffFiles)):
        gff_list = gffFiles[i].strip().split('\t')
        if gff_list[2] != 'exon':
            continue
        transcript_id = re.sub('.*transcript_id |;.*', '', gff_list[8])
        exon_s[transcript_id].append(int(gff_list[3]))
        exon_e[transcript_id].append(int(gff_list[4]))
    trans = {}
    fw = open(("Intron_%s" % (gff)), "w")
    for i in range(parse_line_start, len(gffFiles)):
        gff_list = gffFiles[i].strip().split('\t')
        transcript_id = re.sub('.*transcript_id |;.*', '', gff_list[8])
        gene_id = re.sub('.*gene_id |;.*', '', gff_list[8])
        if transcript_id in trans:
            trans[transcript_id] = 'false'
        elif len(exon_s[transcript_id]) > 1:
            exon_s[transcript_id].sort()
            exon_e[transcript_id].sort()
            for i in range(1, len(exon_s[transcript_id])):
                if gff_list[6] == '+':
                    fw.write(
                        "%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id %s; transcript_id %s; intron_number %s; total_intron_number %s; exon_1 %s_%s; exon_2 %s_%s;\n" % (
                            gff_list[0], gff_list[1], str(int(exon_e[transcript_id][i - 1]) + 1),
                            str(int(exon_s[transcript_id][i]) - 1),
                            gff_list[6], gene_id, transcript_id, i, len(exon_s[transcript_id]) - 1,
                            exon_s[transcript_id][i - 1], exon_e[transcript_id][i - 1], exon_s[transcript_id][i],
                            exon_e[transcript_id][i]))
                else:
                    fw.write(
                        "%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id %s; transcript_id %s; intron_number %s; total_intron_number %s; exon_1 %s_%s; exon_2 %s_%s;\n" % (
                            gff_list[0], gff_list[1], str(int(exon_e[transcript_id][i - 1]) + 1),
                            str(int(exon_s[transcript_id][i]) - 1),
                            gff_list[6], gene_id, transcript_id, len(exon_s[transcript_id]) - i,
                            len(exon_s[transcript_id]) - 1, exon_s[transcript_id][i - 1], exon_e[transcript_id][i - 1],
                            exon_s[transcript_id][i], exon_e[transcript_id][i]))
        trans[transcript_id] = 'true'
    fw.close()


    
# run
gffFiles = read_file('%s' % (gff))
gff2introns(gffFiles)

