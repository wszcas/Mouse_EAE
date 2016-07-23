#!/usr/bin/python

import os,sys


import MySQLdb
import os,sys,re
import string,csv

from Bio import *
from Bio.SeqUtils import *




# build the dictorinaruy storing the substition error count, initalize as 0
dic_error = {k:0 for k in  range(291)}
# build the dictorinaruy storing the substition error count and the length of the reads alignmeng to plasmid
dic_length_and_error = {}


def process_header(header,nseq):
    print 'header is:', header
    
    id,read_length,mutation,mut_loci = header.strip().split(";")
    mut_loci=mut_loci.replace('\n','')
# 	if cdr3 == '':
# 		cdr3 = cdr3_vdj
# 	elif cdr3 in cdr3_vdj:
# 		cdr3 = cdr3_vdj
# 	elif cdr3[1:] in cdr3_vdj:
# 		cdr3 = cdr3_vdj


    id = id.replace('>','')
    mutation = int(mutation)
    return [id,mutation,read_length,mut_loci,nseq]

def processFile(name):


	## initiate a dictionary to store the mutation location and count, for example {'-20':['IgG',5],'-100':['IgG',6] }
    mutation_loci= {}

    lines = file(name).readlines()

    for i in range(0,len(lines)-1,2):


        header = lines[i].strip()
        nseq = lines[i+1].strip()

        if nseq=='' or not (nseq[0] in 'ACTG'):
            continue
        result = process_header(header,nseq)
        insertion_seq = ','.join(["'%s'"%e for e in result])
        print result
        print len(result)
        print insertion_seq
		#sql= "INSERT INTO   EAE_hybridoma (sequence_id,ighv_string,ighj_string,ighd_string,MUTATION,h1_seq,h2_seq,h3_seq,nseq,aasequence,barcode) VALUES (%s)" %insertion_seq
		#c.execute(sql)
		
		
        ##recorde the mutation position
        mut_loci = (header.replace('\n','')).split(';')[-1]
        print mut_loci
        ## fill in the muation dictionary to store the mutation location and count, for example {'-20':['IgG',5],'-100':['IgG',6] }

        if mut_loci != '':
            mut = mut_loci.split('_')
            for m in mut:
                if m not in mutation_loci:
                    mutation_loci[m] = 1
                else:
                    mutation_loci[m] +=1
        
        ## record read length and mutation number
        id,read_length,mutation,t1 = header.strip().split(";")
        id = id.replace('>','')
        
        length = int(read_length)
        mutation = int(mutation)
        if id not in dic_length_and_error.keys():
            dic_length_and_error[id] = [mutation,length]
        if mutation in dic_error.keys():
            dic_error[mutation] +=1
        elif mutation not in dic_error.keys():
            dic_error[mutation] =1

    


    out = 'mutation_position_location_freq.txt'
    f = open(out,'w')

    for key in mutation_loci.keys():
        f.write(key+'\t'+ str(mutation_loci[key])+ '\n')
    f.close()
    
    

input = sys.argv[1]	
processFile(input)



c.close()




print dic_error


out1= open('mutation_count.csv','w')

for mut, count in dic_error.items():    
    out1.write(str(mut) + ',' + str(count) +'\n')
out1.close()

out2= open('mutation_count_and_length.csv','w')

for id, record in dic_length_and_error.items():
    out2.write(id + ','+ str(record[0]) + ',' + str(record[1]) + ','+'\n')
out2.close()

     
            
        
        
          

