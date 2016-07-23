#!/usr/bin/python
from Bio import *
from Bio.SeqUtils import *

# from Indel_Correction import *
# from corrFunctions import *
# from translate_3frame import *
import os,sys



dct = {"TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser", 
"TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp", 
"TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu", 
"CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro", 
"CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg", 
"CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met", 
"ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn", 
"AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg", 
"GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala", 
"GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu", 
"GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}



def translate(seq):
    x = 0
    aaseq = []
    while True:
        try:
            aaseq.append(dct[seq[x:x+3]])
            x += 3
        except (IndexError, KeyError):
            break
    return aaseq

# seq = "TTTCAATACTAGCATGACCAAAGTGGGAACCCCCTTACGTAGCATGACCCATATATATATATATA"
# seq=nseq_output


def translation_3frame(dnaseq):
    aaseq=[]
    for frame in range(3):
        aaseq.append(''.join(item.split('|')[0] for item in translate(dnaseq[frame:])))
    return aaseq

def hammingDistance(a, b):
    # function to caculate the hamming distance between query a and b, notice a is shorter or equal than b
    distance = 0
    for i in xrange(min(len(a),len(b))):
        if a[i]!=b[i]:
            distance +=1
    return distance    
    

def correct_all_Indels(nseq,gm_seq,cdr3_start):
    #deletion_no = 0 
    #correct indels and record the postion of substitutions
    mutation=[]
    result=[char for char in nseq]
    for i in range(len(gm_seq)):

        #deletion
        if nseq[i] == '-':
            if gm_seq[i] != '-':
                result[i] = gm_seq[i]
                #if i < cdr1_start:
                  #  deletion_no += 1
                    
        #insertion
        if gm_seq[i] == '-':
            if nseq[i] != '-':
                result[i] = '#'
        
        ##substitution
        if gm_seq[i]!=nseq[i] and gm_seq[i] in ['A','T','C','G'] and nseq[i] in ['A','T','C','G']:
            mutation.append(i)
            
    
    
    return (''.join(result),[str(i-int(cdr3_start)) for i in mutation])




# def correct_all_Indels(nseq,gm_seq,cdr1_start):
#     cdr1_deletion_no = 0
# #     cdr3_deletion_no = 0
#     result=[char for char in nseq]
#     for i in range(len(gm_seq)):
# 
#         #deletion
#         if nseq[i] == '-':
#             if gm_seq[i] != '-':
#                 result[i] = gm_seq[i]
# #                 cdr3_deletion_no+=1
#                 if i < cdr1_start:
#                     cdr1_deletion_no += 1
#             elif gm_seq[i] == '-':
#                 result[i] = '-'
# #                 cdr3_deletion_no+=1
#                 if i < cdr1_start:
#                     cdr1_deletion_no += 1    
#         #insertion
#         if gm_seq[i] == '-':
#             if nseq[i] != '-':
#                 result[i] = '#'
#     
#     return (''.join(result),cdr1_deletion_no)
# #     return (''.join(result))

def find_aa(dnaseq,gm_seq_aa):
    
    aa_list=[]
    for aa in translation_3frame(dnaseq):
        aa_list.append(aa)
    # get all the possbile translations 
    hamming={}
    for key in aa_list:
        hamming[hammingDistance(gm_seq_aa,key)]=key
        
    #return the one translation with the minial hamming distance compared with germline aa seq 
    return hamming[min(hamming.keys())]
    
def find_real_cdr3_start(nseq,nseq_vgene):
    deletions=0
    
    for i,c in enumerate(nseq_vgene):
        
        j=i+deletions
        if c==nseq[j]:
            print nseq_vgene[i],nseq[j]
            
            continue
        else:
            for k,c2 in enumerate(nseq[j:]):
                if c2==c :
                    deletions +=k
                    break
            
        
    return deletions
    
            
            
# for i,c in enumerate('GAATCT'):
#     print i,c


def process_xGene(xString):

    xList = xString.split(',')
    xgenes = []
    for xgene in xList:
        xgene = xgene.split('*')[0]

        if xgene not in xgenes:
            xgenes.append(xgene)
    ighx = ' '.join(xgenes)

    return ighx

    
#remove indel place holders from sequence to generate output
def remove_holders(seq):

    seq = seq.replace('-','')
    seq = seq.replace('#','')
    return seq


def find_last_occurance(str,query):
    ''' find all occurance of of query in str, and return the end position of last occurance,
    , if not found, then return -1'''
    index = -1
    length= len(query)
    result_list = [n for n,c in enumerate(str) if str[n:n+length]==query]
    if len(result_list)>0:
        index= length + max (result_list)
    return index


def process_file(name):
    inName = "blast_out_A8/%s"%name
#     out = file("extracted_corrected/%s.fa"%name,'w')
    out = file("extract_correct_all_indels_noLocalizedFramShift_A8_new/%s_jan30.fa"%name,'w')
   # out2 = file("extract_human_light_chain3/%s_newCDR1_withaa.fa"%name,'w')

#     inName='redotest'
# 	out2 = file("extracted_locFrame/%s.fa"%name,'w')
#     inName = "kappa_blast_out_test"
#     out = file("%s_errot_test_output.fa"%inName,'w')
    # out2 = file("extracted_locFrame/%s.fa"%name,'w')
##    skipped = 0
    state = 0
    header=''
    
    skipped = 0
    passed = 0
    f= open(inName,'r')
        
    for line in f:
#         break
#         print line
        print 'state is :', state
    
        if state == 0:
            cdr_dict = {}
            cdr_dict2 = {}
            mark = ''
            nseq = ''
            gm_seq = ''
            gm_seq_aa = ''
            cdr3_start = -1
            nseq_vgene = ''
            header =''
            id = ''
            cdr3_vdj = ''
            
                
            if line.find('Query=') > -1:
                header=''
                while line.find('Length')==-1:
                    header=header+line
                    line=f.next()
                header = header.split('Query=')[1].strip()
                print header
            #print line
                    # to grab all the header information which might go to next line marked by \n
                    
                id,cdr3_vdj=header.split(';')
                cdr3_vdj = cdr3_vdj.replace('\n','')
                if cdr3_vdj.find('X')==-1:   ##skip the reads with STOP codon in the CDR3
                    state = 1
                    

                
                   
                
    
        elif state == 1:
            if line.find('(Top V gene match, Top D gene match, Top J gene match') > -1: 
                state = 2
            elif line.find('(Top V gene match, Top J gene match') > -1:
                state = 0
                   
    
        elif state == 2:
    
            list = line.strip().split("\t")
            print list,len(list)
            print list[0]
            print list[1]
            print list[2]
            if list[0].find('IGHV') == -1:
                state = 0
                continue
            
            elif list[0].find('IGHV')> -1:
                vgene,dgene,jgene = process_xGene(list[0]), process_xGene(list[1]),process_xGene(list[2])
                
                print vgene,"\n",dgene,"\n",jgene
                productive = list[-2]
                print productive
                state = 3
#             else:
#                 state = 0
#                 continue
            
            
                
            
    
        elif state == 3:
    
            if line.find('Alignment summary') > -1:
#                 state = 4
                state = 5
#             else:
#                 state = 0
            
    # get the all sequence, so skip the state=4
#         elif state == 4:
#             #only keep the reads with complete CDR1,which means it must have FR1
#             if line.find('FR1') > -1:  
#                 state =5
#             else:
#                 print 'not a complete sequence, go to next entry '
#                 state = 0
        elif state == 5:
            #Extract CDR
            if line.find('CDR') > -1:
                    
                cdr_list = line.split("\t")
                cdr,start,end = cdr_list[:3]
                cdr = cdr.split('-IMGT')[0]
                tup = (int(start)-1,int(end))
                cdr_dict[cdr] = tup
                print cdr_dict
                
            
            #CDR3 isn't alway present. Find the FR3 enpoint in table to determine cdr3 start-point
            elif line.find('FR3') > - 1:				
    
                end =  line.split("\t")[2]
                cdr3_start  = int(end)
    
            #End of sequence stats. Grab mutation count and shift to next state				
            elif line.find('Total') > -1:
                
                mutation = line.split("\t")[-3]
                mutation = int(mutation)
                print mutation
                germline_bool = False
                germline_aa_bool=False
                # make sure the reads have spanning the CDR1 to CDR3 region
#                 if 'cdr1' in cdr_dict and cdr3_start> -1:
                state = 6
    
        elif state == 6:
   
                
                
            #extract input seq subsequence
            if line.find('Query_') > -1:
    
                try:
                    label,start,seq,end = [e for e in line.strip().split(" ") if len(e) > 0]
                except:
                    continue
    
                if len(nseq) == 0:
                    
                    align_start = int(start) - 1
                    #start = align_start-1
                    nseq = align_start*' '+seq
                    
                else:
                    nseq  += seq		
                
                #Entering germline mode. Search for first instance of germline subsequence
                germline_bool = True
    #             print nseq
    
            elif germline_bool:
    
                #Definitely in germline territory. Time to exxtract germline alighnment sequence
    
                if line.find('IGHV') > -1:					
            
                    lineList = [e for e in line.strip().split(" ") if len(e) > 0]
                    seq = lineList[-2]
                        
                    if len(gm_seq) == 0:
                        gm_seq = align_start*' ' + seq
                    else:
                        gm_seq += seq
                    
    #                 print gm_seq
    
                    #Now we look for original sequence again. Exit Germline Mode
                    germline_bool = False
                    germline_aa_bool = True
            elif germline_aa_bool:
                if line.find('IG')== -1:
                    linelist=[e for e in line.strip().split(" ") if len(e) > 0]
                    gm_seq_aa+=''.join(linelist)
                    germline_aa_bool = False
            
                
    
            #We are done extracting sequence information. Now its time to process and output					
            if line.find('Lambda') > -1:
    
                state = 0
                print nseq
                print gm_seq
                #gm_seq = gm_seq[:cdr3_start]
                #Cutt-off sequence at cdr3 end
                nseq_tem = nseq.replace('-','')
                print nseq_tem
                 
                 ## first try find the V_end_query  from nseq using 5bp without dels in between 
                 ## make sure the reads have spanning the CDR1 to CDR3 region
                print 'cdr3 start position is: ', cdr3_start
#                 if cdr3_start ==-1 or not ('CDR1' in cdr_dict):
                if cdr3_start ==-1:
                    # if no cdr3 present, skip the entry
                    state = 0 
                    continue
                elif cdr3_start > -1:
                    nseq_vgene = nseq_tem[:cdr3_start]
                    print nseq_vgene
                    #find the number of shift due to deletions by calling find_reak_cdr3_start function by comparing nseq and nseq_vgene
                    shift= find_real_cdr3_start(nseq,nseq_vgene)
                    print 'cdr3 shift is', shift
                    gm_seq = gm_seq[:cdr3_start+shift]
                    print gm_seq
                    

#                     if 'CDR1' in cdr_dict:
#                         cdr1_start= cdr_dict['CDR1'][0]
#                         
#                     else:
#                         cdr1_start = -1
#                     print 'cdr1_start is',cdr1_start
#                     nseq,cdr1_shift = correct_all_Indels(nseq,gm_seq,cdr1_start)   
                    #record the location of substition postions
                    nseq,mutation_loc  = correct_all_Indels(nseq,gm_seq,cdr3_start+shift )
                         
                    ## find the true CDR1 start location
#                     CDR1_start_query = nseq_tem[cdr1_start:cdr1_start]
#                     if find_last_occurance(nseq,V_end_query) > -1 :
#                         V_end_index = find_last_occurance(nseq,V_end_query)
#                     else :
#                      # if not find, insert dels in between them and find all possible querys from nseq and return the largest one
#                      
#                      #generate all possible V end query sequence
#                         V_end_query_list = [V_end_query[0]+i*'-'+ V_end_query[1]+ j*'-' + V_end_query[2]+ q*'-' + V_end_query[3]+\
#                          p*'-' + V_end_query[4] for i in range(9) for j in range(9) for p in range(9) for q in range(9) ]
#                         try: 
#                             V_end_index = max([find_last_occurance(nseq,element) for element in V_end_query_list])
#                         except:
#                              continue    
                        

    #                 nseq = correct_all_Indels(nseq,gm_seq)
    #                 shift=nseq[:cdr3_start+shift].count('-')+shift
                    print 'nseq after indel correction raw seq is:', nseq
                    print 'subsequence of nseq of V gene part is:', nseq[:cdr3_start+shift]
    #                     print 'shift number is: ', shift
#                     print 'cdr1 shift is:',cdr1_shift
                    print 'cdr3 shift is:',shift
                    nseq_output = remove_holders(nseq.strip())
                    print 'nseq after indel correction is:',nseq
                    print 'fixed nseq is:',nseq_output
                    aaSeq = find_aa(nseq_output,gm_seq_aa)
                    print aaSeq
                    if nseq_output.find('.')> -1:
                        state = 0
                        skipped +=1
                        continue
                
                    if len(aaSeq)>0:
                        
    
        
                        for key in cdr_dict.keys():
            
    ##                        if key.find('CDR3') > -1:
    ##                            continue
            
                            start,end = cdr_dict[key]
#                             if key.find('CDR1') > -1:
#                                 start=start+cdr1_shift
#                                 end=end+cdr1_shift-1
                            
                            if key.find('CDR3') > -1:
                                start=start+shift
                                end=end+shift
    ##        
                            #the sequence begins in the middle of the cdr. No good
    #                         if start == align_start:
    #                             continue
            
                            cdr = nseq[start:end+1]
                            print key,cdr
                            cdr = remove_holders(cdr)
                            aa_cdr = ''
                            trans = six_frame_translations(cdr)
                        
    
                            for i in range(1,4):
    #                             trans = translation_3frame(cdr)
                            
    #                         for i in range(0,3):
                                if len(aaSeq)==0:
                                    break
                                    
                                if aaSeq.find(trans[i]) > -1:
                                    remainder = len(cdr) % i
                                    
                                    #include more of nseq to compensate for left-overs
                                    if remainder > 0:
            
                                        cdr = remove_holders(nseq[start:end+remainder+1])
                                        trans = six_frame_translations(cdr)
                                                                
                                    aa_cdr = trans[i]
                                        
                                    break
                
                                #if aa_cdr.find('*') > -1:
                                #	continue
            
                            cdr_dict2[key] = aa_cdr
        
        
                    try:
                        cdr1 = cdr_dict2['CDR1']
                    except:
                        cdr1 = ''
                    
                    try:
                        cdr2 = cdr_dict2['CDR2']
                    except:
                        cdr2 = ''
                        
                    try:
                        cdr3 = cdr_dict2['CDR3']
                    except:
                        cdr3 = ''   
        
    #                 if len(cdr1) == 0 and len(cdr2) == 0 and len(cdr3) == 0:
    #                     continue
        
                            
                        #print "VDJ_Fasta Results : %s"%id
                    results = [id,vgene,dgene,jgene,mutation,cdr1,cdr2,cdr3,cdr3_vdj,len(gm_seq),'_'.join(mutation_loc)]
    ##                    results = [id,iglv,iglj,mutation,cdr1,cdr2,cdr3,cdr3_vdj]
                    #header = ';'.join(header.replace('\n','').split('_')[:-2]) + '_' + cdr1 + '_' +header.replace('\n','').split('_')[-1]
                    
                    #results = header+'_'+ str(mutation)
                    #print results
                    
                    print results
                    out.write('>' +  ";".join([str(e) for e in results]) + "\n")
                    out.write("%s\n"%nseq_output)
                    out.write("%s\n"%aaSeq)
                    print 'sequence wrote to file'
                    print 10*'**'
                    passed += 1
                    id = ''
                    vgene = ''
                    dgene = ''
                    jgene = ''
                    mutation = ''
                    cdr1 = ''
                    cdr2 = ''
                    cdr3 = ''
                    cdr3_vdj = ''
               
        
#                 out.write('>' + header + "\n" + nseq_output + '\n')
#                 out2.write('>' + header + "\n" + nseq_output + '\n'+ aaSeq + '\n')
	        #out3.write('>' + header.replace('\n','') + "\n" + nseq_output + '\n')


#                     elif triplet:
#                         out3.write('>' + results + "\n" + nseq_output + '\n'+ aaSeq + '\n')
#                         skipped+=1
#                         triplet = False
    # 
    #             if fixedLoc_count > 0:
    #                 out2.write('>' +  ";".join([str(e) for e in results]) + "\n")
    #                 out2.write("%s\n"%nseq_output)
    
    out.close()
    f.close()
    #
    #out3.close()
    print 'number of reads skipped is:',skipped
    print 'number of reads passed is:',passed

input = sys.argv[1]	
process_file(input)

