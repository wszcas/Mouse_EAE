#!/usr/bin/python

import os,sys

input = sys.argv[1]
outdir = sys.argv[2]

n3 = input.replace('.fa','.blastout.txt')



print '\n\n','input file is\n',n3,'\n\n'
print '''###############################################################################################################'''

if n3 not in os.listdir(outdir):
		
	output = 'outdir/%s'%n3
	print 'will perform blast for:\n',output,'\n\n'

##	command = " ../bin/igblastn -germline_db_V imgt_database/mouse_IGHV_new -germline_db_J imgt_database/mouse_IGHJ_new -germline_db_D imgt_database/mouse_IGHD_new -organism mouse -domain_system imgt -query %s -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3 -out %s"%(input,output)
	command = " ../bin/igblastn -germline_db_V imgt_database_818knockinonly/mouse_IGHV_withKnockIn -germline_db_J imgt_database_818knockinonly/mouse_IGHJ_new -germline_db_D imgt_database_818knockinonly/mouse_IGHD_new -organism mouse -domain_system imgt -query %s -auxiliary_data optional_file/mouse_gl.aux -show_translation -outfmt 3 -out %s"%(input,output)

	print command
	print ' finishing writing blasting results to blast_out folder for file:\n',output,'\n\n'

	os.system(command)
else:
	print '''###############################################################################################################'''
	print 'files already exists for:\n', n3,'\n'
	print 'no blast will be performed for:\n',n3
	print '''###############################################################################################################'''
