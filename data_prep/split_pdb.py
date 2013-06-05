#!/opt/local/bin/python2.7
#split_pdb.py
# Splits pdb into its different chain, aka ligand and receptors
#=============== MODULES ======================
##{{{
#-----------------  Suppress stdout --------------------
import sys
from numpy import *
import numpy as np
import scipy
import scipy.stats as ss
import os
import pprint
import math
import re
import time
import atexit
import argparse
##}}}
#==============================================

###############################################
####       Setting Up Parser
###############################################
##{{{
#=====================================
#        Handling Error
#=====================================
##{{{
class MyParser(argparse.ArgumentParser):
    def error(self, message):
		if re.search(r'\.py$',sys.argv[0]): # This mean from shell
			sys.stderr.write('error: %s\n' % message)
			self.print_help()
			sys.exit(2)
		else:
			raise Exception('error: %s\n' % message)
##}}}
#=====================================
#        Setting up
#=====================================
##{{{
parser=MyParser(description='Splitting pdb into separate chains.',usage='%(prog)s 1ycr.pdb -o rec.pdb lig.pdb')
parser.add_argument('input',type=argparse.FileType('r'))
parser.add_argument('-o','--output','--out',nargs='*',default=None,help='If not defined, creates ***_1.pdb, ***_2.pdb etc.')
parser.add_argument('-n',action='store_true',help='Checks for number of chains in input file without output.')
parser.add_argument('-cp',action='store_true',help='Creates different combinations of proteins.')
parser.add_argument('--ow',action='store_true',help='Auto overwrite or not(default).')
p=parser.parse_args()
f=p.input # Assign file
##}}}
##}}}
###########################################
####	  Saving Protein Data
###########################################
##{{{
#====== Extract
print 'Extracting lines from pdb...'
prot=[]
prot_info=[]
chain=[]
for line in f:
	if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
		#Check for chain
		if line.startswith('ATOM'):
			try: #Checks if chain column contains a letter for chain
				chainn=int(line.split()[4])
				raise Exception('No chain letter is assigned.')
			except ValueError:
				chain.append(line.split()[4])
		#Check for chain termination
		try:
			if line.startswith('ATOM') and chain[-1]!=chain[-2]:
				prot.append('TER')
				prot.append(line)
			else:
				prot.append(line)
		except IndexError:
			prot.append(line)
	else:
		prot_info.append(line)

#====== Reporting number of chains
print 'Calculating number of chains...'
num_of_chain=len(np.unique(chain))		
num=num_of_chain
#print num_of_chain
#num=1
#
#for ind,line in enumerate(prot):
#	if line.startswith('TER') and prot[ind+1].startswith('ATOM'):
#		num += 1

if p.n: # Calls to check number of chain
	print 'There are {0:0.0f} separate chains.'.format(num)
	if p.output==None: # If only -n, terminate prematurely
		sys.exit()

#===================================
#     		Saving
#===================================
##{{{
print 'Saving...'
#------- Split input Name
input_name=f.name.split('.') # Splits 1ycr.pdb into ['1ycr','pdb'] for example


#------- Create Filenames
if p.output==None and not p.cp:
	# create filenames eg.1ycr_1.pdb, 1ycr_2.pdb. Name comes from input file
	filenames=[''.join([input_name[0],'_',str(j+1),'.',input_name[1]]) for j in xrange(num)] 
	test=raw_input('Create a default list of proteins '+str(filenames)+'?(y/n)\n')
	if test=='y':
		pass
	elif test=='n':
		raise Exception('Please specify output files or choose -cp.')
	else:
		raise Exception('(Jon)Wrong option selected.')
elif p.output==None and p.cp:
	#For proteins with more than 2 chains, create permutations and combinations
	lis1,lis2=Jon.ntom(range(num))
	pp=[[i,j] for ii,i in enumerate(lis1) for ji,j in enumerate(lis2) if ii==ji]
	# Sorts it so that the 0 is always in the first element, aka receptor
	for ind,i in enumerate(pp):
		if 0 in i[1]:
			pp[ind]=[i[1],i[0]]
		pp[ind][0]=list(pp[ind][0])
		pp[ind][1]=list(pp[ind][1])
		pp[ind][0].sort()
		pp[ind][1].sort()
	filenames=[]
else:
	if len(p.output)!=num:
		sys.stderr.write("ERROR: {0:0d} targets does not match {1:0d} chains in pdb\n".format(len(p.output),num))
		raise Exception("ERROR: {0:0d} targets does not match {1:0d} chains in pdb\n".format(len(p.output),num))
	elif len(p.output)==num:
		filenames=p.output


#------- Check filenames for existence
for i in filenames:
	if os.path.isfile(i):
		if p.ow:
			print i+' exists. Overwriting...'
		else:
			check=raw_input('File {0:s} exist. Overwrite? (y/n) \n'.format(i))
			if check=='y':
				print 'Overwriting...'
			elif check=='n':
				raise Exception('(Jon) File exist, please choose another name.')
			else:
				raise Exception('(Jon) Wrong option')


#------- Save files
if p.cp:
	##{{{
	# For permutation and combination only
	for ind,i in enumerate(pp):
		# Iterate over different sets
		name1=''.join(['rec','_',str(ind),'_',''.join(map(str,i[0])),'.',input_name[1]])
		name2=''.join(['lig','_',str(ind),'_',''.join(map(str,i[1])),'.',input_name[1]])
		Jon.fexist(name1,err='e',option='w')
		Jon.fexist(name2,err='e',option='w')
		check=0 #Check which sub protein we are at
		for ind,line in enumerate(prot):
			# Iterate over lines
			if check in i[0]:
				file_save=open('temp1','a')
				if line.startswith('ATOM'):
					file_save.write(line)
				elif line.startswith('TER'):
					file_save.write('TER')
					file_save.write('\n')
					file_save.close()
					check+=1
				elif line.startswith('END'):
					file_save.write('TER')
					file_save.close()
					check+=1
			elif check in i[1]:
				file_save=open('temp2','a')
				if line.startswith('ATOM'):
					file_save.write(line)
				elif line.startswith('TER'):
					file_save.write('TER')
					file_save.write('\n')
					file_save.close()
					check+=1
				elif line.startswith('END'):
					file_save.write('TER')
					file_save.close()
					check+=1
		# Removes TER at the end and change it to END
		with open(name1,'w') as f:
			ff=open('temp1','r')
			lines=ff.readlines()
			for i in lines[:-1]:
				f.write(i)
			f.write('END')
			ff.close()
		with open(name2,'w') as f:
			ff=open('temp2','r')
			lines=ff.readlines()
			for i in lines[:-1]:
				f.write(i)
			f.write('END')
			ff.close()
		os.remove('temp1')
		os.remove('temp2')
		##}}}
else:
	#for all other options
	filenames.reverse() # Facilitate pop
	file_save=open(filenames.pop(),'w')
	for ind,line in enumerate(prot):
		if line.startswith('ATOM'):
			file_save.write(line)
		elif line.startswith('TER') and prot[ind+1].startswith('ATOM'):
			file_save.write('END')
			file_save.close()
			file_save=open(filenames.pop(),'w') # Save into next file
		elif line.startswith('END'):
			file_save.write(line)
			file_save.close()

##}}}
##}}}

		
		
		

		





