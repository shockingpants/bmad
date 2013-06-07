#!/usr/bin/env python
#!convert_b.py

import numpy as np
import re
import sys
import argparse
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
parser=MyParser(description=\
'Edits the B factor column of the file.\n \
 Assumes that b factor is always the second last column of the pdb file.\n \
 Works regardless of whether there is a chain identifier. \n \
 Outputs a pdb file.\n \
'\
,usage='%(prog)s prot.pdb -d data.dat -o prot_b.pdb')
parser.add_argument('input')
parser.add_argument('-d','--data',help='Please input data file.')
parser.add_argument('-o','--output','--out',help='This creates an output pdb file')
#parser.add_argument('-b','--backup',action='store_true',help='Backs up b factor into file b_factor.BAK')
#parser.add_argument('--ow',action='store_true',help='Auto overwrite or not(default).')
p=parser.parse_args()
pdb=p.input # Assign file
dat_in=p.data
##}}}
##}}}

##############################
###  Checks data type
##############################
with open(dat_in,'r') as f:
	a=f.readline()
	if a[-1:]=='\n':
		a=a[:-1]# This removes '\n'
	a=a.split()
#-----Checking datafile type-----
length=len(a)
if length==1:
	try:
		int(a[0])
		names=0
	except ValueError:
		names=1	
elif length==2:
	try:
		int(a[1])
		names=0
	except ValueError:
		names=1
else:
	print a
	raise Exception('Wrong number of columns.({0:s})'.format(str(len(a))))
#----------------------------
if length==1:
	data=np.genfromtxt(dat_in, skip_header=names)
elif length==2:
	data=np.genfromtxt(dat_in, skip_header=names)[:,1]

##############################
### parse pdb and change stuff
##############################
resid=[] #Unique resid
atm_resid=[] #Resid of each atom
count=0 #Keeps track of number of different residues
atm_count=0 #Keeps track of number of atom lines
prot=[]
with open(pdb,'r') as f:
	for ind,i in enumerate(f):
		if i.startswith('ATOM'):
			atm_count+=1
			prot.append(i)
			try: # Checks to see if chain id is present and pulls out the resid
				res=int(i.split()[5]) #Assumes no chain id, which is usually not the case
				chain=False #Means there is no chain ID
			except ValueError:
				res=int(i.split()[6])
				chain=True
			atm_resid.append(res)
			try: #Adds resid to a list
				if res!=resid[-1]:
					count+=1
					resid.append(res)
			except IndexError: #Accounts for the first residue
				count+=1 
				resid.append(res)

###############################
### Prepare Data
###############################
#Check data length against number of res
#print len(data)
#print data
#print count
if len(data)!=count:
	raise Exception('Data length({0:d}) and number of residues({1:d}) are not the same'.\
	format(len(data),count))


#Generate data based on atom index
data2=[]
which_data=0
for ind,i in enumerate(atm_resid):
	if i!=atm_resid[ind-1] and ind!=0:
			which_data+=1
	data2.append(data[which_data])

#Check data length against number of atm
assert len(data2)==atm_count

###############################
### Apply data to pdb and save
###############################

with open(p.output,'w') as f:
	for ind,i in enumerate(prot):
		if i.startswith('ATOM'):
			p=re.split('(\s+)',i)
			p[-5]='{0:4.2f}'.format(data2[ind]) #B factor field happen to be the 5th last column
			i=''.join(p)
			while len(''.join(p[:-2]))!=78:
				if len(''.join(p[:-2]))>78:
					p[-4]=p[-4][:-1]
				if len(''.join(p[:-2]))<78:
					p[-4]=' '+p[-4]
			i=''.join(p)
		f.write(i)
