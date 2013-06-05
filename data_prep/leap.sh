#!/bin/bash
# leap.sh
# To produce .inpcrd and .parmtop of input .pdb file (Topology and Coordinate)
##########################
####  Get input file
##########################
##{{{
usage()
{
cat << EOF
usage: $0 -i prot.pdb -f leaprc.ff99SB
This produces .inpcrd and .parmtop from input .pdb file at current location.
It uses pdb file as is with no addition or removal or water.

OPTIONS:
   -h      Show this message
   -i      Target protein
   -f      Forcefield
EOF
}

while getopts “h:i:f” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             file=$OPTARG
             ;;
         f)
             FORCEFIELD=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $file ]]; then
     usage
     exit 1
fi

##}}}

###############################
####  Generate Parmtop and Inpcrd
###############################

file1=$(basename "$file")
extension="${file1##*.}"
filename="${file1%.*}"

#---- Check extension
if [[ $extension != "pdb" ]]; then
	echo must be .pdb file
	exit 1
fi

#---- Creating input file for tleap

if [[ -z $FORCEFIELD ]]; then
     FORCEFIELD=leaprc.ff99SB
fi


cat << _EOF > in_leap
source $FORCEFIELD
loadamberparams gaff.dat
loadamberprep AAA.prep
loadamberparams AAA.parm
$filename = loadpdb $file1
set default PBradii mbondi2
saveamberparm $filename $filename.parmtop $filename.inpcrd
quit
_EOF

tleap -f in_leap

#----- Check if files got output successfully
if [[ ! -f $filename.parmtop ]] || [[ ! -f $filename.inpcrd ]]; then
     echo "Unsuccessful. parmtop and inpcrd not produced"
     exit 1
fi

