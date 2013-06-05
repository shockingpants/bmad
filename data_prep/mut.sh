#!/bin/bash
# mut.sh
# Mutates target protein
# Selects the most probably rotamer. EDit accordingly

usage()
{
cat << EOF
usage eg: $0 -i 1ycr.pdb -o prot_mut.pdb -n 24,25 -t TYR

This script mutates target residue in a pdb file.

OPTIONS:
   -h      Show this message
   -i      Input file
   -o      Output file
   -n      Resi number to change
   -t      Target AA to change into
EOF
}
input=
output=
resn=
rest=
while getopts “hi:o:n:t:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             input=$OPTARG
             ;;
         o)
             output=$OPTARG
             ;;
         n)
             resn=$OPTARG
             ;;
         t)
             rest=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $input ]] || [[ -z $output ]] || [[ -z $resn ]] || [[ -z $rest ]]
then
     echo Please use up all the options.
     usage
     exit 1
fi

# Checks for multiple resn to be mutated
if [[ -n $(echo $resn|grep ',') ]]; then 
	IFS=','
	fi

# Writes a script that pymol can understand and run pymol to create mutant. 
cat << _EOF > mut.py
cmd.load('`echo $input`', 'prot')
cmd.wizard("mutagenesis")
_EOF

for i in ${resn[@]}; do
# Iterate over every number
cat << _EOF >> mut.py
print '`echo $i`'
cmd.do("refresh_wizard")
cmd.get_wizard().set_mode('`echo $rest`') # This sets the target AAA we want to change into
cmd.select('sele', 'prot and resi `echo $i`')
cmd.select('sele2', 'prot and resi `echo $i` and name CA')
hihi=cmd.count_atoms('sele2')
if hihi!=1:
	print 'Please select resn via another criteria. There are duplicate resn %s in diff chains.' %(hihi) 
	exit
cmd.get_wizard().do_select("sele")
cmd.frame(1)
cmd.get_wizard().apply()
_EOF
done

cat << _EOF >> mut.py
cmd.set_wizard("done")
cmd.save('`echo $output`', 'prot')
_EOF

/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL mut.py -cq

rm mut.py




