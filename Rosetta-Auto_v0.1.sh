#!/bin/bash	

#author: James Ross
#		 This is a multifunctional script
#		 It accepts (as standard) a single input pdb to initite the process with RELAX (you will need to set up a loop over a list if you want multiple input pdbs)
#		 It can except a list of pdbs for the MUTATE function, supposing this is a pre-relaxed list (essential)
#
#		 USAGE: This script is takes a single pdb and optional relaxes, mutates and calculates an intereaction energy, providing sequence clustering and metrics
#	     USAGE: Fill out the options below.
#		 EXAMPLE: just run the script!
#		 
#		 ADVICE
# 		 It is advised that you provide a 'clean' pdb as input. Solvent removed, no alternative conformations.  Ligands and metal ions may require further support and scripting.
#		 Further information for ligands and metals should be sought from the rosetta help pages
# 		 Have you added "/usr/local/rosetta/main/source/bin/" to your .bashrc i.e.
# 		 export PATH=$PATH:/usr/local/rosetta/main/source/bin/
#		 For CPU management, please make sure that RONST, RONSS and MONST are equal or greater than the number of CPU's, which they should be anyway!!
#
#		 Warnings : This script removes the folders 'rosetta-0-initial', 'rosetta-1-relax', 'rosetta-2-mutate', 'rosetta-4-inter', 'rosetta-5-analysis' from the current directory
#					This script removes all '.tx' files in the current directory
#					Makes use of 'nohup' and deletes nohup.out
#
# 	 	 IMPOVEMENTS TO MAKE:
#					EVERY PART OF IT
#				

########################################################################################################################################################################
# METHOD OPTIONS and OUTPUTS    #
########################################################################################################################################################################

# INPUT OPTIONS
# Option			# Expected			# Comments
CPU=4				# Interger			# Number of CPUs to use  

########################
#       RELAX          #
########################
# It is highly recommended that structures are relaxed prior to any analysis with rosetta
# If you are starting from structures that have not been relaxed in rosetta, you MUST initiate this step
# you should generate a number of stuctures and pick the best to carry forward with the following options
RELAX=True	 		# True False 		# Relax structure : relax into rosetta forcefield, creates output folder 'rosetta-1-relax'
RIPDB=1coi.pdb		# Filename False	# Input pdb file. (if False you can provide a list of files for mutagenesis with MIPDB)
RINST=4 			# Interger			# How many relaxations? : The number relaxations to make
RONSS=4 			# Interger 			# How many results to carry through? : The number of best relaxations to carry through to the next stage 
										#									   Typically do not exceed 10% of the total relaxations if making mutations
ROMET="N/A"			# True False 		# Produce output metrics and graphs for the relaxations. 				!!!!! INCOMPLETE !!!!!


########################
#       MUTATE         #
########################
# Mutagenesis of the protein structure using roestta fast relax and design.
MUTATE=True 		# True False 		# Mutate structure : Use roestta fast relax and design, creates output folder 'rosetta-2-mutate'
MIPDB=False  		# Filename False	# Input filelist : file with list of input pdbs, can contain a single pdb, 
 										#                  only use if no relaxation, otherwise "RONSS" determines input list.
										#				   i.e. 'False' is the default which takes inputs from the relaxation.
MIRES=False			# Filename False	# Res file input : Provide a rosetta res file with mutation options
										# 				   if False provide mutation details below (mires.tx), or provide threading input with MITHD
MITHD=thread.fa     	# Filename False 	# Threading a sequence to a structure : Provide a .pdb or .fa (fasta)
										#                                       Threading sequences MUST have the same number of residues as target structure, must have same numbering.
										#									    must be the same chain, only one chain at present
										#								  	    if Filename, provide threading details below (MITHC)
MONST=4 			# Interger			# How many mutation runs? : The number of fast-relax-design runs to make PER relaxed structure
MOMET="N/A"			# True False 		# Produce output metrics and graphs for the mutagenesis.  				!!!!! INCOMPLETE !!!!!
MOSEQ=True			# True False 		# Produce full output sequences (roestta-full-sequence.fa)
MOMSQ=True			# True False 		# Produce position specific sequence outputs dependent on resfile, required for clustering (roestta-spec-sequence.fa), single chain only
MOSQE=True			# True False 		# For each of the above sequence outputs, name - sequence - total energy

# if "MIRES=False" define mutational parameters here. Alternatively set MITHD=filename.file providing a threading sequence (.fa) or structure (.pdb)
# Provide chain and residue positions for mutagenesis, by default we will mutate to all amino acids except cysteine. Otherwise provide resfile.
# One chain per line, with residue positions, as below
rm mires.tx 2>/dev/null 
echo "
C 7 8 9 10 11
" > mires.tx

# if "MITHD=Filename" define threading parameters here.  
# residue numbers will be ignored, this replaces aminoacids sequenctially from source to target.
MITHC=C  			# Character False	# chain idenifier for replaced sequence. Not currently applicable for mutil chain replacements			

########################
#      INTERFACE       #
########################
# Analysis of the protein interface between chains.
INTER=True			# True False 		# Analyse interface : Use the InterfaceAnalyzer application to calculate ddG in Rosetta energy units (REU), creates output folder 'rosetta-4-inter'
IIFAC="A B"			# chain_names		# qoute chains of single group, if analaysing the interface between chain groups A and B vs C and D then use "A B"
IOSQE=True			# True False 		# For each of the analysed structures output, name - sequence - total energy
ROMET="N/A"			# True False 		# Produce output metrics and graphs for the Interface Energy.			!!!!! INCOMPLETE !!!!!
########################
#       CLUSTER        #
########################
# Sequence clustering and energy analysis.
# 	!!!!! INCOMPLETE !!!!!

########################################################################################################################################################################
# INTERNAL CHECKS AND FUNCTIONS    #
########################################################################################################################################################################

# find the current working directory
oridir=$(pwd)

########################
# PROGRAM DEPENDENCIES #
########################
# rosetta			# 
# pymol 			# This script uses pymol to detect a 10A shell of residues around mutated target residues, which are allow to 
					# move during the relaxation.  Other residues are fixed.
# python			# This script uses python3 to generate graphs for output sequences.

rm nohup.out 2>/dev/null 
nohup pymol -qc 
test=$(grep "Feature:" nohup.out | awk '{print $2}')
if [ $test = "PYMOL_MAIN" ] ; then 
	echo "PyMOL detected"
else	
	echo -en "\rPyMOL ***NOT*** detected"
	if [ MIRES = False ] ; then
		echo "no resfile detected and no PyMOL program detected"
		echo "I cannot make the resfile without PYMOL_MAIN"
		exit 0
	fi
fi


########################
#    OPTIONS CHECK     #head
########################

if [ $RELAX = True ] ; then
	if [ -n $RIPDB ] ; then 
		echo "read initial as pdb file $RIPDB for relaxations"
	else	
		echo "please check input file for relaxation, none found."
		exit 0
	fi
elif [ $MUTATE = True ] ; then 
	if [ -n $MIPDB ] ; then 
		echo "read initial as pdb list file $RIPDB with $(wc -l < $RIPDB) structures for mutations"
	else	
		echo "please check input file for mutations, none found."
		exit 0
	fi
elif [ $INTER = True ] ; then 
	if [ -n $IIPDB ] ; then 
		echo "read initial as pdb list file $IIPDB with $(wc -l < $IIPDB) structures for Interface analysis"
	else	
		echo "please check input file for interface analysis, none found."
		exit 0
	fi
fi
# check threading input.
if [ $MITHD != False ] ; then
	if [ $( echo $MITHD | awk -F'.' '{print $2}' ) = fa ] ; then 
		echo "detected threading with fasta file"
		echo "will thread first chain (line 2) of $MITHD onto chain $MITHC of $RIPDB "
	elif [ $( echo $MITHD | awk -F'.' '{print $2}' ) = pdb ] ; then 
		echo "detected threading with pdb file,"
		echo "will thread chain $MITHC of $MITHD onto chain $MITHC of $RIPDB"
	else
		echo "detected threading but file type not recognised, please use .pdb or .fa"
		exit 0
	fi
	if [ ! -f "$MITHD" ]; then
    echo "$MITHD does not exist, exiting . . ."
	exit 0
	fi
fi

# check res file clasehes.
if [ $MIRES != False ] && [ $MITHD != False ] ; then 
	echo 'We cannnot have both a res file ('$MIRES') and a theading template file ('$MITHD')'
	exit = 0 
fi

########################
#      transpose       #
########################

# pass first argument as $1 which should be any matrix
transpose ()
{
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $1 
}

########################
#       pdb2seq        #
########################

# pass first argument ($1) which should be a pdb file
# pass second argument ($2) which should be a chain character
seqgen ()
{
if [ -n "$2" ] ; then 
	grep 'CA' $1 | grep '^ATOM' | grep "[A-Z][A-Z][A-Z] $2 " | awk '{print $4}' > temp.tx
	sed 's/ALA/A/g; s/CYS/C/g; s/ASP/D/g; s/GLU/E/g; s/PHE/F/g; 
	     s/GLY/G/g; s/HIS/H/g; s/ILE/I/g; s/LYS/K/g; s/LEU/L/g; 
		 s/MET/M/g; s/ASN/N/g; s/PRO/P/g; s/GLN/Q/g; s/ARG/R/g; 
		 s/SER/S/g; s/THR/T/g; s/VAL/V/g; s/TRP/W/g; s/TYR/Y/g' < temp.tx  > transin.tx 
	transpose transin.tx | sed 's/ //g' > temp3.tx
	#sed 's/ //g'  transout.tx > temp3.tx
	echo ">"$1'_'$2 > head.tx
	cat head.tx temp3.tx > temp4.tx
	cat temp4.tx >> seqgen.fa
else
	for i in $(grep 'CA' $1 | grep '^ATOM' | cut -c 22 | sort | uniq) ; do 
		grep 'CA' $1 | grep '^ATOM' | grep "[A-Z][A-Z][A-Z] $i " | awk '{print $4}' > temp.tx
		sed 's/ALA/A/g; s/CYS/C/g; s/ASP/D/g; s/GLU/E/g; s/PHE/F/g; 
		     s/GLY/G/g; s/HIS/H/g; s/ILE/I/g; s/LYS/K/g; s/LEU/L/g; 
			 s/MET/M/g; s/ASN/N/g; s/PRO/P/g; s/GLN/Q/g; s/ARG/R/g; 
			 s/SER/S/g; s/THR/T/g; s/VAL/V/g; s/TRP/W/g; s/TYR/Y/g' < temp.tx  > transin.tx 
		transpose transin.tx | sed 's/ //g' > temp3.tx
		#sed 's/ //g'  transout.tx > temp3.tx
		echo ">"$1'_'$i > head.tx
		cat head.tx temp3.tx > temp4.tx
		cat temp4.tx >> seqgen.fa
	done
fi
}

########################
#    Sequence Matrix   #
########################

seqmatrix ()
{
# here we identify sequence similarities by a naive scoring function.
file=$1
echo "setting up . . ."
# check integrity of input
countnames=$(grep '>' $file | wc -l )
countsequences=$(grep -v '>' $file | wc -l )
if [ $countnames != $countsequences ] ; then
        echo "more than one chain per item, this code is for single chains"
        exit 0
fi

# check sequence lengths
clengthsequences=$(grep -v '>' $file | awk '{print length($0)}' | sort | uniq -c | wc -l )
if [ $clengthsequences != 1 ] ; then
        echo "sequences are of different lengths"
        exit 0
fi
lengthsequences=$( grep -v '>' $file | awk '{print length($0)}' | sort | uniq -c | awk '{print $2}')

# generate working files
echo name > names.tx
grep '>' $file >> names.tx
cp names.tx names2.tx
grep -v '>' $file | sed 's/\(.\{1\}\)/\1,/g' > sequence.tx

for count in $(seq 1 $(grep -v '>' sequence.tx | wc -l) ) ; do
        echo -en "\rprocessing $count"
        initial=$(sed -n ''$count'p' sequence.tx)
        ncount=$(( $count + 1 ))
        sed -n ''$ncount'p' names.tx > current.tx
        rm initial.tx
        for i in $(seq 1 $countsequences) ; do
                echo $initial >> initial.tx
        done

        awk '
BEGIN{
  FS=OFS=","
}
FNR==NR{
  for(i=1;i<=NF;i++){
    array[FNR,i]=$i
  }
  next
}
{
  for(i=1;i<=NF;i++){
    $i=array[FNR,i] $i
  }
}
1
' initial.tx  sequence.tx > acumu.tx

        # MERGE THE FILES
        sed -i 's/AA/0/g;s/CC/0/g;s/DD/0/g;s/EE/0/g;s/FF/0/g;
				s/GG/0/g;s/HH/0/g;s/II/0/g;s/KK/0/g;s/LL/0/g;
				s/MM/0/g;s/NN/0/g;s/PP/0/g;s/QQ/0/g;s/RR/0/g;
				s/SS/0/g;s/TT/0/g;s/VV/0/g;s/WW/0/g;s/YY/0/g' acumu.tx
        sed -i 's/AP/1/g;s/AS/1/g;s/AT/1/g;s/AG/1/g' acumu.tx
        sed -i 's/CV/1/g' acumu.tx
        sed -i 's/DN/1/g;s/DE/1/g' acumu.tx
        sed -i 's/EN/1/g;s/EQ/1/g;s/EH/1/g;s/ED/1/g' acumu.tx
        sed -i 's/FW/1/g;s/FY/1/g;s/FL/1/g;s/FI/1/g;s/FM/1/g' acumu.tx
        sed -i 's/GA/1/g;s/GP/1/g;s/GG/1/g' acumu.tx
        sed -i 's/HK/1/g;s/HR/1/g;s/HE/1/g;s/HQ/1/g' acumu.tx
        sed -i 's/IF/1/g;s/IY/1/g;s/IL/1/g;s/IM/1/g;s/IV/1/g' acumu.tx
        sed -i 's/KR/1/g;s/KQ/1/g;s/KH/1/g' acumu.tx
        sed -i 's/LF/1/g;s/LY/1/g;s/LV/1/g;s/LI/1/g;s/LM/1/g' acumu.tx
        sed -i 's/MF/1/g;s/MY/1/g;s/MV/1/g;s/MI/1/g;s/ML/1/g' acumu.tx
        sed -i 's/NE/1/g;s/ND/1/g;s/NQ/1/g' acumu.tx
        sed -i 's/PA/1/g;s/PS/1/g;s/PT/1/g;s/PG/1/g' acumu.tx
        sed -i 's/QK/1/g;s/QE/1/g;s/QH/1/g;s/QN/1/g' acumu.tx
        sed -i 's/RK/1/g;s/RH/1/g' acumu.tx
        sed -i 's/SA/1/g;s/SP/1/g;s/ST/1/g;s/SG/1/g' acumu.tx
        sed -i 's/TA/1/g;s/TP/1/g;s/TS/1/g' acumu.tx
        sed -i 's/VL/1/g;s/VI/1/g;s/VM/1/g;s/VC/1/g' acumu.tx
        sed -i 's/WF/1/g;s/WY/1/g' acumu.tx
        sed -i 's/YW/1/g;s/YF/1/g;s/YL/1/g;s/YI/1/g;s/YM/1/g' acumu.tx
        sed -i 's/[A-Z][A-Z]/2/g' acumu.tx

        awk -F',' '{ for(i=1; i<=NF;i++) j+=$i; print j; j=0 }' acumu.tx > current2.tx
        cat current.tx current2.tx > current3.tx
        paste -d',' names2.tx current3.tx > temp.tx
        mv temp.tx names2.tx
done
mv names2.tx matrix.csv
}
########################
#   Cluster function   #
########################.

echo "
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from pylab import savefig
from scipy.spatial import distance
from scipy.cluster import hierarchy
import scipy.cluster.hierarchy as sch
df = pd.read_csv('matrix.csv')
df = df.set_index('name')
names=np.array(df.columns.values)
plt.figure(figsize = (12,10))
sns_heat = sns.heatmap(df)
figure = sns_heat.get_figure()    
figure.savefig('heat-native.png', dpi=400)
correlations = df.corr()
correlations_array = np.asarray(df.corr())

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='complete')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='complete')
    # 'single', 'complete', 'average', 'ward'
sns_clus = sns.clustermap(correlations, row_linkage=row_linkage, col_linkage=col_linkage, figsize=(10, 10))
plt.savefig('heat-cluster.png', dpi=400)
# retrieve clusters using fcluster 
d = sch.distance.pdist(df)
L = sch.linkage(d, method='complete')

# 0.2 can be modified to retrieve more stringent or relaxed clusters
cutree = [0.9, 0.5, 0.3, 0.2, 0.15, 0.1]
my_array = []
for i in cutree :
    clusters = sch.fcluster(L, i*d.max(), 'distance')
    my_array.append(clusters)
np.savetxt('heat-clusters.csv', my_array, delimiter=",",  fmt='%s')
np.savetxt('heat-names.csv', names, fmt='%s')
" > cluster.py



########################
#    Graph function    #
########################.

########################################################################################################################################################################
#      INITIAL SCORE      #
########################################################################################################################################################################

rm -r $oridir/rosetta-0-initial score.sc 2>/dev/null 
mkdir $oridir/rosetta-0-initial
cd $oridir/rosetta-0-initial

echo "
-database /usr/local/rosetta/main/database
-nstruct 1
-out:no_nstruct_label
-score:weights ref2015" > sflag.file

if [ $RELAX = True ] ; then
	echo '0.0 Initial Score'
	if [ $RIPDB != False ] ; then
		echo '0.1 - Scoring initial pdb'
		
		cp $oridir/$RIPDB .
		nohup score_jd2.static.linuxgccrelease -in:file:s $RIPDB @sflag.file  2>/dev/null
		wait
	fi
elif [ $MUTATE = True ] ; then 	
	echo '0.0 - Initial Score'
	nuministruc=$(wc -l < $MIPDB)
	echo '0.1 - Scoring initial pdb list'
	if [ $nuministruc -gt $CPU ] ; then 
		for i in $(seq 1 $CPU) ; do
			echo '
			for j in $(seq '$i' '$CPU' $(wc -l < '$MIPDB') ) ; do
			pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$MIPDB')
			cp $oridir/$pdb .
			nohup score_jd2.static.linuxgccrelease -in:file:s $pdb @sflag.file  2>/dev/null 
			done' > run-cpu$i
			chmod 755 run-cpu$i
			./run-cpu$i &
		done
		wait
	else
		for i in $(cat $MIPDB) ; do
			cp $oridir/$i .
			nohup score_jd2.static.linuxgccrelease -in:file:s $i @sflag.file  2>/dev/null
			wait
		done
	fi
fi 
rm sflag.file
cd $oridir

########################################################################################################################################################################
#      BODY RELAX      #
########################################################################################################################################################################

if [ $RELAX = True ] ; then
	rm -r $oridir/rosetta-1-relax 2>/dev/null 
	mkdir $oridir/rosetta-1-relax
	cd $oridir/rosetta-1-relax
	
	echo '1.0 Relaxation'
################################################	
########################  CPU CONTROL
	echo '1.1 - Running relaxations'
	
	# duplicate input files for each CPU
	rm ripdb.tx 2>/dev/null 
	for i in $(seq 1 $CPU) ; do 
		echo $oridir/$RIPDB >> ripdb.tx
	done
	cp $oridir/$RIPDB .
	RIPDBL=ripdb.tx

################################################	
########################  RUN RELAXATION
	# process the input list for the number of CPUs
	for i in $(seq 1 $CPU) ; do
		#make a temp directory per cpu
		mkdir $oridir/rosetta-1-relax/temp$i
		#make a flag file per cpu
		echo "
-out:path:all $oridir/rosetta-1-relax/temp$i
-out:pdb
-relax:fast
-database /usr/local/rosetta/main/database
-nstruct $(( $RINST / $CPU ))
-score:weights ref2015" > rflag$i.file
		#make a rosetta executable per cpu and run them
		echo '
		for j in $(seq '$i' '$CPU' $(wc -l < '$RIPDBL') ) ; do
		pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$RIPDBL')
		nohup relax.static.linuxgccrelease -in:file:s $pdb @rflag'$i'.file  2>/dev/null 
		done' > run-cpu$i
		chmod 755 run-cpu$i
		./run-cpu$i &
	done
	wait
	#for i in $(cat ripdb.tx) ; do
	#	rm $i
	#done
	rm $RIPDB
################################################	
########################  MERGE OUTPUTS	
	echo '1.2 - merging outputs'
	# the output pdbs and score files are in seperate directories, with the same names, lets rename and merge
	for i in $(seq 1 $CPU) ; do
		#count files in destination directory
		#pdbcount=$(ls $oridir/rosetta-1-relax/*.pdb | wc -l )
		if [ $i = 1 ] ; then 
			#move everything from the first temp directory to the destination directory
			mv $oridir/rosetta-1-relax/temp$i/* $oridir/rosetta-1-relax
		else
			#move everything from the other temp directorys to the destination directory
			for j in $(ls temp$i/*.pdb) ; do 
				pdbcount=$(ls *.pdb | wc -l )
				#echo pdb count = $pdbcount
				#define old number (with leading 0's)
				oldnum=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4)
				#define old number (without leading 0's)
				num=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4 | sed 's/0//g')
				#define new number (without leading 0's)
				newnum=$(( 1 + $pdbcount))
				#define new number (with leading 0's)
				newnumname=$(echo 000$newnum | rev | cut -c 1-4 | rev)
				#define file name without number
				nonumname=$(echo $j | rev | cut -c 9- | awk -F'/' '{print $1}' | rev)

				# rename and move pdb and contents of score file
				mv $j $nonumname$newnumname.pdb
				grep $nonumname$oldnum temp$i/score.sc | sed "s/$nonumname$oldnum/$nonumname$newnumname/g" >> score.sc
				#echo check4
			done
			#echo check3
		fi
		#echo check2
		#count files in destination directory
	done
	# clean up
	rm -r temp* run-cpu* rflag*
	#echo check1
	cd $oridir
fi

########################################################################################################################################################################
#      BODY MUTATE     #
########################################################################################################################################################################
#echo checka
if [ $MUTATE = True ] ; then
	rm -r $oridir/rosetta-2-mutate 2>/dev/null 
	mkdir $oridir/rosetta-2-mutate
	cd $oridir/rosetta-2-mutate
	
	echo '2.0 Mutation'
################################################	
########################  INPUT PREPARATION
	echo "2.1 - Preparing mutation inputs"
	
####### 
####### 
####### CHECKING INPUT STURCUTURES
	# if an input pdb is not provided (default option is 'False'), generate an input list from the best relaxed structures
	if [ $MIPDB = False ] ; then 
		cat $oridir/rosetta-1-relax/score.sc | sed '1,2d' | awk '{print $22,$2}' | sort -n -k2 | tail -$RONSS | awk '{print $1}' | sed 's/$/.pdb/g' > mipdb-auto.tx
		for i in $(cat mipdb-auto.tx) ; do 
			cp $oridir/rosetta-1-relax/$i .
		done
	else
		cp $oridir/$MIPDB mipdb-auto.tx
		for i in $(cat mipdb-auto.tx) ; do 
			cp $oridir/$i .
		done
	fi

####### 
####### 	
####### CHECKING RESFILE
	# if an res file is not provided (optional), generate an res file from the MIRES comments
	if [ $MIRES = False ] && [ $MITHD = False ] ; then 
		echo '	generating resfile'
		rm mires-auto.tx 2>/dev/null 
		# print the header of the res file
		echo 'NATRO
start' > mires-auto.tx
		# generate mutable selection
		for chain in $(cat $oridir/mires.tx | sed '/^[[:space:]]*$/d' | awk '{print $1}') ; do
			for res in $(grep $chain $oridir/mires.tx | cut -c 2-) ; do 
				echo "$res $chain ALLAAxc" >> mires-auto.tx
			done
		done
		
####### generate movable selection - using pymol
		respdb=$(head -1 mipdb-auto.tx)
		echo "load $oridir/$RIPDB" > mires.pml
		sed '/^[[:space:]]*$/d' $oridir/mires.tx | sed "s/ /+/g;s/./ \& i. /2;s/^/+ c. /g" | tr '\n' ' ' | sed "s/^+ /select muts, /g
" >> mires.pml
		echo "
select shell, br. all near 10 of muts
save shell.pdb, shell" >> mires.pml
		nohup pymol -qc mires.pml 2>/dev/null 
		wait
		cat shell.pdb | grep ' CA ' | awk '{print $6,$5}' | sed 's/$/ NATAA/g' >> mires-auto.tx
		rm shell.pdb mires.pml
		MIRES=mires-auto.tx
	fi

####### 
####### 
####### CHECKING SEQUENCE THREADING
	if [ $MITHD != False ] ; then
		# print the header of the res file
		echo '	initiating threading'
		echo 'NATRO
start' > mires-auto.tx

####### check if the source file is a fasta or pdb
		if [ $( echo $MITHD | awk -F'.' '{print $2}' ) = fa ] ; then 
			cp $oridir/$MITHD seqgen.fa 
		# if it is a pdb, generate a fasta file of the appropriate chain with seqgen function 
		elif [ $( echo $MITHD | awk -F'.' '{print $2}' ) = pdb ] ; then 
			rm seqgen.fa 2>/dev/null 
			seqgen $oridir/$MITHD $MITHC
		fi
		# transpose horizontal source sequence to vertical sequence
		sed -n '2p' seqgen.fa | sed 's/\(.\{1\}\)/\1 /g' > transin.tx
		transpose transin.tx > source-sequence.tx	
		#mv transout.tx source-sequence.tx	
		
####### generate a fasta file from the first of the relaxed structures as a template
		rm seqgen.fa 
		seqgen $oridir/$RIPDB $MITHC
		# transpose horizontal target sequence to vertical sequence
		sed -n '2p' seqgen.fa | sed 's/\(.\{1\}\)/\1 /g' > transin.tx
		transpose transin.tx > replaced-sequence.tx
		#mv transout.tx replaced-sequence.tx
		
####### generate residue numbering from the target structure
		cat $oridir/$RIPDB | grep ^ATOM | grep ' CA ' | grep "[A-Z][A-Z][A-Z] $MITHC" | awk '{print $6}' > number-sequence.tx
		paste number-sequence.tx replaced-sequence.tx > numrep.tx
		paste number-sequence.tx source-sequence.tx > nupsor.tx
		# use sdiff to compare vertical sequences and print mismatches to the res file
		sdiff numrep.tx nupsor.tx | grep '|' | awk '{print $4,$5}' | sed "s/ / $MITHC PIKAA /g" >> mires-auto.tx
		# check the number of positions and sequences for discontinuity
		ns=$(wc -l < number-sequence.tx)
		nr=$(wc -l < numrep.tx)
		np=$(wc -l < nupsor.tx)

		if [ "$ns" != "$nr" ] || [ "$ns" != "$np" ] ; then 
			echo " please check the source and template for threading. We have found different number of residues between selection and cannot thread.
			number of positions = $ns, target AA = $nr, source AA = $np"
			echo " exiting . . . " ; exit 0
		fi 
####### generate movable selection - using pymol
		respdb=$(sort mipdb-auto.tx | head -1)
		echo "load $oridir/$RIPDB" > mires.pml
		sdiff numrep.tx nupsor.tx | grep '|' | awk '{print $4,"+"}' | tr '\n' ' ' | sed 's/ //g' | sed "s/^/select muts, c. $MITHC \& i. /g" | rev | cut -c2- | rev >> mires.pml
		echo "
select shell, br. all near 10 of muts
save shell.pdb, shell" >> mires.pml
		nohup pymol -qc mires.pml 2>/dev/null 
		wait
		cat shell.pdb | grep ' CA ' | awk '{print $6,$5}' | sed 's/$/ NATAA/g' >> mires-auto.tx
		#cat mires-auto.tx
		rm shell.pdb mires.pml
		MIRES=mires-auto.tx
	fi


################################################
########################  CPU CONTROL
	echo "2.2 - Making mutations"
	# if there are more CPU's than input files, add additional copies of the input to the file list
	# this is a shite way of using multiple CPUs, please improve.
	if [ $CPU -gt $(wc -l < mipdb-auto.tx) ] ; then 
		# we currently only accept single models or more models than the number of CPU's
		if [ $(wc -l < mipdb-auto.tx) != 1 ] ; then
			echo "this script does not currently accept more CPU's than starting models (unless single model)"
			exit 0
		fi
		
		# for a single model input with mutltiple CPU, make a new line in the pdb list per CPU
		rm mipdb.txt
		count=1
		for i in $(seq 1 $CPU) ; do 
			echo $count_$(head -1 $MIPDB) >> mipdb-auto.tx
			count=$(( $count + 1)) 
		done
		
	fi
	MIPDBL=mipdb-auto.tx
################################################	
########################  RUN MUTATION
	# process the input list for the number of CPUs
	for i in $(seq 1 $CPU) ; do
		#make a temp directory per cpu
		mkdir temp$i
		#make a flag file per cpu
		echo "
-out:path:all temp$i
-out:pdb
-relax:fast
-database /usr/local/rosetta/main/database
-nstruct $(( $MONST / $CPU ))
-relax:respect_resfile 
-packing:resfile $MIRES                                                                                                                                                              
-score:weights ref2015" > mflag$i.file
		#make a rosetta executable per cpu and run them
		echo '
		for j in $(seq '$i' '$CPU' $(wc -l < '$MIPDBL') ) ; do
		pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$MIPDBL')
		nohup relax.static.linuxgccrelease -in:file:s $pdb @mflag'$i'.file 2>/dev/null 
		done' > run-cpu$i
		chmod 755 run-cpu$i
		./run-cpu$i &
	done
	wait
################################################	
########################  MERGE OUTPUTS	
	echo '2.3 - consolidating outputs'
	# the output pdbs and score files are in seperate directories, with the same names, lets rename and merge
	
	# as we are not starting from the same initial pdb name, hopefully this will be simple
	rm sfcheck.tx 2>/dev/null 
	for i in $(seq 1 $CPU) ; do
		sed '1,2d' temp$i/score.sc | awk '{print $22}' | rev | cut -c 6- | rev | sort | uniq >> sfcheck.tx
	done
	sfchpre=$(wc -l < sfcheck.tx)
	sfchpost=$(cat sfcheck.tx | sort | uniq | wc -l )
	if [ $sfchpre = $sfchpost ] ; then
		for i in $(seq 1 $CPU) ; do
			mv temp$i/*.pdb .
			if [ $i = 1 ] ; then 
				#move everything from the first temp directory to the destination directory
				mv temp$i/score.sc .
			else
			cat temp$i/score.sc | sed '1,2d' >> score.sc
			fi
		done
	#rm -r temp* run-cpu* mflag*
	else # this is a bit shit but the names should be unique (although naming problems will exist if over 9999 models total cat
		for i in $(seq 1 $CPU) ; do
			#count files in destination directory
			pdbcount=$(ls *.pdb | wc -l )
			if [ $i = 1 ] ; then 
				#move everything from the first temp directory to the destination directory
				mv $oridir/rosetta-2-mutate/temp$i/* $oridir/rosetta-2-mutate
			else
				#move everything from the other temp directorys to the destination directory
				for j in $(ls $oridir/rosetta-2-mutate/temp$i/*.pdb) ; do 
				
					#define old number (with leading 0's)
					oldnum=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4)
					#define old number (without leading 0's)
					num=$(echo $j | rev | cut -c 1-8 | rev | cut -c 1-4 | sed 's/0//g')
					#define new number (without leading 0's)
					newnum=$(( $num + $pdbcount))
					#define new number (with leading 0's)
					newnumname=$(echo 000$newnum | rev | cut -c 1-4 | rev)
					#define file name without number
					nonumname=$(echo $j | rev | cut -c 9- | awk -F'/' '{print $1}' | rev)
					
					# rename and move pdb and contents of score file
					mv $oridir/rosetta-2-mutate/temp$i/$nonumname$oldnum.pdb $oridir/rosetta-2-mutate/$nonumname$newnumname.pdb
					grep $nonumname$oldnum $oridir/rosetta-2-mutate/temp$i/score.sc | sed "s/$nonumname$oldnum/$nonumname$newnumname/g" >> $oridir/rosetta-2-mutate/score.sc
					
				done
			fi
		done
	# clean up
	rm -r temp* run-cpu* mflag* *.tx 2>/dev/null 
	fi
	for i in $(cat $MIPDBL) ; do 
		rm $i 
	done
	ls *.pdb > mopdblist.tx
	MOPDBL=mopdblist.tx
	cd  $oridir
fi


########################################################################################################################################################################
#    BODY INTERFACE    #
########################################################################################################################################################################

if [ $INTER = True ] ; then
	rm -r $oridir/rosetta-4-inter 2>/dev/null 
	mkdir $oridir/rosetta-4-inter
	cd  $oridir/rosetta-4-inter
	
	echo '3.0 Interface Analysis'

	ls -v $oridir/rosetta-1-relax/*.pdb > pdb.list
	ls -v $oridir/rosetta-2-mutate/*.pdb >> pdb.list
	IIPDB=pdb.list

################################################	
########################  Run Interface analysis
	# process the input list for the number of CPUs
	for i in $(seq 1 $CPU) ; do
		#make a temp directory per cpu
		mkdir $oridir/rosetta-4-inter/temp$i
		#make a flag file per cpu
		echo "#specific options for InterfaceAnalyzer
-out:no_nstruct_label
-out:path:all $oridir/rosetta-4-inter/temp$i
-database /usr/local/rosetta/main/database
-fixedchains $IIFAC
-score:weights ref2015
-compute_packstat true
-tracer_data_print false #make a score file with all the important info instead of just printing to the terminal
-out:file:score_only inter_score.sc #This will cause output of all of the info to a file called inter_score.sc
-pack_input false #will not relax the input interface residues (good if you have already minimized/packed your structure)
-pack_separated false #will also not pack the monomers to calculated dG bind.
-add_regular_scores_to_scorefile true #will run the rest of rosettas score function on your complex using score12
#these are some tweeks that we have found helpful
-atomic_burial_cutoff 0.01 #This is set to help rosetta identify buried polar atoms properly
-sasa_calculator_probe_radius 1.4 #This is the default water probe radius for SASA calculations, sometimes lowering the radius helps rosetta more accurately find buried polar atoms
-pose_metrics::interface_cutoff 8.0 # this defines how far away a CBeta atom can be from the other chain to be considered an interface residue" > iflag$i.file

		#make a rosetta executable per cpu and run them
		echo '
		for j in $(seq '$i' '$CPU' $(wc -l < '$IIPDB') ) ; do
		pdb=$(sed -n '"''"'$j'"'"'p'"'"' '$IIPDB')
		nohup InterfaceAnalyzer.static.linuxgccrelease -in:file:s $pdb @iflag'$i'.file 2>/dev/null
		done' > run-cpu$i
		chmod 755 run-cpu$i
		./run-cpu$i &
	done
	wait

	
	for i in $(seq 1 $CPU) ; do
		if [ $i = 1 ] ; then 
			#move everything from the first temp directory to the destination directory
			cat $oridir/rosetta-4-inter/temp$i/inter_score.sc  > $oridir/rosetta-4-inter/inter_score.sc
		else
			cat $oridir/rosetta-4-inter/temp$i/inter_score.sc | sed '1,2d' >> $oridir/rosetta-4-inter/inter_score.sc
		fi
	done
	rm -r temp* iflag* run-* nohup-out pdb.list 2>/dev/null 
fi
cd  $oridir
########################################################################################################################################################################
#    Body Analysis    #
########################################################################################################################################################################
echo "4.0 Analysis"
################################################	
########################  Generate total fasta sequence files. grep
if [ $MOSEQ = True ] || [ $MOMSQ = True ] ; then 

	rm -r $oridir/rosetta-5-analysis 2>/dev/null 
	mkdir $oridir/rosetta-5-analysis
	cd  $oridir/rosetta-5-analysis
fi

if [ $MOSEQ = True ] ; then 
	#echo 'all chain initial'
	rm seqgen.fa  2>/dev/null 
	name=$(echo $RIPDB | sed 's/.pdb//g' )
	#echo 'file = '$RIPDB
	seqgen $oridir/$RIPDB
	#sed -i "s/$RIPDB/$name/g" seqgen.fa

	for i in $(cat $oridir/rosetta-2-mutate/$MOPDBL) ; do 
		#echo 'all chain mutant'
		name=$(echo $i | sed 's/.pdb//g' )
		seqgen $oridir/rosetta-2-mutate/$i 
		#sed -i "s/$i/$name/g" seqgen.fa
	
	done
	mv seqgen.fa roestta-full-sequence.fa
fi 

################################################	
########################  Generate specific fasta sequence files with total energy and interface energy (if appropriate)
if [ $MOMSQ = True ] ; then 
	#echo 'RIPDB = '$RIPDB
	#echo 'spec chain initial'
	rm seqgen.fa  2>/dev/null 
	cat $oridir/rosetta-2-mutate/$MIRES | sed '1,2d' | grep -v 'NATAA' | grep -v 'NATRO' > specseqgen.tx
	rm specseqpdb.tx 2>/dev/null 
	for i in $(seq 1 $(wc -l < specseqgen.tx) ) ; do 
		#sed -n ''$i'p' specseqgen.tx
		specseqpos=$(echo "---$(sed -n ''$i'p' specseqgen.tx | awk '{print $1}')" | rev | cut -c 1-4 | rev)
		specseqcha=$(sed -n ''$i'p' specseqgen.tx | awk '{print $2}')
		specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
		grep ' CA ' $oridir/$RIPDB | grep "$specseqpch" >> specseqpdb.tx
	done
	name=$(echo $RIPDB | sed 's/.pdb//g' )
	mv specseqpdb.tx $name
	#echo 'name ='$name 'chain ='$specseqcha
	seqgen $name $specseqcha
	rm $name

	for j in $(cat $oridir/rosetta-2-mutate/$MOPDBL) ; do 
		#echo 'spec chain mutant'
		for i in $(seq 1 $(wc -l < specseqgen.tx) ) ; do 
			#sed -n ''$i'p' specseqgen.tx
			specseqpos=$(echo "---$(sed -n ''$i'p' specseqgen.tx | awk '{print $1}')" | rev | cut -c 1-4 | rev)
			specseqcha=$(sed -n ''$i'p' specseqgen.tx | awk '{print $2}')
			specseqpch=$(echo $specseqcha$specseqpos | sed 's/-/ /g')
			grep ' CA ' $oridir/rosetta-2-mutate/$j | grep "$specseqpch" >> specseqpdb.tx
		done
		name=$(echo $j | sed 's/\// /g' | rev | cut -d' ' -f 1 | rev | sed 's/.pdb//g' )
		mv specseqpdb.tx $name
		seqgen $name $specseqcha
		rm $name
	done
	mv seqgen.fa roestta-spec-sequence.fa
fi

################################################	
########################  Combine Sequence and scores

if [ $MOSQE = True ] && [ $MOMSQ = True ]; then 
	echo 'name, total REU, Interface REU, Mutations-'$(sed -n '2p' roestta-spec-sequence.fa) > Combined-Output.csv
	for i in $(sed '1,2d' $oridir/rosetta-4-inter/inter_score.sc | awk '{print $42}') ; do 
		scint=$(grep $i$ $oridir/rosetta-4-inter/inter_score.sc | awk '{print $6}')
		sctot=$(grep $i$ $oridir/rosetta-[12]*/score.sc | awk '{print $6}')
		scseq=$(grep -A1 $i'_'[A-Z] $oridir/rosetta-5-analysis/roestta-spec-sequence.fa | sed -n '2p')
		echo $i, $sctot, $scint, $scseq >> Combined-Output.csv
	done

elif [ $MOSQE = True ] ; then 
	echo 'name, total REU, Interface REU' > Combined-Output.csv
	for i in $(sed '1,2d' $oridir/rosetta-4-inter/inter_score.sc | awk '{print $42}') ; do 
		scint=$(grep $i$ $oridir/rosetta-4-inter/inter_score.sc | awk '{print $6}')
		sctot=$(grep $i$ $oridir/rosetta-[12]*/score.sc | awk '{print $6}')
		echo $i, $sctot, $scint >> Combined-Output.csv
	done
fi		
		

rm *.tx  2>/dev/null 
cd $oridir
########################################################################################################################################################################
#    Clean up    #
########################################################################################################################################################################

echo '5.0 Complete!'
#mv score.sc rosetta-0-initial/.
#rm run-cpu* *.tx nohup.out *flag[0-9]* pdb.list cluster.py 2>/dev/null 
