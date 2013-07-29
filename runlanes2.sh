#!/bin/bash

PS3="Select the alignment tool: "
select align in  bwa bowtie
do
if [ "$align" = "" ]; then
	echo "Error in entry."
	exit 1
fi
break
done

echo "Using alignment software: $align."

if [ ! "$(ls -A $HOME/references)" ]; then
    echo "Reference folder is empty!"
    exit 1
fi

PS3="Select the reference from the file list: "
	select reference in `find $HOME/references/ -type f -name "*.fasta"`
	do
		if [ "$reference" = "" ]; then
        		echo "Error in entry."
        		exit 1
		fi

       	break
	done

	echo "Using reference file:  $reference."

cores=`grep -ic ^processor /proc/cpuinfo`
echo "Detected $cores cores!"
echo "How many mismatches shall be allowed?"
read mismatches
echo "How many bases shall be trimmed from the 3' end?"
read threetrim
echo "How many bases shall be trimmed from the 5' end?"
read fivetrim

if [ "$align" = "bowtie" ]
then	

	if [ ! -e $reference.1.ebwt ]; then
		echo "The reference has not yet been indexed by Bowtie! Please index the reference first and try again."
		exit 1
	fi
	mkdir BOWTIE_ANALYSIS--`date +%F__%H-%M`
	FOLDER=$_
	echo "Creating output folder $_"	
	echo `date` >> runlanes_bowtie.log
	echo "Reference: " $reference >> runlanes_bowtie.log
	#for datei in *sequence.txt
	for datei in `find . -maxdepth 1 -name '*.fq' -o -name '*sequence.txt' -o -name '*.fastq'`
	do 
		echo "Aligning $datei"
		echo "Aligning $datei" >> runlanes_bowtie.log
		#bowtie -v $mismatches -y -p $cores -m 1 --norc --best -S --un --an $reference $datei $FOLDER/$(basename $datei .txt).sam 2>> runlanes_bowtie.log
		bowtie -v $mismatches -3 $threetrim -5 $fivetrim -y -p $cores -m 1 --norc --best -S --un $FOLDER/$(basename $datei .txt)_unaligned.fq $reference $datei 1>$FOLDER/$(basename $datei .txt).sam 2>> runlanes_bowtie.log
	done
fi

if [ "$align" = "bwa" ]
then
	if [ ! -e $reference.ann ]; then
		echo "The reference has not yet been indexed by BWA! Please index the reference first and try again."
		exit 1
	fi
	mkdir BWA_ANALYSIS--`date +%F__%H-%M`
	FOLDER=$_
	echo "Creating output folder $_"
	echo `date` >> runlanes_bwa.log
	echo "Reference: " $reference >> runlanes_bwa.log
	for datei in `find . -maxdepth 1 -name '*.fq' -o -name '*sequence.txt' -o -name '*.fastq'`
	do 
		echo "Aligning $datei"
		echo "Aligning $datei" >> runlanes_bwa.log
		bwa aln $reference $datei -t $cores -k $mismatches > $FOLDER/$(basename $datei .txt).sai 2>> runlanes_bwa.log
		bwa samse $reference $FOLDER/$(basename $datei .txt).sai $datei > $FOLDER/$(basename $datei .txt).sam 2>> runlanes_bwa.log

	done
fi
