#!/bin/bash -l

#SBATCH -A b2015285
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 1:00:00
#SBATCH -J pairseqread

FOLDER="."
OUTFILE="seqtable.out"
PYSCRIPT="pairedseqread_cf.py"
STATFILE="readstats.txt"
CFGFILE="confg.yaml"

# change to "bz2" depending on the compression of the sequence files
COMP="gz" 	

pairedread(){

	if [ -z "$1" ]
	then
		echo "No directory given, will look for files in the current dir..."
	else
		echo "Looking for paired read files in \"$1\""
		FOLDER="$1"
	fi

	echo "Dataset_folder"$'\t'"all_seqs"$'\t'"low-Q"$'\t'"SC_seqs"$'\t'"filtered"$'\t'"distinct" $'\t'"promille"> "${STATFILE}"
	for f in "$FOLDER"/*_1.fastq."$COMP"
	do
		echo "Processing \"$f\" "
		
		# Is this pairfile check needed???
		FILENAME=${f%_*}
		pairfile="$FILENAME"_2.fastq."$COMP"

		echo "Looking for the pair file $pairfile"
		if [ -r "$pairfile" ]
		then
			echo "Paired read files found, executing read script..."
			python $PYSCRIPT "$f" "$pairfile" -c $CFGFILE -o $FOLDER/$OUTFILE -a $STATFILE
		else
			echo "Paired read file missing, skipping \"$f\" "
			continue
		fi
	done
	echo
}

TOPDIR="."
if [ $# -gt 0 ]; then
  TOPDIR=$1
fi

while read DIR; do
	test -r "$DIR"/*_1.fastq."$COMP" -a -r "$DIR"/*_2.fastq."$COMP" || continue
	pairedread $DIR &
done < <(find $TOPDIR -type d -print)

echo "Reading sequences..."
wait
echo "Sequence reading done..."

find . -name 'seqdata_paired.json' -print0 | sort -z |xargs -0 python jsonmerge.py -c $CFGFILE

echo "Tables merged..."
