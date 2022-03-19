#!/bin/bash
export WHITE="\e[1;37m"
export LGRAY="\e[0;37m"
export GRAY="\e[1;30m"
export BLACK="\e[0;30m"
export RED="\e[0;31m"
export LRED="\e[1;31m"
export GREEN="\e[0;32m"
export LGREEN="\e[1;32m"
export BROWN="\e[0;33m"
export YELLOW="\e[1;33m"
export MAGENTA="\e[35m"
export BLUE="\e[0;34m"
export LBLUE="\e[1;34m"
export PURPLE="\e[0;35m"
export PINK="\e[1;35m"
export CYAN="\e[0;36m"
export LCYAN="\e[1;36m"
export Z="\e[0m"

input=$(readlink -f $(dirname "$0"))
output="./HBVIF_folder"
TRIM=$(locate trimmomatic.jar | head -n 1 | awk 'BEGIN{FS=OFS="/"}{$NF=""; print $0}' )
adapt=$(readlink -f $TRIM/adapters/TruSeq2-PE.fa)
trimmingwindow=5
trimmingquality=16
trimmingminlen=35
virus=$(NC_003977.2)
skipTrimming=0
threads=1

function helptxt() {
	echo -e "\nHBVIF (Hepatitis B Virus Integration Finder)"
	echo -e "This script take as input raw Illumina reads and an appropriate Human-HBV hybrid
	 genome and search for HBV integration sites.

	 Arguments:
	 	--threads : number of threads (default 1)
	 	--HBV_header : header of the HBV genome used in the Human-HBV hybrid genome (default NC_003977.2)
	 	--input : PATH to reads files (default current PATH)
	 	--genome : PATH to hybrid genome (please insert genome, annotation and other files according to the maunal) [REQUIRED]
	 	--output : PATH to the dir to output all the analysis and results
	 	--skip-trimming : option to diasble trimming during file processing (default off)
	 	--trimmomatic : PATH to the dir containing trimmomatic.jar executable
	 	--adapters : file containing list of adapters used during the trimming (default in $Trimmomatic_PATH/adapters/TruSeq2-PE.fa)
	 	--trimming-window : size of the window used in trimmomatic SLIDINGWINDOW mode (default 5)
	 	--trimming-quality : minimum Phred score value mantained during SLIDINGWINDOW mode (default 16)
	 	--trimming-minlength : minimum read length mantained after trimming (default 35)
	 	--picard : PATH to the dir containing picard.jar 
	 	--help : print this message
	 	"

	exit 1
}

if [ -z "$1" ]; then
    helptxt
else
    while [ $1 ]; do
        case $1 in
            --threads)
                shift
                threads=$1
                ;;
            --HBV_header)
                shift
                virus=$1
                ;;
            --input)
                shift
                input=$(readlink -f $1)
                ;;
            --genome)
                shift
                genome=$(readlink -f $1)
                ;;
            --output)
                shift
                output=$(readlink -f $1)
                ;;
            --skip-trimming)
                skiptrimming=1
                ;;
            --trimmomatic)
				shift
                TRIM=$(readlink -f $1)
                ;;
            --adapters)
				shift
                adapt=$(readlink -f $1)
                ;;
            --trimming-window)
                shift
                trimmingwindow=$1
                ;;
            --trimming-quality)
                shift
                trimmingquaity=$1
                ;;
            --trimming-minlength)
                shift
                trimmingminlen=$1
                ;;
            --picard)
                shift
                PIC=$(readlink -f $1)
                ;;
            --help)
                helptxt
                ;;
            *)
                echo -e "\n${bold}$1 is not a valid parameter.${normal}\n" >&2
                exit 1
                ;;
        esac
        shift
    done
fi

#checking arguments
if [[ -z $genome ]] || [[ ]]; then #add INPUT dir!!!
  helptxt
fi


echo -e $RED "########## Running HBVIF pipeline ##########" $Z
isGzipped=$(ls $input | head -n 1 | grep "gz$")
if [[ -z $isGzipped ]]; then  #if files are not gzipped
	cat $input/*_R1_*q | gzip > R1.fq.gz
	cat $input/*_R2_*q | gzip > R2.fq.gz
else #if files are gzipped
	cat $input/*_R1_*.gz > R1.fq.gz
	cat $input/*_R2_*.gz > R2.fq.gz
fi

echo -e $YELLOW "########## Cleaning reads with Trimmomatic ##########" $Z

if [[ $skiptrimming -eq 0 ]]; then
	java -jar $TRIM/trimmomatic.jar PE -threads $threads R1.fq.gz R2.fq.gz R1_paired_sw_$trimmingwindow_$trimmingquality.fq.gz R1_unpair_sw_$trimmingwindow_$trimmingquality.fq.gz R2_paired_sw_$trimmingwindow_$trimmingquality.fq.gz R2_unpair_sw_$trimmingwindow_$trimmingquality.fq.gz ILLUMINACLIP:$adapt:2:30:10 SLIDINGWINDOW:$trimmingwindow:$trimmingquality MINLEN:$trimmingminlen
	RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "trimmomatic failed at line 139. Exiting..." >&2
		exit $RETURN_CODE
	fi
else
	ln -s R1.fq.gz R1_paired_sw_$trimmingwindow_16.fq.gz
	ln -s R2.fq.gz R2_paired_sw_$trimmingwindow_16.fq.gz
fi

echo -e $YELLOW "## Running FastQC on Raw and Filtered reads ##" $Z

fastqc -t $threads R1.fq.gz R2.fq.gz R1_paired_sw_$trimmingwindow_$trimmingquality.fq.gz R1_unpair_sw_$trimmingwindow_$trimmingquality.fq.gz R2_paired_sw_$trimmingwindow_$trimmingquality.fq.gz R2_unpair_sw_$trimmingwindow_$trimmingquality.fq.gz
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "FASTQC failed at line 152. Exiting..." >&2
		exit $RETURN_CODE
	fi
mkdir -p reads_and_trimming
mv *.html reads_and_trimming/
mv *.zip reads_and_trimming/

echo -e $YELLOW"## Joining unpaired reads ##"$Z

cat R1_unpair_sw_$trimmingwindow_$trimmingquality.fq R2_unpair_sw_$trimmingwindow_$trimmingquality.fq > R12_unpair_sw_$trimmingwindow_$trimmingquality.fq
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "cat of unpair trimmed reads failed at line 164. Exiting..." >&2
		exit $RETURN_CODE
	fi
mv R1_unpair_sw_$trimmingwindow_$trimmingquality.fq reads_and_trimming/
mv R2_unpair_sw_$trimmingwindow_$trimmingquality.fq reads_and_trimming/

echo -e $YELLOW"########## Mapping reads (pair and unpair) on the hybrid genome ##########" $Z

bwa mem -t $threads -V -a -Y $genome/HG_HBV.fa R1_paired_sw_$trimmingwindow_$trimmingquality.fq R2_paired_sw_$trimmingwindow_$trimmingquality.fq | samtools view -h -b -@ $threads > R12_pair.bam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "bwa failed. Exiting..." >&2
	exit $RETURN_CODE
fi

bwa mem -t $threads -V -a -Y $genome/HG_HBV.fa R12_unpair_sw_$trimmingwindow_$trimmingquality.fq | samtools view -h -b -@ $threads > R12_unpair.bam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "bwa failed. Exiting..." >&2
	exit $RETURN_CODE
fi

echo -e $YELLOW"#### Sorting bam files and removing optical duplicates with picard (GATK) ####" $Z

samtools sort -o R12_pair_sort.bam -@ $threads R12_pair.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools paired sort failed. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools sort -o R12_unpair_sort.bam -@ $threads R12_unpair.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools unpaired sort failed. Exiting..." >&2
		exit $RETURN_CODE
fi
java -jar $PIC/picard.jar MarkDuplicates I= R12_pair_sort.bam O= R12_pair_sort_NoDuplicates.bam M= R12_pair_sort_MarkDuplicates.txt REMOVE_DUPLICATES=True 
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "picard on pair failed. Exiting..." >&2
		exit $RETURN_CODE
fi
java -jar $PIC/picard.jar MarkDuplicates I= R12_unpair_sort.bam O= R12_unpair_sort_NoDuplicates.bam M= R12_unpair_sort_MarkDuplicates.txt REMOVE_DUPLICATES=True 
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "picard on unpair failed. Exiting..." >&2
		exit $RETURN_CODE
fi
echo -e $YELLOW "############ Create and moving intermediate bam files into 'intermediate_bam' ##########" $Z

mkdir -p intermediate_bam
mv R12_pair.bam intermediate_bam/
mv R12_unpair.bam intermediate_bam/
mv R12_unpair_sort.bam intermediate_bam/
mv R12_pair_sort.bam intermediate_bam/
mv R12_pair_sort_MarkDuplicates.txt intermediate_bam/
mv R12_unpair_sort_MarkDuplicates.txt intermediate_bam/

echo -e $YELLOW "########## Extracting all CHIMERAs in R12_SA.bam ##########" $Z

samtools view -h -f 2048 R12_pair_sort_NoDuplicates.bam > R12_pair_SA.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools chimeras extraction on paired failed at line 218. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools view -h -f 2048 R12_unpair_sort_NoDuplicates.bam > R12_unpair_SA.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools chimeras extraction on unpaired failed at line 224. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools merge -@ $threads R12_SA_unfiltered.bam R12_pair_SA.bam R12_unpair_SA.bam 
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools merge failed at line 230. Exiting..." >&2
		exit $RETURN_CODE
fi

samtools view -h R12_SA_unfiltered.bam | awk -v vir=$virus -F '\t' '{if($0~/^@/ || $3!="$vir" && $7!="=" && $12~"SA:Z:$vir") {print} else {if($0~/^@/ || $3=="$vir" && $7!="=" && $12!~"SA:Z:$vir") {print}}}' | samtools view -Sb - > R12_SA_to_sort.bam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "samtools chimera extraction in R12_SA_to_sort failed at line 237. Exiting..." >&2
	exit $RETURN_CODE
fi

samtools sort -n -@ $threads R12_SA_to_sort.bam > R12_SA.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools sort R12_SA.bam failed at line 244. Exiting..." >&2
		exit $RETURN_CODE
fi
echo -e $YELLOW "## Extracting CHIMERAs from paired-reads ##" $Z

samtools view -h R12_pair_sort_NoDuplicates.bam | awk -v vir=$virus -F '\t' '{if($0~/^@/ || $3!="$vir" && $7=="$vir") print}' | samtools view -Sb > R12_pair.human.virus.bam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "extracting R12_pair.human.virus.bam failed at line 252. Exiting..." >&2
	exit $RETURN_CODE
fi

samtools view -h R12_pair_sort_NoDuplicates.bam | awk -v vir=$virus -F '\t' '{if($0~/^@/ || $3=="$vir" && $7!="=") print}' | samtools view -Sb - > R12_pair.virus.human.bam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "extracting R12_pair.virus.human.bam failed at line 259. Exiting..." >&2
	exit $RETURN_CODE
fi
samtools sort -n -@ $threads R12_pair.human.virus.bam > R12_pair.human.virus.SortedByReadName.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools sort failed at line 265. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools sort -n -@ $threads R12_pair.virus.human.bam > R12_pair.virus.human.SortedByReadName.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools sort failed at line 271. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools merge -n -@ $threads R12_pair_HG_HBV.bam R12_pair.human.virus.SortedByReadName.bam R12_pair.virus.human.SortedByReadName.bam 
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools merge failed at line 277. Exiting..." >&2
		exit $RETURN_CODE
fi
echo -e $YELLOW "## Remove Chimeras from Hybrid reads to not count it twice ##" $Z

samtools view -h -F 2048 R12_pair_HG_HBV.bam > R12_pair_HG_HBV_no_SA.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools view extracting chimera failed at line 285. Exiting..." >&2
		exit $RETURN_CODE
fi
echo -e $YELLOW "## Merging Hybrid and all chimeric reads ##" $Z

samtools merge -@ $threads R12_SA_HG_HBV.bam R12_pair_HG_HBV_no_SA.bam R12_SA.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools merge on chimeras failedat line 293. Exiting..." >&2
		exit $RETURN_CODE
fi

echo -e $GREEN "### Files produces are:" $Z
echo -e $LRED "# R12_pair_HG_HBV.bam ## It contains all Hybrid reads (Pairs that map one on human and one on virus, or pairs containing in one of those a chimeric sequence #" $Z
echo -e $LCYAN "# R12_SA.bam ## It contains all Chimeric reads derived from both unpair and paired reads #" $Z
echo -e $MAGENTA "# R12_SA_HG_HBV.bam ## It contains all Chimeric reads as well as Hybrid paired-end reads mapping one on human and one on virus genome #" $Z

echo -e $YELLOW"## Moving all intermediate files into the folder 'file_after_filter_useless' and movind all final fastq files in 'final_fastq_files' ## "$Z

mkdir -p file_after_filter_useless
mkdir -p final_fastq_files
mv R12_pair.human.virus.bam file_after_filter_useless/
mv R12_pair.human.virus.SortedByReadName.bam file_after_filter_useless/
mv R12_pair_SA.bam file_after_filter_useless/
mv R12_pair_sort_NoDuplicates.bam file_after_filter_useless/
mv R12_pair.virus.human.bam file_after_filter_useless/
mv R12_pair.virus.human.SortedByReadName.bam file_after_filter_useless/
mv R12_SA_to_sort.bam file_after_filter_useless/
mv R12_SA_unfiltered.bam file_after_filter_useless/
mv R12_unpair_SA.bam file_after_filter_useless/
mv R12_unpair_sort_NoDuplicates.bam file_after_filter_useless/
mv R12_unpair_sw_$trimmingwindow_$trimmingquality.fq final_fastq_files/
mv R1.fq final_fastq_files/
mv R1_paired_sw_$trimmingwindow_$trimmingquality.fq final_fastq_files/
mv R2.fq final_fastq_files/
mv R2_paired_sw_$trimmingwindow_$trimmingquality.fq final_fastq_files/

echo -e $YELLOW"## Generating bed files of viral integrations (primary and secondary alignment) ## "$Z
bedtools bamtobed -cigar -i R12_SA.bam | awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" "PA_"$4 "\t" $5 "\t" $6 "\t" $7}' > SA_primary.bed
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating bed files with primary and secondary alignment failed at line 326. Exiting..." >&2
	exit $RETURN_CODE
fi
samtools view -H R12_SA.bam > SA_secondary_only_header.sam 
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "samtools extracting header failed at line 332. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools view R12_SA.bam | awk -F '\t' '{split($12,a,","); if(a[3] =="+") {print $1 "\t" "2048" "\t" a[1] "\t" a[2] "\t" a[5] "\t" a[4] "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" "SA:Z:" $3 "," $4 "," $2 "," $6 "," $5} else {print $1 "\t" "2064" "\t" a[1] "\t" a[2] "\t" a[5] "\t" a[4] "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" "SA:Z:" $3 "," $4 "," $2 "," $6 "," $5}}' | awk -F '\t' '{gsub(/SA:Z:/,"",$3); print $1 "\t" $2 "\t" $3 "\t" $4 "\t"$5 "\t" $6 "\t"$7 "\t"$8 "\t"$9"\t" $10"\t" $11"\t" $12}' > SA_secondary_no_header.sam
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating SA_secondary_no_header.sam file failed at line 338. Exiting..." >&2
	exit $RETURN_CODE
fi
cat SA_secondary_only_header.sam SA_secondary_no_header.sam > SA_secondary.sam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "Adding header on SA_secondary file failed at line 344. Exiting..." >&2
		exit $RETURN_CODE
fi
samtools view -h -bS SA_secondary.sam > SA_secondary.bam
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "Converting SAM to BAM file failed at line 359. Exiting..." >&2
		exit $RETURN_CODE
fi
bedtools bamtobed -cigar -i SA_secondary.bam | awk -F '\t' '{print $1 "\t" $2 "\t" $3 "\t" "SA_"$4 "\t" $5 "\t" $6 "\t" $7}' > SA_secondary.bed
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating SA_secondary.bed file failed at line 365. Exiting..." >&2
	exit $RETURN_CODE
fi
sortBed -i SA_primary.bed > SA_primary_sorted.bed
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "Sorting and generating SA_primary_sorted.bed file failed at line 371. Exiting..." >&2
		exit $RETURN_CODE
fi
sortBed -i SA_secondary.bed > SA_secondary_sorted.bed
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "Sorting and generating SA_secondary_sorted.bed file failed at line 377. Exiting..." >&2
		exit $RETURN_CODE
fi
cat SA_primary_sorted.bed SA_secondary_sorted.bed > SA_primary_secondary_raw.bed
RETURN_CODE=$(echo $?)
	if [[ $RETURN_CODE -ne 0 ]]; then
		echo "Merging SA_primary_sorted and SA_secondary_sorted files failed at line 383. Exiting..." >&2
		exit $RETURN_CODE
fi
sortBed -i SA_primary_secondary_raw.bed |  bedtools merge -i - | bedtools intersect -wao -a - -b SA_primary_secondary_raw.bed | awk -F '\t' '{print $1"\t" $2 "\t" $3 "\t" $7 "\t" "1" "\t" $10 "\t" $8}' | bedtools groupby -g 1,2,3 -c 4,5 -delim "," -o distinct,sum > SA_primary_secondary_merged.bed
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating SA_primary_secondary_merged.bed file failed at line 389. Exiting..." >&2
	exit $RETURN_CODE
fi
awk 'BEGIN {FS="\t"; OFS="\t"} {if($1 !~ "$vir") print }' SA_primary_secondary_merged.bed | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1, $2, $3, FNR";"$4, $5}' > SA_primary_secondary_no_HBV.bed
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating SA_primary_secondary_no_HBV.bed file failed at line 395. Exiting..." >&2
	exit $RETURN_CODE
fi
## Creating bed file for its use for the MonteCarlo test ##
sortBed -i SA_primary_secondary_raw.bed |  bedtools merge -i - | bedtools intersect -wao -a - -b SA_primary_secondary_raw.bed | awk -F '\t' '{print $1"\t" $2 "\t" $3 "\t" $7 "\t" "1" "\t" $10 "\t" $8 "\t" $9}' | bedtools groupby -g 1,2,3 -c 4,5,8,6 -delim "," -o distinct,sum,distinct,distinct | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"event"NR"_"$7,$5,$6}' > SA_for_MonteCarlo.bed
RETURN_CODE=$(echo ${PIPESTATUS[@]} | awk 'BEGIN{e=0}{split($0,a," ") ; for (i in a) {if (a[i]!=0) {e=1}}}END{print e}')
if [[ $RETURN_CODE -ne 0 ]]; then
	echo "Generating SA_for_MonteCarlo.bed file failed at line 402. Exiting..." >&2
	exit $RETURN_CODE
fi
### adding in the 4 column of the bed file a progressive number. This will be our hash after### 
echo $LBLUE "### Preparing for integration site-specific chimera reconstruction ###" $Z
mkdir -p reads_per_assembly
awk 'BEGIN {FS="\t"; OFS="\t"} {split($4,a,";"); print a[2] > "./reads_per_assembly/"a[1]".txt"}' SA_primary_secondary_no_HBV.bed

echo -e $LRED "########## Collecting for each integration the corresponding supporting reads ##########"$Z

cd reads_per_assembly/
ls | while read line; do sed 's/,/\n/g' $line | cut -f2 -d "_" | cut -f1 -d"/" > $line"_list"; done


ls *_list | while read line; do samtools view ../R12_SA.bam | grep -w -F -f $line - | awk '{print ">"$1 "\n" $10}' > $line".fasta"; done

echo -e $CYAN "### Running the first round of assembly per integration site event (using CAP3) ###"$Z
ls *fasta | while read line; do cap3 $line -x "assembly1_"$line; done

for f in *.ace ; do rm "$f"; done
for f in *.links ; do rm "$f"; done
for f in *.qual ; do rm "$f"; done
for f in *.info ; do rm "$f"; done
ls *.assembly* | paste - - | while read line; do awk '{split(FILENAME,a,"."); print $0 > a[1]".final_assembly_1.fasta"}' $line ; done
for f in *.contigs; do rm "$f"; done
for f in *singlets ; do rm "$f"; done
echo -e $CYAN " Running the second round of assembly per integration site event (using CAP3) ###" $Z

ls *.final_assembly_1.fasta | while read line; do cap3 $line -x "assembly_2"; done

for f in *.ace ; do rm "$f"; done
for f in *.links ; do rm "$f"; done
for f in *.qual ; do rm "$f"; done
for f in *.info ; do rm "$f"; done
ls *.assembly_2* | paste - - | while read line; do awk '{split(FILENAME,a,"."); print $0 > a[1]".final_assembly_2.fasta"}' $line ; done
for f in *.contigs; do rm "$f"; done
for f in *singlets ; do rm "$f"; done

echo -e $RED "###  Running the last round of assembly pr integration site event (using cd-hit) ###" $Z

for i in *assembly_2.fasta*; do echo $i; cd-hit -i $i -o $i".cdhit.fasta"; done

for f in *.clstr ; do rm "$f"; done
ls *.cdhit.fasta | paste - - | while read line; do awk '{split(FILENAME,a,"."); print $0 > a[1]".final_assembly_3.fasta"}' $line ; done

for i in *.final_assembly_3.fasta; do echo $i; awk '{if (/^>/) {f++} ; split(FILENAME,a,".") ; print > a[1]"."f".last.fasta"}' $i ; done

mkdir -p final_assembled_to_blast

mv *.last.fasta final_assembled_to_blast/

cd final_assembled_to_blast/

#echo -e $LGREEN " genero il database del genoma ibrido con makeblastdb"$Z
#cd $genome
#makeblastdb -input_type fasta -dbtype nucl -in HG_HBV.fa -out HG_HBV

#cd $input/reads_per_assembly/final_assembled_to_blast

echo -e $RED "########## Performing BLAST iteration per integration site ##########" $Z
echo -e LCYAN "This can take very long time, depending on the total number of integration sites found. Please be patient" $Z
numberOfSamples=$(ls *.fasta | wc -l)
COUNTER=0
ls *.fasta | while read line; do
	let COUNTER+=1 
	echo $line", which is the integration event "$COUNTER" out of "$numberOfSamples
	blastn -task blastn-short -dust no -soft_masking false -word_size 7 -num_threads $threads -max_target_seqs 2000 -db ../../genome/HG_HBV -query $line -outfmt 6 -penalty -3 -reward 2 -gapopen 5 -gapextend 2 | awk 'BEGIN {FS="\t";OFS = "\t"} {print $2,$9,$10,$7,$8}' > $line.intermediate.blast.out
done
echo -e $LRED "## Collecting and sorting BLAST results ###" $Z

ls *.intermediate.blast.out | while read line; do grep "$virus" $line | head -n 1 > $line.HBV; done

ls *.intermediate.blast.out | while read line; do E=$(echo $line | cut -f1 -d "."); C=$(cat ../../SA_primary_secondary_no_HBV.bed | awk -F '\t' -v E=$E '{split($4,a,";");if(a[1]==E) print $1}'); S=$(cat ../../SA_primary_secondary_no_HBV.bed | awk -F '\t' -v E=$E '{split($4,a,";");if(a[1]==E) print $2}');F=$(cat ../../SA_primary_secondary_no_HBV.bed | awk -F '\t' -v E=$E '{split($4,a,";");if(a[1]==E) print $3}'); awk -F '\t' -v C="$C" -v S="$S" -v F="$F" '{if($1==C && $2>=S-100 && $3<=F+100) print $0}' $line | head -n1 > $line.HUMAN; done

ls *intermediate.blast.out.* | paste - - | while read line; do awk '{split(FILENAME,a,"."); print $0 > a[1]"."a[2]".final.blast.out"}' $line; done

ls *.final.blast.out | while read line; do paste -s $line > $line".one_line"; done

mkdir -p intermediate_results

mv *.blast.out ./intermediate_results/

mv *.blast.out.HBV ./intermediate_results/
mv *.blast.out.HUMAN ./intermediate_results/

echo -e $LCYAN " Collecting blast human and viral coordinates for each integration breakpoint "$Z
ls | paste - - > list_for_cicle

cat list_for_cicle | while read line; do awk '{split(FILENAME,a,"."); print $0 >> a[1]"."a[2]".cicle"}' $line; done

mkdir -p final_cycle
mv *.cicle final_cycle/
cd final_cycle/

ls *.cicle | while read line; do awk 'BEGIN {FS="\t";OFS = "\t"} {if($0 !~/^>/) print}' $line | awk 'BEGIN {FS="\t";OFS= "\t"} {if (NR==1) {print} else {str=str""$0}}END{print str}' | awk 'BEGIN {FS="\t";OFS= "\t"} {if (NR==1) {split($0,a,"\t"); print} else {if(a[10]>a[5] && a[5]>=a[9] && a[9]>a[4]) print substr($0,a[4],a[5]-a[4]+1)"\t"substr($0,a[9],a[10]-a[9]+1)"\t"$0"\t"a[9]"\t"a[5]"\t"substr($0,a[9],a[5]-a[9]+1); else {if(a[5]>a[10] && a[10]>=a[4] && a[4]>a[9]) print substr($0,a[4],a[5]-a[4]+1)"\t"substr($0,a[9],a[10]-a[9]+1)"\t"$0"\t"a[4]"\t"a[10]"\t"substr($0,a[4],a[10]-a[4]+1); else print substr($0,a[4],a[5]-a[4]+1)"\t"substr($0,a[9],a[10]-a[9]+1)"\t"$0"\t""\t""\t""no MH"}}}' | paste -s > $line".final_breakpoints.tab" ; done

mkdir -p integrations_with_sequences/
mv *.final_breakpoints.tab ./integrations_with_sequences/
cd integrations_with_sequences/
ls *.tab | while read line; do awk 'BEGIN {FS=="\t"; OFS=="\t"}{split(FILENAME,a,"."); print a[1]"."a[2]"\t"$0 > a[1]"."a[2]".table_to_append"}' $line; done
echo -e $LBLUE " merging all integrations together " $Z
cat *.table_to_append > final_dataset_to_append.tab
cp ../../../../SA_primary_secondary_no_HBV.bed ./

awk 'BEGIN{FS="\t"; OFS="\t"}{split($4,a,";"); print $1,$2,$3,a[2],$5,a[1]}' SA_primary_secondary_no_HBV.bed > SA_primary_secondary_no_HBV_for_hash.bed 
awk 'BEGIN{FS="\t"; OFS="\t"}{if(FILENAME=="final_dataset_to_append.tab") {split($1,a,".") ; if (var[a[1]]=="") {var[a[1]]=$0} else {var[a[1]]=var[a[1]]"\t"$0}} else {print $1,$2,$3,$4,$5,var[$6]}}' final_dataset_to_append.tab SA_primary_secondary_no_HBV_for_hash.bed | awk 'BEGIN{FS="\t"; OFS="\t"}{split($6,a,"."); cols=$7 ; for (i=8;i<=NF;i++) {cols=cols"\t"$i} ; print $1,$2,$3,a[1]";"$4,$5,cols}' > SA_with_sequences.tab

echo -e $YELLOW "Intersecting results with the genome annotation" 

mkdir -p Final_Results
cp SA_primary_secondary_no_HBV.bed Final_Results/
cp SA_with_sequences.tab Final_Results/

cd Final_Results/

echo -e $LRED "### Collecting Human annotation per integration site ###" $Z
bedtools intersect -wao -a SA_primary_secondary_no_HBV.bed -b $genome/HG_HBV_rm_final.gff > SA_primary_secondary_no_HBV_annot_no_seq.bed

sortBed -i SA_primary_secondary_no_HBV_annot_no_seq.bed | bedtools merge -i - | bedtools intersect -wao -a - -b SA_primary_secondary_no_HBV_annot_no_seq.bed | cut -f1,2,3,7,8,9,10,11,12,13,15,17 | bedtools groupby -g 1,2,3,4,5  -c 8,9,10,11,12 -delim "|" -o distinct,distinct,distinct,distinct,distinct > SA_primary_secondary_no_HBV_annot_no_seq_grouped.bed

bedtools closest -D ref -a SA_primary_secondary_no_HBV_annot_no_seq_grouped.bed -b $genome/HG_HBV_rm_genes_sorted.gff > SA_primary_secondary_no_HBV_annot_no_seq_grouped_closest.bed

awk -F '\t' '{ print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 ":" $14 "-" $15 "\t" $17 "\t" $19 "\t" $20}' SA_primary_secondary_no_HBV_annot_no_seq_grouped_closest.bed | awk -F '\t' '{split($13,a,"Name="); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" a[2] "\t" $14}' | awk -F '\t' '{split($13,b,";"); print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" b[1] "\t" $14}' | awk -F '\t' '{if($1 != "$vir") print $0; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}'  | awk '!uniq[$1 FS $2 FS $3 FS $4 FS $4]++' > SA_primary_secondary_no_HBV_annot_no_seq_grouped_closest_filtered.bed
## Il seguente comando va rivisitato in modo da permettere di avere più righe per lo stesso evento di integrazione, quando ci stanno più risultati di blast per lo stesso evento (quindi quando non riusciamo a ricostruire un singolo consensus) ## modificare SA_primary_secondary_no_HBV_annot_no_seq_grouped_closest_filtered.bed così da avere il numero di righe equivalenti al numero delle vere integrazioni ## appendere 100.1 100.2 etc ##
awk 'BEGIN{FS="\t"; OFS="\t"}{if(FILENAME=="SA_with_sequences.tab") {var[$4]=$0} else {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,var[$4]}}' SA_with_sequences.tab SA_primary_secondary_no_HBV_annot_no_seq_grouped_closest_filtered.bed | awk 'BEGIN {FS="\t"; OFS="\t"}{cols=$20 ; for (i=21;i<=NF;i++) {cols=cols"\t"$i} ; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,cols}' > SA_semi_final.tab

#creo il file di intestazione della tabella finale che contiene le seguenti voci
#chr	start	end	reads	cov	features	start_features	end_features	strand_features	Descriptions	coordinates_closest_gene	strand_closest_gene	closest_gene_name	distance_chimera/closest_gene	HBV_chr	start_HBV	end_HBV	start_on_chimera	end_on_chimera	HG_chr	start_HG	end_HG	start_on_chimera	end_on_chimera	seq HBV	seq HG	chimeras
#nano header_SA_final
cp $genome/header_SA_final ./
cat header_SA_final SA_semi_final.tab > SA_final.tab

echo -e $YELLOW" ### Collecting viral annotation ## "

awk 'BEGIN{FS="\t"; OFS="\t"}{split($4,a,";"); {if($17 >= $16) {print $15,$16,$17,a[1]";cov="$5";seq="$27"\t.\t""+"} else {print $15,$16,$17,a[1]";cov="$5";seq="$27"\t.\t""-"}}}' SA_semi_final.tab > SA_viral_to_annotate.bed
awk 'BEGIN{FS="\t"; OFS="\t"}{if($3>$2) {print} else {print $1,$3,$2,$4,$5,$6}}' SA_viral_to_annotate.bed | sortBed -i - > SA_viral_to_annotate_modified.bed
bedtools intersect -wao -a SA_viral_to_annotate_modified.bed -b $genome/HBV_new_annotation.gff | cut -f1,2,3,4,5,6,9 | awk 'BEGIN{FS="\t"; OFS="\t"}{split($4,a,";"); print $0,a[2]}' | awk 'BEGIN{FS="\t"; OFS="\t"}{gsub(/cov=/,"",$8); print}' > SA_viral_annotate_to_group.bed

bedtools intersect -wao -a SA_viral_to_annotate_modified.bed -b $genome/HBV_new_annotation.gff | cut -f1,2,3,4,5,6,9 | bedtools groupby -g 4 -c 1,2,3,4,7 -o first,distinct,distinct,distinct,distinct -delim ";" | awk 'BEGIN{FS="\t"; OFS="\t"}{split($1,a,";"); print a[1], $2,$3,$4,$5,$6}' | awk 'BEGIN{FS="\t"; OFS="\t"}{split($5,a,";"); print $0, a[2]}' | awk 'BEGIN{FS="\t"; OFS="\t"}{gsub(/cov=/,"",$7); print $2,$3,$4,$1,$5,$6,$7}' | sortBed -i - |  bedtools groupby -g 6 -c 1,2,3,4,5,7 -o first,distinct,distinct,distinct,distinct,sum | sort -k1,1 | bedtools groupby -g 1 -c 2,3,4,5,6,7 -o first,distinct,distinct,distinct,distinct,sum  | sort -k1,1 -k7,7 > SA_viral_annotated_grouped.tab

bedtools intersect -wao -a SA_viral_to_annotate_modified.bed -b $genome/HBV_new_annotation.gff | cut -f1,2,3,4,5,6,9 | awk 'BEGIN{FS="\t"; OFS="\t"}{split($4,a,";"); print $1,$2,$3,$4,$5,$6,$7,a[2],a[1]}' | sort -k9,9n | bedtools groupby -g 9 -c 1,2,3,4,6,7,8 -o first,distinct,distinct,distinct,distinct,distinct,distinct | awk 'BEGIN{FS="\t"; OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$1}' > SA_viral_annotated_to_append_table_human.tab

cd ../
mv Final_Results ../../../../

echo -e $PURPLE "########## Congratulation!!! HBVIF COMPLETED ##########\n########## Thank you for using it ##########" $Z
