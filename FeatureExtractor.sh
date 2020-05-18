#!/bin/bash

args=()
((index=0))
for i in "$@"
do
    args[${index}]="$i"
    ((index++))
done

bamfile="${args[0]}"
peakfile="${args[1]}"
outdir="${args[2]}"
refgenome="${args[3]}"
fasta="${args[4]}"
sepfasta="${args[5]}"
path="${args[6]}/"

keepdups=""
if ${args[7]} == "TRUE"
then
	keepdups="--keepduplicates"
fi

prefix=$(basename ${bamfile} | sed 's/.bam//g')
inDir=$(dirname ${bamfile})
reformattedpeaks=${outdir}${prefix}_peaks.txt


INSERTTHRESH=900
CHARPERLINE=50
HOMERREF=${refgenome}
HOMERMOTIFS=${path}PEAS/humantop_Nov2016_HOMER.motifs
CONSERVATION=${path}PEAS/phastCons46wayPlacental.bed
CTCFMOTIFS=${path}PEAS/CTCF.motifs


cd "${outDir}"

#Step 1: Reformat peaks for 600 bp windows
python ${path}PeakFormatter.py ${peakfile} ${sepfasta} ${reformattedpeaks}


#Step 2: Extract features from BAM
java -jar ${path}BAMExtractor.jar ${bamfile} "${reformattedpeaks}" ${sepfasta} ${CHARPERLINE} ${outdir} ${prefix} ${INSERTTHRESH} ${keepdups}


#Step 3: Extract PEAS features
${path}DeepLearningPEAS.sh ${inDir} ${prefix} ${outdir} ${fasta} ${HOMERREF} ${HOMERMOTIFS} ${CONSERVATION} ${CTCFMOTIFS} "${reformattedpeaks}_original.txt" ${path}PEAS ${keepdups}


#Step 4: Move PEAS features to directory
#mv ${outDir}/peak_features/${prefix}_features.txt ${outDir}

#Done!

