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
fasta="${args[3]}"
sepfasta="${args[4]}"
path="${args[5]}/"

prefix=$(basename ${bamfile} | sed 's/.bam//g')
inDir=$(dirname ${bamfile})
reformattedpeaks=${outdir}${prefix}_peaks.txt


INSERTTHRESH=900
CHARPERLINE=50
HOMERREF="hg38"
HOMERMOTIFS=${path}PEAS/humantop_Nov2016_HOMER.motifs
CONSERVATION=${path}PEAS/phastCons30way_hg38.bed
CTCFMOTIFS=${path}PEAS/CTCF.motifs


cd "${outDir}"

#Step 1: Reformat peaks for 600 bp windows
python ${path}PeakFormatter.py ${peakfile} ${sepfasta} ${reformattedpeaks}


#Step 2: Extract features from BAM
java -jar ${path}BAMExtractor.jar ${bamfile} "${reformattedpeaks}" ${sepfasta} ${CHARPERLINE} ${outdir} ${prefix} ${INSERTTHRESH} --keepduplicates #!!!TODO!!!#


#Step 3: Extract PEAS features
${path}DeepLearningPEAS.sh ${inDir} ${prefix} ${outdir} ${fasta} ${HOMERREF} ${HOMERMOTIFS} ${CONSERVATION} ${CTCFMOTIFS} "${reformattedpeaks}_original.txt" ${path}PEAS


#Step 4: Move PEAS features to directory
#mv ${outDir}/peak_features/${prefix}_features.txt ${outDir}

#Done!

