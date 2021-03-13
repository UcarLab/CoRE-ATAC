#!/bin/bash

ARGUMENTS=() 
KEEPDUPS=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --keep-duplicates) KEEPDUPS=" --keep-duplicates"; shift ;;
        --*) echo "Unknown option: $1"; exit 1 ;;
        *)  ARGUMENTS+=("$1"); shift ;;
    esac
done

if [[ ${#ARGUMENTS[@]} -ne 6 ]]; then
    echo "Usage: CoRE-ATACFeatureExtraction-singularity.sh bamfile peakfile genome fastapath fastasplitdirectory outdir"
    echo "bamfile: Path to the bam file."
    echo "peakfile: Path to the peak file."
    echo "genome: Genome used by HOMER and conservation scores (i.e., hg19, hg38, mm10)."
    echo "Note: For other genomes, please add conservation and filter files for the genome to /PEAS/extraction_files/ in the singularity directory."
    echo "fastapath: Path to the .fa file for the genome."
    echo "fastasplitdirectory: Path to directory containing chromosome fasta (.fa) files (i.e., chr1.fa, chr2.fa etc)."
    echo "outdir: The output directory to save output files."
    echo "--keep-duplicates: Provide this option to keep duplicates."
    exit 0
fi

BAMPATH=${ARGUMENTS[0]}
PEAKPATH=${ARGUMENTS[1]}
GENOME=${ARGUMENTS[2]}
FASTA=${ARGUMENTS[3]}
FASTASPLIT=${ARGUMENTS[4]}
OUTDIR=${ARGUMENTS[5]}

/CoRE-ATAC/FeatureExtractor.sh ${BAMPATH} ${PEAKPATH} ${OUTDIR} ${GENOME} ${FASTA} ${FASTASPLIT} /CoRE-ATAC/${KEEPDUPS} 

