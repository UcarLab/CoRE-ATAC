#!/bin/bash
args=()
((index=0))
for i in "$@"
do
    args[${index}]="$i"
    ((index++))
done

inDir="${args[0]}"
prefix="${args[1]}"
outDir="${args[2]}"
fasta="${args[3]}"
homerref="${args[4]}"
homermotifs="${args[5]}"
conservation="${args[6]}"
ctcfmotifs="${args[7]}"

inputpeaks="${args[8]}"

cd "${outDir}"

mkdir peak_features
cd peak_features

jarpath="${args[9]}/"

keepdups="${args[10]}"

##Start with peaks and extract features.
##Sort BAM  #TODO add option to skip this step if already sorted
echo "--- Sorting bam file. ---"
#echo "${inDir}/${prefix}.bam"
samtools sort  -T PEASEXTRACT_${prefix} -o ${prefix}_sorted.bam "${inDir}/${prefix}.bam"
samtools index ${prefix}_sorted.bam

#cp "${inDir}/${prefix}.bam" ${prefix}_sorted.bam

#

echo "--- Calling annotations & known motifs. ---"
#HOMER Annotations
annotatePeaks.pl "${inputpeaks}" "${homerref}" -m "${homermotifs}" -nmotifs > ${prefix}_peaks_annotated.bed

#call denovo motifs
echo "--- Calling denovo motifs. ---"
findMotifsGenome.pl "${inputpeaks}" "${fasta}" "${outDir}/denovo"

mkdir "${outDir}/denovo/merge"
cp "${outDir}/denovo/homerResults/"*.motif "${outDir}/denovo/merge"
rm "${outDir}/denovo/merge/"*.similar*
rm "${outDir}/denovo/merge/"*RV.motif
cat "${outDir}/denovo/merge/"*.motif >> "${outDir}/denovo/merge/merged.motifs"
#call motifs with homer again using denovo motifs file homerMotifs.all.motifs
annotatePeaks.pl "${inputpeaks}" "${homerref}" -m "${outDir}/denovo/merge/merged.motifs" -nmotifs > ${prefix}_peaks_denovo.bed

echo "--- Calling CTCF motifs. ---"
annotatePeaks.pl "${inputpeaks}" "${homerref}" -m "${ctcfmotifs}" -nmotifs > ${prefix}_peaks_ctcf.bed

#Get the insert size threshold to remove outlier inserts
echo "--- Getting insert size threshold. ---"
java -jar "${jarpath}PEASTools.jar" insertsizethresh "${prefix}_sorted.bam" "${outDir}/peak_features" ${keepdups}
thresh=$(cat "thresh.txt")


#Get Insert features
echo "--- Getting insert features. ---"
for i in {1..22}
do
    chr=chr$i
    java -jar "${jarpath}PEASTools.jar" insertmetrics "${chr}" "${chr}.bam" "${inputpeaks}" "${prefix}_${chr}_insertmetrics.txt" "$thresh" ${keepdups}
rm ${chr}.bam
    cat ${prefix}_${chr}_insertmetrics.txt >> ${prefix}_insertmetrics.txt
rm "${prefix}_${chr}_insertmetrics.txt"
done

echo "--- Getting conservation scores. ---"
#Get Conservation Scores
java -jar "${jarpath}PEASTools.jar" conservation "${inputpeaks}" "${conservation}" "${prefix}_conservation.txt"

echo "--- Merging features. ---"
java -jar "${jarpath}PEASTools.jar" mergedl "${inputpeaks}" "${prefix}_peaks_annotated.bed" "${prefix}_insertmetrics.txt" "${prefix}_conservation.txt" "${prefix}_peaks_denovo.bed" "${prefix}_peaks_ctcf.bed" "${prefix}_features.txt" "MERGED"