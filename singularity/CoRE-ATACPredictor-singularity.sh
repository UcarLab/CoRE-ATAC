#!/bin/bash

#ARGUMENTS=() 
#while [[ "$#" -gt 0 ]]; do
#    ARGUMENTS+=("$1"); shift ;
#done

#if [[ ${#ARGUMENTS[@]} -ne 4 ]]; then
#    echo "Usage: CoRE-ATACPredictor-singularity.sh featuredirectory model outputfile"
#    echo "featuredirectory: The feature extraction output directory."
#    echo "basname: The prefix used for the input files in the feature directory."
#    echo "model: The trained model (.h5 file)."
#    echo "outputfile: The output file for predictions."
#    echo ""
#    echo "Optional Arguments:"
#    echo "--channelslast - Set the channels to last column. (Default = first)"
#    echo "Optional Arguments:"
#    echo "Optional Arguments:"
#
#    exit 0
#fi

#FEATUREDIR=${ARGUMENTS[0]}
#BASENAME=${ARGUMENTS[1]}
#MODEL=${ARGUMENTS[2]}
#OUTPUTFILE=${ARGUMENTS[3]}

python3 /CoRE-ATAC/CoREATAC_PredictionTool.py "$@"

