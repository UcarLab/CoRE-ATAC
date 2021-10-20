# CoRE-ATAC #

Classification of Regulatory Elements with ATAC-seq (CoRE-ATAC).

CoRE-ATAC is split into 

1. Feature Extraction/Prediction
2. Model Training. 

CoRE-ATAC makes use of features extracted by PEAS[1], using a modified version of the original source code available from https://github.com/UcarLab/PEAS

[1] Thibodeau, A., Uyar, A., Khetan, S. et al. A neural network based model effectively predicts enhancers from clinical ATAC-seq samples. Sci Rep **8**, 16048 (2018). https://doi.org/10.1038/s41598-018-34420-9

# Feature Extraction/Prediction #

For feature extraction, we have provided a singularity image file in releases that users can immediately load into their systems with minimal set up. 

### System Requirements/Recommendations: ###

1. System with Singularity installed: https://sylabs.io/
2. GPU required?: No GPU is required. Both feature extraction and predictions use CPU. 
3. Memory: We recommend at least 16gb of memory.

Note: For GPU use, a modified singularity image will need to be created, pulling from the corresponding [NVIDIA docker image](https://ngc.nvidia.com/catalog/containers/nvidia:tensorflow/tags).

### Installation ###

To install the feature extraction and prediction framework, download the latest FeaturePredictor singularity image file (.sif) from [releases](https://github.com/UcarLab/CoRE-ATAC/releases).

Once the image is downloaded, following the following steps to complete the installation:

**Step 1: Convert the image file into a sandbox directory:**

`singularity build --sandbox ./CoRE-ATAC-FeaturePredictor-hg19/ CoRE-ATAC-FeaturePredictor.sif`

**Step 2: Complete the HOMER installation:**

Enter the sandbox directory:

`singularity shell --writable ./CoRE-ATAC-FeaturePredictor-hg19/`

and install the hg19 reference for HOMER

`perl /HOMER/configureHomer.pl -install hg19`


**Step 3: Download the hg19 .fa files from UCSC.**

1. Download the hg19 reference: `ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz`
2. Download the chromosome files: `ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/`

Only the chromosomes that will be used for feature extraction are required 

**Step 4: Create a list of chromosome fasta files for future reference.**

The format of this file is a list for the path of each chromosome reference, where each line contains:

1. The chromosome name
2. The absolute path to the zipped fasta file.

separated by tabs

Example:

line1: `chr1    /hg19/chr1.fa.gz`

line2: `chr2    /hg19/chr2.fa.gz`

---

After these steps, you should now have a directory that can be used with singularity and the necessary files to run both feature extraction and model predictions!


### Running Feature Extraction ###

To extract feaures, use the singularity sandbox with the following code, replacing arguments and paths as necessary:

`singularity exec ./CoRE-ATAC-FeaturePredictor-hg19/ /CoRE-ATAC/CoRE-ATACFeatureExtraction-singularity.sh <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>`

**Arg1:** The path to the ATAC-seq alignment file (.bam)

**Arg2:** The path to the ATAC-seq peak file (3 column tab delimited file for chr, start, end)

**Arg3:** The reference genome. (e.g., hg19). This is used to inform HOMER.

**Arg4:** The path to the reference genome fasta file (e.g., hg19.fa)

**Arg5:** The path to the file specifying the reference genome chromosomes as noted in the installation section.

**Arg6:** The output directory. Note: This directory should already be created.

**IMPORTANT: The full path must be used for the output directory!**

After running, multiple feature files will be generated in the specific output directory.


### Predicting *cis*-REs with CoRE-ATAC ###

To predict Promoters, Enhancers, and Inuslators with CoRE-ATAC, you will need:

1. The singularity image/sandbox directory
2. The extracted features after running feature extraction.
3. The model file.

For the extracted features, we refer to the previous section.

The CoRE-ATAC model file can be downloaded from: 

[https://github.com/UcarLab/CoRE-ATAC/releases/download/CoRE-ATAC_Model_CPU/CoRE-ATAC_CPUCompatible.h5.zip](https://github.com/UcarLab/CoRE-ATAC/releases/download/CoRE-ATAC_Model_CPU/CoRE-ATAC_CPUCompatible.h5.zip)

Unzip this file to make it readable by the prediction tool.

With the singularity image or sandbox directory, directory of features extracted from the previous step and the model file, all of the inputs are now available to predict *cis*-REs!

To predict *cis*-REs (i.e., Pormoters, Insulators, Enhancers) run singularity with the following code, replacing the arguments and paths as necessary:

`singularity exec ./CoRE-ATAC-FeaturePredictor-hg19/ /CoRE-ATAC/CoRE-ATACPredictor-singularity.sh <arg1> <arg2> <arg3> <arg4>`

**Arg1:** The path to the feature extraction directory.

**Arg2:** The prefix used the feature extraction files. This can be identified by looking at the filenames before files such as "x_cutmatrices.txt". For example if the name is MCF7_cutmatrices.txt, the prefix will be MCF7.

**Arg3:** The model file (i.e., CoRE-ATAC_CPUCompatible.h5)

**Arg4:** The output file (e.g., Predictions.txt)

After running this code, predictions will be available in the output specified by Arg4.

The output file is a tab delimited text file with the following columns:

1. The chromosome
2. The start position
3. The end position
4. Promoter prediction probability
5. Enhancer prediction probability
6. Insulator prediction probability
7. Other prediction probability

# Model Training #

We do not provide image files for training models. However, we do provide a [definition file](https://github.com/UcarLab/CoRE-ATAC/blob/master/singularity/CoRE-ATAC-ModelTrainer.def) template and the necessary releases for building a singularity image. Due to different architectures and CUDA installations, it is important that the correct tensorflow GPU image is used.

### Requirements ###

**GPU:** A tensorflow compatible GPU (i.e., NVIDIA) is required to train models. Due to the amount of data used to train these models, we do not recommend training models via CPU. 

**Memory:** For training a model from scratch, it is necessary to have enough memory for all of the training data. For reference: CoRE-ATAC's model was trained using 128gb of memory on a high performance computing cluster.


### Building the singularity image file for training models ###

CoRE-ATAC has been tested using the (https://github.com/UcarLab/CoRE-ATAC/blob/master/singularity/CoRE-ATAC-ModelTrainer.def)[definition file] provided in the singularity directory. Due to size limits on release, we do not provide the built singularity image file on github.

However, it may be necessary to build a new singularity image that better fits the computing environment. For example, CUDA installations require earlier versions of tensorflow.

To build the singularity image file, first identify the NVIDIA docker image](https://ngc.nvidia.com/catalog/containers/nvidia:tensorflow/tags) that best suits your system.

Next, update the singularity definition file, available from: (https://github.com/UcarLab/CoRE-ATAC/blob/master/singularity/CoRE-ATAC-ModelTrainer.def)[https://github.com/UcarLab/CoRE-ATAC/blob/master/singularity/CoRE-ATAC-ModelTrainer.def]

Edit this file, replacing: 

`nvcr.io/nvidia/tensorflow:19.01-py3`

with the appropriate tensorflow docker image


Next, build the singularity image:

`sudo singularity build CoRE-ATACModelTrainer.sif CoRE-ATAC-ModelTrainer.def` 

The singularity image should now be built and ready to use to train models!


### Training Models ###

Once the image file is downloaded/created, we can now use this image to train models. For this follow the following steps:

**Step 1: Create a file listing the feature directories.**

Create a list of feature directories for each sample to be included in model training.

For example:

line 1: `/path/to/MCF7/`
line 2: `/path/to/K562/`

**Step 2: Create a file listing the base names.**

Create a list of base names for features corresponding to the directories listed in step 1. These can be identified by looking at the filenames before files such as "x_cutmatrices.txt". For example if the name is MCF7_cutmatrices.txt, the base name prefix will be MCF7.

For example:

line 1: `MCF7`
line 2: `K562`

**Step 3: Create a list of PEAS features.** These are the PEAS feature extracted during feature extraction. These are located in the `peak_features` directory of the feature extraction directory with the `_features.txt` suffix.

For example:

line 1: `/path/to/MCF7/peak_features/MCF7_features.txt`
line 2: `/path/to/K562/peak_features/K562_features.txt`

**Step 4: Specify the training, validation, and test chromosomes**
Create a file specifying the train, validation, and text chromosomes. This is 2 column tab delimited file where the first column specifies the category (train, val, or test) and the second specifies the chromosomes (separated by commas)

line 1: `train	chr1,chr4,chr5,chr6,chr7,chr8,chr9,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22`
line 2: `val	chr2,chr10`
line 3: `test	chr3,chr11`
 
 
With all 4 of these files, we are now ready to train models!

Training CoRE-ATAC is a 3 step process which first trains CoRE-ATAC and PEAS features separately and then merges them for a final round of training.

**Step 5: Train Signal and Sequence**
To train the signal and sequence model run the following:

`singularity exec --nv ./CoRE-ATAC-ModelTrainer.sif /CoRE-ATAC/CoRE-ATACFeatureExtraction-singularity.sh <arg1> <arg2> <arg3> <arg4>`

**Arg1:** The path of the file listing the base names from Step 2.

**Arg2:** The path of the file listing the feature directories from Step 1.

**Arg3:** The path of the file listing the train, val, and test chromosomes from Step 4.

**Arg4:** The file path of the output model (e.g., "MyModel-SigSeq.h5")




**Step 6: Train PEAS** To train the PEAS model run the following:

`singularity exec --nv ./CoRE-ATAC-ModelTrainer.sif /CoRE-ATAC/CoRE-ATACFeatureExtraction-singularity.sh <arg1> <arg2> <arg3> <arg4> <arg5>`

**Arg1:** The path of the file listing the base names from Step 2.

**Arg2:** The path of the file listing the feature directories from Step 1.

**Arg3:** The path of the file listing the PEAS features from Step 3.

**Arg4:** The path of the file listing the train, val, and test chromosomes from Step 4.

**Arg5:** The file path of the output model (e.g., "MyModel-PEAS.h5")



**Step 7: Train the combined model**

To train the PEAS model run the following:

`singularity exec --nv ./CoRE-ATAC-ModelTrainer.sif /CoRE-ATAC/CoRE-ATACFeatureExtraction-singularity.sh <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7>`

**Arg1:** The path of the file listing the base names from Step 2.

**Arg2:** The path of the file listing the feature directories from Step 1.

**Arg3:** The path of the file listing the PEAS features from Step 3.

**Arg4:** The path of the file listing the train, val, and test chromosomes from Step 4.

**Arg5:** The file path of sig-seq model previously trained (e.g., "MyModel-SigSeq.h5")

**Arg6:** The file path of PEAS model previously trained (e.g., "MyModel-PEAS.h5")

**Arg7:** The file path of the output model (e.g., "MyModel-Merged.h5")


This method is used to prevent component overfitting. For example, signal and sequence models may need more time for training, meanwhile PEAS components need less time. To prevent overfitting on PEAS features, these components are trained independently. After merging, both models should be at a point close to overfitting. Therefore, we only want to train on 1-5 epochs so that the model can learn how to integrate the features, but not begin overfitting them.


# Feature Extraction Output Files #

Note: 600bp is the default used. These tools can be used for arbitrary window sizes.

**\_peaks.txt:** A table delimited file containing the original and 600bp window peak locations. The 600bp window peak location is defined by the center of the original peak provided to CoRE-ATAC. The columns are defined as:

1. Chromosome
2. 600bp start
3. 600bp end
4. Peak index
5. Original peak start
6. Original peak end

**\_cutpileups.txt:** Cut pileups for each 600bp peak window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated integer values corresponding to the number of cut sites (5' or 3' read locations) within the 600bp window.

**\_forwardcutpileups.txt:** Forward strand cut pileups for each 600bp peak window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated integer values corresponding to the number of forward strand cut sites (5' read locations) within the 600bp window.

**\_reversecutpileups.txt:** Reverse strand cut pileups for each 600bp peak window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated integer values corresponding to the number of reverse strand cut sites (3' read locations) within the 600bp window.

**\_insertsizes.txt:** The insert sizes of reads overlapping the 600bp window. Each peak is defined and follwed by a list of integers.

Peak definition line: chromosome, 600bp start, 600bp end, peak index
each peak definition line is follwed by a single integer per line until reaching the next peak definition or the end of the file.

**\_insertpileups.txt:** Insert pileups for each 600bp peak window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated integer values corresponding to the number of reads observed at each base in the 600bp window.

**\_forwardmedianinsert.txt:** The median insert size of reads with forward cut sites at the current bp in the 600bp window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated float values corresponding to the median insert size of forward cut sites.

**\_reversemedianinsert.txt:** The median insert size of reads with reverse cut sites at the current bp in the 600bp window. Every two lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Comma separated float values corresponding to the median insert size of reverse cut sites.

**\_sequencefreq.txt:** The DNA sequence frequency. Every 6 lines contains the following information:

1st line: Tab delimited information regarding the peak location: chromosome, 600bp start, 600bp end, peak index
2nd line: Reference DNA sequence.
3rd line: Frequency of adenine (A) for each bp in the 600bp window.
4th line: Frequency of cytosine (C) for each bp in the 600bp window.
5th line: Frequency of guanine (G) for each bp in the 600bp window.
6th line: Frequency of thymine (T) for each bp in the 600bp window.

Note: One hot encoding of the reference is used when < 10 read sequences overlap the position.

**\_peaks.txt_original.txt:** Tab delimited file of the original peak locations with index information. This is used for extracting features from PEAS.
1. Chromosome
2. Original peak start
3. Original peak end
4. Peak index

**\_cutmatrices.txt:** Unused. This file encodes the cut sites of paired end reads as a sparse matrix. This was tested, but did not produce a valuable model on its own. This file is formatted to first specify the peak location on the first line, specify the number of sparse matrix entries on the second line, and finally specify wher each cut pair is located in the subsequent lines equal to the number specified by the 2nd line. This process repeeats until the end of the file.

All other files correspond to outputs defined by PEAS.

