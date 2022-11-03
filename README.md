# Backend code for https://amr.poppunk.net/

## Overview

https://amr.poppunk.net/ is a tool to predict antimicrobial resistance in __*S.pneumoniae*__. It applies linear models on genetic sequence data to return a probability of resistance to five antibiotics. The code for the website can be found [here](https://github.com/bacpop/AMR_ReactApp).

The backend code takes a FASTA file as input, creates an FM-index of the genetic sequence and queries specific nucleotide sequences ([unitigs](https://pubmed.ncbi.nlm.nih.gov/30419019/#&gid=article-figures&pid=fig-1-uid-0)), which are relevant for prediction. 

Linear models are applied to calculate probabilities of resistance to the respective antibiotic. 

Coefficients for the models have been derived from logistic [ElasticNet](https://en.wikipedia.org/wiki/Elastic_net_regularization) models trained on data from the USA and South Africa from the [GPS](https://www.pneumogen.net/gps/) database.

## Testing

Example FASTA-files are provided in `files/fa_examples`.

The GitHub workflow C/C++ CI will run a test comparing the generated results to the expected results stored in `test_result.txt`.

To run tests locally, download the repository with
```
git clone git@github.com:bacpop/AMR_prediction.git
```
From the `/src` folder run the following command to compile the script and run the tests:
```
make
make check
```

## Running on the command line

If you want to run predictions on a set of fasta files from the command line, follow these steps:

Download the repository with
```
git clone git@github.com:bacpop/AMR_prediction.git
```
From the `/src` folder run the following command to compile the script:
```
make amr_for_CLI
```
You'll need a .txt file listing all files you want to process including their path from your current directory (src). An example can be seen in `file_list.txt`, which lists all .fasta files that are stored in the folder `/files/fa-examples`. 

To run the prediction model, run
```
./amr_for_CLI file_list.txt 
```

Results will be stored in a csv file.

## Compiling into WebAssembly
To compile into WebAssembly, you need [emscripten] (https://emscripten.org/docs/getting_started/downloads.html) installed on your machine.

Download the repository with
```
git clone git@github.com:bacpop/AMR_prediction.git
```
Create a new empty folder in the root of the directory called `build`. This is the place where you will find the output files.
From `/src` run the following command to compile the script into WebAssembly:
```
make web
```