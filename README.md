# Backend code for https://amr.poppunk.net/

## Overview

https://amr.poppunk.net/ is a tool to predict antimicrobial resistance in __*S.pneumoniae*__. It applies linear models on genetic sequence data to return a probability of resistance to five antibiotics. The code for the website can be found [here](https://github.com/bacpop/AMR_ReactApp).

The backend code takes a FASTA file as input, creates an FM-index of the genetic sequence and queries specific nucleotide sequences ([unitigs](https://pubmed.ncbi.nlm.nih.gov/30419019/#&gid=article-figures&pid=fig-1-uid-0)), which are relevant for prediction. 

Linear models are applied to calculate probabilities of resistance to the respective antibiotic. 

Coefficients for the models have been derived from logistic [ElasticNet](https://en.wikipedia.org/wiki/Elastic_net_regularization) models trained on data from the USA and South Africa from the [GPS](https://www.pneumogen.net/gps/) database.

## Testing

Two example FASTA-files are provided in `files/fa_examples`.

The GitHub workflow C/C++ CI will run a test comparing the generated results to the expected results stored in `test_result.txt`.
