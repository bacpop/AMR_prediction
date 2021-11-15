# AMR_prediction

## Structure of amr_pred.cpp:

void **read_model** (string antibiotic, vector<string\>& unitigs, vector<double\>& coefs)
{ 
  read .txt file and store unitigs & coefficients in 2 vectors
}
  
sds::csa_wt<\> **create_index**(string filepath)
{
  open and read fasta file
  generate fm-index
}
  
char **complement**(char n)
{
  change bases to complementary base
}
 
string **invert**(string seq)
{
  **complement()** nucleotides, then reverse order
}
 
vector<double\> **lookup_unitigs**(sdsl::csa_wt<\> fm_index, vector<string\> unitigs)
{
  for-loop through vector of unitigs:
      add 1 for intercept and present unitigs, 0 otherwise
      check also reverse sequence with **invert()**
}
 
void **main**()
{
  **create_index()** 
  for-loop through vector of antibiotics:
      **read_model()**
      **lookup_unitigs()**
      calculate dot_product & apply sigmoid function
      print results
}
 
