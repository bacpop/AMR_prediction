#include "../Matrix.h"
#include "../kseq.h"
#include "../zlib.h"
#include <sdsl/suffix_arrays.hpp>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>                                                                                                                                                                                         
#include<algorithm>
#include<cassert>
#include <chrono>

KSEQ_INIT(gzFile, gzread)

using namespace std;

void read_model(string antibiotic, vector<string>& unitigs, vector<double>& coefs)
{
    ostringstream filename; //set filepath and open as istream
    filename<<"/mnt/c/Users/Marie/Documents/Imperial/ResearchAssistant/Model_coefficients/"<<antibiotic<<"_coefficients.txt";
    ifstream ist {filename.str()};
    if(!ist) perror("Can't open model input file ");

    string unitig;  
    float coef;
    while (ist>>unitig>>coef)   //fill data in 2 vectors
    {   
        unitigs.push_back(unitig);
        coefs.push_back(coef);
    }
}

sdsl::csa_wt<> create_index(string filepath) 
{ 
    string reference_seq;

    // open the file handler
    gzFile fp = gzopen(filepath.c_str(), "r");
    if(fp == 0) {
        perror("Cannot open .fasta file");
        exit(1);
    }

    // initialize seq
    kseq_t *seq = kseq_init(fp);

    // read sequence
    int l;
    while ((l = kseq_read(seq)) >= 0)
    {
        reference_seq += seq->seq.s;
    }

    // destroy seq and fp objects
    kseq_destroy(seq);
    gzclose(fp);

    //initialise & generate fm-index
    sdsl::csa_wt<> ref_index;

    auto start_fm = std::chrono::steady_clock::now();

    sdsl::construct_im(ref_index, reference_seq, 1);
    
    auto end_fm = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds_fm = end_fm-start_fm;
    std::cout << "elapsed time for creating fm-index: " << elapsed_seconds_fm.count() << "s\n";

    

    return ref_index;
}

char complement(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
    return ' ';
}

string invert(string seq){
    
    transform(
        begin(seq),
        end(seq),
        begin(seq),
        complement);
    reverse(seq.begin(), seq.end());
    return seq;
}

vector<double> lookup_unitigs(sdsl::csa_wt<> fm_index, vector<string> unitigs)
{
    vector<double> pa; //vector to store presence/absence
    for(string i:unitigs){
        if (i=="Intercept"){pa.push_back(1.0);} //start each vector with 1 for the Intercept
        else{
            int c = count(fm_index,i);
            if (c==0) {
                int c2 = count(fm_index,invert(i));//check reverse sequence
                if (c2==0){pa.push_back(0.0);}
                else {pa.push_back(1.0);}
            }
            else {pa.push_back(1.0);}
        }
        
    }
    return pa;
}



void main2(string assembly_filename)
{
    auto start = std::chrono::steady_clock::now();

    sdsl::csa_wt<> fm_index = create_index(assembly_filename); // create fm-Index from .fasta file

    vector<string> antibiotics ={"Penicillin","Chloramphenicol","Erythromycin","Tetracycline", "Trim_sulfa"};

    for(string antibiotic:antibiotics){

        vector<string> unitigs; //store data from model
        vector<double> coefs;
        read_model(antibiotic, unitigs, coefs);  //fill vectors
        Numeric_lib::Matrix<double,1>coefs_mat{&coefs[0], int(coefs.size())}; // translate into Matrix

        vector<double> pa = lookup_unitigs(fm_index, unitigs); // get vector of presence/absence of unitigs
        Numeric_lib::Matrix<double,1>pa_mat{&pa[0], int(pa.size())}; // translate into Matrix for calculations

        double pred = Numeric_lib::dot_product(pa_mat,coefs_mat); // calculate dot product of both Matrices
        double prob =1/(1+exp(-pred)); // apply sigmoid function for probability of resistance

        cout<<"Probability of resistance to "<<antibiotic<<": "<<prob<<'\n';

    }  

    auto fm_time = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = fm_time-start;
    std::cout << "elapsed total time: " << elapsed_seconds.count() << "s\n";

}

int main()
{
    string id_file = "/mnt/c/Users/Marie/Documents/Imperial/ResearchAssistant/mass_assemblies_ids_short.txt";
    ifstream ist {id_file};
    if(!ist) perror("Can't open file with strain IDs");   

    string id;
    while(ist>>id)
    {
        string assembly_filename = "/mnt/c/Users/Marie/Documents/Imperial/ResearchAssistant/assemblies_mass/"+id+".fa";
        cout<<'\n'<<assembly_filename<<'\n';
        main2(assembly_filename);
    }
    
   
}


// g++ -I ~/include -L ~/lib -g amr_pred.cpp -lz -o amr_pred -lsdsl -ldivsufsort -ldivsufsort64