#include "./Matrix.h"
#include "./Model.h"
#include "./kseq.h"
#include <zlib.h>
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
#include <emscripten/bind.h> 

#include <nlohmann/json.hpp> // to store results as JSON string that can be passed to web worker
using json = nlohmann::ordered_json;

KSEQ_INIT(gzFile, gzread)

sdsl::csa_wt<> create_index(std::string filepath) 
{ 
    gzFile fp = gzopen(filepath.c_str(), "r"); // open the file handler
    if(fp == 0) {
        perror("Cannot open .fasta file");
        exit(1);
    }
    
    int l;
    std::string reference_seq;
    kseq_t *seq = kseq_init(fp);    // initialize seq
    while ((l = kseq_read(seq)) >= 0) // read sequence
    {
        reference_seq += seq->seq.s;
    }
    kseq_destroy(seq);  // destroy seq and fp objects
    gzclose(fp);

    if(reference_seq.length()<1500000||reference_seq.length()>3000000){abort();} //check if length of sequence is reasonable

    sdsl::csa_wt<> ref_index;   //initialise fm-index

    auto start_fm = std::chrono::steady_clock::now(); // start timer

    sdsl::construct_im(ref_index, reference_seq, 1); // generate fm-index
    
    auto end_fm = std::chrono::steady_clock::now(); //end timer
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

std::string invert(std::string seq)
{
    transform(
        begin(seq),
        end(seq),
        begin(seq),
        complement);
    reverse(seq.begin(), seq.end());
    return seq;
}

Numeric_lib::Matrix<double,1> lookup_unitigs(const sdsl::csa_wt<>& fm_index, const std::vector<std::string>& unitigs)
{
    Numeric_lib::Matrix<double,1> pa_mat{int(unitigs.size())}; // Matrix to store presence/absence
    for(int i=0; i<unitigs.size(); ++i){
        if (unitigs[i]=="Intercept"){pa_mat[i]=1.0;} //start each vector with 1 for the Intercept
        else{
            int c = count(fm_index,unitigs[i]);
            if (c==0) {
                int c2 = count(fm_index,invert(unitigs[i]));//check reverse sequence
                if (c2==0){pa_mat[i]=0.0;}
                else {pa_mat[i]=1.0;}
            }
            else {pa_mat[i]=1.0;}
        } 
    }
    return pa_mat;
}


std::string make_prediction(std::string assembly_filename)
{
    auto start = std::chrono::steady_clock::now(); // start timer

    sdsl::csa_wt<> fm_index = create_index(assembly_filename); // create fm-Index from .fasta file

    std::vector<std::string> antibiotics ={"Penicillin","Chloramphenicol","Erythromycin","Tetracycline", "Trim_sulfa"};
    
    json results_json;
    results_json["filename"] = assembly_filename.erase (0,9);

    for(std::string const &antibiotic:antibiotics){ 

        Model thismodel;
        thismodel.read_model(antibiotic);  // create Model containing coefs and unitigs

        Numeric_lib::Matrix<double,1> pa_mat = lookup_unitigs(fm_index, thismodel.getUnitigs()); // get vector of presence/absence of unitigs

        double probability = round(thismodel.get_prob(pa_mat)*1000.0)/1000.0; // calculate prob of resistance

        std::cout<<"Probability of resistance to "<<antibiotic<<": "<<probability<<'\n'; 

        results_json[antibiotic] = probability;
        
    }  
    auto fm_time = std::chrono::steady_clock::now();    //end timer
    std::chrono::duration<double> elapsed_seconds = fm_time-start;
    std::cout << "elapsed total time: " << elapsed_seconds.count() << "s\n";

    std::string result = results_json.dump();
    return result;
}

/*
int main()
{
    std::string id_file = "files/mass_assemblies_ids_short.txt";
    std::ifstream ist {id_file};
    if(!ist) perror("Can't open file with strain IDs");   
    std::string id;
    while(ist>>id)      // loop through all IDs in .txt file to make predictions for each
    {
        std::string assembly_filename = "files/assemblies_mass_short/"+id+".fa";
        std::cout<<'\n'<<assembly_filename<<'\n';
        make_prediction(assembly_filename);
    }
}
*/

EMSCRIPTEN_BINDINGS(my_module) {        // include bindings for emscripten
    emscripten::function("make_prediction", &make_prediction);
    emscripten::register_map<std::string,double>("map<string,double>");
}