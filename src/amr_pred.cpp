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

#ifdef WEB
#include <emscripten/bind.h> 
#endif

/// to store results as JSON string that can be passed to web worker
#include "./json.hpp" 
using json = nlohmann::ordered_json;

KSEQ_INIT(gzFile, gzread)

/**
 * Generate an fm-index from a .fasta file sequence.
 * The fm-index provides a quick way to search for unitigs (shorter sequences)
 * within the overall sequence.
 * Also, a length check is performed at this stage to check whether the overall
 * sequence has a sensible length for S. pneumoniae (currently very generously,
 * with an accepted range between 1.5 mio and 3 mio basepairs, usually S. pneumoniae
 * genomes have a length ~2.1 mio basepairs.)
 * 
 * @param filepath the path to the .fasta file holding the sequence
 * 
 * @returns pair holding the index and a bool value stating the lenght check outcome
 */
std::pair<sdsl::csa_wt<>,bool> create_index(std::string filepath) 
{ 
    /// open the file handler
    gzFile fp = gzopen(filepath.c_str(), "r");
    if(fp == 0) {
        perror("Cannot open .fasta file");
        exit(1);
    }
    int l;
    std::string reference_seq;
    /// initialize seq
    kseq_t *seq = kseq_init(fp);
    /// read sequence
    while ((l = kseq_read(seq)) >= 0)
    {
        reference_seq += seq->seq.s;
    }
    // destroy seq and fp objects
    kseq_destroy(seq);
    gzclose(fp);
    /// perform length check
    bool lengthcheck;
    if(reference_seq.length()<1500000||reference_seq.length()>3000000){
        lengthcheck=false;
    } else {
        lengthcheck=true;
    }
    /// translate sequence into all uppercase letters
    std::transform(reference_seq.begin(), reference_seq.end(),reference_seq.begin(), ::toupper);
    /// initialise fm-index
    sdsl::csa_wt<> ref_index;   
    /// generate fm-index
    sdsl::construct_im(ref_index, reference_seq, 1); 
    /// add index and lengthcheck result to pair 
    std::pair<sdsl::csa_wt<>,bool> result;
    result.first=ref_index;
    result.second=lengthcheck;

    return result;
}

/**
 * This function translates a base into its complementary base.
 * A pairs with T, G pairs with C
 * 
 * @param n character that encodes a DNA base (A, T, G or C)
 * @return char complementary base
 */
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

/**
 * Since DNA is double stranded, a sequence will always have a complementary
 * sequence. To check whether a unitig (small sequence) is present in the overall DNA,
 * we can either search for the unitig in both the given and its inverted sequence, or
 * alternatively invert the unitig and search for this in the given sequence.
 * Generating the inverted sequence is done by exchanging As and Ts, as well as Cs and Gs.
 * 
 * @param seq sequence from .fasta file
 * @return std::string complementary sequence
 */
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

/**
 * We perform a very simple species check here.
 * The function checks in a first step whether unitigs, that we expect to be present
 * in all S.pneumoniae samples, are indeed present in the sequence. If not present we check
 * for the inverted unitig.
 * We derived these unitigs from the GPS database samples, selecting a random subset of
 * unitigs that were present in all isolates.
 * 
 * The second step is implemented only rudimentary yet, and should include checking for
 * sequences, that should not appear in S.pneumoniae isolates. So far the sequences in 
 * 'unitigs_non_pneumo.txt' are just arbitrary dummy sequences, but these could be replaced
 * in the future for more sensible sequences, e.g. from other species.
 * 
 * @param fm_index the index generated from the .fasta file
 * @return true if species check was successfull, false if not
 */
bool check_species(const sdsl::csa_wt<>& fm_index){
    /// 1st check: pneumo sequences present
    std::string filename_pneumo= "./files/species_check/unitigs_pneumo.txt"; 
    std::ifstream ist_p {filename_pneumo};
    if(!ist_p) perror("Can't open species check input file ");
    /// fill unitigs from file in vector
    std::string unitig;
    std::vector<std::string> unitigs;
    while (ist_p>>unitig)
    {   
        unitigs.push_back(unitig);
    }
    /// check if all unitigs are present in sequence
    bool result = true;
    for(int i=0; i<unitigs.size(); ++i){
        int c = count(fm_index,unitigs[i]);
        /// if not present, check inverted unitig
        if (c==0) {
            /// invert unitig and check if this is present in sequence
            int c2 = count(fm_index,invert(unitigs[i]));
            /// exit loop if any unitig is not found
            if (c2==0){result=false; break;} 
        }
    }

    //2nd check: non-pneumo sequences absent
    if(result==true){
        std::string filename_non_pneumo= "./files/species_check/unitigs_non_pneumo.txt"; 
        std::ifstream ist_np {filename_non_pneumo};
        if(!ist_np) perror("Can't open species check input file ");
        /// fill unitigs in vector
        std::string kmer;
        std::vector<std::string> kmers;
        while (ist_np>>kmer)
        {   
            kmers.push_back(kmer);
        }
        /// check if none of the unitigs is present in sequence
        for(int i=0; i<kmers.size(); ++i){
            int c = count(fm_index,kmers[i]);
            if (c==0) {
                //check inverted unitig
                int c2 = count(fm_index,invert(kmers[i]));
                //exit loop if any unitig is found
                if (c2!=0){result=false; break;} 
            } else {result=false; break;}
        }
    }
    return result;
}

/**
 * Generate a Vector that holds the information whether unitigs are
 * present (=1) or absent (=0) in a given sequence.
 * 
 * @param fm_index holding infomation about overall sequence
 * @param unitigs a vector holding multiple unitigs (smaller sequences)
 * @return Numeric_lib::Matrix<double,1> vector with 0s and 1s
 */
Numeric_lib::Matrix<double,1> lookup_unitigs(const sdsl::csa_wt<>& fm_index, const std::vector<std::string>& unitigs)
{
    /// Matrix to store presence/absence
    Numeric_lib::Matrix<double,1> pa_mat{int(unitigs.size())};
    /// looping through all unitigs in vector
    for(int i=0; i<unitigs.size(); ++i){
        /// start each vector with 1 for the Intercept
        if (unitigs[i]=="Intercept"){pa_mat[i]=1.0;} 
        else{
            /// check for presence
            int c = count(fm_index,unitigs[i]);
            if (c==0) {
                /// if not present, check inverted sequence
                int c2 = count(fm_index,invert(unitigs[i]));
                if (c2==0){pa_mat[i]=0.0;}
                else {pa_mat[i]=1.0;}
            }
            else {pa_mat[i]=1.0;}
        } 
    }
    return pa_mat;
}

/**
 * Make prediction for AMR to 5 antibiotics: 
 * "Penicillin","Chloramphenicol","Erythromycin","Tetracycline", "Trim_sulfa"
 * Results are stored in json string to be used in a web browser.
 * 
 * @param assembly_filename filename of the .fasta file
 * @return std::string json string with prediction results
 */
std::string make_prediction(std::string assembly_filename)
{
    /// generate fm-index from .fasta file to look up unitigs
    std::pair<sdsl::csa_wt<>,bool> fm_result = create_index(assembly_filename);
    sdsl::csa_wt<> fm_index = fm_result.first;

    /// perform species check
    bool check = check_species(fm_index);

    /// define set of antibiotics to be tested
    std::vector<std::string> antibiotics ={"Penicillin","Chloramphenicol","Erythromycin","Tetracycline", "Trim_sulfa"};
    
    /// initialise json object
    json results_json;

    /// since the webworker stores files in '/working/', we delete this added bit of text from the filename
    results_json["filename"] = assembly_filename.erase (0,9);

    /// if species check was successful, make predictions forr all antibiotics
    if(check==true){
        for(std::string const &antibiotic:antibiotics){ 
            /// create Model containing unitigs and their coefficients from the Elastic Net models
            Model thismodel;
            thismodel.read_model(antibiotic);  
            /// get vector of presence/absence of unitigs
            Numeric_lib::Matrix<double,1> pa_mat = lookup_unitigs(fm_index, thismodel.getUnitigs());
            /// calculate prob of resistance
            double probability = round(thismodel.get_prob(pa_mat)*1000.0)/1000.0; 
            results_json[antibiotic] = probability;
        }  
        results_json["length"]=fm_result.second;
        results_json["species"]=true;
    } else {
        results_json["species"]=false;
    }
    std::string result = results_json.dump();
    return result;
}

#ifndef WEB
/**
 * To run tests, the main function runs predictions on 2 example fasta files.
 * For the WebAssembly code, we do not need this main function
 * 
 * @return int 
 */
int main(){
    std::string result1 = make_prediction("./files/fa_examples/6999_3#3.fa.gz");
    std::string result2 = make_prediction("./files/fa_examples/6999_3#5.fa.gz");
    std::string test_result = result1+result2;

    std::ifstream testfile("./files/fa_examples/test_result.txt");
    std::string true_result((std::istreambuf_iterator<char>(testfile)),
                            (std::istreambuf_iterator<char>() ) );

    if(test_result==true_result){
        std::cout<<"Test successful!\n";
        return 0;
    }
    else{
        std::cout<<"Test failed!\n";
        return 1;
    }

} 
#else
/**
 * When compiled with emscripten to generate a WebAssembly Script
 * we need these bindings so the webworker can access the make_prediction()
 * function.
 */
EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::function("make_prediction", &make_prediction);
    emscripten::register_map<std::string,double>("map<string,double>");
}
#endif
