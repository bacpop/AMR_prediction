#include "amr_pred_functions.hpp"

/**
 * This script allows running the AMR prediction from the command line.
 * It expects a filename as input from a .txt file that lists all samples
 * to be tested (including their path from the current working directory)
 * 
 * @return int 
 */
int main(int argc, char** argv){
    /// read in txt file that lists all .fasta files to be checked
    std::string id_file = argv[1];
    std::ifstream ist {id_file};
    if(!ist) perror("Can't open file with filenames");  
    /// get total line count of file to display progress
    std::string line;
    int total_count=0;
    while(ist.peek()!=EOF)
    {
        getline(ist, line);
        total_count++;
    }
    ist.clear();
    ist.seekg(0);
    /// define set of antibiotics to be tested
    std::vector<std::string> antibiotics ={"Penicillin","Chloramphenicol","Erythromycin","Tetracycline", "Trim_sulfa"};
    /// generate header row for csv file
    std::string header_row;
    header_row.append("Filename,");
    for(std::string const &antibiotic:antibiotics){
        header_row.append(antibiotic+",");
    }
    header_row.append("length,");
    header_row.append("species,");
    header_row.append("\n");
    /// generate csv file
    std::ofstream myfile;
    myfile.open ("AMR_prediction_"+id_file+".csv");
    myfile << header_row;
    /// loop through all files in .txt file to make predictions for each
    int current_count=0;
    std::string filepath;
    while(ist>>filepath)
    {
        current_count++;
        std::cout<<"processing sample "<<current_count<<" of "<<total_count<<'\n';
        std::string result_row = make_prediction_csv(filepath, antibiotics);
        /// add prediction results to csv
        myfile << result_row;
    }
    myfile.close();
    std::cout<<"Results saved in AMR_prediction_"<<id_file<<".csv\n";
}
