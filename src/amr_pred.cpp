#include "amr_pred_functions.hpp"

/**
 * To run tests, the main function runs predictions on 2 example fasta files and
 * compares the output with the expected results.
 * 
 * @return int 
 */
int main(){
    std::string result1 = make_prediction_json("./files/fa_examples/6999_3#3.fa.gz");
    std::string result2 = make_prediction_json("./files/fa_examples/6999_3#5.fa.gz");
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
