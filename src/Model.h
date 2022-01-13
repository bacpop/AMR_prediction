#include<vector>
#include<string>
#include<fstream>
#include<cmath>

class Model{
        std::vector<std::string> unitigs;
        std::vector<double> coefs;
    public:
        std::vector<std::string> getUnitigs(){return unitigs;}
        
        void read_model(const std::string& antibiotic)
        {
            //std::ostringstream filename; //set filepath and open as istream
            //filename<<"files/model_coefficients/"<<antibiotic<<"_coefficients.txt";
            std::string filename= "../files/model_coefficients/"+antibiotic+"_coefficients.txt"; // more or less efficient?? other alternative: .append()
            std::ifstream ist {filename};
            if(!ist) perror("Can't open model input file ");

            std::string unitig;  
            double coef;
            while (ist>>unitig>>coef)   //fill data in 2 vectors
            {   
                unitigs.push_back(unitig);
                coefs.push_back(coef);
            }
        }

        double get_prob(const Numeric_lib::Matrix<double,1>& pa_mat)
        {
            Numeric_lib::Matrix<double,1>coefs_mat{&coefs[0], int(coefs.size())}; // translate coefs vector into Matrix
            double pred = Numeric_lib::dot_product(pa_mat,coefs_mat); // calculate dot product
            double prob =1/(1+exp(-pred)); // apply sigmoid function for probability of resistance
            return prob;

        }

};
