/******************************************************************************

File        : Protein.cpp

Date        : 21/01/2021

Author      : Agne Bieliajevaite

Description :  Overridden input operator definition

References  : Programming 2 - C++ lecture slides and lab exercises

History     : 12/01/2021 - v1.00 Initial class

******************************************************************************/

#include "Protein.h"

/******************************************************************************
    Overridden Input operator definition
******************************************************************************/

inline std::istream &operator>>(std::istream& is, Protein &protein) {

    std::vector<std::string> v;

    std::string line;

    // check if there's input to read
    if (getline(is, line)){

        // for tesing the string captured
        //cout << line << endl;

        // store retrieved line in string stream for further processing
        std::stringstream s(line);

        // string to store each item in the stringstream s
        std::string w;

        // get items separated by a comma
        while(getline(s, w, ','))
        {
            // check if item read has quotes on it
            if ((w.back() == '"') && (w[0] == '"'))
            {
                // For testing if quotes are removed correctly
                //cout << w.substr(1,w.length()-2) << endl;

                // insert item to the vector
                v.push_back(w.substr(1,w.length()-2));
            }
                // if item does not have quotes insert  to vector v as it is
            else{

                v.push_back(w);

                // For testing if everything works correctly (prints only the read foldIndex)
                //cout << w << endl;
            }
        }

        // if all 5 values were found, create new Protein object
        if (v.size() == 4){
            protein  = Protein(v[0], v[1], v[2],stod(v[3]));
        }
        else{
            std::cerr << "Input file scan failed, incorrect number of columns." << std::endl;
        }
    }
    else {
        is.clear(std::ios_base::failbit);
    }

    return is;
}
