/******************************************************************************

File        : ProteinTest.cpp

Date        : 21/01/2021

Author      : Agne Bieliajevaite

Description : Unit testing: testing accessors and constructors,
              setFoldIndex function and full set of relational operators.

******************************************************************************/

#include <iomanip>
#include "Protein.h"

using namespace std;

class ProteinTest{

public:
    static bool testProtein()
    {
        bool result = false;
        if (test1() &&
            test2() &&
            test3())
        {
            result = true;
        }
        return result;
    }

private:
    // testing accessors and constructors
    static bool test1()
    {
        // default constructor test
        Protein p;

        // ordinary constructor
        Protein p1 = Protein("Q9HFQ6",
                             "MSTEASVSYA",
                             "#########",
                             7.0);

        if (!(p.getId().empty() &&
            p.getSequence().empty() &&
            p.getLabel().empty() &&
            p.getFoldIndex() == 0.0)){

            cerr << "Error testing default constructor and its accessors" << endl;
            return false;
        } else if(!(p1.getId() == "Q9HFQ6" &&
           p1.getSequence() == "MSTEASVSYA" &&
           p1.getLabel() == "#########" &&
           p1.getFoldIndex() == 7.0
           ))
        {
            cerr << "Error testing ordinary constructor and its accessors" << endl;
            return false;
        }

        return true;
    }

    // function testing (setFoldIndex)
    static bool test2 (){
        Protein p;
        Protein protein = Protein("Q9HFQ6",
                                  "MSTEASVSYA",
                                  "#########",
                                  7.0);

        // assign any double numbers
        protein.setFoldIndex(-8.0);
        p.setFoldIndex(27.9);

        if(p.getFoldIndex() != 27.9)
        {
            cerr << "Error testing setFoldIndex function with default constructor" << endl;
            return false;
        }
        else if (protein.getFoldIndex() != -8.0)
        {

            cerr << "Error testing setFoldIndex function with ordinary constructor" << endl;
            return false;
        }

        return true;
    }

    // relational operator testing
    static bool test3()
    {
        Protein protein = Protein("Q9HFQ6",
                                  "MSTEASVSYA",
                                  "##########",
                                  7.0);

        Protein protein1 = Protein("PQ6Q5H",
                                  "VSEIFQSSKEKLQGDA",
                                  "################",
                                  99.6);

        if(!(protein < protein1))
        {
            cerr << "Error testing less than operator" << endl;
            return false;
        }
        else if (!(protein1 > protein))
        {
            cerr << "Error testing more than operator" << endl;
            return false;
        }
        else if (protein1 == protein)
        {
            cerr << "Error testing equal operator" << endl;
            return false;
        }
        else if (!(protein1 != protein))
        {
            cerr << "Error testing not equal operator" << endl;
            return false;
        }
        else if (!(protein <= protein1)&&!(protein <= protein))
        {
            cerr << "Error testing less than or equal to operator" << endl;
            return false;
        }
        else if (!(protein1 >= protein)&&!(protein1 <= protein1))
        {
            cerr << "Error testing more than or equal to operator" << endl;
            return false;
        }

        return true;
    }




};

