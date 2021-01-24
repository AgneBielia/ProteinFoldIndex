/******************************************************************************

File        : Protein.h

Date        : 21/01/2021

Author      : Agne Bieliajevaite

Description : Header file for Protein class, constructors, all appropriate
              accessors, method setFoldIndex for setting a new foldIndex value,
              a full set of relational operators (compared by foldIndex value)
              Stream-based I/O operator declarations, Output operator definition

References  : Programming 2 - C++ lecture slides and lab exercises

History     : 12/01/2021 - v1.00 Initial class header
              12/01/2021 - v1.10 Refractored Output operator definition to accommodate
                                 Protein object pointers for efficiency

******************************************************************************/

#ifndef PROTEIN_H
#define PROTEIN_H

#include<iostream>
#include <iomanip>
#include <utility>
#include <vector>

class Protein {

public:

    /*** Constructor declarations ***/

    /* Default Constructor */
    Protein() = default;

    /* Ordinary constructor */
    Protein(std::string ID, std::string sequence, std::string label , double index);

    /*** Accessor declarations ***/

    [[nodiscard]] std::string getId() const;

    [[nodiscard]] std::string getSequence() const;

    [[nodiscard]] std::string getLabel() const;

    [[nodiscard]] double getFoldIndex() const;

    /*** Method declarations ***/

    /* Set a new protein fold index by passing a double as an argument */
    void setFoldIndex(double index);

private:
    std::string ID, sequence, label;
    double foldIndex{};

};

/******************************************************************************
    Overridden relational and I/O operators declarations
******************************************************************************/

inline std::ostream& operator<<(std::ostream &os, const Protein &protein);

inline std::istream& operator>>(std::istream &is, const Protein &protein);

inline bool operator==(const Protein &p1, const  Protein &p2);

inline bool operator!=(const Protein &p1, const Protein &p2);

inline bool operator<(const Protein &p1, const Protein &p2);

inline bool operator>(const Protein &p1, const Protein &p2);

inline bool operator<=(const Protein &p1, const Protein &p2);

inline bool operator>=(const Protein &p1, const Protein &p2);

/******************************************************************************
    Constructor definition
******************************************************************************/

inline Protein::Protein(std::string ID, std::string sequence, std::string label, double index) {
    this->ID        = std::move(ID);
    this->sequence  = std::move(sequence);
    this->label     = std::move(label);
    this->foldIndex = index;
}

/******************************************************************************
    Accessor definitions
******************************************************************************/

/* Returns protein sequence as a string */
inline std::string Protein::getSequence() const{
    return this->sequence;
}

/* Returns protein label as a string */
inline std::string Protein::getLabel() const{
    return this->label;
}

/* Returns protein fold index as a double */
inline double Protein::getFoldIndex() const {
    return this->foldIndex;
}

/* Returns protein id as a string */
inline std::string Protein::getId() const {
    return this->ID;
}

/******************************************************************************

Function    : setFoldIndex

Parameters  : index - a double representing Fold Index

Description : Changes value of Protein fold index to the double index passed
              to the function

******************************************************************************/

inline void Protein::setFoldIndex(double index) {
    this->foldIndex = index;
}

/******************************************************************************
Overridden Output operator definition
******************************************************************************/

inline std::ostream& operator<<(std::ostream &os, const Protein *protein)
{
    return os << std::quoted(protein->getId()) << "," << std::quoted(protein->getSequence()) << "," << std::quoted(protein->getLabel())
          << "," << protein->getFoldIndex();
}

/******************************************************************************

Overridden relational operator definitions

Description : Compare Proteins by Fold Index value

******************************************************************************/

bool operator==(const Protein &p1, const Protein &p2) {
    return p1.getFoldIndex() == p2.getFoldIndex();
}

bool operator!=(const Protein &p1, const Protein &p2) {
    return !(p1 == p2);
}

bool operator<(const Protein &p1, const Protein &p2) {
    return p1.getFoldIndex() < p2.getFoldIndex();
}

bool operator>(const Protein &p1, const Protein &p2) {
    return p1.getFoldIndex() > p2.getFoldIndex();
}

bool operator<=(const Protein &p1, const Protein &p2) {
    return ((p1==p2)||(p1<p2));
}

bool operator>=(const Protein &p1, const Protein &p2) {
    return ((p1==p2)||(p1>p2));
}

#endif //PROTEIN_H
