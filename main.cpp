#include "berlecamp.hpp"
#include <vector>
//--------------------------------------------------------------------------------------------------------------------------------------------
int modul = 2;                                // Sets the GF(moduls) field
//--------------------------------------------------------------------------------------------------------------------------------------------
int n = 512;                                      // Specifies the degree of the polynomial deg f(x) = n
//--------------------------------------------------------------------------------------------------------------------------------------------
bool printToFile = false;                       // Output results to a file (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
std::string filename = "new.txt";               // File name
//--------------------------------------------------------------------------------------------------------------------------------------------
bool printMatrix = false;                        // Matrix output (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
bool printPolynomialValues = false;              // Output of intermediate polynomials (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
bool factorizePolynomial = true;                // Polynomial factorisation (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
bool generateRandomPoly = true;                // Random polynomial generation (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
bool printFirstIrreducibleOnly = true;         // Output only the first irreducible polynomial (TRUE - on, FALSE - off)
//--------------------------------------------------------------------------------------------------------------------------------------------
int maxNonzeroTerms = 10;                        // Maximum number of summands (leave 0 for large degrees)
//--------------------------------------------------------------------------------------------------------------------------------------------

int main()
{
    Berlekamp berlekampObj;
    berlekampObj.printToFile = printToFile;                              
    berlekampObj.printMatrix_global = printMatrix;                       
    berlekampObj.filename = filename;
    berlekampObj.printPolynomialValues = printPolynomialValues;                      
    berlekampObj.factorizePolynomial = factorizePolynomial;                       
    berlekampObj.randomPoly = generateRandomPoly;                                
    berlekampObj.printFirstIrreducibleOnly = printFirstIrreducibleOnly;
    berlekampObj.max_nonzero_terms = maxNonzeroTerms;                          
    berlekampObj.n = n;
    berlekampObj.modulus = modul;


    if (!berlekampObj.randomPoly)
    {

        ZZ p = conv<ZZ>(modul);
        ZZ_p::init(p);
        ZZ_pX f;

        //x^5+3x^4+7x^3+x^2+5x+6
//        SetCoeff(f, 0, conv<ZZ_p>(6));
//        SetCoeff(f, 1, conv<ZZ_p>(5));
//        SetCoeff(f, 2, conv<ZZ_p>(1));
//        SetCoeff(f, 3, conv<ZZ_p>(7));
//        SetCoeff(f, 4, conv<ZZ_p>(3));
//        SetCoeff(f, 5, conv<ZZ_p>(1));

        //x^2+x+1
//        SetCoeff(f, 0, conv<ZZ_p>(1));
//        SetCoeff(f, 1, conv<ZZ_p>(1));
//        SetCoeff(f, 2, conv<ZZ_p>(1));
//        x^7+6x^5+5x^4+4x^3+3x^2+2x^1+7
        SetCoeff(f, 0, conv<ZZ_p>(7));
        SetCoeff(f, 1, conv<ZZ_p>(2));
        SetCoeff(f, 2, conv<ZZ_p>(3));
        SetCoeff(f, 3, conv<ZZ_p>(4));
        SetCoeff(f, 4, conv<ZZ_p>(5));
        SetCoeff(f, 5, conv<ZZ_p>(6));
        SetCoeff(f, 6, conv<ZZ_p>(0));
        SetCoeff(f, 7, conv<ZZ_p>(1));


        berlekampObj.processIrreduciblePoly(f);

    }
    else
    {
        berlekampObj.generatePolynomials(n,modul);
    }

    return 0;
}
    

    

