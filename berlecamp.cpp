#include "berlecamp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>
using namespace std;
using namespace NTL;

Berlekamp::Berlekamp()
{
    printFirstIrreducibleOnly = true;
    int n;
    int modulus;
}


void Berlekamp::writeToFile(const string& content, const string& filename)
{
    ofstream file;
    file.open(filename, ios::out | ios::app);
    if (file.is_open())
    {
        file << content;
        file.close();
    }
    else
    {
        cout << "Unable to open file " << filename << endl;
    }
}

 void Berlekamp::BerlekampMatrix(Mat<ZZ_p>& B, const ZZ_pX& f, bool printMatrix_global,bool printPolynomialValues,const string& filename,bool printToFile)
{
    long n = deg(f);
    B.SetDims(n, n);
    clear(B);

    ZZ p = ZZ_p::modulus();

    for (long i = 1; i < n; i++)
    {
        ZZ_pX g = PowerXMod(p * conv<ZZ>(i), f) - PowerXMod(conv<ZZ>(i), f);
        if (printPolynomialValues)
        {
            if (printToFile)
            {
                stringstream ss;
                ss << "At i = " << i << ":" << endl;
                ss << " g(x)=x^(" << conv<long>(p) << "*" << i << ") - x^" << i << ":" << endl;
                
                g %= f;
                ss << "Polynomial after taking modulo g(x): " << PolynomialToString(g) << endl << "\n";
                writeToFile(ss.str(), filename);
            }
            else
            {
                cout << "At i = " << i << ":" << endl;
                printf(" g(x)=x^(%li*%li) - x^%li:\n", conv<long>(p), i, i);
                g %= f;
                if (printPolynomialValues)
                {
                    cout << "Polynomial after taking modulo g(x): " << PolynomialToString(g) << endl ;
                    printf("\n");
                }
            }

        }
    
        for (long j = 0; j < n; j++)
        {
            B[j][i] = coeff(g, j);
        }
    }

    if (printMatrix_global)
    {
        if (!printToFile)
        {
            cout << "Berlekamp Matrix:" << endl;
            printMatrix(B);
        }
        else if (printToFile)
        {
            writeMatrixToFile(B, filename);
        }
    }
    
}

void Berlekamp::printMatrix(const Mat<ZZ_p>& B)
{
    for (long i = 0; i < B.NumRows(); i++)
    {
        cout << "[";
        for (long j = 0; j < B.NumCols(); j++)
        {
            cout << B[i][j];
            if (j < B.NumCols() - 1)
            {
                cout << " ";
            }
        }
        cout << "]" << endl;
    }
}

void Berlekamp::runBerlekampAlgorithm(const Mat<ZZ_p>& B, int n, int mod, const ZZ_pX& f, bool factorizePolynomial, const string& filename, bool printToFile)
{
    vec_pair_ZZ_pX_long factors;
    berlekamp(factors, f);
    if (factorizePolynomial)
    {
        stringstream ss;
        ss << "\nThe factors of the polynomial are: " << endl;
        for (long i = 0; i < factors.length(); i++)
        {
            ss << "(" << PolynomialToString(factors[i].a) << ") (multiplicity " << factors[i].b << ")";
            if (i < factors.length() - 1)
            {
                ss << " * ";
            }
        }
        ss << endl;

        if (printToFile)
        {
            writeToFile(ss.str(), filename);
        }
        else
        {
            cout << ss.str();
        }
    }
}

long Berlekamp::MatrixRank(const Mat<ZZ_p>& B, int n, int mod, const string& filename, bool printToFile)
{
    Mat<ZZ_p> temp(B);
    long rank = gauss(temp);
    stringstream ss;
    ss << "\nMatrix rank = " << rank << " => polynomial ";
    if (rank==(n-1))
    {
        ss << "Polynomial is irreducible" << endl;
        if (printToFile)
        {
            writeToFile(ss.str(), filename);
        }
        else
        {
            cout << ss.str();
        }
        return 0;
    }
    else
    {
        ss << "Polynomial is reduced" << endl;
        if (printToFile)
        {
            writeToFile(ss.str(), filename);
        }
        else
        {
            cout << ss.str();
        }
        return 1;
    }
}

string Berlekamp::PolynomialToString(const ZZ_pX& poly)
{
    stringstream ss;
    long degree = deg(poly);
    for (long i = degree; i >= 0; i--)
    {
        ZZ_p coef = coeff(poly, i);
        if (coef != 0)
        {
            if (i < degree && conv<int>(coef) > 0)
            {
                ss << " + ";
            }
            else if (i < degree && conv<int>(coef) < 0)
            {
                ss << " - ";
                coef = -coef;
            }
            if (i == 0 || coef != 1)
            {
                ss << coef;
            }
            if (i > 0)
            {
                ss << "x";
                if (i > 1)
                {
                    ss << "^" << i;
                }
            }
        }
    }

    return ss.str();
}

void Berlekamp::writeMatrixToFile(const Mat<ZZ_p>& matrix, const string& filename)
{
    ofstream file;
    file.open(filename, ios::out | ios::app);
    if (file.is_open())
    {
        for (long i = 0; i < matrix.NumRows(); ++i)
        {
            for (long j = 0; j < matrix.NumCols(); ++j)
            {
                file << matrix[i][j] << " ";
            }
            file << endl;
        }
        file.close();
    }
    else
    {
        cout << "Unable to open file " << filename << endl;
    }
}

void Berlekamp::berlekamp(vec_pair_ZZ_pX_long& factors, const ZZ_pX& f)
{
    CanZass(factors, f);
}

void Berlekamp::processIrreduciblePoly(const ZZ_pX& f)
{
    Mat<ZZ_p> B;
    BerlekampMatrix(B, f, printMatrix_global, printPolynomialValues,filename,printToFile);
    long rank = MatrixRank(B, n, modulus, filename, printToFile);
    if (factorizePolynomial)
    {
        runBerlekampAlgorithm(B, n, modulus, f, factorizePolynomial,filename,printToFile);
    }
}

void Berlekamp::generatePolynomials(int n, int modulus) {
    bool found = false;
    ZZ_pX f; 
    ZZ p = conv<ZZ>(modulus);
    ZZ_p::init(p);
    while (!found)
    {

        f.SetLength(n + 1);
        SetCoeff(f, n, conv<ZZ_p>(1));

        for (int i = 0; i < n; ++i)
        {
            if (RandomBnd(p) % (n/2) == 0)
            {
                SetCoeff(f, i, random_ZZ_p());
            }
            else
            {
                SetCoeff(f, i, conv<ZZ_p>(0));  
            }
        }
        if (max_nonzero_terms != 0)
        {
            int nonzero_terms = 0;
            for (int i = 0; i <= n; ++i)
            {
                if (!IsZero(coeff(f, i)))
                {
                    nonzero_terms++;
                }
            }
            
            if (nonzero_terms > max_nonzero_terms)
            {
                continue;
            }
        }
        found = IterIrredTest(f);

    }
    string poly_str = "Random irreducible polynomial: " + PolynomialToString(f) + "\n\n";
    if (printToFile)
    {
        writeToFile(poly_str, filename);
    }
    else
    {
        cout << poly_str;
    }
    processIrreduciblePoly(f);
}
