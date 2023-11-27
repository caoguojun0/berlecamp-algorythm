#ifndef BERLEKAMP_HPP
#define BERLEKAMP_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <vector>
#include <thread>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;

class Berlekamp
{
public:
    Berlekamp();
    void generatePolynomials(int n, int modulus);
    void processIrreduciblePoly(const ZZ_pX& f);
    bool printMatrix_global;
    bool printToFile;
    string filename;
    bool printPolynomialValues;
    bool randomPoly;
    bool printFirstIrreducibleOnly;
    bool factorizePolynomial;
    long n;
    int max_nonzero_terms;
    int modulus;

private:
    static void writeToFile(const string& text, const string& filename);
    void BerlekampMatrix(Mat<ZZ_p>& B, const ZZ_pX& f, bool printMatrix_global,bool printPolynomialValues,const string& filename,bool printToFile);
    string PolynomialToString(const ZZ_pX& poly);
    void writeMatrixToFile(const Mat<ZZ_p>& B, const string& filename);
    void berlekamp(vec_pair_ZZ_pX_long& factors, const ZZ_pX& f);
    void runBerlekampAlgorithm(const Mat<ZZ_p>& B, int n, int mod, const ZZ_pX& f, bool factorizePolynomial, const string& filename, bool printToFile);
    long MatrixRank(const Mat<ZZ_p>& B, int n, int mod, const string& filename, bool printToFile);
    void printMatrix(const Mat<ZZ_p>& B);
};

#endif

