/*
    [TV] marks functions taken more or less directly from Torrance and Vannieuwenhoven
*/

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define __FFLASFFPACK_USE_OPENMP

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Eigen>
#include <fflas-ffpack/fflas-ffpack.h>
#include <chrono>
#include <omp.h>
#include <string.h>

// General purpose integers, vectors and matrices using 4-byte integers
#define int_t int32_t
#define vector_t Eigen::Matrix<int_t,Eigen::Dynamic,1>  // vectors are columns
#define matrix_t Eigen::Matrix<int_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>

// Integers, vectors and matrices for representing the elements of a finite field, using 2-byte integers
#define ffint_t int16_t // int where we want to save memory
#define ffvector_t Eigen::Matrix<ffint_t,Eigen::Dynamic,1>
#define ffmatrix_t Eigen::Matrix<ffint_t,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>

#define time_point std::chrono::time_point<std::chrono::system_clock>
#define time_duration std::chrono::duration<double>
#define time_now std::chrono::system_clock::now

// Mersenne prime modulo which we will work.
// This needs to be chosen in such a way that MPRIME^2 doesn't overflow an ffint_t
#define BITS 7L
#define MPRIME ((1L<<BITS)-1L)

////////////////////
// BASE UTILITIES //
////////////////////

// [TV]
// Reduce a number to its canonical representation in Z_p.
inline ffint_t modReduce(ffint_t a) {
	assert(a >= 0);
	a = (a & MPRIME) + (a >> BITS);
	while( a > MPRIME )
		a = (a & MPRIME) + (a >> BITS);
	if(a == MPRIME) a = 0;
	return a;
}

inline ffint_t myrand() {
    ffint_t result = std::rand();
    // if RAND_MAX is bigger than the maximum of ffint_t, this could be negative and we don't want that
    if (result < 0) {
        result = -result-1;
    }
    return modReduce(result);
}

// Write random integers as elements of Z_p into the given positions of a preexisting vector_t
// (Assume the given positions are already in the range of indices in the ffvector_t)
void randomizePositions(ffvector_t &v, std::vector<int_t> &positions) {
    std::vector<int_t>::iterator it;
    for (it = positions.begin(); it != positions.end(); ++it) {
        v[*it] = myrand();
    }
}


// [TV]
// Compute the rank of A over the finite field Z_p.
int_t ffRank( ffmatrix_t& A ) {
	Givaro::Modular<ffint_t> FMPRIME(MPRIME);
	ffint_t* ptrA = A.data();
	const int_t m = A.rows();
	const int_t n = A.cols();
	const int_t rank = FFPACK::Rank(FMPRIME, m, n, ptrA, A.outerStride() );
	return rank;
}


inline double timeElapsed(const time_point &start, const time_point &end) {
    time_duration dur = end - start;
    return dur.count();
}


//////////////////////////////////
// MONOMIALS AND THEIR HANDLING //
//////////////////////////////////

// Change this vector into the next one in lexicographic order amongst non-decreasing sequences of non-negative integers at most max
// Return whether or not overflow occurred
// Presume that seq is already a non-decreasing sequence with entries at most max
bool incrementRepeatingSeq(const int_t &length, vector_t &seq, const int_t &max) {
    int_t i = length-1;
    bool overflow = false;
    while (i >= 0) {
        if (seq[i] < max) {
            ++seq[i];
            break;
        }
        --i;
    }
    int_t val;
    if (i >= 0) {
        val = seq[i];
    }
    else {  // if i = -1, all entries were m, so overflow back to all 0s
        val = 0;
        overflow = true;
    }
    ++i;
    while (i < length) {
        seq[i] = val;
        ++i;
    }
    return overflow;
}

class Monomial {
    public:
        int_t a, b; // bidegree
        vector_t xVars, yVars;   // list of variables

    Monomial(int_t myA, int_t myB) {    // default initialization
        a = myA;    xVars.setZero(a);
        b = myB;    yVars.setZero(b);
    }
    Monomial(int_t myA, int_t myB, vector_t myX, vector_t myY) : Monomial(myA, myB) {
        assert(myX.size() == a);    xVars = myX;
        assert(myY.size() == b);    yVars = myY;
    }

    // Print this monomial to a stream
    void print(std::ostream &os) {
        bool first = true;
        for (int_t i = 0; i < a; ++i) {
            if (not first) {
                os << " ";
            }
            os << "x_" << xVars[i];
        }
        for (int_t i = 0; i < b; ++i) {
            os << " y_" << yVars[i];
        }
    }

    bool incrementX(int_t m) {
        return incrementRepeatingSeq(a, xVars, m);
    }

    bool incrementY(int_t n) {
        return incrementRepeatingSeq(b, yVars, n);
    }

    // Take dm/dv evaluated at (x,y), where v means y_i for i >= 0 or x_(-i-1) for i < 0
    int_t diff(int_t varIndex, ffvector_t x, ffvector_t y) {
        ffint_t result = 1;
        int_t exponent = 0;
        if (varIndex < 0) {
            for (int_t i = 0; i < a; ++i) {
                if (xVars[i] == -varIndex - 1) {
                    if (exponent > 0) {
                        result *= x[xVars[i]];
                    }
                    ++exponent;
                }
                else {
                    result *= x[xVars[i]];
                }
                result = modReduce(result);
            }
            if (exponent == 0) {
                return 0;
            }
            else {
                result *= exponent;
                result = modReduce(result);
                for (int_t i = 0; i < b; ++i) {
                    result *= y[yVars[i]];
                    result = modReduce(result);
                }
                return result;
            }
        }
        else {
            for (int_t i = 0; i < b; ++i) {
                if (yVars[i] == varIndex) {
                    if (exponent > 0) {
                        result *= y[yVars[i]];
                    }
                    ++exponent;
                }
                else {
                    result *= y[yVars[i]];
                }
                result = modReduce(result);
            }
            if (exponent == 0) {
                return 0;
            }
            else {
                result *= exponent;
                result = modReduce(result);
                for (int_t i = 0; i < a; ++i) {
                    result *= x[xVars[i]];
                    result = modReduce(result);
                }
                return result;
            }
        }
    }
};
std::ostream& operator<<(std::ostream &os, Monomial &monom) {
    monom.print(os);
    return os;
}

///////////////////////////////
// COORDINATE CONFIGURATIONS //
///////////////////////////////

// Return if overflow happened
bool incrementBinaryString(vector_t &b) {
    bool overflow = true;
    vector_t::iterator it;
    for (it = b.begin(); it != b.end(); ++it) {
        if (*it) {
            *it = 0;
        }
        else {
            *it = 1;
            overflow = false;
            break;
        }
    }
    return overflow;
}

// Object to handle a particular configuration in a particular P^m x P^n,
// with some particular subvarieties specified by a collection of linear monomials
// and some amounts of general double points on various intersections of these varieties
// (e.g. points on the whole space, points on one of the varieties, points on an intersection of two etc.)
class CoordConfig {
    public:
        int_t m, n; // The ambient P^m x P^n
        int_t a, b; // The bidegree considered
        int_t numberOfSubvars;
        std::vector<vector_t> subvarietiesX, subvarietiesY;  // list ofboolean vectors indicating with 1's the monomial equations in x-, resp. y-variables that define the subvariety (so the matrix should have m+1, resp. n+1 rows)
        vector_t numbersOfPoints;   // Integer at index i indicates the number of points that should be (randomly) considered at the intersection of the subvarieties given by 1's in i considered as a base-2 number

    CoordConfig(int_t myM, int_t myN, int_t myA, int_t myB, std::vector<vector_t> &xSubvars, std::vector<vector_t> &ySubvars, vector_t &points) {
        m = myM;    n = myN;
        a = myA;    b = myB;
        subvarietiesX = xSubvars;
        subvarietiesY = ySubvars;
        numbersOfPoints = points;
        numberOfSubvars = subvarietiesX.size();
        assert(subvarietiesY.size() == numberOfSubvars);

        vector_t::iterator it;
        for (it = numbersOfPoints.begin(); it != numbersOfPoints.end(); ++it) {
            assert(*it >= 0);
        }

        monomialsListed = false;
        variablesListed = false;
        pointsAllocated = false;
        matrixBuilt = false;
        hasVirtDim = false;
    }

    // Print this coordinate configuration to a stream
    void print(std::ostream &os) {
        os << "----------------------------------------" << std::endl;
        os << "A coordinate configuration in P^" << m << " x P^" << n  << " with O(" << a << "," << b << ")" << std::endl;
        os << "with groups of ";
        vector_t::iterator it;
        for (it = numbersOfPoints.begin(); it != numbersOfPoints.end(); ++it) {
            os << *it << " ";
        }
        os << "points" << std::endl;
        os << "and " << numberOfSubvars << " subvarieties:" << std::endl;
        for (int_t k = 0; k < numberOfSubvars; ++k) {
            os << subvarietiesX[k].transpose() << " | " << subvarietiesY[k].transpose() << std::endl;
        }
        os << "----------------------------------------" << std::endl;
    }

    // Find all relevant monomials. If save = true, remember them, if save = false, just count them
    int_t monomialCount;
    std::vector<Monomial> monomials;
    bool monomialsListed;
    void getMonomials(bool save) {
        monomialCount = 0;
        if (save) {
            monomials = {};
        }

        Monomial mon(a, b);

        bool xloop, yloop, fits, accepted;
        xloop = false;
        while (not xloop) {
            yloop = false;
            while (not yloop) {
                // does this monomial lie in each of the ideals of the subvarieties?
                accepted = true;
                for (int_t k = 0; k < numberOfSubvars; ++k) {
                    fits = false;
                    for (int_t i = 0; i < a; ++i) {
                        if (subvarietiesX[k][mon.xVars[i]]) {
                            fits = true;
                            break;
                        }
                    }
                    if (not fits) {
                        for (int_t i = 0; i < b; ++i) {
                            if (subvarietiesY[k][mon.yVars[i]]) {
                                fits = true;
                                break;
                            }
                        }
                    }
                    if (not fits) {
                        accepted = false;
                        break;
                    }
                }
                if (accepted) {
                    ++monomialCount;
                    if (save) {
                        monomials.push_back(mon);
                    }
                }
                yloop = mon.incrementY(n);
            }
            xloop = mon.incrementX(m);
        }

        if (save) {
            monomialsListed = true;
        }
    }

    bool variablesListed;
    bool pointsAllocated;
    std::vector<ffvector_t> generatedPointsX, generatedPointsY;   // these sequences will potentially hold the generated points
    std::vector<std::vector<int_t>> relevantVariables;  // for each relevant intersection of subvarieties, list of variables with respect to which we will take derivatives
    std::vector<int_t> pointNums;  // rewrite the numbers of points into a new vector whose indexation agrees with that of relevantVariables, i.e. one that skips intersections with no relevant points. This will only be used for virtual dimension calculation.
    std::vector<std::vector<int_t>> constraintsX, constraintsY;   // for each relevant intersection of subvarieties, vector where 0s denote coordinates that are to be contrained to zero
    std::vector<int_t> variableCount;   // lengths of entries in relevantVariables
    int_t columnCount;
    std::vector<int_t> startColumn; // At which column do derivatives pertaining to each point start?
    std::vector<int_t> group;   // The index in variableCount, constraintsX/Y and relevantVariables of the group to which a point at the given index should belong
    // To sum up, these vectors should have length = # of relevant intersections:
    //      relevantVariables, pointNums, constraintsX, constraintsY, variableCount
    // And the following should have length = total # of points:
    //      generatedPointsX, generatedPointsY, startColumn, groups

    // Prepare relevant variables and other infrastructure relating to the would-be columns of the matrix to be built
    // The argument denotes whether or not vectors for the points that will be randomly generated should be allocated during this. This is more or less equivalent to whether or not we are expecting to write out a certificate after generating these points.
    void getVariables(bool allocateForPoints) {
        if (allocateForPoints) {
            generatedPointsX = {};
            generatedPointsY = {};
        }
        relevantVariables = {};
        pointNums = {};
        constraintsX = {};
        constraintsY = {};
        variableCount = {};
        columnCount = 0;
        startColumn = {};
        group = {};


        vector_t subset(numberOfSubvars);   // boolean vector to iterate through possible intersections of subvarieties
        // vectors for the calculation of constraints and variables
        vector_t constrX(m+1), varX(m+1);
        vector_t constrY(n+1), varY(n+1);
        ffvector_t allocX(m+1), allocY(n+1);
        bool overflow = false;
        int_t z = 0;
        int_t groupIndex = 0;
        while (not overflow) {
            if (numbersOfPoints[z] > 0) {
                relevantVariables.push_back({});
                pointNums.push_back(numbersOfPoints[z]);
                constraintsX.push_back({});
                constraintsY.push_back({});
                constrX.setOnes();
                constrY.setOnes();
                varX.setOnes();
                varY.setOnes();

                // iterate through all subvarieties and apply relevant ones to constraints and variables
                for (int_t k = 0; k < numberOfSubvars; ++k) {
                    if (subset[k]) {
                        for (int_t i = 0; i <= m; ++i) {
                            constrX[i] *= 1-subvarietiesX[k][i];
                            varX[i] *= subvarietiesX[k][i];
                        }
                        for (int_t i = 0; i <= n; ++i) {
                            constrY[i] *= 1-subvarietiesY[k][i];
                            varY[i] *= subvarietiesY[k][i];
                        }
                    }
                }

                // push to vectors
                int_t varNum = 0;
                for (int_t i = 0; i <= m; ++i) {
                    if (varX[i]) {
                        relevantVariables.back().push_back(-i-1);
                        ++varNum;
                    }
                    if (constrX[i]) {
                        constraintsX.back().push_back(i);
                    }
                }
                for (int_t i = 0; i <= n; ++i) {
                    if (varY[i]) {
                        relevantVariables.back().push_back(i);
                        ++varNum;
                    }
                    if (constrY[i]) {
                        constraintsY.back().push_back(i);
                    }
                }
                variableCount.push_back(varNum);

                // deal with information for individual vectors
                for (int_t i = 0; i < numbersOfPoints[z]; ++i) {
                    startColumn.push_back(columnCount);
                    columnCount += varNum;
                    group.push_back(groupIndex);
                    if (allocateForPoints) {
                        generatedPointsX.push_back(allocX);
                        generatedPointsY.push_back(allocY);
                    }
                }
                ++groupIndex;
            }
            overflow = incrementBinaryString(subset);
            ++z;
        }

        variablesListed = true;
        if (allocateForPoints) {
            pointsAllocated = true;
        }
    }
    void getVariables() {
        getVariables(false);
    }

    // Construct the main matrix
    ffmatrix_t calcMatrix;
    bool matrixBuilt;
    // std::vector<vector_t> generatedX, generatedY;
    void buildMatrix(bool save) {   // save = should the randomly generated points be saved?
        if (not monomialsListed) {
            getMonomials(true);
        }
        if (not variablesListed or (save and not pointsAllocated)) {
            getVariables(save);
        }

        calcMatrix.resize(monomialCount, columnCount);

        #pragma omp parallel for
        for (int_t p = 0; p < group.size(); ++p) {
            int_t g = group[p];
            ffvector_t x(m+1);    x.setZero();
            randomizePositions(x, constraintsX[g]);
            ffvector_t y(n+1);    y.setZero();
            randomizePositions(y, constraintsY[g]);

            if (save) {
                generatedPointsX[p] = x;
                generatedPointsY[p] = y;
            }

            for (int_t j = 0; j < monomialCount; ++j) {
                for (int_t i = 0; i < variableCount[g]; ++i) {
                    calcMatrix(j, startColumn[p]+i) = monomials[j].diff(relevantVariables[g][i], x, y);
                }
            }
        }

        matrixBuilt = true;
    }
    void buildMatrix() {
        buildMatrix(false);
    }

    // Get the dimension of the ideal of this coordinate configuration
    int_t getDim() {
        if (not matrixBuilt) {
            buildMatrix();
        }
        const int_t rank = ffRank(calcMatrix);
        matrixBuilt = false;    // rank computation messes up the matrix, pretend we've lost it
        return monomialCount - rank;
    }

    int_t virtDim;
    bool hasVirtDim;
    int_t getVirtDim() {
        if (not hasVirtDim) {
            virtDim = makeVirtDim();
            hasVirtDim = true;
        }
        return virtDim;
    }
    // Get the virtual dimension of this coordinate configuration
    int_t makeVirtDim() {
        if (not monomialsListed) {
            getMonomials(false);
        }
        if (not variablesListed) {
            getVariables();
        }

        int_t expDirs = 0;

        for (int_t z = 0; z < relevantVariables.size(); ++z) {
            expDirs += pointNums[z] * std::min(n+m+1, variableCount[z]);
        }

        return monomialCount - expDirs;
    }
    int_t expectedDim() {
        return std::max(0, getVirtDim());
    }
    // Is this coordinate configuration (sub/super/equi)abundant?
    bool subabundant() {
        return (getVirtDim() >= 0);
    }
    bool superabundant() {
        return (getVirtDim() <= 0);
    }
    bool equiabundant() {
        return (getVirtDim() == 0);
    }
    std::string abundancy() {
        if (equiabundant()) {
            return "EQUIABUNDANT";
        }
        else if (subabundant()) {
            return "SUBABUNDANT";
        }
        else {
            return "SUPERABUNDANT";
        }
    }

    // Certificate printing
    std::ofstream certificateDest;
    // Set where the certificate is to be written
    void openCertificate(const std::string filename) {
        certificateDest.open(filename);
    }
    void closeCertificate() {
        certificateDest.close();
    }
    // Make the determination of (non)-defectivity of this CoordConfig
    // Return whether or not the coordinate configuration was found non-defective
    bool computeWithCertificate() {
        if (certificateDest.is_open()) {
            time_point absoluteStart = time_now();
            auto seed = int(time(0));
            srand(seed);
            certificateDest << "Using random seed " << seed << " and the finite field with " << MPRIME << " elements" << std::endl;
            certificateDest << "In the ambient space P^" << m << " x P^" << n << ", computing the dimension of the O(" << a << "," << b << ") part of the ideal of a collection of";
            if (numberOfSubvars > 0) {
                certificateDest << " subvarieties given by" << std::endl;
            }
            // write individual subvarieties
            bool first;
            for (int_t k = 0; k < numberOfSubvars; ++k) {
                first = true;
                certificateDest << "{";
                for (int_t i = 0; i <= m; ++i) {
                    if (subvarietiesX[k][i]) {
                        if (not first) {
                            certificateDest << " = ";
                        }
                        certificateDest << "x_" << i;
                        first = false;
                    }
                }
                for (int_t j = 0; j <= n; ++j) {
                    if (subvarietiesY[k][j]) {
                        if (not first) {
                            certificateDest << " = ";
                        }
                        certificateDest << "y_" << j;
                        first = false;
                    }
                }
                certificateDest << " = 0}" << std::endl;
            }
            if (numberOfSubvars > 0) {
                certificateDest << "and";
            }
            certificateDest << " double points supported at" << std::endl;

            time_point timeBeforeMatrix = time_now();
            buildMatrix(true);
            time_point timeAfterMatrix = time_now();

            // write out the points
            for (int_t k = 0; k < generatedPointsX.size(); ++k) {
                certificateDest << "([" << generatedPointsX[k].transpose() << "], [" << generatedPointsY[k].transpose() << "])" << std::endl;
            }

            certificateDest << "The " << monomialCount << " x " << columnCount << " matrix was built in " << timeElapsed(timeBeforeMatrix, timeAfterMatrix) << " s." << std::endl;

            time_point timeBeforeDim = time_now();
            int_t dim = getDim();
            time_point timeAfterDim = time_now();

            certificateDest << "Computing the rank of the " << monomialCount << " x " << columnCount << " matrix took " << timeElapsed(timeBeforeDim, timeAfterDim) << " s." << std::endl;
            certificateDest << "The dimension is " << dim << " vs. expected " << expectedDim() << std::endl;

            certificateDest << "The configuration is ";
            if (dim == expectedDim()) {
                certificateDest << "NON-DEFECTIVE";
            }
            else {
                certificateDest << "probably DEFECTIVE";
            }
            certificateDest << " (" << abundancy() << ")." << std::endl;

            time_point absoluteEnd = time_now();
            certificateDest << "The entire computation took " << timeElapsed(absoluteStart, absoluteEnd) << " s." << std::endl;

            return (dim == expectedDim());
        }
        else {
            std::cout << "Certificate file is not open. Abort." << std::endl;
            return false;
        }
    }
    bool computeWithCertificate(const std::string filename) {
        openCertificate(filename);
        bool result = computeWithCertificate();
        closeCertificate();
        return result;
    }


};
std::ostream& operator<<(std::ostream &os, CoordConfig &cc) {
    cc.print(os);
    return os;
}


///////////////////////////
// BIVARIATE POLYNOMIALS //
///////////////////////////

#define MAXDEGREE 4

// [TV]
// Computes the binomial coefficient n choose d.
inline int_t binomial(const int_t n, const int_t d) {
	if (d > n) {
        return 0;
    }
	int_t r = 1;
	for (int_t k = n; k > n-d; --k) {
		r *= k;
    }
	for (int_t k = 2; k <= d; ++k) {
		r /= k;
    }
	return r;
}

// Return the vector [1 x ... x^MAXDEGREE]^T
// (for polynomial evaluation)
inline vector_t moment(const int_t x) {
    vector_t result(MAXDEGREE+1);
    int_t power = 1;
    for (int_t i = 0; i <= MAXDEGREE; ++i){
        result[i] = power;
        power *= x;
    }
    return result;
}

// atoi that takes a string and trims all whitespace
int_t myatoi(std::string str) {
    std::string newstr = "";
    std::string::iterator it;
    for (it = str.begin(); it != str.end(); ++it) {
        if (not isspace(*it)) {
            newstr += *it;
        }
    }
    return atoi(newstr.c_str());
}

#define DEFAULTMCHAR 'm'
#define DEFAULTNCHAR 'n'

// Object to handle matrices denoting bivariate polynomials.
// The two variables will be implicitly called m and n.
// In evaluation, the moment vectors of m and n will be multiplied to the matrix from left and right respectively
class Bipol {
    public:
        matrix_t coefs;
        int_t denominator;

    Bipol() {   // default initialization
        coefs.setZero(MAXDEGREE+1, MAXDEGREE+1);
        denominator = 1;
        hasDegs = false;
        hasString = false;
    }
    Bipol(const Bipol &other) : Bipol() {
        coefs = other.coefs;
        denominator = other.denominator;
        getDegs();
    }
    Bipol(matrix_t myCoefs, int_t myDenom) : Bipol() {
        coefs = myCoefs;
        denominator = myDenom;
        getDegs();
    }
    Bipol(const int_t constValue) : Bipol() {
        coefs(0,0) = constValue;
        degM = 0;   degN = 0;   hasDegs = true;
    }

    // Read data into self from a string.
    // This should be of the form (...+ c m^k n^j +...)/denominator.
    // Parentheses may be omitted, whitespace should be ignored. No nested parentheses inside are allowed.
    void readFromString(const std::string &src, char mVar, char nVar){
        hasString = false;
        // first check for the denominator
        int_t stop = src.find('/');
        if (stop != std::string::npos) {
            denominator = myatoi(src.substr(stop+1));
        }
        // if '/' is not found, we just leave the denominator as the default 1
        else {
            stop = src.length();
        }

        int_t start, end, mPosition, nPosition;
        // move our indices so as to disregard possible parentheses
        start = src.find('(');
        if (start == std::string::npos) {
            start = 0;
        }
        else {
            start += 1;
            stop = src.find(')');
        }

        int_t coef, mExp, nExp;
        int_t caret, limit;
        bool sign, numeralEncountered;
        // go until you reach the stopping point
        while (start < stop) {
            end = std::min(src.find('+', start+1), src.find('-', start+1));
            if (end == std::string::npos) {
                end = stop;
            }

            // find u and v
            mPosition = src.find(mVar, start);
            if (mPosition == std::string::npos or mPosition > end) {
                mPosition = end;
                mExp = 0;
            }
            nPosition = src.find(nVar, start);
            if (nPosition == std::string::npos or nPosition > end) {
                nPosition = end;
                nExp = 0;
            }

            // find the coefficient, handling a special case of no numerals
            limit = std::min(mPosition, nPosition);
            sign = false;
            numeralEncountered = false;
            for (int_t i = start; i < limit; ++i) {
                if ('9' >= src[i] and src[i] >= '0') {
                    numeralEncountered = true;
                    break;
                }
                else if (src[i] == '-') {
                    sign = not sign;
                }
            }
            if (numeralEncountered) {
                coef = myatoi(src.substr(start, limit-start));
            }
            else {
                if (sign) {
                    coef = -1;
                }
                else {
                    coef = 1;
                }
            }

            // find the exponents (unless the variable was already not present)
            if (mPosition < end) {
                limit = nPosition;
                if (nPosition < mPosition) {
                    limit = end;
                }
                caret = src.find('^', mPosition+1);
                if (caret == std::string::npos or caret >= limit) {
                    mExp = 1;
                }
                else {
                    mExp = myatoi(src.substr(caret+1, limit-caret-1));
                }
            }
            if (nPosition < end) {
                limit = mPosition;
                if (mPosition < nPosition) {
                    limit = end;
                }
                caret = src.find('^', nPosition+1);
                if (caret == std::string::npos or caret >= limit) {
                    nExp = 1;
                }
                else {
                    nExp = myatoi(src.substr(caret+1, limit-caret-1));
                }
            }

            // Write the coefficient. If the exponents were higher than MAXDEGREE, this will throw an error and that is correct
            coefs(mExp, nExp) = coef;

            start = end;
        }
        reduce();
    }
    Bipol(const std::string &src, char mVar, char nVar) : Bipol() {
        readFromString(src, mVar, nVar);
        getDegs();
    }
    Bipol(const std::string &src) : Bipol(src, DEFAULTMCHAR, DEFAULTNCHAR) {}
    Bipol(const char src[], char mVar, char nVar) : Bipol(std::string(src), mVar, nVar) {}
    Bipol(const char src[]) : Bipol(std::string(src)) {}


    bool hasDegs;
    int_t degM, degN, deg;
    // Construct the degree information
    void getDegs() {
        degM = -1;
        degN = -1;
        deg = -1;
        for (int_t i = 0; i <= MAXDEGREE; ++i) {
            for (int_t j = 0; j <= MAXDEGREE; ++j) {
                if (coefs(i,j) != 0) {
                    if (i > degM) {
                        degM = i;
                    }
                    if (j > degN) {
                        degN = j;
                    }
                    if (i+j > deg) {
                        deg = i+j;
                    }
                }
            }
        }
        hasDegs = true;
    }

    // Evaluate the polynomial at a point.
    int_t eval(int_t m, int_t n) const {
        int_t beforeDivision = (moment(m)).dot(coefs * moment(n));
        assert(beforeDivision % denominator == 0);
        return beforeDivision / denominator;
    }

    // Is this zero?
    bool isZero() const {
        for (int_t i = 0; i <= MAXDEGREE; ++i) {
            for (int_t j = 0; j <= MAXDEGREE; ++j) {
                if (coefs(i,j)) {
                    return false;
                }
            }
        }
        return true;
    }

    // Convert to a string compatible with readFromString
    std::string toString(char mVar, char nVar) const {
        std::string result = "";
        int_t termCount = 0;
        std::string term;
        int_t coef;
        for (int_t i = MAXDEGREE; i >= 0; --i) {
            for (int_t j = MAXDEGREE; j >= 0; --j) {
                coef = coefs(i,j);
                if (not coef) {
                    continue;
                }

                if (i == 0 and j == 0) {
                    term = std::to_string(coef);
                    if (termCount > 0 and coef > 0) {
                        term = "+" + term;
                    }
                }
                else {
                    term = "";
                    if (i > 0) {
                        term += mVar;
                        if (i > 1) {
                            term += "^" + std::to_string(i);
                        }
                    }
                    if (j > 0) {
                        term += nVar;
                        if (j > 1) {
                            term += "^" + std::to_string(j);
                        }
                    }
                    if (coef > 0) {
                        if (coef != 1) {
                            term = std::to_string(coef) + term;
                        }
                        if (termCount > 0) {
                            term = "+" + term;
                        }
                    }
                    else {
                        if (coef == -1) {
                            term = "-" + term;
                        }
                        else {
                            term = std::to_string(coef) + term;
                        }
                    }
                }
                // append the term
                termCount += 1;
                result += term;
            }
        }
        if (termCount == 0) {
            result = "0";
        }
        else if (denominator != 1) {
            if (termCount >= 2) {
                result = "(" + result + ")";
            }
            result += " / " + std::to_string(denominator);
        }
        return result;
    }
    std::string toString() const {
        return toString(DEFAULTMCHAR, DEFAULTNCHAR);
    }
    bool hasString;
    std::string string;
    std::string getString() {
        if (not hasString) {
            string = toString();
            hasString = true;
        }
        return string;
    }

    // Remove any common (integer) factors from the numerator and denominator
    void reduce() {
        // this procedure might change how the polynomial is written, so throw away any string we might have
        hasString = false;
        int_t commonFactor = std::abs(denominator);
        for (int_t i = 0; i <= MAXDEGREE; ++i) {
            for (int_t j = 0; j <= MAXDEGREE; ++j) {
                if (coefs(i,j)) {
                    commonFactor = std::gcd(commonFactor, coefs(i,j));
                }
            }
        }
        // make sure the new denominator will be positive
        if (denominator < 0) {
            commonFactor *= -1;
        }
        denominator /= commonFactor;
        coefs /= commonFactor;
    }

    // addition in place
    Bipol& operator+=(const Bipol &other) {
        hasString = false;
        int_t newDenom = std::lcm(denominator, other.denominator);
        coefs = coefs * (newDenom/denominator) + other.coefs * (newDenom/other.denominator);
        denominator = newDenom;
        return *this;
    }
    // addition
    Bipol operator+(const Bipol &other) const {
        Bipol result(*this);
        result += other;
        return result;
    }
    // subtraction in place
    Bipol& operator-=(const Bipol &other) {
        hasString = false;
        int_t newDenom = std::lcm(denominator, other.denominator);
        coefs = coefs * (newDenom/denominator) - other.coefs * (newDenom/other.denominator);
        denominator = newDenom;
        return *this;
    }
    // subtraction
    Bipol operator-(const Bipol &other) const {
        Bipol result(*this);
        result -= other;
        return result;
    }

    // scalar multiplication in place
    Bipol& operator*=(const int_t other) {
        hasString = false;
        int_t commonFactor = std::gcd(denominator, other);
        denominator /= commonFactor;
        coefs *= (other / commonFactor);
        return *this;
    }
    // scalar multiplication
    Bipol operator*(const int_t other) const {
        Bipol result(*this);
        result *= other;
        return result;
    }
    // multiplication
    Bipol operator*(const Bipol &other) const {
        if (hasDegs and other.hasDegs) {
            assert(std::max(degM+other.degM, degN+other.degN) <= MAXDEGREE);
        }

        Bipol result;
        result.denominator = denominator * other.denominator;
        int_t highestSingleDegree = 0;
        for (int_t i = 0; i <= MAXDEGREE; ++i) {
            for (int_t j = 0; j <= MAXDEGREE; ++j) {
                result.coefs(i,j) = 0;
                for (int_t ii = 0; ii <= i; ++ii) {
                    for (int_t jj = 0; jj <= j; ++jj) {
                        result.coefs(i,j) += coefs(ii,jj) * other.coefs(i-ii,j-jj);
                    }
                }
            }
        }
        return result;
    }

    // substitution
    Bipol substitute(const Bipol &m, const Bipol &n) const {
        Bipol result, term;
        for (int_t i = 0; i <= MAXDEGREE; ++i) {
            for (int_t j = 0; j <= MAXDEGREE; ++j) {
                if (coefs(i,j)) {
                    term = coefs(i,j);
                    for (int_t k = 0; k < i; ++k) {
                        term = term * m;
                    }
                    for (int_t k = 0; k < j; ++k) {
                        term = term * n;
                    }
                    result += term;
                }
            }
        }
        result.hasString = false;
        result.denominator *= denominator;
        result.reduce();
        return result;
    }

    bool operator==(const Bipol &other) const {
        return ((coefs * other.denominator) == (other.coefs * denominator));
    }

    void matrixPrint(std::ostream &os) {
        os << coefs << std::endl << "/ " << denominator;
    }
};
std::ostream& operator<<(std::ostream &os, Bipol &f) {
    os << f.getString();
    return os;
}


//////////////////////////////////////////
// FAMILIES OF COORDINATE CONFIGURATION //
//////////////////////////////////////////

class CoordConfigFamily {
    public:
        int_t a, b;
        int_t subvarCount;
        Bipol polyM, polyN; // polynomial functions to describe the ambient space
        std::vector<Bipol> eqsX, eqsY;   // Polynomial functions describing how many x- and y-equations constitute exactly the intersection (and no other subvariety) indicated by the index as a binary string.
        std::vector<Bipol> points;  // similarly, functions for the numbers of points at each intersection

    CoordConfigFamily(int_t myA, int_t myB, const Bipol &myM, const Bipol &myN, int_t subCount, const std::vector<Bipol> &myEqsX, const std::vector<Bipol> &myEqsY, const std::vector<Bipol> &myPts) {
        a = myA;    b = myB;
        polyM = myM;    polyN = myN;
        subvarCount = subCount;
        eqsX = myEqsX;   eqsY = myEqsY;
        points = myPts;
        assert(eqsX.size() == std::pow(2, subvarCount) and eqsX.size() == eqsY.size() and eqsX.size() == points.size());
    }

    // Constructor where we (possibly) automatically compute eqs[0], i.e. the number of equations that go unused
    CoordConfigFamily(int_t myA, int_t myB, const Bipol &myM, const Bipol &myN, int_t subCount, const std::vector<Bipol> &myEqsX, const std::vector<Bipol> &myEqsY, const std::vector<Bipol> &myPts, bool recheckUnusedEqs) : CoordConfigFamily(myA, myB, myM, myN, subCount, myEqsX, myEqsY, myPts) {
        eqsX[0] = polyM+1;
        eqsY[0] = polyN+1;
        for (int_t z = 1; z < eqsX.size(); ++z) {
            eqsX[0] -= eqsX[z];
            eqsY[0] -= eqsY[z];
        }
    }

    // Remove points from intersections that are irrelevant
    // Return whether any change actually occurred
    bool clean() {
        bool result = false;
        Bipol polyZero;
        std::vector<bool> noEqs;
        std::vector<bool> irrelevantGroup;
        for (int_t z = 0; z < points.size(); ++z) {
            noEqs.push_back(eqsX[z].isZero() and eqsY[z].isZero());
            irrelevantGroup.push_back(true);   // just initializing this to length 2^subvarCount
        }

        // iteratively go through each group to see if it and all finer groups have no eqs
        int_t binpow;
        vector_t negative(subvarCount); // we really imagine this to be bitwise negatives, it's just that I want to reuse the increment function instead of coding a decrement
        negative.setZero();
        for (int_t z = points.size()-1; z >= 0; --z) {
            // if this has no eqs
            if (noEqs[z]) {
                // check that also all finer groups have no eqs
                binpow = 1;
                 // irrelevantGroup[z] is initially true
                for (int_t i = 0; i < subvarCount; ++i) {
                    if (z+binpow >= points.size()) {    // we don't want to look at garbage data outside the bounds of our vectors
                        break;
                    }
                    if (negative[i]) {
                        if (not irrelevantGroup[z+binpow]) {
                            irrelevantGroup[z] = false;
                            break;
                        }
                    }
                    binpow <<= 1;
                }
                if (irrelevantGroup[z]) {
                    points[z] = polyZero;
                    result = true;
                }
            }
            else {
                irrelevantGroup[z] = false;
            }
            incrementBinaryString(negative);
        }
        return result;
    }

    // Return an (mSource, nSource)-inductant of the current family.
    // For implementational reasons, we actually allows the "smaller instance" to come from a different family.
    // (This is because we can't have quasipolynomials, so some families are actually implemented as two separate ones.)
    CoordConfigFamily inductant(const Bipol &mSource, const Bipol &nSource, const CoordConfigFamily &source) {
        std::vector<Bipol> newEqsX, newEqsY, newPoints, tmpPoints;

        // deal with patches of the Venn diagram that don't contain the new subvariety
        for (int_t z = 0; z < points.size(); ++z) {
            // equations outside need to look like a smaller instance
            newEqsX.push_back(source.eqsX[z].substitute(mSource, nSource));    newEqsX.back().reduce();
            newEqsY.push_back(source.eqsY[z].substitute(mSource, nSource));    newEqsY.back().reduce();
        }
        // deal with ones that do
        for (int_t z = 0; z < points.size(); ++z) {
            // just put whatever number of equations is left
            newEqsX.push_back(eqsX[z] - newEqsX[z]);    newEqsX.back().reduce();
            newEqsY.push_back(eqsY[z] - newEqsY[z]);    newEqsY.back().reduce();
            // for points, we need the inner numbers first, so we put them to a temporary vector
            tmpPoints.push_back(source.points[z].substitute(mSource, nSource));
        }
        // calculate outer point numbers
        for (int_t z = 0; z < points.size(); ++z) {
            newPoints.push_back(points[z] - tmpPoints[z]);
            newPoints.back().reduce();
        }
        // put the inner ones back
        for (int_t z = 0; z < points.size(); ++z) {
            newPoints.push_back(tmpPoints[z]);
        }

        CoordConfigFamily result(a, b, polyM, polyN, subvarCount+1, newEqsX, newEqsY, newPoints);
        result.clean();
        return result;
    }
    CoordConfigFamily inductant(const Bipol &mSource, const Bipol &nSource) {
        return inductant(mSource, nSource, *this);
    }

    // Evaluate at integers m, n to get a concrete coordinate configuration.
    CoordConfig eval(int_t m, int_t n) {
        const int_t mDim = polyM.eval(m, n);
        const int_t nDim = polyN.eval(m, n);

        vector_t pts(points.size());
        // deal with point
        for (int_t z = 0; z < pts.size(); ++z) {
            pts[z] = points[z].eval(m, n);
        }

        vector_t subset(subvarCount);
        std::vector<vector_t> xSubvars;
        std::vector<vector_t> ySubvars;
        vector_t xdud(m+1); // initialize these two to zero, just for
        vector_t ydud(n+1); // the purpose of initializing the above two vectors
        for (int_t k = 0; k < subvarCount; ++k) {
            xSubvars.push_back(xdud);
            ySubvars.push_back(ydud);
        }
        bool stop = false;
        int_t xdir, ydir;
        int_t xfree = 0, yfree = 0; // first available x-, resp. y-monomial
        int_t z = 0;
        while (not stop) {
            xdir = eqsX[z].eval(m,n);   assert(xdir >= 0);
            ydir = eqsY[z].eval(m,n);   assert(ydir >= 0);
            // allocate these equations to the indicated subvarieties
            for (int_t i = 0; i < xdir; ++i) {
                for (int_t k = 0; k < subvarCount; ++k) {
                    if (subset[k]) {
                        xSubvars[k][xfree+i] = 1;
                    }
                }
            }
            xfree += xdir;
            for (int_t i = 0; i < ydir; ++i) {
                for (int_t k = 0; k < subvarCount; ++k) {
                    if (subset[k]) {
                        ySubvars[k][yfree+i] = 1;
                    }
                }
            }
            yfree += ydir;

            ++z;
            stop = incrementBinaryString(subset);
        }

        assert(xfree == mDim+1 and yfree == nDim+1);  // check that we've successfully used up all variables

        return CoordConfig (mDim, nDim, a, b, xSubvars, ySubvars, pts);
    }

    void print(std::ostream &os) {
        os << "----------------------------------------" << std::endl;
        os << "A family of coordinate configurations in P^{" << polyM << "} x P^{" << polyN << "} with O(" << a << "," << b << ")" << std::endl;
        os << "with groups of " << std::endl;
        std::vector<Bipol>::iterator it;
        for (it = points.begin(); it != points.end(); ++it) {
            os << *it << std::endl;
        }
        os << "points and " << subvarCount << " subvarieties specified by:" << std::endl;
        for (int_t k = 0; k < eqsX.size(); ++k) {
            os << eqsX[k] << " and "<< eqsY[k] << std::endl;
        }
        os << "----------------------------------------";
    }
};
std::ostream& operator<<(std::ostream &os, CoordConfigFamily &ccFam) {
    ccFam.print(os);
    return os;
}
