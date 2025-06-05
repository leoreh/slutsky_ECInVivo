#include "mex.h"
#include "matrix.h"
#include <random>
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check inputs
    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumInputs",
            "Must have 1 or 2 inputs: rasterMat and optional nSwaps");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumOutputs",
            "Must have exactly 1 output");
    }

    // Get input matrix
    const mxArray* rasterMat = prhs[0];
    if (!mxIsDouble(rasterMat) && !mxIsLogical(rasterMat)) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidInputType",
            "Input must be double or logical matrix");
    }

    // Get dimensions
    mwSize nUnits = mxGetM(rasterMat);
    mwSize nBins = mxGetN(rasterMat);

    // Get number of swaps (default to nBins if not provided)
    mwSize nSwaps = nBins;
    if (nrhs > 1) {
        nSwaps = static_cast<mwSize>(mxGetScalar(prhs[1]));
    }

    // Create output matrix (copy of input)
    plhs[0] = mxDuplicateArray(rasterMat);
    double* shuffledMat = mxGetPr(plhs[0]);

    // Initialize random number generator
    std::mt19937 gen(12345); // Fixed seed for reproducibility

    // Pre-allocate vectors for random indices
    std::vector<mwSize> units(nUnits);
    std::vector<mwSize> bins(nBins);
    for (mwSize i = 0; i < nUnits; i++) units[i] = i;
    for (mwSize i = 0; i < nBins; i++) bins[i] = i;

    // Perform swaps
    for (mwSize iSwap = 0; iSwap < nSwaps; iSwap++) {
        // Generate random indices using Fisher-Yates shuffle
        std::shuffle(units.begin(), units.begin() + 2, gen);
        std::shuffle(bins.begin(), bins.begin() + 2, gen);

        // Get 2x2 submatrix values using MATLAB-style indexing
        double v11 = shuffledMat[units[0] + bins[0]*nUnits];
        double v12 = shuffledMat[units[0] + bins[1]*nUnits];
        double v21 = shuffledMat[units[1] + bins[0]*nUnits];
        double v22 = shuffledMat[units[1] + bins[1]*nUnits];

        // Check patterns using direct comparison (like MATLAB's isequal)
        bool pattern1 = (v11 == 1 && v12 == 0 && v21 == 0 && v22 == 1);
        bool pattern2 = (v11 == 0 && v12 == 1 && v21 == 1 && v22 == 0);

        if (pattern1) {
            shuffledMat[units[0] + bins[0]*nUnits] = 0;
            shuffledMat[units[0] + bins[1]*nUnits] = 1;
            shuffledMat[units[1] + bins[0]*nUnits] = 1;
            shuffledMat[units[1] + bins[1]*nUnits] = 0;
        }
        else if (pattern2) {
            shuffledMat[units[0] + bins[0]*nUnits] = 1;
            shuffledMat[units[0] + bins[1]*nUnits] = 0;
            shuffledMat[units[1] + bins[0]*nUnits] = 0;
            shuffledMat[units[1] + bins[1]*nUnits] = 1;
        }
    }
} 