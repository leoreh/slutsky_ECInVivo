#include "mex.h"
#include "matrix.h" // Usually included by mex.h
#include <random>    // For std::mt19937, std::uniform_int_distribution
#include <vector>    // Though not strictly needed for this revised version

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --- Input and Output Checks ---
    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumInputs",
                          "Must have 1 or 2 inputs: rasterMat and optional nSwaps");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumOutputs",
                          "Must have exactly 1 output");
    }

    // --- Get Input Matrix ---
    const mxArray* rasterMat = prhs[0];
    if (!mxIsDouble(rasterMat) && !mxIsLogical(rasterMat)) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidInputType",
                          "Input must be double or logical matrix");
    }

    // --- Get Dimensions ---
    mwSize nUnits = mxGetM(rasterMat);
    mwSize nBins = mxGetN(rasterMat);

    // --- Get Number of Swaps ---
    mwSize nSwaps = nBins; // Default
    if (nrhs > 1) {
        if (!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) {
             mexErrMsgIdAndTxt("shuffle_raster_mex:nSwapsInvalid", "nSwaps must be a numeric scalar.");
        }
        double nSwaps_double = mxGetScalar(prhs[1]);
        if (nSwaps_double < 0) {
            mexErrMsgIdAndTxt("shuffle_raster_mex:nSwapsInvalid", "nSwaps must be non-negative.");
        }
        nSwaps = static_cast<mwSize>(nSwaps_double);
    }

    // --- Create Output Matrix (Copy of Input) ---
    plhs[0] = mxDuplicateArray(rasterMat);
    double* shuffledMat = mxGetPr(plhs[0]); // mxGetDoubles is fine too

    // --- Early Exit for Trivial Cases ---
    if (nUnits < 2 || nBins < 2 || nSwaps == 0) {
        return; // Output is already a copy of input, nothing to shuffle
    }

    // --- Initialize Random Number Generator ---
    // Use random_device for a non-deterministic seed for general use
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // Or a fixed seed for reproducible testing:
    std::mt19937 gen(12345); 

    std::uniform_int_distribution<mwSize> distUnits(0, nUnits - 1);
    std::uniform_int_distribution<mwSize> distBins(0, nBins - 1);

    // --- Perform Swaps ---
    for (mwSize iSwap = 0; iSwap < nSwaps; ++iSwap) {
        mwSize r1 = distUnits(gen);
        mwSize r2;
        do {
            r2 = distUnits(gen);
        } while (r1 == r2);

        mwSize c1 = distBins(gen);
        mwSize c2;
        do {
            c2 = distBins(gen);
        } while (c1 == c2);

        // Get 2x2 submatrix values
        double v11 = shuffledMat[r1 + c1 * nUnits];
        double v12 = shuffledMat[r1 + c2 * nUnits];
        double v21 = shuffledMat[r2 + c1 * nUnits];
        double v22 = shuffledMat[r2 + c2 * nUnits];

        // Check patterns and swap
        // Using 1.0 and 0.0 for explicit double comparison
        bool pattern1 = (v11 == 1.0 && v12 == 0.0 && v21 == 0.0 && v22 == 1.0);
        bool pattern2 = (v11 == 0.0 && v12 == 1.0 && v21 == 1.0 && v22 == 0.0);

        if (pattern1 || pattern2) {
            shuffledMat[r1 + c1 * nUnits] = 1.0 - v11;
            shuffledMat[r1 + c2 * nUnits] = 1.0 - v12;
            shuffledMat[r2 + c1 * nUnits] = 1.0 - v21;
            shuffledMat[r2 + c2 * nUnits] = 1.0 - v22;
        }
    }
}