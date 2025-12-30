#include "mex.h"
#include "matrix.h"
#include <random>
#include <vector>
#include <utility>   // For std::pair, std::swap
#include <algorithm> // For std::swap
#include <chrono>    // For high_resolution_clock

// Fast RNG implementation (xoshiro256++)
class Xoshiro256PP {
    uint64_t s[4];
public:
    Xoshiro256PP(uint64_t seed = 12345) {
        std::mt19937_64 gen(seed);
        for(int i = 0; i < 4; i++) s[i] = gen();
    }
    
    uint64_t next() {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }
    
private:
    static uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};

// Bit-packed matrix representation for efficient memory usage and bitwise operations
class BitPackedMatrix {
    std::vector<uint64_t> data;
    size_t rows, cols;
    size_t words_per_row;
    
public:
    BitPackedMatrix(const mxArray* mat) {
        rows = mxGetM(mat);
        cols = mxGetN(mat);
        words_per_row = (cols + 63) / 64;
        data.resize(rows * words_per_row, 0);
        
        // Convert input MATLAB matrix (logical or double) to bit-packed format
        if (mxIsLogical(mat)) {
            const mxLogical* input = mxGetLogicals(mat);
            for(size_t i = 0; i < rows; i++) {
                for(size_t j = 0; j < cols; j++) {
                    if(input[i + j * rows]) {
                        data[i * words_per_row + j/64] |= (1ULL << (j % 64));
                    }
                }
            }
        } else {
            const double* input = mxGetPr(mat);
            for(size_t i = 0; i < rows; i++) {
                for(size_t j = 0; j < cols; j++) {
                    if(input[i + j * rows] != 0) {
                        data[i * words_per_row + j/64] |= (1ULL << (j % 64));
                    }
                }
            }
        }
    }
    
    // Convert back to a MATLAB logical matrix
    mxArray* toMatlab() const {
        mxArray* result = mxCreateLogicalMatrix(rows, cols);
        mxLogical* output = mxGetLogicals(result);
        
        for(size_t i = 0; i < rows; i++) {
            for(size_t j = 0; j < cols; j++) {
                output[i + j * rows] = (data[i * words_per_row + j/64] & (1ULL << (j % 64))) != 0;
            }
        }
        return result;
    }
    
    // Get bit value at a given position
    bool get(size_t r, size_t c) const {
        return (data[r * words_per_row + c/64] & (1ULL << (c % 64))) != 0;
    }
    
    // Flip all four corners of a 2x2 rectangle defined by (r1,c1) and (r2,c2)
    void flipPattern(size_t r1, size_t r2, size_t c1, size_t c2) {
        size_t word_idx1 = c1 / 64;
        uint64_t bit_mask1 = 1ULL << (c1 % 64);
        size_t word_idx2 = c2 / 64;
        uint64_t bit_mask2 = 1ULL << (c2 % 64);

        if (word_idx1 == word_idx2) {
            uint64_t combined_mask = bit_mask1 | bit_mask2;
            data[r1 * words_per_row + word_idx1] ^= combined_mask;
            data[r2 * words_per_row + word_idx1] ^= combined_mask;
        } else {
            data[r1 * words_per_row + word_idx1] ^= bit_mask1;
            data[r1 * words_per_row + word_idx2] ^= bit_mask2;
            data[r2 * words_per_row + word_idx1] ^= bit_mask1;
            data[r2 * words_per_row + word_idx2] ^= bit_mask2;
        }
    }
    
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // --- Input and Output Checks ---
    if (nrhs < 1 || nrhs > 3) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumInputs",
                          "Must have 1 to 3 inputs: rasterMat, optional nSwaps, optional seed");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidNumOutputs",
                          "Must have exactly 1 output");
    }

    // --- Get Input Matrix ---
    const mxArray* rasterMat = prhs[0];
    if (!mxIsDouble(rasterMat) && !mxIsLogical(rasterMat)) {
        mexErrMsgIdAndTxt("shuffle_raster_mex:invalidInputType",
                          "Input must be a double or logical matrix.");
    }

    // --- Get Dimensions ---
    const mwSize nUnits = mxGetM(rasterMat);
    const mwSize nBins = mxGetN(rasterMat);

    // --- Get Number of Swaps ---
    mwSize nSwaps = nBins * nUnits;  // Default swap count
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

    // --- Get Seed ---
    bool use_explicit_seed = false;
    uint64_t explicit_seed = 0;
    if (nrhs > 2) {
        if (!mxIsNumeric(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1) {
             mexErrMsgIdAndTxt("shuffle_raster_mex:seedInvalid", "Seed must be a numeric scalar.");
        }
        explicit_seed = static_cast<uint64_t>(mxGetScalar(prhs[2]));
        use_explicit_seed = true;
    }

    // --- Early Exit for Trivial Cases ---
    if (nUnits < 2 || nBins < 2 || nSwaps == 0) {
        plhs[0] = mxDuplicateArray(rasterMat);
        return;
    }

    // --- Convert to Bit-packed Format ---
    BitPackedMatrix matrix(rasterMat);

    // --- Find all non-zero elements (spikes) ---
    std::vector<std::pair<mwSize, mwSize>> ones_coords;
    ones_coords.reserve(nUnits * nBins * 0.01); // Pre-allocate assuming sparsity
    for (mwSize r = 0; r < nUnits; ++r) {
        for (mwSize c = 0; c < nBins; ++c) {
            if (matrix.get(r, c)) {
                ones_coords.emplace_back(r, c);
            }
        }
    }

    const size_t nOnes = ones_coords.size();
    if (nOnes < 2) {
        plhs[0] = matrix.toMatlab(); // No swaps possible, return original
        return;
    }

    // --- Initialize Random Number Generator ---
    uint64_t seed;
    if (use_explicit_seed) {
        seed = explicit_seed;
    } else {
        // Use random_device combined with high_resolution_clock for robust seeding
        std::random_device rd;
        seed = (static_cast<uint64_t>(rd()) << 32) | rd();
        
        // Mix in time-based entropy to ensure MinGW/Windows non-determinism
        auto now = std::chrono::high_resolution_clock::now();
        uint64_t time_entropy = now.time_since_epoch().count();
        seed ^= time_entropy;
    }
    
    Xoshiro256PP rng(seed);

    mwSize successful_swaps = 0;
    // Set a limit on attempts to avoid infinite loops in pathological cases
    const mwSize max_attempts = nSwaps * 5; 
    mwSize attempts = 0;

    // --- Perform Swaps ---
    while (successful_swaps < nSwaps && attempts < max_attempts) {
        attempts++;
        
        // Pick two distinct random indices from the list of 'ones'
        size_t idx1 = rng.next() % nOnes;
        size_t idx2;
        do {
            idx2 = rng.next() % nOnes;
        } while (idx1 == idx2);

        auto& coord1 = ones_coords[idx1];
        auto& coord2 = ones_coords[idx2];
        
        // Check if the two spikes are in different rows and columns
        if (coord1.first == coord2.first || coord1.second == coord2.second) {
            continue;
        }

        // Check if the other two corners of the 2x2 rectangle are zeros
        if (!matrix.get(coord1.first, coord2.second) && !matrix.get(coord2.first, coord1.second)) {
            // Perform the swap by flipping all 4 corners of the rectangle
            matrix.flipPattern(coord1.first, coord2.first, coord1.second, coord2.second);

            // The 'ones' have moved, so update their column coordinates in our list
            std::swap(coord1.second, coord2.second);
            
            successful_swaps++;
        }
    }
    
    if (successful_swaps < nSwaps && nSwaps > 0) {
        mexWarnMsgIdAndTxt("shuffle_raster:fewerSwaps", 
            "Performed fewer swaps than requested. The matrix may be too sparse or dense for effective shuffling.");
    }

    // --- Convert Back to MATLAB Matrix ---
    plhs[0] = matrix.toMatlab();
}