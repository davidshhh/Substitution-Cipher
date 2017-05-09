#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <regex.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <queue>
#include <time.h>
#include <sstream>

using namespace std;

/* --------------------------------- Global Definitions ------------------------------------ */

// In order:  the number of plaintext letters, the number of ciphertext letters, and the
//            length of the cipher.

int plain_num;
int cipher_num;
long cipher_length;

// Character vectors to represent the plaintext and ciphertext letters.

// These are only used in setting up: in the actual algorithm, the first ciphertext letter
// is given the index 0, the second the index 1, and so on.  So if the ciphertext letters are
// '_', 'a', and 'b', then 'ab ba' would be transformed into '1, 2, 0, 2, 1'.
// This allows us to just use the ciphertext letters and plaintext letters as indices, without
// an extra read into memory.  It has proven to speed up the program.

vector<char> cipher_alpha(0);
vector<char> plain_alpha(0);

// These hold the cipher string (as said before, as numbers, not as characters).
// The vector is used in setting up, since we don't know the size of the cipher beforehand.
// Once the string is read and the vector is full, we allocate the array for all future use.
// This is done for speed.

vector<int> cipher_string_vec(0);
int * cipher_string;

// Arrays to hold the unigram, bigram, and trigram probabilities.
// They are all one - dimensional.  Probabilities are accessed as follows:
// Unigram: P(x) = *(unigram + x)
// Bigram: P(x | y) = *(bigram + x + y * plain_num)
// Trigram: P(x | zy) = *(trigram + x + y * plain_num + z * plain_num * plain_num)

float * unigram;
float * bigram;
float * trigram;

// Maps characters to the integers representations

map<char, int> inv_cipher;
map<char, int> inv_plain;

// These are CUDA variables (defined here so that both the main program and the
// CUDA device can see them).
// greenhouse and backpointers are the tables - greenhouse holds the probabilities, 
// and backpointers the back pointers.
// cuda(Uni|Bi|Tri)grams are the CUDA versions of the above unigram, bigram, and trigram
// arrays.

__device__ float * greenhouse;
__device__ int * backpointers;
__device__ float * cudaUnigrams;
__device__ float * cudaBigrams;
__device__ float * cudaTrigrams;

// This is just a bunch of stuff put in to make the CUDA part of the program work right.

texture<float, 1, cudaReadModeNormalizedFloat> texRef;
cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
const textureReference* texRefPtr;

/* ------------------------------End Global Definitions ------------------------------------ */

/* ---------------------------------------- Classes ---------------------------------------- */

// This class is just a comparison function for the A* heap.  It looks as pairs of probailities
// and partial solutions, and orders them according to the probability.  On a tie, the larger 
// solution is given precedence.
// (Note: the probabilities are in the negative log domain, so will seem backwards.)

class pair_compare
{

public:
  bool operator() (const pair<float, map<int, int> > &lhs, const pair<float, map<int, int> > &rhs) const
  {
    return (lhs.first > rhs.first) || ((lhs.first == rhs.first) && (lhs.second.size() < rhs.second.size()));
  }
};

/* ------------------------------------ End Classes ---------------------------------------- */

/* ----------------------------------- CUDA Functions -------------------------------------- */

// These parts are the actual function calls to the CUDA architecture.  Each one fills one
// plaintext cell (k) of one row (l) of the table. The variables k and l are called in
// parallel, so they aren't part of the arguments: instead, they're taken from the number
// of the thread created. This takes up the first two lines of all of the functions.

// This fills the first column of the table.  It just gives every cell its corresponding 
// unigram probability.  It gives zero probability if the cell is inconsistant with a solution
// (i.e., if 'a' maps to 'b' in the partial solution, the cell for 'a' maps to 'c' goes to -1),
// or if it is in conflict with the path through that cell (i.e., if we are lookong at a table 
// indicting the probabilities for 'g' mapping to 'h' *at that point*, we cannot map 'g' to 'j',
// as well.)
// Note: the partial solution used is the given solution plus one extra letter, given in the 
// arguments and runin parallel.
// 
// Parameters:
//   - greenhouse: the array holding the probabilities for the paths.  
//                   Note (1): To save space, we only store 3 values on the 'i' dimension.
//                             We wrap around using i % 3.
//                   Note (2): It is one - dimensional: greenhouse[i, l, k] is taken as
//                             *(greenhouse + (((i) % 3) * plain_num + l) * plain_num + k).
//   -plain_num: the number of plaintext letters.
//   -b: the next ciphertext letter to be added.
//   -part_soln: an array holding the current guess of the partial solution.  It always holds
//               all of the plaintext letters, but the ones that are not yet set all point to -1.
//               (indices are plaintext letters, values are ciphertext)
//   -part_inv: same as above, but ciphertext to plain text.
//   -c: the cihertext letter at the beginning of the cipher (i.e., it the ciphertext is '5, 2, 5, 6'
//       after being tranformed into ints, c here would be '5').
//   -unigram: the unigram array.
//   -countFnCalls: a variable used to determine l and k after breaking the thread blocks into 
//                  manageable sizes.  Basically, it just lets the program be used on CUDA.
//                  we can only call a certain number of threads at once, so we call blocks of
//                  fixed sizes.
//                  This is the number of blocks that have been called.
//   -numK: same as above, but this is the size of the blocks.
//
// Returns: Nothing, but changes greenhouse (the CUDA version) in place.

__global__ void rowZero(float* greenhouse, int plain_num, int b, int * part_soln, int * part_inv, int c, float * unigram, int countFnCalls, int numK){
  int l = threadIdx.x;
  int k = __fadd_rn(threadIdx.y, __fmul_rn(countFnCalls, numK));
  if (k < plain_num){
    int temp = __fadd_rn(k, __fmul_rn(plain_num, l));
    // The following check is common to all CUDA functions.
    // (l == k) == (c == b) is true iff either both ciphertexts
    // and both plaintexts are the same, of if both are different.
    // The other four checks determine compatibility with the partial solution.

    if(((l == k) == (c == b)) 
       && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
       && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
       && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
       && ((*(part_inv + b) == -1) || (*(part_inv + b) == k)) ){
      *(greenhouse + temp) = *(unigram + l);
    } else {
      // In the log linear domain, -1 represents a zero probability.

      *(greenhouse + temp) = -1;
    }
  }
}

// This function is the same as above, but it fills the second column (index one).  It runs
// the same checks as before, and adds checks of the bigram probabilities.
//
// Parameters: 
//   same as rowZero, except for the following:
//   - backpointers: the array holding the backpointers.  The same notes apply as for the
//                   greenhouse array.
//   -bigram: the CUDA version of the bigram probabilities.
//
// Returns: Nothing, but changes the greenhouse and backpointers arrays in place.

__global__ void rowOne(float* greenhouse, int* backpointers, int plain_num, int b, int * part_soln, int * part_inv, int c, float * unigram, float * bigram, int countFnCalls, int numK){
  int l = threadIdx.x;
  int k = __fadd_rn(threadIdx.y, __fmul_rn(countFnCalls, numK));
  if (k < plain_num){
    int temp = __fadd_rn(k, __fmul_rn(plain_num, __fadd_rn(l, plain_num)));
    // Same as rowZero.
    // The last check is to ensure that there is a unigram count for the chosen 
    // plaintext (likely enough that the check may be unneeded).

    if(((l == k) == (c == b)) 
       && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
       && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
       && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
       && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
       && (*(unigram + l) > 0) ){
      int j;
      float best_prob = -1;
      int back = -1;
      int gHptr = k;
      int biPtr = l;
      // Loop over possible previous letters.

      for(j = 0; j < plain_num; j++){
        if((*(bigram + biPtr) > 0) && (*(greenhouse + gHptr) > 0)
           && ((best_prob < 0) 
               || (best_prob > *(bigram + biPtr) + *(greenhouse + gHptr)))){
              best_prob = *(bigram + biPtr) + *(greenhouse + gHptr);
              back = j;
        }
        gHptr = __fadd_rn(plain_num, gHptr);
        biPtr = __fadd_rn(plain_num, biPtr);
      }
        *(greenhouse + temp) = best_prob;
        *(backpointers + temp) = back;
      } else {
        *(greenhouse + temp) = -1;
        *(backpointers + temp) = -1;
    }
  }
}

// This function is the same as above, but it fills the third column (index two).  It runs
// the same checks as before, but adds trigram probability checks.  Also, it adds consistency
// checks on the letters: i.e., if the ciphertest is ' ..., 7, 9, 7, ... ', then both sevens
// are forced to be the same plaintext letters.  This window only lasts for two characters.
//
// Note: the loops and if statements are commented for readability.  I should probably have 
//       seperated this into smaller functions, but wanted to reduce the number of calls.
//       The loops are over indices j and j2: j is the plaintext of the previous letter
//       (index 1), and j2 is the second last plaintext letter (index 0).
//       The branches of the if statements are as follows:
//         if # 1: same as in rowZero.  Checks consistency of path and with partial solutions.
//                 If this doesn't pass, we assign zero probability.
//           if # 1.1 (1st branch of #1): Check if c1 (previous letter) is fixed in the solution.
//                                        If so, we only need to check the plaintext letter that
//                                        maps to it.
//             if # 1.1.1: Checks that the the value at the cell greenhouse(1, j, k) is nonzero
//                         (i.e., there is any path through that cell).  If this fails, we give zero
//                         probability, since this is the only path.
//               if # 1.1.1.1: Check to see if the ciphertext at c2 is in the partial solution.
//                             If it is, we only have to check one plaintext letter.  This is similar to
//                             check # 1.1.
//                 if # 1.1.1.1.1: A few checks here: Checks that there is a path through
//                                 greenhouse(0, j2, k), like in check 1.1.1, and checks the
//                                 bigram and trigram constraints.  Doesn't need to check the
//                                 '..., 7, 9, 7, ...' constraint, since this is covered by the partial 
//                                 solution.
//                                 If this fails, give zero probability.
//               else # 1.1.1.1: c2 is not known, so check all possible j2:
//                 if # 1.1.1.2.1 (second branch of fourth if): We're looping over j2s here, so do the same
//                                                              checks as in 1.1.1.1.1, but add the 
//                                                              '..., 7, 9, 7, ...' constraint, and check to
//                                                              see if the current path is both possible and
//                                                              beats anything we've checked before.  If it fails,
//                                                              we just go to the next j2.
//           else # 1.1: c1 is not known, so check all possible j:
//             if # 1.1.2: Check that there is a path through greenhouse(1, j, k), as in check #1.1.1.  Also, check
//                         the '..., 7, 9, 7, ...' constraint, and that the probability at greenhouse(1, j, k) is high enough
//                         to actually give us a better solution.  If this fails, go to the next j.
//               if # 1.1.2.1: Same as 1.1.1.1, different branch.
//                  if # 1.1.2.1.1: Same as 1.1.1.1.1, different branch.
//               else # 1.1.2.1: c2 is not known, just like else # 1.1.1.1.
//                  if # 1.1.2.2.1: Same as 1.1.1.2.1, different branch.
//
// Parameters: 
//   same as rowOne, except for the following:
//   - trigram: the CUDA array for the trigram probabilities.
//   - c1, c2: the cipher letters at the previous two indices (c1 is index 1, c2 is index 0).
//
// Returns: Nothing, but changes the greenhouse and backpointers arrays in place.

__global__ void rowTwo(float* greenhouse, int* backpointers, int plain_num, int b, int c, int c1, int c2, int * part_soln, int * part_inv, float * unigram, float * bigram, float * trigram, int countFnCalls, int numK){
  int l = threadIdx.x;
  int k = __fadd_rn(threadIdx.y, __fmul_rn(countFnCalls, numK));
   if (k < plain_num){
     int temp = __fadd_rn(k, __fmul_rn(plain_num, __fadd_rn(l, __fmul_rn(plain_num, 2))));
     if(((l == k) == (c == b)) // if #1
        && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
        && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
        && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
        && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
        && (*(unigram + l) > 0) ){
        int j;
        int j2;
        float best_prob = -1;
        int back = -1;
        if(*(part_inv + c1) >= 0){ // if #1.1
          j = *(part_inv + c1);
          if(*(greenhouse + (plain_num + j) * plain_num + k) > 0){ // if #1.1.1
            if(*(part_inv + c2) >= 0){ // if # 1.1.1.1
              j2 = *(part_inv + c2);
              if((*(greenhouse + j2 * plain_num + k) >= 0)
                 && (*(bigram + j + plain_num * j2) >= 0)
                 && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)){ // if # 1.1.1.1.1
                best_prob = *(greenhouse + (plain_num + j) * plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                back = j;
              } // if # 1.1.1.1.1
            } else { // if # 1.1.1.1
              for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                if((*(greenhouse + (plain_num + j2) * plain_num + k) >= 0)
                   && ((j2 == l) == (c == c2))
                   && ((j == j2) == (c1 == c2))
                   && (*(bigram + j + plain_num * j2) >= 0)
                   && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.1.2.1
                   && ((best_prob < 0) 
                       || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + j2 * plain_num + k)))){
                  best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + j2 * plain_num + k);
                  back = j;
                } // if #1.1.1.2.1 
              } // for: index j2 
            } // if #1.1.1.1
          } // if # 1.1.1
        } else { // if # 1.1
           for(j = 0; j < plain_num; j++){ // for: index j
             if((*(greenhouse + (plain_num + j) * plain_num + k) >= 0)
                && ((j == l) == (c == c1))
                && ((best_prob < 0) || (*(greenhouse + (plain_num + j) * plain_num + k) < best_prob))){ // if #1.1.2
               if(*(part_inv + c2) >= 0){ // if # 1.1.2.1
                 j2 = *(part_inv + c2);
                 if((*(greenhouse + j2 * plain_num + k) >= 0)
                    && (*(bigram + j + plain_num * j2) >= 0)
                    && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)){ // if # 1.1.2.1.1
                   best_prob = *(greenhouse + (plain_num + j) * plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                   back = j;
                 } // if # 1.1.2.1.1
               } else { // if # 1.1.2.1
                 for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                   if((*(greenhouse + j2 * plain_num + k) >= 0)
                     && ((j2 == l) == (c == c2))
                     && ((j == j2) == (c1 == c2))
                     && (*(bigram + j + plain_num * j2) >= 0)
                     && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.2.2.1
                     && ((best_prob < 0) 
                         || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + j2 * plain_num + k)))){
                    best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + j2 * plain_num + k);
                    back = j;
                  } // if #1.1.2.2.1 
                } // for: index j2 
              } // if #1.1.2.1
            } // if # 1.1.2
          } // for: index j
        } // if #1.1
      *(greenhouse + temp) = best_prob;
      *(backpointers + temp) = back;
    } else { // if # 1
      *(greenhouse + temp) = -1;
      *(backpointers + temp) = -1;
    } // if # 1
  }
}

// This function is the same as above, but it fills the every column column after the third (index > 2).  
// If is exactly the same as rowTwo, except now it starts to use the backpointer information.  Otherwise,
// the loops, checks, and arguments are almost exactly the same.  There are a few places where the pointers
// are fixed to save a bit of time, too.  Those have been done only on this function since the last three 
// functions are called once per run, while this one is called ~cipher_length times.
//
// Note: the loops and if statements are commented for readability.  I should probably have 
//       seperated this into smaller functions, but wanted to reduce the number of calls.
//       The loops are over indices j and j2: j is the plaintext of the previous letter
//       (index 1), and j2 is the second last plaintext letter (index 0).
//       The branches of the if statements are as follows:
//         if # 1: same as in rowZero.  Checks consistency of path and with partial solutions.
//                 If this doesn't pass, we assign zero probability.
//           if # 1.1 (1st branch of #1): Check if c1 (previous letter) is fixed in the solution.
//                                        If so, we only need to check the plaintext letter that
//                                        maps to it.
//             if # 1.1.1: Checks that the the value at the cell greenhouse(i - 1, j, k) is nonzero
//                         (i.e., there is any path through that cell).  If this fails, we give zero
//                         probability, since this is the only path.
//               if # 1.1.1.1: Check to see if the ciphertext at c2 is in the partial solution.
//                             If it is, we only have to check one plaintext letter.  This is similar to
//                             check # 1.1.
//                 if # 1.1.1.1.1: A few checks here: Checks that there is a path through
//                                 greenhouse(i - 2, j2, k), like in check 1.1.1, and checks the
//                                 bigram and trigram constraints.  Doesn't need to check the
//                                 '..., 7, 9, 7, ...' constraint, since this is covered by the partial 
//                                 solution.
//                                 If this fails, give zero probability.
//               else # 1.1.1.1: c2 is not known, so check all possible j2:
//                 if # 1.1.1.2.1 (second branch of fourth if): We're looping over j2s here, so do the same
//                                                              checks as in 1.1.1.1.1, but add the 
//                                                              '..., 7, 9, 7, ...' constraint, and check to
//                                                              see if the current path is both possible and
//                                                              beats anything we've checked before.  If it fails,
//                                                              we just go to the next j2.
//           else # 1.1: c1 is not known, so check all possible j:
//             if # 1.1.2: Check that there is a path through greenhouse(i - 1, j, k), as in check #1.1.1.  Also, check
//                         the '..., 7, 9, 7, ...' constraint, and that the probability at greenhouse(i - 1, j, k) is high 
//                         enough to actually give us a better solution.  If this fails, go to the next j.
//               if # 1.1.2.1: Same as 1.1.1.1, different branch.
//                  if # 1.1.2.1.1: Same as 1.1.1.1.1, different branch.
//               else # 1.1.2.1: c2 is not known, just like else # 1.1.1.1.
//                  if # 1.1.2.2.1: Same as 1.1.1.2.1, different branch.
//
// Parameters: 
//   same as rowTwo, except:
//   -i: the index of the current column.
//
// Returns: same as rowTwo.

__global__ void rowThreePlus(float* greenhouse, int* backpointers, int plain_num, int b, int c, int c1, int c2, int i, long cipher_length, int * part_soln, int * part_inv, float * unigram, float * bigram, float * trigram, int countFnCalls, int numK){
  int l = threadIdx.x;
  int k = __fadd_rn(threadIdx.y, __fmul_rn(countFnCalls, numK));
  int j;
  int j2;
  int back = -1;
  int step1;
  int triPtr;
  int lastPtr;
  float best_prob = -1;
  if (k < plain_num){
    int temp = __fadd_rn(k, __fmul_rn(plain_num, __fadd_rn(l, __fmul_rn(plain_num, i % 3))));
            if(((l == k) == (c == b)) // if #1
               && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
               && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
               && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
               && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
               && (*(unigram + l) > 0) ){
            step1 = plain_num * ((i - 1) % 3);
            if(*(part_inv + c1) >= 0){ // if #1.1
              j = *(part_inv + c1);
              lastPtr = __fadd_rn(__fmul_rn(__fadd_rn(step1, j), plain_num), k);
              if(*(greenhouse + lastPtr) > 0){ // if #1.1.1
                if(*(part_inv + c2) >= 0){ // if # 1.1.1.1
                  j2 = *(part_inv + c2);
                  triPtr = __fadd_rn(l, __fmul_rn(plain_num, (__fadd_rn(j, __fmul_rn(plain_num, j2)))));
                  if((*(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0)
                       && (*(trigram + triPtr) >= 0)){ // if # 1.1.1.1.1
                    best_prob = *(greenhouse + lastPtr) + *(trigram + triPtr);
                    back = j;
                  } // if # 1.1.1.1.1
                } else { // if # 1.1.1.1
                  for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                    if((*(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0)
                        && (*(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0)
                        && ((j2 == l) == (c == c2))
                        && ((j == j2) == (c1 == c2))
                        && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.1.2.1
                        && ((best_prob < 0) 
                              || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k))) + *(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k)))){
                        best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k))) + *(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k);
                        back = j;
                    } // if #1.1.1.2.1 
                  } // for: index j2 
                } // if #1.1.1.1
              } // if # 1.1.1
            } else { // if # 1.1
              for(j = 0; j < plain_num; j++){ // for: index j
                if((*(greenhouse + (((i - 1) % 3) * plain_num + j) * plain_num + k) >= 0)
                   && ((j == l) == (c == c1))
                   && ((best_prob < 0) || (*(greenhouse + (((i - 1) % 3) * plain_num + j) * plain_num + k) < best_prob))){ // if #1.1.2
                  if(*(part_inv + c2) >= 0){ // if # 1.1.2.1
                    j2 = *(part_inv + c2);
                    if((*(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0) // if # 1.1.2.1.1
                       && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)
                       && ((best_prob < 0) || (best_prob > *(greenhouse + (((i - 1) % 3) * plain_num + j) * plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2))))){ 
                      best_prob = *(greenhouse + (((i - 1) % 3) * plain_num + j) * plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                      back = j;
                    } // if # 1.1.2.1.1
                  } else { // if # 1.1.2.1
                    for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                      if((*(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0)
                        && (*(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k) >= 0)
                        && ((j2 == l) == (c == c2))
                        && ((j == j2) == (c1 == c2))
                        && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.2.2.1
                        && ((best_prob < 0) 
                              || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k))) + *(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k)))){
                        best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + (((i - 2) % 3) * plain_num + j2) * plain_num + k))) + *(greenhouse + (((i - 2) % 3) * plain_num + j2) * plain_num + k);
                        back = j;
                      } // if #1.1.2.2.1 
                    } // for: index j2 
                  } // if #1.1.2.1
                } // if # 1.1.2
              } // for: index j
            } // if #1.1
            *(greenhouse + temp) = best_prob;
            *(backpointers + temp) = back;
          } else { // if # 1
            *(greenhouse + temp) = -1;
            *(backpointers + temp) = -1;
          } // if # 1
  }
}

/* --------------------------------End CUDA Functions -------------------------------------- */

/* -------------------------------- General Functions -------------------------------------- */

// This is the actual generalized Viterbi algorithm.  The conceptual details should already be known, so the details are:
//   (1) set up the block size to determine the number of threads per block.
//   (2) run the first column in parallel.
//   (3) run the second column in parallel.
//   (4) run the third column in parallel.
//   (4) succesively run the every remaining column in parallel.
//   (5) once finished, find the best probability over every remaining plaintext letter k,
//       where k points to the next ciphertext letter to be fixed.  This is done by looking at 
//       greenhouse(cipher_length - 1, l, k) for every possible l, and taking the best result.
//       put these into the result array.
//
//   The columns are treated differently for these reasons:
//     Column one can only use the unigram statistics.
//     Column two can only use bigram statistics.
//     Column three can't use back pointers.
//     The remaining columns can use everything.
//     I've seperated them into different functions, rather than different branches, so that the time taken to 
//     decide which branch to use is mininmized.
//
// Parameters:
//   -b: The next ciphertext letter to be fixed.
//   -plain_num: the number of plaintext letters.
//   -cipher_num: the number of ciphertext letters.
//   -part_soln: an array holding the current guess of the partial solution.  It always holds
//               all of the plaintext letters, but the ones that are not yet set all point to -1.
//               (indices are plaintext letters, values are ciphertext).  This is a CUDA array, but isn't
//               listed as such.
//   -cuda_inv_soln_arr: An inverse solution - like the above parameter, but maps ciphertext to plaintext.
//   -cipher_string: The array of integers representing the cipher.
//   -cipher_length: the length of the cipher.
//   -result: an array of length plain_num giving, for each plaintext letter, the best possible probability
//            of a path that assigns that plaintext letter to b, the next ciphertext letter to be fixed.
//    This is basically an array to be changed in place for the return values.
//   -(uni|bi|tri)gram: CUDA arrays representing the unigrams, bigrams, and trigrams.
//   -greenhouse, backpointers: CUDA arrays holding the greenhouse table and the associated backpointers.
//    the greenhouse table holds the probabilities, and the backpointer table holds the back pointers.
//   -tempResult: a table meant to hold the final greenhouse table so it can be refined into the results array.
//    (recall: (1) We can't read directly from a CUDA array, so we have to transfer it. (2) If we pass it in this way,
//    we can reuse it, and so save a hassle with memory management.)
//    -totalSize: The size of the greenhouse table.  We need to have it to transfer between the greenhouse and the 
//               tempResult table.
//
// Returns: Nothing, but alters the results array in place.

void genviterbi(int b, int plain_num, int cipher_num, int * part_soln, int * cuda_inv_soln_arr, int * cipher_string, long cipher_length, float * result, float * unigram, float * bigram, float * trigram, float * greenhouse, int * backpointers, float * tempResult, int totalSize){

// Variables:
//   -i: the index of the column of the greenhouse and backpointer tables that are being filled.
//   -(c|c1|c2): The current / previous / second previous letter in the cipher at index i.
//   -(l|j): loop variables used to iterate over all plaintext letters when filling the result
//           array.
//   -countFnCalls: Used when breaking the CUDA threads into blacks.  Counts the number of blocks
//    used.

  int i;
  int c = 0;
  int c1 = 0;
  int c2 = 0;
  int l;
  int j;
  int countFnCalls = 0;

  // CUDA architecture drops function calls if there more than 512 of them.  
  // This section determines how large the chunks we can use are.
  // We're assuming that the size of the alpohabet (plain_num) in much
  // less than the size of the blocks (threshold).  This works for English,
  // where we're dealing with ~54 at max, but we will likely need to change how
  // the blocks are called when we switch to Hindi (hundreds of characters).

  int threshold = 512;
  int numK = threshold / plain_num;

  // A thread block used to call the CUDA threads.

  dim3 db(plain_num, numK);

  // Fill the first column.

  c = *cipher_string;
  while  (countFnCalls * numK < plain_num){
    rowZero<<< 1, db >>>(greenhouse, plain_num, b, part_soln, cuda_inv_soln_arr, c, unigram, countFnCalls, numK);
    countFnCalls++;
  }

  // Fill the second column.

  countFnCalls = 0;
  c = *(cipher_string + 1);
  while  (countFnCalls * numK < plain_num){
    rowOne<<< 1, db >>>(greenhouse, backpointers, plain_num, b, part_soln, cuda_inv_soln_arr, c, unigram, bigram, countFnCalls, numK);
    countFnCalls++;
  }
  
  // Fill the third column.

  countFnCalls = 0;
  c = *(cipher_string + 2);
  c1 = *(cipher_string + 1);
  c2 = *(cipher_string);
  while  (countFnCalls * numK < plain_num){
    rowTwo<<< 1, db >>>(greenhouse, backpointers, plain_num, b, c, c1, c2, part_soln, cuda_inv_soln_arr, unigram, bigram, trigram, countFnCalls, numK);
    countFnCalls++;
  }  
  
  // Fill the remaining columns.

  for(i = 3; i < cipher_length; i++){
    countFnCalls = 0;
    c = *(cipher_string + i);
    c1 = *(cipher_string + i-1);
    c2 = *(cipher_string + i-2);
    while  (countFnCalls * numK < plain_num){
      rowThreePlus<<< 1, db >>>(greenhouse, backpointers, plain_num, b, c, c1, c2, i, cipher_length, part_soln, cuda_inv_soln_arr, unigram, bigram, trigram, countFnCalls, numK);
      countFnCalls++;
    }  
  }

  // Copy the greenhouse to the tempResult table.

  cudaMemcpy(tempResult, greenhouse, (totalSize * sizeof(float)), cudaMemcpyDeviceToHost);

  // Fill the result array.

  for(l = 0; l < plain_num; l++){
    *(result + l) = -1;
  }
  for(l = 0; l < plain_num; l++){
    for(j = 0; j < plain_num; j++){
      if(*(tempResult + (((cipher_length - 1) % 3) * plain_num + l) * plain_num + j) > 0){
        if(*(result + j) <= 0){
          *(result + j) = *(tempResult + (((cipher_length - 1) % 3) * plain_num + l) * plain_num + j);
        } else {
          if(*(result + j) > *(tempResult + (((cipher_length - 1) % 3) * plain_num + l) * plain_num + j)){
            *(result + j) = *(tempResult + (((cipher_length - 1) % 3) * plain_num + l) * plain_num + j);
          }
        }
      }
    }
  }
}

/* ---------------------------- End General Functions -------------------------------------- */


/* ------------------------------------ Main------------------------------------------------ */

// The main program.
//
// Calling procedure: If this program is compiled with the name "greenhouse", type
//   >> greenhouse profile
// in the command line, where profile is the name of the profile to be used.
//
//   profiles have the format:
//   line 1: all ciphertext letters, seperated by spaces.
//   line 2: all plaintext letters, seperated by spaces.
//   line 3: the name of the unigram file.
//   line 4: the name of the bigram file.
//   line 5: the name of the trigram file.
//   line 6: the name of the ciphertext file.
//   Current profiles also have a seventh line that gives the
//   location of the plaintext file, which can be used for setup and evaluation
//   purposes.  It won't be read by this program, though.
//
// Note: this program returns a set of stats for each run of the viterbi algorithm, and
//       also the final solution, number of runs per solution size, and the total running time.
//
// e.g., the final lines are of the format:
// 
//Time taken for process: 56:52:38
//
//  Solution: {0:0,1:1,2:3,3:4,4:5,5:6,6:7,7:8,8:9,9:10,10:11,11:12,12:13,13:14,14:15,15:16,16:17,17:18,18:19,19:20,20:21,21:22,22:23,23:24,24:25,25:26,26:27,27:2}
//
//  Solution sizes: {1:1,2:24,3:340,4:938,5:6090,6:27485,7:13285,8:19922,9:10783,10:6764,11:3143,12:202,13:15,14:4,15:1,16:4,17:3,18:2,19:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1}
// 
// (this was a large run.)
//
// while the information given in a particular run of the algorithm is of the fromat:
//
//Starting Viterbi Algorithm: 
//
//  Queue size: 418
//
//  Solution size: 3
//
//  Number of passes: 61
//
//  Solution sizes: {1:1,2:24,3:340,4:113}
// 
// Information on the number of zeroed solutions can be obtained by looking at the information given in the different runs of the algorithm.
// Correctness of the solution is not checked by this algorithm.
// Both of these jobs are done by simple python scripts.

int main(int argc, char * argv[]){

  // A generic loop variable.

  int i;

  // This will be the total size of the greenhouse table.

  int totalSize;

  // The array used to read the argument in the commandline function call.

  char * s = "profile";

  // make sure that the calling procedure is correct.

  if (argc != 2){
    //printf("\nUsage: ");
    //printf("\t<insert usage description here>\n\n");
    //exit(0);
  } else {
    s = argv[1];
  }

  // read profile
  ifstream f_profile(s);

  // Lines from the profile.
  // profiles have the format:
  // line 1: all ciphertext letters, seperated by spaces.
  // line 2: all plaintext letters, seperated by spaces.
  // line 3: the name of the unigram file.
  // line 4: the name of the bigram file.
  // line 5: the name of the trigram file.
  // line 6: the name of the ciphertext file.
  // Current profiles also have a seventh line that gives the
  // location of the plaintext file, which can be used for setup and evaluation
  // purposes.  It won't be read by this program, though.

  string cipher_line;
  string plain_line;
  string unigram_name;
  string bigram_name;
  string trigram_name;
  string cipher_name;

  getline(f_profile, cipher_line, '\n');
  getline(f_profile, plain_line, '\n');
  getline(f_profile, unigram_name, '\n');
  getline(f_profile, bigram_name, '\n');
  getline(f_profile, trigram_name, '\n');
  getline(f_profile, cipher_name, '\n');
  
  f_profile.close();

  // We don't know the sizes of the alphabets beforehand, so we're reading them into 
  // vectors.  Recall that they're only being used for setup: in the actual algorithm, 
  // we'll only use numbers.
  //
  // Use the cipher and plaintext lines to create the corresponding lists:

  for(i = 0; i < cipher_line.size(); i++){
    if(cipher_line[i] != ' '){
      cipher_alpha.push_back(cipher_line[i]);
    }
  }
  cipher_num = cipher_alpha.size();
  for(i = 0; i < cipher_num; i++){
    inv_cipher[cipher_alpha[i]] = i;
  }
  for(i = 0; i < plain_line.size(); i++){
    if(plain_line[i] != ' '){
      plain_alpha.push_back(plain_line[i]);
    }
  }
  plain_num = plain_alpha.size();
  for(i = 0; i < plain_num; i++){
    inv_plain[plain_alpha[i]] = i;
  }

  // Read the unigrams, bigrams, and trigrams

  // for the unigrams, read each line, process it into its parts,
  // and add each entry into a unigram map.
  // The basic map has size plain_num, with a default entry of -1.
  // unigram is an array whose ith index is *(unigram + i)

  unigram = (float *) malloc(plain_num * sizeof(float));
  for(i = 0; i < plain_num; i++){
    *(unigram + i) = -1;
  }
  ifstream unigrams_file(unigram_name.c_str());
  string temp;
  float uni_total = 0;
  while(getline(unigrams_file, temp, '\n')){
    char char_1;
    i = 0;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_1 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    string str = "";
    while((temp[i] != ' ') && (temp[i] != '\t')){
      str.push_back(temp[i]);
      i++;
    }
    float uni_prob = atof(str.c_str());
    *(unigram + inv_plain[char_1]) = -1 * log(uni_prob);
    uni_total += uni_prob;
  }
  for(i = 0; i < plain_num; i++){
     if((*(unigram + i) != -1) && isfinite(*(unigram + i))){
       *(unigram + i) += log(uni_total);
     } else {
       *(unigram + i) = -1;
     }
  }
  unigrams_file.close();

  // Read the bigram file.
  // This is an array of size plain_num^2, with default values -1.
  // the index for i following j is *(bigram + plain_num * j + i).

  int p_2 = plain_num * plain_num;
  bigram = (float *) malloc(p_2 * sizeof(float));
  for(i = 0; i < p_2; i++){
    *(bigram + i) = -1;
  }
  ifstream bigrams_file(bigram_name.c_str());
  float bi_total = 0;
  while(getline(bigrams_file, temp, '\n')){
    char char_1;
    char char_2;
    i = 0;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_1 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_2 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    string str = "";
    while((temp[i] != ' ') && (temp[i] != '\t')){
      str.push_back(temp[i]);
      i++;
    }
    float bi_prob = atof(str.c_str());
    *(bigram + inv_plain[char_2] + (inv_plain[char_1] * plain_num)) = -1 * log(bi_prob);
    bi_total += bi_prob;
  }
  for(i = 0; i < p_2; i++){
     if((*(bigram + i) != -1) && isfinite(*(bigram + i))){
       // TODO: I think it shouldn't normalize like this, but experiments showed better results with this
       *(bigram + i) += log(bi_total);
     } else {
       *(bigram + i) = -1;
     }
  }
  bigrams_file.close();

  // Read the trigram file.
  // This is an array of size plain_num^3, with default values -1.
  // the index for i following j following k is *(trigram + plain_num * (plain_num * k + j) + i).

  int p_3 = plain_num * plain_num * plain_num;
  trigram = (float *) malloc(p_3 * sizeof(float));
  for(i = 0; i < p_3; i++){
    *(trigram + i) = -1;
  }
  ifstream trigrams_file(trigram_name.c_str());
  float tri_total = 0;
  while(getline(trigrams_file, temp, '\n')){
    char char_1;
    char char_2;
    char char_3;
    i = 0;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_1 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_2 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    char_3 = temp[i];
    i++;
    while((temp[i] == ' ') || (temp[i] == '\t')){
      i++;
    }
    string str = "";
    while((temp[i] != ' ') && (temp[i] != '\t')){
      str.push_back(temp[i]);
      i++;
    }
    float tri_prob = atof(str.c_str());
    *(trigram + inv_plain[char_3] + (inv_plain[char_2] * plain_num) + (inv_plain[char_1] * p_2)) = -1 * log(tri_prob);
    tri_total += tri_prob;
  }
  for(i = 0; i < p_3; i++){
     if((*(trigram + i) != -1) && isfinite(*(trigram + i))){
       // TODO: I think it shouldn't normalize like this, but experiments showed better results with this
       *(trigram + i) += log(tri_total);
     } else {
       *(trigram + i) = -1;
     }
  }
  trigrams_file.close();

  // Create the cipher string.
  // Recall that we're going to turn everything into numbers, and first read everything into a vector,
  // the transfer the vector to an array.

  ifstream cipher_file(cipher_name.c_str());
  while( getline(cipher_file, temp, '\n')){
    for(i = 0; i < temp.size(); i++){
      if(temp[i] != ' '){
        cipher_string_vec.push_back(inv_cipher[temp[i]]);
      }
    }
  }
  cipher_length = cipher_string_vec.size();
  cipher_string = (int *) malloc(cipher_length * sizeof(int));
  for(i = 0; i < cipher_length; i++){
    *(cipher_string + i) = cipher_string_vec[i];
  }
  cipher_file.close();

  // get the locations of the last characters in the cipher.
  // This will be used to determine the order in which letters are 
  // added to the solution.

  int num_last = 0;
  map<int, int> last;
  map<int, int> lastcount;
  for(i = cipher_length - 1; i >= 0; i--){
    if(last.count(*(cipher_string + i)) == 0){
      last[*(cipher_string + i)] = i;
      lastcount[num_last] = i;
      num_last++;
    }
  }

  // Set up the table in the GPU, and copy out the unigrams, bigrams, and trigrams into their own tables here.
  // Put the backpointers table there, too.

  totalSize = cipher_num * plain_num * plain_num * 3;

  int * test_array;
  cudaMalloc((void**)&test_array, (plain_num * sizeof(int)));

  float * result = (float *) malloc(plain_num * sizeof(float));
  float * tempResult = (float *) malloc(totalSize * sizeof(float));

  // Create the greenhouse array.

  cudaMalloc((void**)&greenhouse, (totalSize * sizeof(float)));
  cudaMalloc((void**)&backpointers, (totalSize * sizeof(int)));

  cudaMemcpy(greenhouse, tempResult, (totalSize * sizeof(float)), cudaMemcpyHostToDevice);

  // Copy the unigrams, bigrams, and trigrams into the GPU.

  cudaMalloc((void**)&cudaUnigrams, (plain_num * sizeof(float)));
  cudaMemcpy(cudaUnigrams, unigram, (plain_num * sizeof(float)), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&cudaBigrams, (plain_num * plain_num * sizeof(float)));
  cudaMemcpy(cudaBigrams, bigram, (plain_num * plain_num * sizeof(float)), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&cudaTrigrams, (plain_num * plain_num * plain_num * sizeof(float)));
  cudaMemcpy(cudaTrigrams, trigram, (plain_num * plain_num * plain_num * sizeof(float)), cudaMemcpyHostToDevice);

  // Make sure the part_soln and part_inv are in the GPU here.

  int * inv_soln_arr = (int *) malloc(cipher_num*sizeof(int));  
  int * curr_soln_arr = (int *) malloc(plain_num * sizeof(int));
  int * cuda_curr_soln_arr;
  cudaMalloc((void**)&cuda_curr_soln_arr, (plain_num * sizeof(int)));
  int * cuda_inv_soln_arr;
  cudaMalloc((void**)&cuda_inv_soln_arr, (cipher_num * sizeof(int)));

  // Some stuff to get the GPU to work the way we want.

  cudaGetTextureReference(&texRefPtr, "texRef");
  texRef.addressMode[0] = cudaAddressModeWrap;
  texRef.addressMode[1] = cudaAddressModeWrap;
  texRef.filterMode     = cudaFilterModeLinear;
  texRef.normalized     = true;

  // set up priority heap.
  // set up a map of sizes.
  // Set up a starting solution and add it to the heap.

  map<int, int> start_soln;

  // Normally, the solution is empty. here.
  // The following line adds the restriction that spaces map to spaces
  // to the solution (the ciphers we're using always end in a space).
  // Comment or uncomment this line as needed.

  start_soln[0] = 0;
  map<int, int> soln_sizes;
  soln_sizes[start_soln.size()] = 1;
  long pass_num = 1;
  
  // A boolean variable that will tell us when we're finished looking for the solution.

  bool found = false;

  // prime the priority queue with the first solution.

  int l;
  time_t start_time = time (NULL);   
  priority_queue< pair<float, map<int, int> >, vector<pair<float, map<int, int> > >, pair_compare > aStar;
  aStar.push(pair<float, map<int, int> >(-1, start_soln));

  // The main part of the program: pop solutions and run the genviterbi algorithm to grow solutions until
  // the final solution is found.

  while(!(aStar.empty()) && !found){

    // pop the best solution, and set up the next run.

    map<int, int> curr_soln = aStar.top().second;
    float curr_soln_prob = aStar.top().first;
    for(i = 0; i < plain_num; i++){  // Copy all curr soln to curr_soln_arr
      if(curr_soln.count(i) > 0){
        *(curr_soln_arr + i) = curr_soln[i];
      } else {
        *(curr_soln_arr + i) = -1;
      }
    }
    aStar.pop();

    // print some stats.

    cout << "Starting Viterbi Algorithm: \n" << endl;
    cout << "  Queue size: " << (aStar.size() + 1) << "\n" << endl;
    cout << "  Solution size: " << curr_soln.size() << "\n" << endl;
    cout << "  Solution probability " << curr_soln_prob << "\n" << endl;
    cout << "  Number of passes: " << pass_num++ << "\n" << endl;
    
    stringstream temp_soln;
    map<int, int>::iterator soln_iter;
    for(soln_iter = soln_sizes.begin(); soln_iter != soln_sizes.end(); ++soln_iter){
      temp_soln << "," <<  soln_iter->first << ":" << soln_iter->second;
    }
    string soln_string = temp_soln.str();
    soln_string.replace(0, 1, "{");
    soln_string += "}";
    cout << "  Solution sizes: " << soln_string << "\n" << endl;

    // Check to see if we're done (the solution is large enough to cover all ciphertext letters).
    // If so, set the solution to the current solution and exit the loop.
    // Otherwise, find the next letter to be added.

    if(curr_soln.size() >= cipher_num){
      found = true;
      start_soln = curr_soln;
    } else {
      int curr_endpoint = *(cipher_string + lastcount[curr_soln.size()]);

      // This part will create the inverse partial solution.

      for(i = 0; i < cipher_num; i++){
        *(inv_soln_arr + i) = -1;
      }
      for(i = 0; i < plain_num; i++){
        if(*(curr_soln_arr + i) >= 0){
          *(inv_soln_arr + *(curr_soln_arr + i)) = i;
        }
      }

      // Copy the partial solution and the inverse partial solution to the CUDA arrays.

      cudaMemcpy(cuda_curr_soln_arr, curr_soln_arr, (plain_num * sizeof(int)), cudaMemcpyHostToDevice);
      cudaMemcpy(cuda_inv_soln_arr, inv_soln_arr, (plain_num * sizeof(int)), cudaMemcpyHostToDevice);

      // run the actual genviterbi algorithm.

      genviterbi(curr_endpoint, plain_num, cipher_num, cuda_curr_soln_arr, cuda_inv_soln_arr, cipher_string, cipher_length, result, cudaUnigrams, cudaBigrams, cudaTrigrams, greenhouse, backpointers, tempResult, totalSize);

      // For each possible extension to the solution, format the next partial solution
      // and add it to the heap.

      for(l = 0; l < plain_num; l++){
        if(*(result + l) >= 0){
          map<int, int> next_soln;
          map<int, int>::iterator curr_iter;
          for(curr_iter = curr_soln.begin(); curr_iter != curr_soln.end(); ++curr_iter){
            next_soln[curr_iter->first] = curr_iter->second;
          }
          next_soln[l] = curr_endpoint;
          float next_prob = *(result + l);
          if(next_prob > 0){
            aStar.push(pair<float, map<int, int> >(next_prob, next_soln));
            int size_count = 1;
            if(soln_sizes.count(curr_soln.size() + 1) > 0){
              size_count = soln_sizes[curr_soln.size() + 1] + 1;
            }
            soln_sizes[curr_soln.size() + 1] = size_count;
          }
        }
      }
    }

  }

  // return the results.

  time_t stop_time = time (NULL);
  long total_time = stop_time - start_time;
  int seconds = total_time % 60;
  int minutes = ((total_time - seconds) / 60) % 60;
  long hours = (((total_time - seconds) / 60) - minutes) / 60;
  cout << "Time taken for process: " << hours << ":" << minutes << ":" << seconds << "\n"<< endl;
  if(start_soln.size() == cipher_num){
    stringstream soln_stream_1;
    map<int, int>::iterator soln_iter;
    for(soln_iter = start_soln.begin(); soln_iter != start_soln.end(); ++soln_iter){
      soln_stream_1 << "," <<  soln_iter->first << ":" << soln_iter->second;
    }
    string soln_string = soln_stream_1.str();
    soln_string.replace(0, 1, "{");
    soln_string += "}";
    cout << "  Solution: " << soln_string << "\n" << endl;

    stringstream soln_stream_2;
    for(soln_iter = soln_sizes.begin(); soln_iter != soln_sizes.end(); ++soln_iter){
      soln_stream_2 << "," <<  soln_iter->first << ":" << soln_iter->second;
    }
    soln_string = soln_stream_2.str();
    soln_string.replace(0, 1, "{");
    soln_string += "}";
    cout << "  Solution sizes: " << soln_string << "\n" << endl;
  } else {
    cout << "Warning: could not complete solution." << endl;
  }

  // free up everything.
  free(curr_soln_arr);
  free(result);
  free(tempResult);
  free(unigram);
  free(bigram);
  free(trigram);
  free(cipher_string);
  free(inv_soln_arr);

  cudaFree(greenhouse);
  cudaFree(backpointers);
  cudaFree(cudaUnigrams);
  cudaFree(cudaBigrams);
  cudaFree(cudaTrigrams);
  cudaFree(cuda_curr_soln_arr);
  cudaFree(cuda_inv_soln_arr);
  cudaFree(test_array);

  return 0;
}

/* -------------------------------- End Main------------------------------------------------ */
