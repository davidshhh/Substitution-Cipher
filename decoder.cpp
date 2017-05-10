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

double * unigram;
double * bigram;
double * trigram;

// Pattern list to map from pattern string "abb" to words "all, odd, ..."
map<string, vector<string> > patlist;

// Word frequencies.
map<string, double> word_unigram;

// Maps characters to the integers representations

map<char, int> inv_cipher;
map<char, int> inv_plain;

// greenhouse and backpointers are the tables - greenhouse holds the probabilities, 
// and backpointers the back pointers.

double * greenhouse;
int * backpointers;

/* ------------------------------End Global Definitions ------------------------------------ */

/* ---------------------------------------- Classes ---------------------------------------- */

// This class is just a comparison function for the A* heap.  It looks as pairs of probailities
// and partial solutions, and orders them according to the probability.  On a tie, the larger 
// solution is given precedence.
// (Note: the probabilities are in the negative log domain, so will seem backwards.)

class pair_compare
{

public:
  bool operator() (const pair<double, pair<map<int, int>, vector<double> > > &lhs, const pair<double, pair<map<int, int>, vector<double> > > &rhs) const
  {
    return (lhs.first > rhs.first) || ((lhs.first == rhs.first) && (lhs.second.first.size() < rhs.second.first.size()));
  }
};

/* ------------------------------------ End Classes ---------------------------------------- */








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
// TODO: This documentation is from the CUDA version decoder.cu, we may want to re-read the code if we use this function.

void genviterbi(int len, int plain_num, int cipher_num, int * part_soln, int * cipher_string, long cipher_length, double * result, double * unigram, double * bigram, double * trigram, double * greenhouse, int * backpointers, double uselessvariableoriginallynamedthreshold){

  // loop variables
  //   -i: the index of the column of the greenhouse and backpointer tables that are being filled.
  int i;
  int l;
  int b;
  int k;

  // the value of the cipher at the psotion i, i - 1, and i - 2 in the code.
  //   -(c|c1|c2): The current / previous / second previous letter in the cipher at index i.
  int c;
  int c1;
  int c2;
  
  // The value of the ciphertext at the endpoint of the code.
  int c_1 = len;


  // This part will create a reverse map.
  int * part_inv = (int *) malloc(cipher_num*sizeof(int));
  for(l = 0; l < cipher_num; l++){
    *(part_inv + l) = -1;
  }
  for(l = 0; l < plain_num; l++){
    if(*(part_soln + l) >= 0){
      *(part_inv + *(part_soln + l)) = l;
    }
  }

  // Set up the zeroth column:
  c = *cipher_string;
  for(l = 0; l < plain_num; l++){
    for(b = 0; b < cipher_num + 1; b++){
      for(k = 0; k < plain_num; k++){
        if(((l == k) == (c == b)) 
              && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
              && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
              && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
              && ((*(part_inv + b) == -1) || (*(part_inv + b) == k)) 
              && (*(result + b * plain_num + k) >= 0) 
              && (*(result + c * plain_num + l) >= 0) ){
          *(greenhouse + ((l * cipher_num) + b)* plain_num + k) = *(unigram + c);
        } else {
          *(greenhouse + ((l * cipher_num) + b)* plain_num + k) = -1;
        }
      }
    }
  }

  // Set up the first column:
  c = *(cipher_string + 1);
  for(l = 0; l < plain_num; l++){
    for(b = 0; b < cipher_num + 1; b++){
      for(k = 0; k < plain_num; k++){
        if(((l == k) == (c == b)) 
              && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
              && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
              && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
              && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
              && (*(unigram + l) > 0)  
              && (*(result + c * plain_num + l) >= 0) ){
          int j;
          double best_prob = -1;
          int back = -1;
          for(j = 0; j < plain_num; j++){
            if((*(bigram + l + j * plain_num) > 0) && (*(greenhouse + ((j * cipher_num) + b)* plain_num + k) > 0)
                && ((best_prob < 0) 
                     || (best_prob > *(bigram + l + j * plain_num) + *(greenhouse + ((j * cipher_num) + b)* plain_num + k)))){
              best_prob = *(bigram + l + j * plain_num) + *(greenhouse + ((j * cipher_num) + b)* plain_num + k);
              back = j;
            }
          }
          *(greenhouse + ((plain_num + l) * cipher_num + b)* plain_num + k) = best_prob;
          *(backpointers + (((1) * plain_num + l) * cipher_num + b)* plain_num + k) = back;
        } else {
          *(greenhouse + ((plain_num + l) * cipher_num + b)* plain_num + k) = -1;
          *(backpointers + (((1) * plain_num + l) * cipher_num + b)* plain_num + k) = -1;
        }
      }
    }
  }
  
  // Set up the second column:
  c = *(cipher_string + 2);
  c1 = *(cipher_string + 1);
  c2 = *(cipher_string);
  b = c_1;
  for(l = 0; l < plain_num; l++){ // for: index l
     for(b = 0; b < cipher_num + 1; b++){ // for: index b
      for(k = 0; k < plain_num; k++){ // for: index k
        if(((l == k) == (c == b)) // if #1
             && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
             && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
             && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
             && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
             && (*(unigram + l) > 0)  
             && (*(result + c * plain_num + l) >= 0) ){
          int j;
          int j2;
          double best_prob = -1;
          int back = -1;
          if(*(part_inv + c1) >= 0){ // if #1.1
            j = *(part_inv + c1);
            if(*(greenhouse + (((1) * plain_num + j) * cipher_num + b)* plain_num + k) > 0){ // if #1.1.1
              if(*(part_inv + c2) >= 0){ // if # 1.1.1.1
                j2 = *(part_inv + c2);
                if((*(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                     && (*(bigram + j + plain_num * j2) >= 0)
                     && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)
                     && ((best_prob < 0) 
                            || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) 
                                + *(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k)))){ // if # 1.1.1.1.1
                  best_prob = *(greenhouse + (((1) * plain_num + j) * cipher_num + b)* plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                  back = j;
                } // if # 1.1.1.1.1
              } else { // if # 1.1.1.1
                for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                  if((*(greenhouse + (((1) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                      && ((j2 == l) == (c == c2))
                      && ((j == j2) == (c1 == c2))
                      && (*(bigram + j + plain_num * j2) >= 0)
                      && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.1.1.1
                      && ((best_prob < 0) 
                            || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k)))){
                      best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k);
                      back = j;
                  } // if #1.1.1.1.1 
                } // for: index j2 
              } // if #1.1.1.1
            } // if # 1.1.1
          } else { // if # 1.1
            for(j = 0; j < plain_num; j++){ // for: index j
              if((*(greenhouse + (((1) * plain_num + j) * cipher_num + b)* plain_num + k) >= 0)
                 && ((j == l) == (c == c1)) //(changed this since the July 10 run)
                 && ((best_prob < 0) || (*(greenhouse + (((1) * plain_num + j) * cipher_num + b)* plain_num + k) < best_prob))){ // if #1.1.2
                if(*(part_inv + c2) >= 0){ // if # 1.1.2.1
                  j2 = *(part_inv + c2);
                  if((*(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                     && (*(bigram + j + plain_num * j2) >= 0)
                     && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)){ // if # 1.1.2.1.1
                    best_prob = *(greenhouse + (((1) * plain_num + j) * cipher_num + b)* plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                    back = j;
                  } // if # 1.1.2.1.1
                } else { // if # 1.1.2.1
                  for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                    if((*(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                      && ((j2 == l) == (c == c2))
                      && ((j == j2) == (c1 == c2))
                      && (*(bigram + j + plain_num * j2) >= 0)
                      && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.2.1.1
                        && ((best_prob < 0) 
                            || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k)))){
                      best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(bigram + j + plain_num * j2) + *(greenhouse + (((0) * plain_num + j2) * cipher_num + b)* plain_num + k);
                      back = j;
                    } // if #1.1.2.1.1 
                  } // for: index j2 
                } // if #1.1.2.1
              } // if # 1.1.2
            } // for: index j
          } // if #1.1
          *(greenhouse + (((2) * plain_num + l) * cipher_num + b)* plain_num + k) = best_prob;
          *(backpointers + (((2) * plain_num + l) * cipher_num + b)* plain_num + k) = back;
        } else { // if # 1
          *(greenhouse + (((2) * plain_num + l) * cipher_num + b)* plain_num + k) = -1;
          *(backpointers + (((2) * plain_num + l) * cipher_num + b)* plain_num + k) = -1;
        } // if # 1
      } // for: index k
    } // for: index b
  } // for: index l



  // Then fill in the rest up to the end point:
  for(i = 3; i < cipher_length; i++){ // for: index i
    c = *(cipher_string + i);
    c1 = *(cipher_string + i-1);
    c2 = *(cipher_string + i-2);
    for(l = 0; l < plain_num; l++){ // for: index l
      for(b = 0; b < cipher_num + 1; b++){ // for: index b
        for(k = 0; k < plain_num; k++){ // for: index k
          if(((l == k) == (c == b)) // if #1
               && ((*(part_soln + l) == -1) || (*(part_soln + l) == c))
               && ((*(part_inv + c) == -1) || (*(part_inv + c) == l)) 
               && ((*(part_soln + k) == -1) || (*(part_soln + k) == b))
               && ((*(part_inv + b) == -1) || (*(part_inv + b) == k))
               && (*(unigram + l) > 0)  
               && (*(result + c * plain_num + l) >= 0) ){
            int j;
            int j2;
            
            double best_prob = -1;
            int back = -1;
            if(*(part_inv + c1) >= 0){ // if #1.1
              
              j = *(part_inv + c1);
              if(*(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) > 0){ // if #1.1.1
                if(*(part_inv + c2) >= 0){ // if # 1.1.1.1
                  j2 = *(part_inv + c2);
                  if((*(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                       && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)){ // if # 1.1.1.1.1
                    best_prob = *(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                    back = j;
                  } // if # 1.1.1.1.1
                } else { // if # 1.1.1.1
                  for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                    if((*(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                        && (*(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                        && ((j2 == l) == (c == c2))
                        && ((j == j2) == (c1 == c2))
                        && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.1.1.1
                        && ((best_prob < 0) 
                              || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k))) + *(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k)))){
                        best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k))) + *(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k);
                        back = j;
                    } // if #1.1.1.1.1 
                  } // for: index j2 
                } // if #1.1.1.1
              } // if # 1.1.1
            } else { // if # 1.1
              for(j = 0; j < plain_num; j++){ // for: index j
                if((*(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) >= 0)
                   && ((j == l) == (c == c1)) //(changed this since the July 10 run)
                   && ((best_prob < 0) || (*(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) < best_prob))){ // if #1.1.2
                  if(*(part_inv + c2) >= 0){ // if # 1.1.2.1
                    j2 = *(part_inv + c2);
                    if((*(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0) // if # 1.1.2.1.1
                       && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0)
                       && ((best_prob < 0) || (best_prob > *(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2))))){ 
                      best_prob = *(greenhouse + ((((i - 1) % 3) * plain_num + j) * cipher_num + b)* plain_num + k) + *(trigram + l + plain_num * (j + plain_num * j2));
                      back = j;
                    } // if # 1.1.2.1.1
                  } else { // if # 1.1.2.1
                    for(j2 = 0; j2 < plain_num; j2++){ // for: index j2
                      if((*(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                        && (*(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k) >= 0)
                        && ((j2 == l) == (c == c2))
                        && ((j == j2) == (c1 == c2))
                        && (*(trigram + l + plain_num * (j + plain_num * j2)) >= 0) // if # 1.1.2.1.1
                        && ((best_prob < 0) 
                              || (best_prob > *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k))) + *(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k)))){
                        best_prob = *(trigram + l + plain_num * (j + plain_num * j2)) + *(trigram + j + plain_num * (j2 + plain_num * *(backpointers + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k))) + *(greenhouse + ((((i - 2) % 3) * plain_num + j2) * cipher_num + b)* plain_num + k);
                        back = j;
                      } // if #1.1.2.1.1 
                    } // for: index j2 
                  } // if #1.1.2.1
                } // if # 1.1.2
              } // for: index j
            } // if #1.1
            *(greenhouse + (((i % 3) * plain_num + l) * cipher_num + b)* plain_num + k) = best_prob;
            *(backpointers + ((((i) % 3) * plain_num + l) * cipher_num + b)* plain_num + k) = back;
          } else { // if # 1
            *(greenhouse + (((i % 3) * plain_num + l) * cipher_num + b)* plain_num + k) = -1;
            *(backpointers + ((((i) % 3) * plain_num + l) * cipher_num + b)* plain_num + k) = -1;
          } // if # 1
        } // for: index k
      } // for: index b
    } // for: index l
  } // for: index i

  // Find and return the solutions:
  for(l = 0; l < plain_num * cipher_num; l++){
    *(result + l) = -1;
  }
  for(l = 0; l < plain_num; l++){
    for(b = 0; b < cipher_num; b++){
      int j;
      for(j = 0; j < plain_num; j++){
        if(*(greenhouse + ((((cipher_length - 1) % 3) * plain_num + l) * cipher_num + b)* plain_num + j) > 0){
          if(*(result + b * plain_num + j) <= 0){
            *(result + b * plain_num + j) = *(greenhouse + ((((cipher_length - 1) % 3) * plain_num + l) * cipher_num + b)* plain_num + j);
          } else {
            *(result + b * plain_num + j) = fmin(*(result + b * plain_num + j), *(greenhouse + ((((cipher_length - 1) % 3) * plain_num + l) * cipher_num + b)* plain_num + j));
          }
        }
      }
    }
  }

  free(part_inv);
}


// Returns a pattern string like "abc" if cipher_string.substring(i,j) is "the"
string GetPattern(int *cipher_string, int i, int j) {

    int count = 0;
    map<char, int> dict;
    string pattern = "";

    for( int ind = i; ind <= j; ind++ ) {
      int cipher_char = *(cipher_string + ind);
      if (dict.count(cipher_char) < 1)
        dict[cipher_char] = count++;
      pattern.push_back(char(dict[cipher_char] + 'a'));
    }
   return pattern;
}


void readchrmodel(
  double *unigram, double *bigram, double *trigram,
  string unigram_name, string bigram_name, string trigram_name,
  int p_2, int p_3
) {

  // loop variable.
  int i;

  // for the unigrams, read each line, process it into its parts,
  // and add each entry into a unigram map.
  // The basic map has size plain_num, with a default entry of -1.
  // unigram is an array whose ith index is *(unigram + i)

  for(i = 0; i < plain_num; i++){
    *(unigram + i) = -1;
  }
  ifstream unigrams_file(unigram_name.c_str());
  string temp;
  double uni_total = 0;
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
    double uni_prob = atof(str.c_str());
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
  for(i = 0; i < p_2; i++){
    *(bigram + i) = -1;
  }
  ifstream bigrams_file(bigram_name.c_str());
  double bi_total = 0;
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
    double bi_prob = atof(str.c_str());
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

  for(i = 0; i < p_3; i++){
    *(trigram + i) = -1;
  }
  ifstream trigrams_file(trigram_name.c_str());
  double tri_total = 0;
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
    double tri_prob = atof(str.c_str());
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
}


void ReadWordFrequencies(string word_unigram_name) {

  // loop variable.
  int i;
  map<string, double>::iterator word_iter;
  string temp;

  double total_prob = 0.0;
  cerr << "Reading word frequencies...";
  ifstream word_file(word_unigram_name.c_str());
  while (getline(word_file, temp, '\n')) {
    i = 0;

    // Read word
    string word = "";
    while((temp[i] == ' ') || (temp[i] == '\t')) { i++; }
    while((temp[i] != ' ') && (temp[i] != '\t')) {
      word.push_back(temp[i]);
      i++;
    }
    // Done reading word, now read probability
    string str = "";
    while((temp[i] == ' ') || (temp[i] == '\t')) { i++; }
    while((temp[i] != ' ') && (temp[i] != '\t')) {
      str.push_back(temp[i]);
      i++;
    }
    double prob = atof(str.c_str());
    word_unigram[word] = -1 * log(prob);
    total_prob += prob;
  }
  for(word_iter = word_unigram.begin(); word_iter != word_unigram.end(); ++word_iter) {
    word_iter->second += log(total_prob);
  }
  word_file.close();
  cerr << " Done." << endl;
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
  // line 3: the name of the word unigram file.
  // line 4: the name of the pattern list file.
  // line 5: the name of the unigram file.
  // line 6: the name of the bigram file.
  // line 7: the name of the trigram file.
  // line 8: the name of the ciphertext file.
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
// Both of these jobs are done by simple python scripts that were not given by the original author.

int main(int argc, char * argv[]){

  // A generic loop variable.

  int i;

  // Generic temp string.
  string temp;

  double uselessvariableoriginallynamedthreshold;

  // make sure that the calling procedure is correct.

  if (argc != 2){
    printf("\nUsage: ");
    printf("\t./decoder <profile>\n\n");
    exit(0);
  }

  /* read profile */
  ifstream f_profile(argv[1]);

  // Lines from the profile.
  // profiles have the format:
  // line 1: all ciphertext letters, seperated by spaces.
  // line 2: all plaintext letters, seperated by spaces.
  // line 3: the name of the word unigram file.
  // line 4: the name of the pattern list file.
  // line 5: the name of the unigram file.
  // line 6: the name of the bigram file.
  // line 7: the name of the trigram file.
  // line 8: the name of the ciphertext file.
  // Current profiles also have a seventh line that gives the
  // location of the plaintext file, which can be used for setup and evaluation
  // purposes.  It won't be read by this program, though.

  string cipher_line;
  string plain_line;
  string word_name;
  string pat_name;
  string unigram_name;
  string bigram_name;
  string trigram_name;
  string cipher_name;

  getline(f_profile, cipher_line, '\n');
  getline(f_profile, plain_line, '\n');
  getline(f_profile, word_name, '\n');
  getline(f_profile, pat_name, '\n');
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

  // Read the word frequencies (stored with negative log scale).
  ReadWordFrequencies(word_name);

  // Read the unigrams, bigrams, and trigrams

  unigram = (double *) malloc(plain_num * sizeof(double));
  int p_2 = plain_num * plain_num;
  bigram = (double *) malloc(p_2 * sizeof(double));
  int p_3 = plain_num * plain_num * plain_num;
  trigram = (double *) malloc(p_3 * sizeof(double));
  readchrmodel(unigram, bigram, trigram, unigram_name, bigram_name, trigram_name, p_2, p_3);

  // Create the cipher string.
  // Recall that we're going to turn everything into numbers, and first read everything into a vector,
  // then transfer the vector to an array.

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
  map<int, int> invlastcount;  // From cipher letter to the position of its order from last
  for(i = cipher_length - 1; i >= 0; i--){
    if(last.count(*(cipher_string + i)) == 0){
      last[*(cipher_string + i)] = i;
      lastcount[num_last] = i;
      invlastcount[*(cipher_string + i)] = num_last;
      num_last++;
    }
  }

  // Create the greenhouse array.
  greenhouse = (double *) malloc(plain_num * plain_num * cipher_num * 3 * sizeof(double));
  backpointers = (int *) malloc(plain_num * plain_num * cipher_num * 3 * sizeof(int));

  // set up priority heap.
  // set up a map of sizes.
  // Set up a starting solution and add it to the heap.
  double * result = (double *) malloc((plain_num * cipher_num) * sizeof(double));
  map<int, int> start_soln;

  // Normally, the solution is empty. here.
  // The following line adds the restriction that spaces map to spaces
  // to the solution (the ciphers we're using always end in a space).
  // In english, spaces are easy to find, so we can set them here.  We'll try to avoid it, though.
  // (to set the spaces, just uncomment here, and set all of the codes in the test to space -> space)
  start_soln[0] = 0;
  map<int, int> soln_sizes;
  soln_sizes[start_soln.size()] = 1;
  long pass_num = 0;  // TODO: the CUDA version initializes to 1 here
  int * letter_order = (int *) malloc(cipher_num * sizeof(int));

  // A boolean variable that will tell us when we're finished looking for the solution.

  bool found = false;

  // prime the priority queue with the first solution.

  int l, b, k;

  int * curr_soln_arr = (int *) malloc((plain_num + plain_num * cipher_num) * sizeof(int));
  time_t start_time = time (NULL);   
  priority_queue< pair<double, pair<map<int, int>, vector<double> > >, vector<pair<double, pair<map<int, int>, vector<double> > > >, pair_compare > aStar;
  vector<double> full_curr_soln = vector<double>(plain_num * cipher_num);
  for(i = 0; i < plain_num * cipher_num; i++){
    full_curr_soln[i] = 0;
  }
  aStar.push(pair<double, pair<map<int, int>, vector<double> > >(-1, pair<map<int, int>, vector<double> >(start_soln, full_curr_soln ) ) );

  // The main part of the program: pop solutions and run the genviterbi algorithm to grow solutions until
  // the final solution is found.

  while(!(aStar.empty()) && !found){

    // pop the best solution, and set up the next run.

    map<int, int> curr_soln = aStar.top().second.first;
    full_curr_soln = aStar.top().second.second;
    double curr_soln_prob = aStar.top().first;
    for(i = 0; i < plain_num; i++){  // Copy all curr soln to curr_soln_arr
      if(curr_soln.count(i) > 0){
        *(curr_soln_arr + i) = curr_soln[i];
      } else {
        *(curr_soln_arr + i) = -1;
      }
    }
    for(i = 0; i < plain_num * cipher_num; i++){
      *(curr_soln_arr + plain_num + i) = full_curr_soln[i];
    }
    aStar.pop();

    // print some stats.

    cout << "Starting Viterbi Algorithm: \n" << endl;
    cout << "  Queue size: " << (aStar.size() + 1) << "\n" << endl;
    cout << "  Solution size: " << curr_soln.size() << "\n" << endl;
    cout << "  Solution probability " << curr_soln_prob << "\n" << endl;
    cout << "  Number of passes: " << ++pass_num << "\n" << endl;
    
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
      // Uncomment the preferred order for adding solutions.
      // If using the most constrained first setting, uncomment the last first and loop lines.
      // ********** last first ********** //
      int curr_endpoint = *(cipher_string + lastcount[curr_soln.size()]);
      // ********** max freq ********** //
      //int curr_endpoint = *(cipher_uni_inv + curr_soln.size());
      // ********** min freq ********** //
      //int curr_endpoint = *(cipher_uni_inv + ((cipher_num-curr_soln.size()) % cipher_num));

      for(i = 0; i < plain_num * cipher_num; i++){    
          *(result + i) = full_curr_soln[i];
      }

      // use the gen_viterbi algorithm to grow larger solutions.
      genviterbi(curr_endpoint, plain_num, cipher_num, curr_soln_arr, cipher_string, cipher_length, result, unigram, bigram, trigram, greenhouse, backpointers, uselessvariableoriginallynamedthreshold);

      if(pass_num == 1){
        // order the solutions.

        // ********** Last First ********** //
          //for(l = 0; l < cipher_num; l++){
          //  *(letter_order + l) = *(cipher_string + lastcount[l]);
          //}

        // ********** Most Constrained First ********** //
          vector<int> numconstrained = vector<int>(cipher_num, -1);
          int most_constrained_index;
          int most_constrained_count;

          // Set up the vector numconstrained to tell how many solutions each cipher letter has.
          for(i = 0; i < cipher_num; i++){
            most_constrained_count = 0;
            for(l = 0; l < plain_num; l++){
              if(*(result + i * plain_num + l) >= 0){
                most_constrained_count++;
              }
            }
            numconstrained[i] = most_constrained_count;
          }

          // Go through the numconstrained vector and put the most constrained into the letter order array.
          for(i = 0; i < cipher_num; i++){
            most_constrained_index = 0;
            for (l = 0; l < cipher_num; l++){
              if((numconstrained[most_constrained_index] < 0)
                  || ( (numconstrained[l] > 0)
                       &&(numconstrained[most_constrained_index] > numconstrained[l])) 
                  || ( (numconstrained[l] > 0)
                       &&(numconstrained[most_constrained_index] == numconstrained[l])
                       && (invlastcount[l] < invlastcount[most_constrained_index])) ){
                most_constrained_index = l;
              }
            }
            *(letter_order + i) = most_constrained_index;
            numconstrained[most_constrained_index] = -1;
          }
      }
      curr_endpoint = *(letter_order + curr_soln.size());

      // About here, put in a filtering algorithm for the partial solutions.
        // (1) a. For each remaining cipher letter, find the possible plaintext letters.
        //        Count the number r of these letters
        //     b. Go over all of the remaining cipher letters, and see if their possiblilities
        //        are contained in this list.  It so, increment a counter cr by 1 (start it at 0).
        //     c. If cr > r, set solution to 0 probability.
        //        If cr = r, set any cells poining to the possible plain text letters other than the 
        //        ones in the cr count to 0, and proceed to the next letter.
        //     d. Repeat until there are no changed in one round of ciphertext letters.
        //     Complexity:  if there are K rounds, each round can take at most c^2*p steps, for Kpc^2 steps.
        //                  But every round, the program will continue unless at least one cell is set to 0.
        //                  So K <= pc + 1, and the complexity is at most (cp + 1)pc^2 in O(p^2c^3) steps.
      vector<int> constraints = vector<int>(plain_num, 0);
      int finished = 0;
      int diff;
      int i7;
      int constraint_count;
      int match_count;
      while(!finished){
        finished = 1;
        for(i = 0; i < cipher_num; i++){
          // for each cipher letter
          // find the constraints.
          constraint_count = 0;
          for(l = 0; l < plain_num; l++){
            if(*(result + i * plain_num + l) >= 0){
              constraints[l] = 1;
              constraint_count++;
            } else {
              constraints[l] = 0;
            }
          }

          // find other letters whose constraints match (or, that don't match)
          match_count = cipher_num - constraint_count;
          for(i7 = 0; i7 < cipher_num; i7++){
            diff = 0;
            for(l = 0; l < plain_num; l++){
              if((constraints[l] == 0) && (*(result + i7 * plain_num + l) >= 0)){
                diff = 1;
              }
            }
            if(diff == 1){
              match_count--;
            }
          }

          // If the number of matches is high enough, change the rest of the table.
          if(match_count >= 0){
            for(i7 = 0; i7 < cipher_num; i7++){
              diff = 0;
              for(l = 0; l < plain_num; l++){
                if((constraints[l] == 0) && (*(result + i7 * plain_num + l) >= 0)){
                  diff = 1;
                }
              }
              if(diff == 1){
                for(l = 0; l < plain_num; l++){
                  if((constraints[l] > 0) && (*(result + i7 * plain_num + l) >= 0)){
                    *(result + i7 * plain_num + l) = -1;
                    finished = 0;
                  }
                }
              }
            }
          } // if(match_count >= 0)

          // done for each cipher letter
        } // for(i = 0; i < cipher_num; i++)
      } // while(!finished)

      // For each possible extension to the solution, format the next partial solution
      // and add it to the heap.

      for(l = 0; l < plain_num; l++){
        if(*(result + curr_endpoint * plain_num + l) >= 0){
          map<int, int> next_soln;
          map<int, int>::iterator curr_iter;
          for(curr_iter = curr_soln.begin(); curr_iter != curr_soln.end(); ++curr_iter){
            next_soln[curr_iter->first] = curr_iter->second;
          }
          next_soln[l] = curr_endpoint;
          double next_prob = *(result + curr_endpoint * plain_num + l);
          if(next_prob > 0){
            vector<double> temp_full = vector<double>(plain_num * cipher_num);
            for(i = 0; i < plain_num * cipher_num; i++){
              temp_full[i] = *(result + i);
            }
            aStar.push(pair<double, pair<map<int, int>, vector<double> > >(next_prob, pair<map<int, int>, vector<double> >(next_soln, temp_full ) ) );
            if(soln_sizes.count(curr_soln.size() + 1) > 0){
              soln_sizes[curr_soln.size() + 1] += 1;
            } else {
              soln_sizes[curr_soln.size() + 1] = 1;
            };
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

  free(curr_soln_arr);
  free(result);
  free(greenhouse);
  free(unigram);
  free(bigram);
  free(trigram);
  free(cipher_string);
  return 0;
}

/* -------------------------------- End Main------------------------------------------------ */
