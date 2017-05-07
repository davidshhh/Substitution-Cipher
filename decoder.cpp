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
//#include "genviterbi.h"

using namespace std;

// Global definitions

int plain_num;
int cipher_num;
long cipher_length;

double * greenhouse;
int * backpointers;

vector<char> cipher_alpha(0);
vector<char> plain_alpha(0);
vector<int> cipher_string_vec(0);
int * cipher_string;

double * unigram;
double * bigram;
double * trigram;

map<char, int> inv_cipher;
map<char, int> inv_plain;

// classes
class pair_compare
{

public:
  bool operator() (const pair<double, pair<map<int, int>, vector<double> > > &lhs, const pair<double, pair<map<int, int>, vector<double> > > &rhs) const
  {
    return (lhs.first > rhs.first) || ((lhs.first == rhs.first) && (lhs.second.first.size() < rhs.second.first.size()));
  }
};


class partsoln
{
};













void genviterbi(int len, int plain_num, int cipher_num, int * part_soln, int * cipher_string, long cipher_length, double * result, double * unigram, double * bigram, double * trigram, double * greenhouse, int * backpointers, double threshold){

  // loop variables
  int i;
  int l;
  int b;
  int k;

  // the value of the cipher at the psotion i, i - 1, and i - 2 in the code.
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











int main(int argc, char * argv[]){

  int i;
  double threshold;

  int * LM_uni_sorted;
  int * cipher_uni_sorted;

  // make sure that the calling procedure is correct.
  if (argc != 2){
    printf("\nUsage: ");
    printf("\t<insert usage description here>\n\n");
    exit(0);
  }

  /* read profile */
  ifstream f_profile(argv[1]);

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

  /*Use the cipher and plaintext lines to create the corresponding lists:*/
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
  
  // Create the greenhouse array.
  greenhouse = (double *) malloc(plain_num * plain_num * cipher_num * 3 * sizeof(double));
  backpointers = (int *) malloc(plain_num * plain_num * cipher_num * 3 * sizeof(int));

  /*Read the unigrams, bigrams, and trigrams*/
  // for the unigrams, read each line, process it into its parts,
  // and add each entry into a unigram map.
  // unigram is array such that index for i = base(bigram) + i
  unigram = (double *) malloc(plain_num * sizeof(double));
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
     if(*(unigram + i) != -1){
       *(unigram + i) += log(uni_total);
     } else {
       *(unigram + i) = -1;
     }
  }
  unigrams_file.close();

  // for the bigrams, read each line, process it into its parts,
  // and add each entry into a unigram map of maps, adding intermediate
  // entries where needed.
  // bigram is array such that index for ij = base(bigram) + plain_num * i + j
  int p_2 = plain_num * plain_num;
  bigram = (double *) malloc(p_2 * sizeof(double));
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
     if(*(bigram + i) != -1){
       *(bigram + i) += log(bi_total);
     } else {
       *(bigram + i) = -1;
     }
  }
  bigrams_file.close();
  
  threshold = 0;  
  // for the trigrams, read each line, process it into its parts,
  // and add each entry into a unigram map of maps, adding intermediate
  // entries where needed.
  // trigram is array such that index for ijk = base(bigram) + plain_num^2 * i + j * plain_num + k
  int p_3 = plain_num * plain_num * plain_num;
  trigram = (double *) malloc(p_3 * sizeof(double));
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
     if(*(trigram + i) != -1){
       *(trigram + i) += log(tri_total);
       threshold = fmax(threshold, *(trigram + i));
     } else {
       *(trigram + i) = -1;
     }
  }
  trigrams_file.close();

  /* Create the cipher string*/
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
  int num_last = 0;
  map<int, int> last;
  map<int, int> lastcount;
  map<int, int> invlastcount;
  for(i = cipher_length - 1; i >= 0; i--){
    if(last.count(*(cipher_string + i)) == 0){
      last[*(cipher_string + i)] = i;
      lastcount[num_last] = i;
      invlastcount[*(cipher_string + i)] = num_last;
      num_last++;
    }
  }
  
// This part sets the array LM_uni_sorted to give the order of the most frequent letters in the LM (i.e., the most frequent at LM_uni_sorted[0]).
// It'll be used to change the order in which letters are processed in the viterbi algorithm: algorithmically, it's the same, since everything is parallel, 
// but computationally it'll do just a little bit better timewise on CUDE, because of the way the memory calls work: it keeps recent calls in cache, 
// and this increases the likelihood a bit of repeating a call before it goes out of cache (aat least, it gives better performance, and I think that's what's happening).
      // Set up the unigram mapping for the LM:
        double * temp_uni = (double *) malloc(plain_num * sizeof(double));
        for(i = 0; i < plain_num; i++){
          *(temp_uni + i) = *(unigram + i);
        }
        LM_uni_sorted = (int *) malloc(plain_num * sizeof(int));
        int best_loc = 0;
        for(i = 0; i < plain_num; i++){
          int best = 0;
          int i2;
          for(i2 = 0; i2 < plain_num; i2++){
            if(*(temp_uni + best) < *(temp_uni + i2)){
              best = i2;
            }
          }
          *(LM_uni_sorted + i) = best;
          *(temp_uni + best) = -2;
        }
      // Set up the unigram model for the cipher:
        for(i = 0; i < plain_num; i++){
          if(i < cipher_num){
            *(temp_uni + i) = 0;
          } else {
            *(temp_uni + i) = -2;
          }
        }
        for(i = 0; i < cipher_length; i++){
          *(temp_uni + *(cipher_string + i)) += 1;
        }
        for(i = 0; i < cipher_num; i++){
          *(temp_uni + i) = *(temp_uni + i) / cipher_length;
        }
        cipher_uni_sorted = (int *) malloc(plain_num * sizeof(int));
        for(i = 0; i < cipher_num; i++){
          int best = 0;
          int i2;
          for(i2 = 0; i2 < cipher_num; i2++){
            if(*(temp_uni + best) < *(temp_uni + i2)){
              best = i2;
            }
          }
          *(cipher_uni_sorted + i) = best;
          *(temp_uni + best) = -2;
        }
      // Compare:
        int *cipher_uni_inv = (int *) malloc(cipher_num * sizeof(int));
        for(i = 0; i < cipher_num; i++){
          *(cipher_uni_inv + *(cipher_uni_sorted + i)) = i;
        }
        int c_1 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string + 2)));
        int c_2 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string + 1)));
        int c_3 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string)));
        double total = 0;
        if((*(unigram + c_3) > 0) && (*(bigram + c_2 + plain_num * c_3) > 0) 
            && (*(trigram + c_1 + plain_num *( c_2 + plain_num * c_3)) > 0)  ){
          total = *(unigram + c_3) + *(bigram + c_2 + plain_num * c_3) + *(trigram + c_1 + plain_num *( c_2 + plain_num * c_3));
          for(i = 3; i < cipher_length; i++){
            c_1 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string + i)));
            c_2 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string + i - 1)));
            c_3 = *(LM_uni_sorted + *(cipher_uni_inv + *(cipher_string + i - 2)));
            if(*(trigram + c_1 + plain_num *( c_2 + plain_num * c_3)) > 0){
              total += *(trigram + c_1 + plain_num *( c_2 + plain_num * c_3));
            } else {
              total += threshold*cipher_length;
            }
          }
          threshold = total;
        } else {
          threshold = threshold * cipher_length;
        }

  // set up priority heap
  // set up a map of sizes
  double * result = (double *) malloc((plain_num * cipher_num) * sizeof(double));
  map<int, int> start_soln;
  // In english, spaces are easy to find, so we can set them here.  We'll try to avoid it, though.
  // (to set the spaces, just uncomment here, and set all of the codes in the test to space -> space)
  start_soln[0] = 0;

  map<int, int> soln_sizes;
  soln_sizes[start_soln.size()] = 1;
  long pass_num = 0;
  int * letter_order = (int *) malloc(cipher_num * sizeof(int));

  bool found = false;

  int l, b, k;

  int * curr_soln_arr = (int *) malloc((plain_num + plain_num * cipher_num) * sizeof(int));
  time_t start_time = time (NULL);   
  priority_queue< pair<double, pair<map<int, int>, vector<double> > >, vector<pair<double, pair<map<int, int>, vector<double> > > >, pair_compare > aStar;
  vector<double> full_curr_soln = vector<double>(plain_num * cipher_num);
  for(i = 0; i < plain_num * cipher_num; i++){
    full_curr_soln[i] = 0;
  }
  aStar.push(pair<double, pair<map<int, int>, vector<double> > >(-1, pair<map<int, int>, vector<double> >(start_soln, full_curr_soln ) ) );
  // while the priority heap is not empty and while the solution size
  // is less than the plaintext size:
  while(!(aStar.empty()) && !found){
    map<int, int> curr_soln = aStar.top().second.first;
    full_curr_soln = aStar.top().second.second;
    double curr_soln_prob = aStar.top().first;
    // change here
    int curr_soln_size = curr_soln.size();
    for(i = 0; i < plain_num; i++){
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

    cout << "Starting Viterbi Algorithm: \n" << endl;
    cout << "  Queue size: " << (aStar.size() + 1) << "\n" << endl;
    cout << "  Solution size: " << curr_soln_size << "\n" << endl;
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

    if(curr_soln.size() >= cipher_num){
      found = true;
      start_soln = curr_soln;
    } else {
      // Uncomment the preferred order for adding solutions.  If using the most constrained first setting, uncomment the last first and loop lines.
      // ********** last first ********** //
      int curr_endpoint = *(cipher_string + lastcount[curr_soln_size]);
      // ********** max freq ********** //
      //int curr_endpoint = *(cipher_uni_inv + curr_soln.size());  
      // ********** min freq ********** //
      //int curr_endpoint = *(cipher_uni_inv + ((cipher_num-curr_soln.size()) % cipher_num)); 

      for(i = 0; i < plain_num * cipher_num; i++){    
          *(result + i) = full_curr_soln[i];
      }

      // use the gen_viterbi algorithm to grow larger solutions.
      genviterbi(curr_endpoint, plain_num, cipher_num, curr_soln_arr, cipher_string, cipher_length, result, unigram, bigram, trigram, greenhouse, backpointers, threshold);

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
      curr_endpoint = *(letter_order + curr_soln_size);

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

          // for each cipher letter, find the constraints.
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
          }
        }
      }

      // Create and apply the next solutions.
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
    // put the valid ones on the heap.

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
