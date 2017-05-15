# Substitution-Cipher
Jiawei Huang jhujhuang@gmail.com
David Hsieh  thsieh10@gmail.com

### Credits

This is our final project for Declarative Methods at the Johns Hopkins University.

Code and the original paper were from Eric Corlett at University of Toronto.

Word frequency data originally from https://github.com/bmhauer/hhk-decipherment who derived these from the New York Times corpus.
We modified the data for compatibility.

## Files

#### Main program

`decoder.cpp` is our main solver program.

To run the solver, first compile it
```
g++ -o solve decoder.cpp
```
Then run it with a profile file specified
```
./solve profile
```
where `profile` is an input file with the format as described below.

#### Input files

`profile` has the format the solver `decoder.cpp` is looking for.

  profiles have the format:
-  line 1: all ciphertext letters, seperated by spaces.
-  line 2: all plaintext letters, seperated by spaces.
-  line 3: the name of the word unigram file.
-  line 4: the name of the pattern list file.
-  line 5: the name of the character unigram file.
-  line 6: the name of the character bigram file.
-  line 7: the name of the character trigram file.
-  line 8: the name of the ciphertext file.

`lmtrain_nyt_word.unk.wlm1` is the word unigram file we used.

`lmtrain_nyt_word.patlist.wlm1` is the word pattern list file we used.

`lmtrain_nyt_word.chr.clm1.noEOS` is the character unigram file we used.

`lmtrain_nyt_word.chr.clm2.noEOS` is the character bigram file we used.

`lmtrain_nyt_word.chr.clm3.noEOS.nomiss` is the character trigram file we used.

Our example ciphertext files include:

`smallcipher.txt`

`multiwordcipher.txt`

`3000cipher.txt`

The plaintext files are not actual inputs but were used to generate the ciphertext:

`A.txt`
TODO: cite source of A.txt

`smallplain.txt`

`multiwordplain.txt`

Some accompanying `.key` files generated with ciphertexts

`multiword.key`

`3000.key`


#### Scripts

`cipher.py` A script to generate ciphertext of certain length from a plain text.

Run it like:
```
python cipher.py A.txt 50cipher.txt --length 50
```
where `A.txt` is the plain text file with a space at the beginning.
It will write the ciphertext to `50cipher.txt` and the random key mapping to standard output.

`decode.py`

`most_complex_pattern.py`

`nomiss.py`

`alphabet.py`
Because the orignal code didn't support alphabet input that contains letters unused in the cipher text,
we wrote this script to extract the actual alphabet in a cipher text to feed it into their system.
We do not need to do this for our input since we made a simple fix for this problem.

#### Original code from Eric Corlett and compatible input files

The author Eric Corlett of the paper _An Exact A* Method for Solving Letter Substitution Ciphers_
provided these two programs on his website:

`decoder_ut.cpp`

`decoder.cu`

We cannot find accompanying data files and scripts so we created the following profile file as an example that is compatible with his programs:

`500profileut`
