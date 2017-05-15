"""
This script takes in a file and outputs the alphabet it uses.
"""
import sys

f = open(sys.argv[1], 'r')
alphabet = []
for line in f:
    for c in line:
        if c not in alphabet:
            alphabet.append(c)
            print c,

