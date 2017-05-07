import argparse
from random import randint
from pprint import pprint

''' Creates a random mapping (key) of letters to letters.
'''


def get_key():
    key = {}
    mapping = [chr(i) for i in range(ord('a'), ord('z') + 1)]
    for i in range(ord('a'), ord('z') + 1):
        r = randint(0, len(mapping) - 1)
        key[chr(i)] = mapping.pop(r)
    # ' ' is assumed to be known
    key[' '] = '_'
    return key


''' Encodes plaintext file with a random mappign into cipher text.
'''


def cipher(plaintext_filename, outfile, length):
    key = get_key()
    f = open(plaintext_filename, 'r')
    out = open(outfile, 'w')
    count = 0
    for l in f:
        for c in l:
            out.write(key[c])
            count += 1
            if count >= length and c == ' ':
                pprint(key)
                out.close()
                print count
                return
    pprint(key)
    out.close()


parser = argparse.ArgumentParser()
parser.add_argument("plaintext_filename")
parser.add_argument("output_filename")
parser.add_argument("--length", type=int, default=6000)

args = parser.parse_args()

pf = args.plaintext_filename
o = args.output_filename
l = args.length

cipher(pf, o, l)
