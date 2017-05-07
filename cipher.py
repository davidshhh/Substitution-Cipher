import argparse
from random import randint
from pprint import pprint

''' Creates a random mapping (key) of letters to letters.
'''
def getKey():
	key = {}
	mapping = [chr(i) for i in range(ord('a'), ord('z')+1)]
	for i in range(ord('a'), ord('z')+1):
		r = randint(0, len(mapping)-1)
		key[chr(i)] = mapping.pop(r)
	# ' ' is assumed to be known
	key[' '] = ' '
	return key

''' Encodes plaintext file with a random mappign into cipher text.
'''
def cipher(plaintext_filename):
	key = getKey();
	f = open(plaintext_filename, 'r')
	s = ""
	for l in f:
		for c in l:
			s += key[c]
	pprint(key)
	print s

parser = argparse.ArgumentParser()
parser.add_argument("plaintext_filename");

args = parser.parse_args()

plaintext_filename = args.plaintext_filename

cipher(plaintext_filename)









