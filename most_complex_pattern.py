import argparse
from pprint import pprint

def count(pattern_file):
    pf = open(pattern_file, 'r')

    outfile = open("unigram_lm", 'w')
    check = set()
    total_count = 0.0
    total_lines = 0.0
    max_count = 0.0

    for l in pf:
        items = l.split()

        if items[0] == '2':
            outfile.close()
            break

        outfile.write(l)
        pattern = items[1]
        num_pattern = len(items) - 2
        max_count = len(pattern) if len(pattern) > max_count else max_count
        total_count += num_pattern
        total_lines += 1
        for c in pattern:
            if c not in check:
                check.add(c)
                print pattern
    
    print 'average words per pattern is {}'.format(total_count / total_lines)
    print 'max number of characters per pattern is {}'.format(max_count)
    pprint(check)
    print '{} characters were used in pattern'.format(len(check))


parser = argparse.ArgumentParser()
parser.add_argument("pattern_file", help='file containing the pattern list')


args = parser.parse_args()

pattern_file = args.pattern_file
count(pattern_file)