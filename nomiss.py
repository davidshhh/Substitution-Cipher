import string

f = open('lmtrain_nyt_word.chr.clm3.noEOS')
alphabet = list(string.ascii_lowercase)
alphabet = ['_'] + alphabet
i = 0
j = 0
k = 0
SIZE = 27
for line in f:
    parts = line.split()
    pi = alphabet.index(parts[0])
    pj = alphabet.index(parts[1])
    pk = alphabet.index(parts[2])
    while not (i == pi and j == pj and k == pk):
        print '%s\t%s\t%s\t0.000000001' % (alphabet[i], alphabet[j], alphabet[k])
        if k == SIZE - 1:
            k = 0
            if j == SIZE - 1:
                j = 0
                if i == SIZE - 1:
                    # Done
                    pass
                else:
                    i += 1
            else:
                j += 1
        else:
            k += 1
    print '%s\t%s\t%s\t%d' % (parts[0], parts[1], parts[2], float(parts[3]))

    if k == SIZE - 1:
        k = 0
        if j == SIZE - 1:
            j = 0
            if i == SIZE - 1:
                # Done
                pass
            else:
                i += 1
        else:
            j += 1
    else:
        k += 1
