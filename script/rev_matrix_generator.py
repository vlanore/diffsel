# thanks to http://www.petercollingridge.co.uk/book/export/html/474
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

short_aa_table = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

def cno_to_nuc(cno):
    if (cno <10):
        return codons[cno]
    elif (cno == 10 or cno == 11):
        return codons[cno+1]
    else:
        return codons[cno+2]

def synonym(c1, c2):
    return codon_table[cno_to_nuc(c1)] == codon_table[cno_to_nuc(c2)]

def diff(cno1, cno2):
    c1 = cno_to_nuc(cno1)
    c2 = cno_to_nuc(cno2)
    r = []
    for i in range(len(c1)):
        if c1[i] != c2[i]:
            r.append((c1[i], c2[i]))
    return r


def case(c1, c2):
    d = diff(c1, c2)
    if (len(d) != 1):
        print("R[%d][%d] = 0" % (c1, c2))
    elif synonym(c1, c2):
        print("R[%d][%d] = Q[%s][%s]" % (c1, c2, d[0][0], d[0][1]))
    else:
        print("R[%d][%d] = Q[%s][%s] * sqrt(F[%s]/F[%s])" % (c1, c2, d[0][0], d[0][1], short_aa_table[c2], short_aa_table[c1]))


for cno in range(61):
    for cno2 in range(61):
        case(cno, cno2)
