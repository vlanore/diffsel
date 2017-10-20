# thanks to http://www.petercollingridge.co.uk/book/export/html/474
bases = ['a', 'c', 'g', 't']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

reverse_bases = dict(zip(bases, range(4)))
short_aa_table = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# short_aa_table = "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
aas = "".join(sorted("FSYCWLPHQRIMTNKVADEG"))
reverse_aas = dict(zip(aas, range(20)))

def cno_to_nuc(cno):
    return codons[cno]

# def cno_to_nuc(cno):
#     if (cno <10):
#         return codons[cno]
#     elif (cno == 10 or cno == 11):
#         return codons[cno+1]
#     else:
#         return codons[cno+2]

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
<<<<<<< b10459504dc9d172d322b53d09b6fcef0db5da23
        a = 1
        #print("\tQcodons[cod][%d][%d] := 0" % (c1, c2))
    elif synonym(c1, c2):
        print("\tQcodons[cod][%d][%d] := abs(Q[%d][%d])" % (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1))
    else:
        print("\tQcodons[cod][%d][%d] := abs(Q[%d][%d] * sqrt(fitness[cod][%d]/fitness[cod][%d]))" %
              (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1, reverse_aas[short_aa_table[c2]]+1, reverse_aas[short_aa_table[c1]]+1))

def print_all():
    print("for (cod in 1:nsites) {")
    print("  for (i in 1:64) {")
    print("    for (j in 1:64) {")
    print("      Qcodons[cod][i][j] := abs(0)")
    print("    }")
    print("  }")
    print("}")
    print("for (cod in 1:nsites) {")
    for cno in range(61):
        for cno2 in range(61):
=======
        print("\tQcodons[cod][%d][%d] := 0" % (c1, c2))
    elif synonym(c1, c2):
        print("\tQcodons[cod][%d][%d] := Q[%d][%d]" % (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1))
    elif (amino_acids[c1] == '*' or amino_acids[c2] == '*'):
        print("\tQcodons[cod][%d][%d] := 0" % (c1, c2))
    else:
        print("\tQcodons[cod][%d][%d] := Q[%d][%d] * sqrt(fitness[cod][%d]/fitness[cod][%d])" %
              (c1+1, c2+1, reverse_bases[d[0][0]]+1, reverse_bases[d[0][1]]+1, reverse_aas[amino_acids[c2]]+1, reverse_aas[amino_acids[c1]]+1))

def print_all():
    print("for (cod in 1:nsites) {")
    for cno in range(64):
        for cno2 in range(64):
>>>>>>> script: rev gen now uses base 64
            case(cno, cno2)
    print("\tR[cod] := fnFreeK(Qcodons[cod])")
    print("}")

print_all()
