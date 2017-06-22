import sys
import random

if len(sys.argv)<3 or int(sys.argv[1])<=0:
    print "usage:\n\n\tpython generate_toy.py <n: stricly positive integer> <fn: filename>\n\nThis will generate a tree with 2^n leaves and a corresponding alignment. Alignment and tree will be written to fn.ali and fn.tree respectively"
else:
    def tree(k, i, mystr):
        if k == 0:
            return mystr+"(S"+str(i)+":0,S"+str(i+1)+":1):0", i+2
        else:
            a, b = tree(k-1, i, mystr)
            c, d = tree(k-1, b, "("+a+",")
            return c+"):0", d

    bases = ["A", "C", "G", "T"]

    def rand_codon():
        return random.choice(bases)+random.choice(bases)+random.choice(bases)

    def selected_codon():
        aa1 = ["AAT", "AAC"]
        return random.choice(aa1)

    def mutate(codon):
        result = list(codon)
        result[random.randint(0,2)] = random.choice(bases)
        return "".join(result)

    def align(length):
        mystr = ""
        mystr += str(2**n)+" "+str(length)+'\n'
        lastcodon = ""
        for i in range(2**(n-1)):
            codon = selected_codon()
            codon2 = mutate(codon)
            # codon2 = rand_codon()
            mystr+="S"+str(2*i)+"\t"+codon2+'\n'
            mystr+="S"+str(2*i+1)+"\t"+codon+'\n'
        return mystr

    n = int(sys.argv[1])-1
    fn = sys.argv[2]
    treefile = open(fn+".tree", "w+")
    alifile = open(fn+".ali", "w+")

    a,b = tree(n, 0, "")
    treefile.write(a)
    alifile.write(align(3))
