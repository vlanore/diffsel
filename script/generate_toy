# Copyright or Copr. Centre National de la Recherche Scientifique (CNRS) (2017/11/27)
# Contributors:
# - Vincent Lanore <vincent.lanore@gmail.com>

# This software is a computer program whose purpose is to provide small tools and scripts related to phylogeny and bayesian
# inference.

# This software is governed by the CeCILL-B license under French law and abiding by the rules of distribution of free software.
# You can use, modify and/ or redistribute the software under the terms of the CeCILL-B license as circulated by CEA, CNRS and
# INRIA at the following URL "http://www.cecill.info".

# As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users
# are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive
# licensors have only limited liability.

# In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or
# reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated
# to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth
# computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements
# in conditions enabling the security of their systems and/or data to be ensured and, more generally, to use and operate it in
# the same conditions as regards security.

# The fact that you are presently reading this means that you have had knowledge of the CeCILL-B license and that you accept
# its terms.

from diffsel_script_utils import *

if len(sys.argv)<3 or int(sys.argv[1])<=0:
    print("usage:\n\n\generate_toy <n: stricly positive integer> <fn: filename>\n\nThis will generate a tree with 2^n leaves and a corresponding alignment. Alignment and tree will be written to fn.ali and fn.tree respectively")
else:
    def tree(k, i, mystr):
        if k == 0:
            return mystr+"(S"+str(i)+":0,S"+str(i+1)+":1):0", i+2
        else:
            a, b = tree(k-1, i, mystr)
            c, d = tree(k-1, b, "("+a+",")
            return c+"):0", d

    def align(length):
        mystr = ""
        mystr += str(2**n)+" "+str(length)+'\n'
        lastcodon = ""
        for i in range(2**(n-1)):
            codon = mutate(selected_codon(), 20)
            codon2 = mutate(codon, 80)
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
