#!/usr/bin/python3.5

from diffsel_script_utils import *

#===================================================================================================
print(step("Parsing command line arguments"))

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='Removes indels from a sequence file.')
parser.add_argument('inputFile', metavar="input", type=FileType('r'), nargs=1, help='the sequence file (phylip format)')

args = parser.parse_args()

seq_file = args.inputFile[0]
print("-- Sequence file is "+param(seq_file.name))
out_file = seq_file.name+".woindels"
print("-- Output file is "+param(out_file))

#===================================================================================================
print(step("Extracting and verifying data from file"))

nb_taxa, length = [int(x) for x in seq_file.readline().split()]
print("-- Number of taxa: "+data(nb_taxa))
print("-- Sequence length before indel removal: "+data(length))

sequences = [l.split() for l in seq_file]
if (nb_taxa == len(sequences)):
    print(good("Number of sequences matches number of taxa"))
else:
    print(failure("number of sequences is "+str(len(sequences))+" instead of "+str(nb_taxa)))
    exit(1)
for i in range(nb_taxa):
    if (len(sequences[i][1]) != length):
        print(failure("sequence "+str(i)+" ("+sequences[i][0]+") has length "+str(len(sequences[i][1]))+" instead of "+str(length)))
        exit(1)
print(good("Lengths of sequences all match expected sequence length"))

#===================================================================================================
print(step("Removing indels!"))
indel_pos = []
for s in sequences:
    for pos in range(length):
        if s[1][pos] == '-' or s[1][pos] == '?':
            if not (pos in indel_pos):
                indel_pos.append(pos)
print("-- Found "+data(len(indel_pos))+" positions with indels")
print("-- Building resulting alignment... ", end="")
non_indel_pos = set(range(length)).difference(indel_pos)
result = [[s[0], [s[1][i] for i in non_indel_pos]] for s in sequences]
print("Done")
print("-- Resulting alignment length is "+data(len(result[1][1]))+" ("+data(round(100.*len(indel_pos)/length))+"% of sites were removed)")

#===================================================================================================
print(step("Writing result to file"))
out = open(out_file, 'w')
out.write(str(nb_taxa)+'\t'+str(len(non_indel_pos))+'\n')
for s in result:
    out.write(s[0]+'\t')
    out.write(str(''.join(s[1])))
    out.write('\n')
out.close()
