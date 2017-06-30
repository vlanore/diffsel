#include "CodonSequenceAlignment.hpp"
#include <cstdlib>
#include <iostream>
#include "Random.hpp"
using namespace std;

CodonSequenceAlignment::CodonSequenceAlignment(SequenceAlignment *from, bool force_stops,
                                               GeneticCodeType type) {
    try {
        DNAsource = from;

        if (from->Nsite % 3 != 0) {
            cerr << "not multiple of three\n";
            exit(1);
        }
        Nsite = from->Nsite / 3;
        Ntaxa = from->Ntaxa;
        auto tempstatespace = new CodonStateSpace(type);
        statespace = tempstatespace;
        // DNAStateSpace* nucspace = tempstatespace->GetDNAStateSpace();

        taxset = DNAsource->GetTaxonSet();

        // make my own arrays
        // make translation
        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
            for (int j = 0; j < Nsite; j++) {
                try {
                    Data[i][j] = GetCodonStateSpace()->GetCodonFromDNA(
                        DNAsource->GetState(i, 3 * j), DNAsource->GetState(i, 3 * j + 1),
                        DNAsource->GetState(i, 3 * j + 2));
                    if (Data[i][j] == -1) {
                        if ((DNAsource->GetState(i, 3 * j) != -1) &&
                            (DNAsource->GetState(i, 3 * j + 1) != -1) &&
                            (DNAsource->GetState(i, 3 * j + 2) != -1)) {
                            // cerr << "in CodonSequenceAlignment: taxon " <<
                            // taxset->GetTaxon(i) <<
                            // " and codon " << j+1 << " (site " << 3*j+1 << ") :";
                            // cerr << nucspace->GetState(DNAsource->GetState(i, 3*j)) <<
                            // nucspace->GetState(DNAsource->GetState(i, 3*j+1)) <<
                            // nucspace->GetState(DNAsource->GetState(i, 3*j+2)) << '\n';
                        }
                    }
                } catch (...) {
                    // catch(Exception e)	{
                    // cerr << "in CodonSequenceAlignment: taxon " << i << " and codon "
                    // << j << "
                    // (site " << 3*j << ")\n";
                    // cerr << "taxon : " << taxset->GetTaxon(i) << '\n';
                    if (force_stops) {
                        // Data[i][j] = -2;
                        Data[i][j] = -1;
                    } else {
                        throw;
                    }
                }
            }
        }

    } catch (...) {
        cerr << "Codon Sequence Alignment: failed to read the datafile\n";
        exit(1);
    }
}
