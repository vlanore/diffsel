/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

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
