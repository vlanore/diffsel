#include "CodonStateSpace.hpp"
#include <cstdlib>
#include <iostream>
#include <sstream>
using namespace std;

CodonStateSpace::CodonStateSpace(GeneticCodeType type) {
    nucstatespace = new DNAStateSpace;
    protstatespace = new ProteinStateSpace;

    code = type;
    if (code == Universal) {
        Nstate = Ncodon - UniNStopCodons;

        CodonCodeWithStops = new int[Ncodon];  // stops included
        CodonCode = new int[Nstate];           // stops excluded

        CodonPos = new int *[Npos];
        for (int pos = 0; pos < Npos; pos++) {
            CodonPos[pos] = new int[Nstate];  // stops excluded
        }

        int k = 0;
        for (int i = 0; i < Ncodon; i++) {
            CodonCodeWithStops[i] = UniCodonCode[i];
            if (CodonCodeWithStops[i] != -1) {
                CodonCode[k] = CodonCodeWithStops[i];
                for (int pos = 0; pos < Npos; pos++) {
                    CodonPos[pos][k] = codonpos[pos][i];
                }
                k++;
            }
        }

        Nstop = UniNStopCodons;
        StopCodons = new int[Nstop];
        StopPos1 = new int[Nstop];
        StopPos2 = new int[Nstop];
        StopPos3 = new int[Nstop];
        for (int i = 0; i < Nstop; i++) {
            StopCodons[i] = UniStopCodons[i];
            StopPos1[i] = UniStopPos1[i];
            StopPos2[i] = UniStopPos2[i];
            StopPos3[i] = UniStopPos3[i];
        }

    } else if (code == MtInv) {
        Nstate = Ncodon - MtInvNStopCodons;

        CodonCodeWithStops = new int[Ncodon];  // stops included
        CodonCode = new int[Nstate];           // stops excluded

        CodonPos = new int *[Npos];
        for (int pos = 0; pos < Npos; pos++) {
            CodonPos[pos] = new int[Nstate];  // stops excluded
        }

        int k = 0;
        for (int i = 0; i < Ncodon; i++) {
            CodonCodeWithStops[i] = MtInvCodonCode[i];
            if (CodonCodeWithStops[i] != -1) {
                CodonCode[k] = CodonCodeWithStops[i];
                for (int pos = 0; pos < Npos; pos++) {
                    CodonPos[pos][k] = codonpos[pos][i];
                }
                k++;
            }
        }

        Nstop = MtInvNStopCodons;
        StopCodons = new int[Nstop];
        StopPos1 = new int[Nstop];
        StopPos2 = new int[Nstop];
        StopPos3 = new int[Nstop];
        for (int i = 0; i < Nstop; i++) {
            StopCodons[i] = MtInvStopCodons[i];
            StopPos1[i] = MtInvStopPos1[i];
            StopPos2[i] = MtInvStopPos2[i];
            StopPos3[i] = MtInvStopPos3[i];
        }
    } else if (code == MtMam) {
        Nstate = Ncodon - MtMamNStopCodons;

        CodonCodeWithStops = new int[Ncodon];  // stops included
        CodonCode = new int[Nstate];           // stops excluded

        CodonPos = new int *[Npos];
        for (int pos = 0; pos < Npos; pos++) {
            CodonPos[pos] = new int[Nstate];  // stops excluded
        }

        int k = 0;
        for (int i = 0; i < Ncodon; i++) {
            CodonCodeWithStops[i] = MtMamCodonCode[i];
            if (CodonCodeWithStops[i] != -1) {
                CodonCode[k] = CodonCodeWithStops[i];
                for (int pos = 0; pos < Npos; pos++) {
                    CodonPos[pos][k] = codonpos[pos][i];
                }
                k++;
            }
        }

        Nstop = MtMamNStopCodons;
        StopCodons = new int[Nstop];
        StopPos1 = new int[Nstop];
        StopPos2 = new int[Nstop];
        StopPos3 = new int[Nstop];
        for (int i = 0; i < Nstop; i++) {
            StopCodons[i] = MtMamStopCodons[i];
            StopPos1[i] = MtMamStopPos1[i];
            StopPos2[i] = MtMamStopPos2[i];
            StopPos3[i] = MtMamStopPos3[i];
        }
    } else {
        cerr << "genetic code not recognised\n";
        cerr << type << '\n';
        exit(1);
    }
}

CodonStateSpace::~CodonStateSpace() {
    delete[] CodonCode;
    delete[] CodonCodeWithStops;
    for (int pos = 0; pos < Npos; pos++) {
        delete[] CodonPos[pos];
    }
    delete[] CodonPos;

    delete nucstatespace;
    delete protstatespace;
}

string CodonStateSpace::GetState(int codon) const {
    ostringstream s;
    if (codon == -1) {
        s << "---";
    } else {
        s << DNAletters[GetCodonPosition(0, codon)] << DNAletters[GetCodonPosition(1, codon)]
          << DNAletters[GetCodonPosition(2, codon)];
    }
    if (s.str().length() != 3) {
        cerr << "error in translation\n";
        exit(1);
    }
    return s.str();
}

int CodonStateSpace::GetState(string word) const {
    return GetCodonFromDNA(GetDNAStateSpace()->GetState(word.substr(0, 1)),
                           GetDNAStateSpace()->GetState(word.substr(1, 1)),
                           GetDNAStateSpace()->GetState(word.substr(2, 1)));
}

bool CodonStateSpace::CheckStop(int pos1, int pos2, int pos3) const {
    if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown)) {
        return false;
    }
    int l = 0;
    while ((l < Nstop) &&
           ((pos1 != StopPos1[l]) || (pos2 != StopPos2[l]) || (pos3 != StopPos3[l]))) {
        l++;
    }
    return (l < Nstop);
}

int CodonStateSpace::GetCodonFromDNA(int pos1, int pos2, int pos3) const {
    if ((pos1 == unknown) || (pos2 == unknown) || (pos3 == unknown)) {
        return unknown;
    }
    int l = 0;
    while ((l < GetNstate()) &&
           ((pos1 != GetCodonPosition(0, l)) || (pos2 != GetCodonPosition(1, l)) ||
            (pos3 != GetCodonPosition(2, l)))) {
        l++;
    }
    if (l == GetNstate()) {
        cerr << "warning in CodonStateSpace::GetCodonFromDNA : out of bound : "
             << GetDNAStateSpace()->GetState(pos1) << GetDNAStateSpace()->GetState(pos2)
             << GetDNAStateSpace()->GetState(pos3) << '\n';
        if (code == Universal) {
            cerr << "universal\n";
        } else if (code == MtMam) {
            cerr << "mt mam\n";
        } else if (code == MtInv) {
            cerr << "mt inv\n";
        }
        throw 1;
    }
    return l;
}

int CodonStateSpace::GetDifferingPosition(int i, int j) const {
    // identical
    if ((GetCodonPosition(0, i) == GetCodonPosition(0, j)) &&
        (GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
        (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
        return -1;
    }
    if (GetCodonPosition(0, i) != GetCodonPosition(0, j)) {
        if ((GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
            (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
            return 0;
        }
        return 3;
    }
    if (GetCodonPosition(1, i) != GetCodonPosition(1, j)) {
        if ((GetCodonPosition(0, i) == GetCodonPosition(0, j)) &&
            (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
            return 1;
        }
        return 3;
    }
    if (GetCodonPosition(2, i) != GetCodonPosition(2, j)) {
        if ((GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
            (GetCodonPosition(0, i) == GetCodonPosition(0, j))) {
            return 2;
        }
        return 3;
    }
    return 3;
}
