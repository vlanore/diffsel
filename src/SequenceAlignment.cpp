#include "SequenceAlignment.hpp"
#include <fstream>
#include <sstream>
#include "BiologicalSequences.hpp"
#include "Random.hpp"
using namespace std;

const char digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

inline int IsInt(std::string s) {
    int returnValue = 1;
    unsigned int i = 0;
    if ((s[0] == '+') || (s[0] == '-')) {
        i++;
    }
    if (i == s.length()) {
        returnValue = 0;
    }

    while ((returnValue != 0) && (i < s.length())) {
        int j = 0;
        while ((j < 10) && (digit[j] != s[i])) {
            j++;
        }
        if (j == 10) {
            returnValue = 0;
        }
        i++;
    }
    return returnValue;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     SequenceAlignment
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
int Int(string s) { return atoi(s.c_str()); }

void SequenceAlignment::ToStream(ostream &os) const {
    os << Ntaxa << '\t' << Nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) {
            max = l;
        }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) {
            os << ' ';
        }
        for (int j = 0; j < Nsite; j++) {
            os << statespace->GetState(GetState(i, j));
        }
        os << '\n';
    }
    os << '\n';
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     FileSequenceAlignment
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
FileSequenceAlignment::FileSequenceAlignment(string filename) {
    SpeciesNames.clear();
    ReadDataFromFile(filename, 0);
    taxset = new TaxonSet(SpeciesNames, Ntaxa);
}

int FileSequenceAlignment::ReadDataFromFile(string filespec, int forceinterleaved) {
    string tmp;
    ifstream is((Path + filespec).c_str());
    if (!is) {
        cerr << "error : cannot find data file " << filespec << '\n';
        cerr << "\n";
        exit(1);
    }
    is >> tmp;
    try {
        cerr << "-- [SequenceAlignment] Alignment file uses Phylip format" << endl;
        if (forceinterleaved == 0) {
            int returnvalue = TestPhylipSequential(filespec);
            if (returnvalue != 0) {
                cerr << "-- [SequenceAlignment] Alignment file is sequential" << endl;
                ReadPhylipSequential(filespec);
                return 1;
            }
        }
        cerr << "-- [SequenceAlignment] Alignment file is interleaved";
        int returnvalue = TestPhylip(filespec, 1);
        if (returnvalue != 0) {
            cerr << ", taxon names repeated" << endl;
            ReadPhylip(filespec, 1);
            return 1;
        }
        cerr << ", taxon names not repeated" << endl;
        TestPhylip(filespec, 0);
        ReadPhylip(filespec, 0);
        return 1;
    } catch (...) {
        exit(1);
    }
    return 1;
}

void FileSequenceAlignment::FailFormat() {
    cerr << "error when reading data\n"
         << "data should be formatted as follows:\n"
         << "#taxa #sites\n"
         << "name1 seq1.....\n"
         << "name2 seq2.....\n"
         << "...\n\n";
    exit(1);
}

// ---------------------------------------------------------------------------
//     ReadPhylip()
// ---------------------------------------------------------------------------
int FileSequenceAlignment::TestPhylipSequential(string filespec) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Ntaxa = atoi(temp.c_str());
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Nsite = atoi(temp.c_str());

        SpeciesNames = std::vector<std::string>(Ntaxa, "");

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            do {
                char c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        while (c != ')') {
                            theStream >> c;
                        }
                    } else if (c == '{') {
                        while (c != '}') {
                            theStream >> c;
                        }
                    } else {
                        int p = 0;
                        if (DNAcomp != 0) {
                            while ((p < DNAN) && (c != DNAset[p])) {
                                p++;
                            }
                            if (p == DNAN) {
                                DNAcomp = 0;
                            }
                        }
                        p = 0;
                        if (RNAcomp != 0) {
                            while ((p < RNAN) && (c != RNAset[p])) {
                                p++;
                            }
                            if (p == RNAN) {
                                RNAcomp = 0;
                            }
                        }
                        p = 0;
                        if (AAcomp != 0) {
                            while ((p < AAN) && (c != AAset[p])) {
                                p++;
                            }
                            if (p == AAN) {
                                AAcomp = 0;
                            }
                        }
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            if (theStream.eof()) {
                if (nsite < Nsite) {
                    return 0;
                }
            }
            ntaxa++;
        }
        if (theStream.eof()) {
            if (ntaxa < Ntaxa) {
                return 0;
            }
        }
        if (DNAcomp != 0) {
            statespace = new DNAStateSpace;
        } else if (RNAcomp != 0) {
            statespace = new RNAStateSpace;
        } else if (AAcomp != 0) {
            statespace = new ProteinStateSpace;
        } else {
            return 0;
        }
    } catch (...) {
        return 0;
    }
    return 1;
}

void FileSequenceAlignment::ReadPhylipSequential(string filespec) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Nsite = Int(temp);

        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }
        SpeciesNames = std::vector<std::string>(Ntaxa, "");

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            do {
                char c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != ')') {
                            theStream >> c;
                        }
                    } else if (c == '{') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != '}') {
                            theStream >> c;
                        }
                    } else {
                        ostringstream s;
                        s << c;
                        Data[ntaxa][nsite] = statespace->GetState(s.str());
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        exit(1);
    }
}

int FileSequenceAlignment::TestPhylip(string filespec, int repeattaxa) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Nsite = Int(temp);

        SpeciesNames = std::vector<std::string>(Ntaxa, "");

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((l == 0) || (repeattaxa != 0)) {
                    string temp;
                    theStream >> temp;
                    if (l == 0) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) {
                            return 0;
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            while (c != ')') {
                                theStream >> c;
                            }
                        } else if (c == '{') {
                            while (c != '}') {
                                theStream >> c;
                            }
                        } else {
                            int p = 0;
                            if (DNAcomp != 0) {
                                while ((p < DNAN) && (c != DNAset[p])) {
                                    p++;
                                }
                                if (p == DNAN) {
                                    DNAcomp = 0;
                                }
                            }
                            p = 0;
                            if (RNAcomp != 0) {
                                while ((p < RNAN) && (c != RNAset[p])) {
                                    p++;
                                }
                                if (p == RNAN) {
                                    RNAcomp = 0;
                                }
                            }
                            p = 0;
                            if (AAcomp != 0) {
                                while ((p < AAN) && (c != AAset[p])) {
                                    p++;
                                }
                                if (p == AAN) {
                                    AAcomp = 0;
                                }
                            }
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13) && (c != 10));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    theStream.get();
                    c = theStream.peek();
                }
                if (m == 0) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "in test phylip\n";
                        cerr << "error when reading data non matching number of sequences "
                                "in block "
                                "number "
                             << block << " for taxon " << i + 1 << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        if (DNAcomp != 0) {
            statespace = new DNAStateSpace;
            // cerr << "dna sequences\n";
        } else if (RNAcomp != 0) {
            statespace = new DNAStateSpace;
            // cerr << "rna sequences\n";
        } else if (AAcomp != 0) {
            statespace = new ProteinStateSpace;
            // cerr << "protein sequences\n";
        } else {
            // cerr << "format not recognised\n";
            return 0;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return 1;
}

void FileSequenceAlignment::ReadPhylip(string filespec, int repeattaxa) {
    ifstream theStream((Path + filespec).c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            FailFormat();
        }
        Nsite = Int(temp);
        // cerr << Ntaxa << '\t' << Nsite << '\n';

        Data = new int *[Ntaxa];
        for (int i = 0; i < Ntaxa; i++) {
            Data[i] = new int[Nsite];
        }
        SpeciesNames = std::vector<std::string>(Ntaxa, "");

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((l == 0) || (repeattaxa != 0)) {
                    string temp;
                    theStream >> temp;
                    if (l == 0) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) {
                            cerr << "error when reading data: read " << temp << " instead of "
                                 << SpeciesNames[i] << '\n';
                            exit(1);
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            Data[i][k] = unknown;
                            while (c != ')') {
                                theStream >> c;
                            }
                        } else if (c == '{') {
                            Data[i][k] = unknown;
                            while (c != '}') {
                                theStream >> c;
                            }
                        } else {
                            ostringstream s;
                            s << c;
                            Data[i][k] = statespace->GetState(s.str());
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    theStream.get();
                    c = theStream.peek();
                }

                if (m == 0) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "error when reading data non matching number of sequences "
                                "in block "
                                "number "
                             << block << " for taxon " << i << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
    } catch (...) {
        cerr << "error while reading data file\n";
    }
}
