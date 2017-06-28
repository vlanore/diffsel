#ifndef STRINGSTREAMUTILS_H
#define STRINGSTREAMUTILS_H

#include <sstream>

const char digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

inline int GoPastNext(std::istream &is, const char inChar) {
    if (!is.eof()) {
        unsigned char c;
        do {
            is >> c;
        } while (c != inChar && !is.eof());
    }
    return static_cast<int>(!is.eof());
}

inline void GoPastNextWord(std::istream &is, const std::string &inWord) {
    unsigned int k = 0;
    char c;
    while ((!is.eof()) && (k < inWord.length())) {
        is.get(c);
        if ((c >= 65) && (c <= 90)) {
            c += 32;
        }
        char ca = inWord[k];
        if ((ca >= 65) && (ca <= 90)) {
            ca += 32;
        }
        if (c == ca) {
            k++;
        } else {
            k = 0;
        }
    }
}

inline int EquivalentStrings(std::string a, std::string b) {
    if (a.length() != b.length()) {
        return 0;
    }
    unsigned int k = 0;
    int cont = 1;
    while ((k < a.length()) && ((cont) != 0)) {
        char ca = a[k];
        char cb = b[k];
        if ((ca >= 65) && (ca <= 90)) {
            ca += 32;
        }
        if ((cb >= 65) && (cb <= 90)) {
            cb += 32;
        }
        if (ca != cb) {
            cont = 0;
        }
        k++;
    }
    return cont;
}

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

#endif
