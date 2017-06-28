#include "Tree.hpp"
#include <fstream>
#include <list>
#include <sstream>
#include "TaxonSet.hpp"
using namespace std;

void Tree::ToStream(ostream &os) const {
    ToStream(os, GetRoot());
    os << ";\n";
}

void Tree::ToStream(ostream &os, const Link *from) const {
    if (!from->isLeaf()) {
        os << '(';
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            ToStream(os, link->Out());
            if (link->Next() != from) {
                os << ',';
            }
        }
        os << ')';
    }
    os << GetNodeName(from);
    if (!from->isRoot()) {
        string brval = GetBranchName(from);
        if (brval != "") {
            os << ':' << brval;
        }
    }
}

void Tree::DeleteNextLeaf(Link *previous) {
    Link *link = previous->Next();
    if (!link->Out()->isLeaf()) {
        cout << "Bad call of DeleteNextLeaf, it must be call on a link pointing on "
                "a leaf\n";
        exit(1);
    }
    previous->Next()->Next()->AppendTo(previous);
    delete link->Out()->GetNode();
    delete link->Out();
    delete link->GetBranch();
    delete link;
}

void Tree::DeleteUnaryNode(Link *from) {
    if (!from->isUnary()) {
        cout << "Bad call of DeleteUnaryNode, node is not unary\n";
        exit(1);
    }
    if (from->isRoot()) {
        Link *newroot = from->Next()->Out();
        newroot->SetBranch(from->GetBranch());
        delete from->Next()->GetBranch();
        delete from->GetNode();
        delete from->Next();
        delete from;
        root = newroot;
        root->InsertOut(root);
    } else {
        ostringstream sum;
        sum << atof((from->GetBranch()->GetName()).c_str()) +
                   atof((from->Next()->GetBranch()->GetName()).c_str());
        from->GetBranch()->SetName(sum.str());
        from->Next()->Out()->SetBranch(from->GetBranch());
        from->Out()->InsertOut(from->Next()->Out());
        delete from->Next()->GetBranch();
        delete from->GetNode();
        delete from->Next();
        delete from;
    }
}

void Tree::RootAt(Link *from) {
    Link *link = from->Next();
    root->AppendTo(from);
    link->AppendTo(root);
    root->SetNode(from->GetNode());
}

void Tree::EraseInternalNodeName() { EraseInternalNodeName(GetRoot()); }

void Tree::EraseInternalNodeName(Link *from) {
    if (!from->isLeaf()) {
        from->GetNode()->SetName("");
    }
    for (Link *link = from->Next(); link != from; link = link->Next()) {
        EraseInternalNodeName(link->Out());
    }
}

void Tree::RegisterWith(const TaxonSet *taxset) {
    int tot = 0;
    if (!RegisterWith(taxset, GetRoot(), tot)) {
        cout << "There is no match between the tree and the sequences.\n";
        exit(1);
    }
    if (tot != taxset->GetNtaxa()) {
        cerr << "error : non matching number of taxa : " << tot << '\t' << taxset->GetNtaxa()
             << '\n';
        cerr << "some taxa in the dataset are not present in the tree\n";
        exit(1);
    }
}

bool Tree::RegisterWith(const TaxonSet *taxset, Link *from, int &tot) {
    if (from->isLeaf()) {
        int i = taxset->GetTaxonIndex(from->GetNode()->GetName());
        if (i != -1) {
            from->GetNode()->SetIndex(i);
            tot++;
        }
        return (i != -1);
    }
    Link *previous = from;
    while (previous->Next() != from) {
        if (RegisterWith(taxset, previous->Next()->Out(), tot)) {
            previous = previous->Next();
        } else {
            // cout << "delete !!\n";
            DeleteNextLeaf(previous);
        }
    }
    if (from == nullptr) {
        cerr << "from is null\n";
        exit(1);
    }
    if (from->isUnary()) {
        // cerr << "DELETE UNARY NODE\n";
        DeleteUnaryNode(from);
    }
    return (!from->isLeaf());
}

Tree::Tree(string filename) {
    ifstream is(filename.c_str());
    if (!is) {
        cerr << "cannot find file : " << filename << '\n';
        exit(1);
    }
    ReadFromStream(is);
}

void Tree::ReadFromStream(istream &is) {
    string expr = "";
    int cont = 1;
    while (cont != 0) {
        string s;
        is >> s;
        unsigned int k = 0;
        while ((k < s.length()) && (s[k] != ';')) {
            k++;
        }
        expr += s.substr(0, k);
        cont = static_cast<int>((!is.eof()) && (k == s.length()));
    }
    SetRoot(ParseGroup(expr, nullptr));
}

Link *Tree::ParseList(string input, Node *node) {
    try {
        // parse input as a list of strings separated by ','
        list<string> lst;
        int n = input.size();
        int k = 0;
        int brack = 0;
        int b = 0;
        while (k < n) {
            char c = input[k];
            if (c == '(') {
                brack++;
            }
            if (c == ')') {
                brack--;
            }
            if ((brack == 0) && (c == ',')) {
                lst.push_back((string)(input.substr(b, k - b)));
                b = k + 1;
            }

            if (brack < 0) {
                cout << "in parse list : too many )\n";
                cout << input.substr(0, k) << '\n';
                cout << input << '\n';
                throw;
            }
            k++;
        }
        if (brack != 0) {
            cout << "in parse list : too many (\n";
            cout << input << '\n';
            throw;
        }
        lst.push_back(input.substr(b, k - b));

        // make a circular single link chain around the node
        // with one link for each term of the list
        // and call parse group on each term
        auto firstlink = new Link;
        Link *prevlink = firstlink;
        firstlink->SetNode(node);
        for (auto &i : lst) {
            auto link = new Link;
            link->SetNode(node);
            link->AppendTo(prevlink);
            ParseGroup(i, link);
            prevlink = link;
        }
        firstlink->AppendTo(prevlink);
        return firstlink;
    } catch (...) {
        cout << "exit in parse list\n";
        exit(1);
    }
}

Link *Tree::ParseGroup(string input, Link *from) {
    try {
        // parse input as (body)nodeval:branchval
        string body = "";
        unsigned int k = 0;
        if (input[0] == '(') {
            int brack = 1;
            k = 1;
            while ((k < input.length()) && (brack != 0)) {
                char c = input[k];
                if (c == '(') {
                    brack++;
                }
                if (c == ')') {
                    brack--;
                }
                k++;
            }
            if (brack != 0) {
                cout << "in parse group: too many (\n";
                cout << input << '\n';
                throw;
            }
            body = input.substr(1, k - 2);
        }

        int b = k;
        while ((k < input.length()) && (input[k] != ':')) {
            k++;
        }
        string nodeval = input.substr(b, k - b);

        string branchval = "";
        if (k < input.length()) {
            branchval = input.substr(k + 1, input.length() - k);
        }

        // make a new node and a new branch
        Node *node = new Node(nodeval);

        // call parse body
        Link *link = nullptr;
        if (body != "") {
            link = ParseList(body, node);
        } else {
            link = new Link;
            link->SetNode(node);
        }
        if (from != nullptr) {
            Branch *branch = new Branch(branchval);
            link->SetBranch(branch);
            from->SetBranch(branch);
            link->InsertOut(from);
        }
        return link;
    } catch (...) {
        cout << "exit in parse group\n";
        exit(1);
    }
}

TEST_CASE("test") { printf("Hello world!\n"); }
