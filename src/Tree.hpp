#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <map>
#include <string>
#include "Random.hpp"
#include "StringStreamUtils.hpp"
using namespace std;

class TaxonSet;  // forward decl

class Node {
  private:
    int index;
    std::string name;

  public:
    Node() : index(0), name("") {}
    Node(std::string s) : index(0), name(std::move(s)) {}
    Node(const Node *from) : index(from->index), name(from->name) {}

    virtual ~Node() = default;

    virtual std::string GetName() const { return name; }
    virtual void SetName(std::string inname) { name = inname; }
    int GetIndex() const { return index; }
    void SetIndex(int i) { index = i; }
};

class Branch {
  private:
    int index;
    std::string name;

  public:
    Branch() : index(0), name("") {}
    Branch(std::string s) : index(0), name(std::move(s)) {}
    Branch(const Branch *from) : index(from->index), name(from->name) {}

    virtual ~Branch() = default;

    virtual std::string GetName() const { return name; }
    virtual void SetName(std::string inname) { name = inname; }
    int GetIndex() const { return index; }
    void SetIndex(int i) { index = i; }
};

class Link {
  private:
    Link *next;
    Link *out;
    Branch *branch;
    Node *node;
    int index;

  public:
    Link() {
        next = out = this;
        branch = nullptr;
        node = nullptr;
    }

    Link(const Link * /*unused*/) {
        next = out = this;
        node = nullptr;
        branch = nullptr;
    }

    Link *Next() const { return next; }
    Link *Out() const { return out; }
    Branch *GetBranch() const { return branch; }
    Node *GetNode() const { return node; }

    void SetBranch(Branch *inbranch) { branch = inbranch; }
    void SetNode(Node *innode) { node = innode; }

    void SetIndex(int i) { index = i; }

    int GetIndex() const { return index; }

    void SetOut(Link *inout) { out = inout; }

    void SetNext(Link *innext) { next = innext; }

    void AppendTo(Link *link) {
        if (link != nullptr) {
            link->next = this;
        }
    }

    void Insert(Link *link) {  // insert link after this
        link->next = next;
        next = link;
    }

    void InsertOut(Link *link) {  // insert link as out
        link->out = this;
        out = link;
    }

    bool isLeaf() const { return (next == this); }

    bool isUnary() const { return (next->Next() == this && !isLeaf()); }

    bool isRoot() const { return (out == this); }

    // degree : number of branches connecting to the node associated to this link
    int GetDegree() const {
        int d = 1;
        const Link *link = next;
        while (link != this) {
            d++;
            link = link->next;
        }
        return d;
    }

    const Link *GetUp(int &d) const {
        std::cerr << "in getup\n";
        exit(1);
        d = 1;
        const Link *link = Out();
        // const Link* link = link->Out();
        while (link->GetDegree() == 2) {
            link = link->Next()->Out();
            d++;
        }
        return link;
    }
};

class NewickTree {
  public:
    virtual ~NewickTree() = default;
    virtual Link *GetRoot() const = 0;

    void ToStream(std::ostream &os) const;
    void ToStream(std::ostream &os, const Link *from) const;
    double ToStreamSimplified(std::ostream &os, const Link *from) const;

    const Link *GetLeftMostLink(const Link *from) const {
        if (from->isLeaf()) {
            return from;
        }
        return GetLeftMostLink(from->Next()->Out());
    }

    const Link *GetRightMostLink(const Link *from) const {
        if (from->isLeaf()) {
            return from;
        }
        const Link *link = from->Next();
        while (link->Next() != from) {
            link = link->Next();
        }
        return GetRightMostLink(link->Out());
    }

    std::string GetLeftMost(const Link *from) const {
        if (from->isLeaf()) {
            return GetNodeName(from);
        }
        return GetLeftMost(from->Next()->Out());
    }

    std::string GetRightMost(const Link *from) const {
        if (from->isLeaf()) {
            return GetNodeName(from);
        }
        const Link *link = from->Next();
        while (link->Next() != from) {
            link = link->Next();
        }
        return GetRightMost(link->Out());
    }

    static void Simplify() { simplify = true; }

    void PrintTab(std::ostream &os) const { RecursivePrintTab(os, GetRoot()); }

    void RecursivePrintTab(std::ostream &os, const Link *from) const {
        os << GetLeftMost(from) << '\t' << GetRightMost(from) << '\t' << GetNodeName(from) << '\n';
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            RecursivePrintTab(os, link->Out());
        }
    }

    int GetNnode() const { return RecursiveGetNnode(GetRoot()); }

    int GetNinternalNode() const { return RecursiveGetNinternalNode(GetRoot()); }

    int RecursiveGetNinternalNode(const Link *from) const {
        int n = 0;
        if (!from->isLeaf()) {
            n++;
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            n += RecursiveGetNinternalNode(link->Out());
        }
        return n;
    }

    int RecursiveGetNnode(const Link *from) const {
        int n = 1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            n += RecursiveGetNnode(link->Out());
        }
        return n;
    }

  protected:
    virtual std::string GetNodeName(const Link *link) const = 0;
    virtual std::string GetBranchName(const Link *link) const = 0;

    virtual std::string GetLeafNodeName(const Link *link) const { return GetNodeName(link); }

    static bool simplify;
};

class Tree : public NewickTree {
  public:
    Tree();
    // default constructor: set member pointers to 0

    Tree(const Tree *from);
    // clones the entire Link structure
    // but does NOT clone the Nodes and Branches
    // calls RecursiveClone

    Tree(std::string filename);
    // create a tree by reading into a file (netwick format expected)
    // calls ReadFromStream

    ~Tree() override;
    // calls RecursiveDelete
    // but does NOT delete the Nodes and Branches

    // Delete the leaf pointing by the next link and set everithing right.
    void DeleteNextLeaf(Link *previous);

    // Delete the unary Node wich from is paart of and set everithing right.
    void DeleteUnaryNode(Link *from);

    Link *GetRoot() const override { return root; }
    const TaxonSet *GetTaxonSet() const { return taxset; }

    void RootAt(Link *from);

    void RegisterWith(const TaxonSet *taxset);
    // Registers all leaves of the tree with an external TaxonSet
    // the taxon set defines a map between taxon names and indices (between 0 and
    // P-1)
    // the tree is recursively traversed
    // each leaf's name is looked for in the map of the taxon set
    // if not found : an error is emitted
    // otherwise, the leaf's index is set equal to the index of the corresponding
    // taxon

    bool RegisterWith(const TaxonSet *taxset, Link *from, int &tot);
    // recursive function called by RegisterWith

    double GetBranchLength(const Link *link) const { return atof(GetBranchName(link).c_str()); }

    double GetMaxHeight(const Link *from) const {
        double max = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = GetMaxHeight(link->Out());
            if (max < tmp) {
                max = tmp;
            }
        }
        if (!from->isRoot()) {
            max += GetBranchLength(from);
        }
        return max;
    }

    double GetMinHeight(const Link *from) const {
        double min = -1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double tmp = GetMinHeight(link->Out());
            if ((min == -1) || (min > tmp)) {
                min = tmp;
            }
        }
        if (!from->isRoot()) {
            min += GetBranchLength(from);
        }
        return min;
    }

    void ToStreamRenorm(const Link *from, std::ostream &os, double normfactor) const {
        if (from->isLeaf()) {
            os << GetNodeName(from);
        } else {
            os << "(";
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                ToStreamRenorm(link->Out(), os, normfactor);
                if (link->Next() != from) {
                    os << ",";
                }
            }
            os << ")";
            os << GetNodeName(from);
        }
        if (from->isRoot()) {
            os << ";\n";
        } else {
            os << ":";
            os << GetBranchLength(from) * normfactor;
        }
    }

    std::string GetBranchName(const Link *link) const override {
        return link->GetBranch()->GetName();
    }

    std::string GetNodeName(const Link *link) const override {
        return link->GetNode()->GetName();
        /*
          if (! link->isLeaf())	{
          return link->GetNode()->GetName();
          }
          std::string s = link->GetNode()->GetName();
          unsigned int l = s.length();
          unsigned int i = 0;
          while ((i < l) && (s[i] != '_')) i++;
          if (i == l)	{
                  std::cerr << "error in get name\n";
          exit(1);
          }
          i++;
          return s.substr(i,l-i);
        */
    }
    // trivial accessors
    // they can be useful to override, so as to bypass Branch::GetName() and
    // Node::GetName()

    void EraseInternalNodeName();
    void EraseInternalNodeName(Link *from);

    // void Print(std::ostream& os,const Link* from) const ;
    // void Print(std::ostream& os) const;
    // printing int netwick format

    unsigned int GetSize() const { return GetSize(GetRoot()); }

    int GetSize(const Link *from) const {
        if (from->isLeaf()) {
            return 1;
        }
        int total = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetSize(link->Out());
        }
        return total;

        return 0;
    }

    int GetFullSize(const Link *from) const {
        if (from->isLeaf()) {
            return 1;
        }
        int total = 1;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            total += GetFullSize(link->Out());
        }
        return total;

        return 0;
    }

    virtual const Link *GetLCA(std::string tax1, std::string tax2) const {
        bool found1 = false;
        bool found2 = false;
        const Link *link = RecursiveGetLCA(GetRoot(), tax1, tax2, found1, found2);
        return link;
    }

    virtual const Link *GetLCA(const Link *from1, const Link *from2) const {
        bool found1 = false;
        bool found2 = false;
        const Link *link = RecursiveGetLCA(GetRoot(), from1, from2, found1, found2);
        return link;
    }

    void Subdivide(Link *from, int Ninterpol);

    std::string Reduce(const Link *from = nullptr) {
        if (from == nullptr) {
            from = GetRoot();
        }
        if (from->isLeaf()) {
            std::cerr << from->GetNode()->GetName() << '\n';
            ;
            return from->GetNode()->GetName();
        }
        std::string name = "None";
        for (Link *link = from->Next(); link != from; link = link->Next()) {
            std::string tmp = Reduce(link->Out());
            if (tmp == "diff") {
                name = "diff";
            } else if (name == "None") {
                name = tmp;
            } else if (name != tmp) {
                name = "diff";
            }
        }
        std::cerr << '\t' << name << '\n';
        from->GetNode()->SetName(name);
        return name;

        return "";
    }

    void PrintReduced(std::ostream &os, const Link *from = nullptr) {
        if (from == nullptr) {
            from = GetRoot();
        }
        if (from->GetNode()->GetName() != "diff") {
            os << from->GetNode()->GetName();
        } else {
            os << '(';
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                PrintReduced(os, link->Out());
                if (link->Next() != from) {
                    os << ',';
                }
            }
            os << ')';
        }
        if (from->isRoot()) {
            os << ";\n";
        }
    }

    const Link *ChooseInternalNode() const {
        int n = CountInternalNodes(GetRoot());
        int m = (int)(n * Random::Uniform());
        const Link *tmp;
        const Link *chosen = ChooseInternalNode(GetRoot(), tmp, m);
        if (chosen == nullptr) {
            std::cerr << "error in choose internal node: null pointer\n";
            exit(1);
        }
        return chosen;
    }

    int CountInternalNodes(const Link *from) const;
    const Link *ChooseInternalNode(const Link *from, const Link *&fromup, int &n) const;
    int CountNodes(const Link *from) const;
    const Link *ChooseNode(const Link *from, const Link *&fromup, int &n) const;

    // index links, nodes and branches through a recursive traversal of the tree
    // nodes: tip nodes are assumed already numbered between 0 and Ntaxa-1 based on their
    // correspondance with data
    // branches: start at 0
    // alternative (currently inactivated): branches: start at 1 (there is a null branch (null
    // pointer) behind the root, of index 0)
    void SetIndices() {
        Nlink = 0;
        Nnode = GetSize();
        Nbranch = 0;
        linkmap.clear();
        nodemap.clear();
        branchmap.clear();
        SetIndices(GetRoot(), Nlink, Nnode, Nbranch);
    }

    int GetNlink() const { return Nlink; }

    int GetNbranch() const { return Nbranch; }

    int GetNnode() const { return Nnode; }

    const Node *GetNode(int index) const { return nodemap[index]; }
    const Branch *GetBranch(int index) const { return branchmap[index]; }

    Link *GetLink(int index) const { return linkmap[index]; }

  protected:
    mutable map<int, const Node *> nodemap;
    mutable map<int, const Branch *> branchmap;
    mutable map<int, Link *> linkmap;

    void CheckIndices(Link *from) const {
        if (!from->isRoot()) {
            if (from->GetBranch() != branchmap[from->GetBranch()->GetIndex()]) {
                cerr << "branch index : " << from->GetBranch()->GetIndex() << '\n';
                exit(1);
            }
        } else {
            if (branchmap[0] != 0) {
                cerr << "root branch index\n";
                exit(1);
            }
        }

        if (from->GetNode() != nodemap[from->GetNode()->GetIndex()]) {
            cerr << "node index: " << from->GetNode()->GetIndex() << '\n';
            exit(1);
        }

        if (!from->isRoot()) {
            if (from->Out() != linkmap[from->Out()->GetIndex()]) {
                cerr << "link index : " << from->Out()->GetIndex() << '\n';
            }
        }
        if (from != linkmap[from->GetIndex()]) {
            cerr << "link index : " << from->GetIndex() << '\n';
        }


        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            CheckIndices(link->Out());
        }
    }

    void SetIndices(Link *from, int &linkindex, int &nodeindex, int &branchindex) {
        if (!from->isRoot()) {
            from->GetBranch()->SetIndex(branchindex);
            branchmap[branchindex] = from->GetBranch();
            branchindex++;
        }

        if (!from->isLeaf()) {
            from->GetNode()->SetIndex(nodeindex);
            nodemap[nodeindex] = from->GetNode();
            nodeindex++;
        } else {
            nodemap[from->GetNode()->GetIndex()] = from->GetNode();
        }

        if (!from->isRoot()) {
            from->Out()->SetIndex(linkindex);
            linkmap[linkindex] = from->Out();
            linkindex++;
        }
        from->SetIndex(linkindex);
        linkmap[linkindex] = from;
        linkindex++;

        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            SetIndices(link->Out(), linkindex, nodeindex, branchindex);
        }
    }

  protected:
    // returns 0 if not found
    // returns link if found (then found1 and found2 must
    const Link *RecursiveGetLCA(const Link *from, std::string tax1, std::string tax2, bool &found1,
                                bool &found2) const {
        const Link *ret = nullptr;
        if (from->isLeaf()) {
            // found1 |= (from->GetNode()->GetName() == tax1);
            // found2 |= (from->GetNode()->GetName() == tax2);
            std::string name1 = GetLeafNodeName(from).substr(0, tax1.size());
            std::string name2 = GetLeafNodeName(from).substr(0, tax2.size());
            found1 |= static_cast<int>(name1 == tax1);
            found2 |= static_cast<int>(name2 == tax2);
            /*
              found1 |= (GetLeafNodeName(from) == tax1);
              found2 |= (GetLeafNodeName(from) == tax2);
            */
            if (ret == nullptr) {
                if (found1 && found2) {
                    ret = from;
                }
            }
        } else {
            for (const Link *link = from->Next(); link != from; link = link->Next()) {
                bool tmp1 = false;
                bool tmp2 = false;
                const Link *ret2 = RecursiveGetLCA(link->Out(), tax1, tax2, tmp1, tmp2);
                found1 |= static_cast<int>(tmp1);
                found2 |= static_cast<int>(tmp2);
                if (ret2 != nullptr) {
                    if (ret != nullptr) {
                        std::cerr << "error : found node twice\n";
                        std::cerr << tax1 << '\t' << tax2 << '\n';
                        ToStream(std::cerr, ret2->Out());
                        std::cerr << '\n';
                        ToStream(std::cerr, ret->Out());
                        std::cerr << '\n';
                        exit(1);
                    }
                    ret = ret2;
                }
            }
            if (ret == nullptr) {
                if (found1 && found2) {
                    ret = from;
                }
            }
        }
        return ret;
    }

    const Link *RecursiveGetLCA(const Link *from, const Link *from1, const Link *from2,
                                bool &found1, bool &found2) const {
        const Link *ret = nullptr;
        found1 |= static_cast<int>(from == from1);
        found2 |= static_cast<int>(from == from2);
        if (ret == nullptr) {
            if (found1 && found2) {
                ret = from;
            }
        }
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            bool tmp1 = false;
            bool tmp2 = false;
            const Link *ret2 = RecursiveGetLCA(link->Out(), from1, from2, tmp1, tmp2);
            found1 |= static_cast<int>(tmp1);
            found2 |= static_cast<int>(tmp2);
            if (ret2 != nullptr) {
                if (ret != nullptr) {
                    std::cerr << "error : found node twice\n";
                    std::cerr << from1 << '\t' << from2 << '\n';
                    ToStream(std::cerr, ret2->Out());
                    std::cerr << '\n';
                    ToStream(std::cerr, ret->Out());
                    std::cerr << '\n';
                    exit(1);
                }
                ret = ret2;
            }
        }
        if (ret == nullptr) {
            if (found1 && found2) {
                ret = from;
            }
        }
        return ret;
    }

    void ReadFromStream(std::istream &is);
    // reading a tree from a stream:
    // recursively invokes the two following functions

    Link *ParseGroup(std::string input, Link *from);
    // a group is an expression of one of the two following forms:
    //
    //  (Body)Node_name
    //  (Body)Node_name:Branch_name
    //
    // where Body is a list, and Node_name and Branch_name are 2     std::strings
    // Node_name may be an empty     std::string
    //
    // Node_name cannot contain the ':' character, but Branch_name can
    // thus, if the group reads "(BODY)A:B:C"
    // then Node_name = "A" and Branch_name = "B:C"

    Link *ParseList(std::string input, Node *node);
    // a list is an expression of the form X1,X2,...Xn
    // where Xi is a group

    void RecursiveClone(const Link *from, Link *to);
    // Used by Tree(const Tree* from)
    // only clone the Links, and their mutual relations
    // does not copy the Node or Branch objects

    // deletes the link structure
    // does not delete the Node or Branch objects
    void RecursiveDelete(Link *from);

    void SetRoot(Link *link) { root = link; }

    // data fields
    // just 2 pointers, to the root and to a list of taxa
    Link *root;
    const TaxonSet *taxset;
    mutable int Nlink;
    mutable int Nbranch;
    mutable int Nnode;
};

#endif  // TREE_H
