#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <map>
#include <string>

class TaxonSet;  // forward decl

class Node {
  private:
    int index{0};
    std::string name{""};

  public:
    Node() = default;
    explicit Node(const std::string &s) : name(s) {}
    explicit Node(const Node *from) : index(from->index), name(from->name) {}

    std::string GetName() const { return name; }
    void SetName(std::string inname) { name = inname; }
    int GetIndex() const { return index; }
    void SetIndex(int i) { index = i; }
};

using Branch = Node;  // the two classes were completely identical

class Link {
    Link *next;
    Link *out;
    Branch *branch{nullptr};
    Node *node{nullptr};
    int index{-1};

  public:
    Link() : next(this), out(this) {}

    Link *Next() const { return next; }
    Link *Out() const { return out; }
    Branch *GetBranch() const { return branch; }
    Node *GetNode() const { return node; }
    int GetIndex() const { return index; }

    void SetBranch(Branch *inbranch) { branch = inbranch; }
    void SetNode(Node *innode) { node = innode; }
    void SetOut(Link *inout) { out = inout; }
    void SetIndex(int i) { index = i; }
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
};

class NewickTree {
  public:
    virtual ~NewickTree() = default;

    virtual Link *GetRoot() const = 0;

    void ToStream(std::ostream &os) const;
    void ToStream(std::ostream &os, const Link *from) const;
    double ToStreamSimplified(std::ostream &os, const Link *from) const;

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

    explicit Tree(const Tree *from);
    // clones the entire Link structure
    // but does NOT clone the Nodes and Branches
    // calls RecursiveClone

    explicit Tree(std::string filename);
    // create a tree by reading into a file (netwick format expected)
    // calls ReadFromStream

    ~Tree();
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

    std::string GetBranchName(const Link *link) const override {
        return link->GetBranch()->GetName();
    }

    std::string GetNodeName(const Link *link) const override { return link->GetNode()->GetName(); }

    void EraseInternalNodeName();
    void EraseInternalNodeName(Link *from);

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
    mutable std::map<int, const Node *> nodemap;
    mutable std::map<int, const Branch *> branchmap;
    mutable std::map<int, Link *> linkmap;

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
