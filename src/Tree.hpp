/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

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

#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <map>
#include <sstream>  // only needed for tests
#include <string>

#include <doctest.h>

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

class Tree {
  public:
    // create a tree by reading into a file (netwick format expected) ; calls ReadFromStream
    explicit Tree(std::string filename);

    Tree(const Tree &) = delete;  // forbidding copy construction
    Tree() = delete;

    void ToStream(std::ostream &os) const;  // calls recursive private ToStream method

    // Delete the leaf pointing by the next link and set everithing right.
    void DeleteNextLeaf(Link *previous);  // TODO

    // Delete the unary Node wich from is paart of and set everithing right.
    void DeleteUnaryNode(Link *from);  // TODO

    Link *GetRoot() const { return root; }
    const TaxonSet *GetTaxonSet() const { return taxset; }

    void RootAt(Link *from);

    // Registers all leaves of the tree with an external TaxonSet
    // the taxon set defines a map between taxon names and indices (between 0 and
    // P-1)
    // the tree is recursively traversed
    // each leaf's name is looked for in the map of the taxon set
    // if not found : an error is emitted
    // otherwise, the leaf's index is set equal to the index of the corresponding
    // taxon
    void RegisterWith(const TaxonSet *taxset);

    bool RegisterWith(const TaxonSet *taxset, Link *from,
                      int &tot);  // recursive function called by RegisterWith
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
    }

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

    std::string GetBranchName(const Link *link) const { return link->GetBranch()->GetName(); }
    std::string GetNodeName(const Link *link) const { return link->GetNode()->GetName(); }
    int GetNlink() const { return Nlink; }
    int GetNbranch() const { return Nbranch; }
    int GetNnode() const { return Nnode; }
    const Node *GetNode(int index) const { return nodemap.at(index); }
    const Branch *GetBranch(int index) const { return branchmap.at(index); }
    Link *GetLink(int index) const { return linkmap.at(index); }

  private:
    // data fields
    // just 2 pointers, to the root and to a list of taxa
    Link *root{nullptr};
    const TaxonSet *taxset{nullptr};
    int Nlink;
    int Nbranch;
    int Nnode;
    std::map<int, const Node *> nodemap;
    std::map<int, const Branch *> branchmap;
    std::map<int, Link *> linkmap;

    void ToStream(std::ostream &os, const Link *from) const;

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

    // reading a tree from a stream:
    // recursively invokes the two following functions
    void ReadFromStream(std::istream &is);

    // a group is an expression of one of the two following forms:
    //  (Body)Node_name
    //  (Body)Node_name:Branch_name
    // where Body is a list, and Node_name and Branch_name are 2     std::strings
    // Node_name may be an empty     std::string
    // Node_name cannot contain the ':' character, but Branch_name can
    // thus, if the group reads "(BODY)A:B:C"
    // then Node_name = "A" and Branch_name = "B:C"
    Link *ParseGroup(std::string input, Link *from);

    // a list is an expression of the form X1,X2,...Xn
    // where Xi is a group
    Link *ParseList(std::string input, Node *node);

    void SetRoot(Link *link) { root = link; }
};

/*
#### TESTS #########################################################################################
*/
TEST_CASE("Tree test") {
    Tree mytree{"data/toy.tree"};
    mytree.SetIndices();
    CHECK(mytree.GetNnode() == 7);  // 4 leaves + 3 internal nodes
    CHECK(mytree.GetNbranch() == 6);
    CHECK(mytree.GetNlink() == 13);  // assuming 2 by branch + 1 for the root

    std::stringstream ss;
    mytree.ToStream(ss);
    CHECK(ss.str() ==
          "((S0:0,S1:1):0,(S2:0,S3:1):0);\n");  // apparently it removes the ':0' at the root

    // for (int i = 0; i < 50; i++) {
    //     try {
    //         printf("%d: %s\n", i, mytree.GetNode(i)->GetName().c_str());
    //     } catch (...) {
    //     }
    // }
}

#endif  // TREE_H
