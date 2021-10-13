/*
  Copyright (c) 2003 by Stefan Kurtz and The Institute for
  Genomic Research.  This is OSI Certified Open Source Software.
  Please see the file LICENSE for licensing information and
  the file ACKNOWLEDGEMENTS for names of contributors to the
  code base.
*/

#ifndef STREETYP_H
#define STREETYP_H

#include "st_types.h"
#include "st_arraydef.h"

/*
  A \texttt{Reference} consists of an \texttt{address} pointing a leaf,
  or to a branching node. The boolean \texttt{toleaf} is \texttt{True} if
  and only if \texttt{address} points to a leaf.
*/

typedef struct
{
  BOOL toleaf;
  Uint *address;
} Reference;            // \Typedef{Reference}

/*
  The following types are used for references to leaves and
  branching nodes, respectively. We will always identify a leaf and
  and branching node with their references.
*/

typedef Uint * Bref;    // \Typedef{Bref}
typedef Uint * Lref;    // \Typedef{Lref}

/*
  For each branching node we store five values, as described in Section
  \ref{Representation}. These values comprise the following structure.
*/

typedef struct
{
  Uint headposition,        // the head position of the branching node
       depth;               // the depth of the branching node
  Bref suffixlink;          // the suffix link is always to a branching node
  Reference firstchild,     // the reference to the first child
            branchbrother;  // the reference to the right brother; 
                            // if this doesn't exist then it's \texttt{NULL}
} Branchinfo;               // \Typedef{Branchinfo}

/*
  For each leaf, we store a reference to its right brother, which is 
  \texttt{NULL}, if the right brother does not exist. This is
  expressed in the type synonym \texttt{Leafinfo}. 
*/

typedef Reference Leafinfo;  // \Typedef{Leafinfo}

/*
  A suffix tree is implemented by the type \texttt{Suffixtree}. 
  This structure contains several components which are mostly only used 
  during the suffix tree construction. For applications, assume the 
  following definition. Note that the input sequence is represented
  as an array of elements of type \texttt{SYMBOL}. The latter
  is a synonym for \texttt{Uchar} by default.

@typedef struct
@{ 
@  SYMBOL *text;     // points to the input string
@  Uint textlen;     // the length of the input string
@  Uint *branchtab;  // stores the infos for the branching nodes
@  Uint *leaftab;    // stores the brother-references of the leaves
@} Suffixtree;       // \Typedef{Suffixtree}

*/

//\Ignore{

// Added by TLB
typedef struct {
  Uchar *head;                  // the head of the string.
  Uint min_width;               // minimum width of valid words ending on node's in-branch
  Uint max_width;               // maximum width of valid words ending on node's in-branch
  Uint pos_count;		// number of positive sequences at leaves of subtree
  Uint neg_count;		// number of negative sequences at leaves of subtree
  Uint *pos_seqno_bittab;   	// the bit string of pos. sequences at leaves in subtree
  Uint *neg_seqno_bittab;   	// the bit string of neg. sequences at leaves in subtree
  Uint *pos_seqno_list;   	// the list pos. sequences at leaves in subtree
  Uint *neg_seqno_list;   	// the list neg. sequences at leaves in subtree
  double dtc_sum;		// sum of distances of (left edge of) sites to the sequence centers
  double log_pvalue;		// hypergeometric (log) p-value of the counts
  Uint seqno;			// sequence number of the head.
  Uint seqpos;			// position within the sequence of the head.
  // The next value is only set if scanning both strands.
  Uchar *rc_head;		// the head of the reverse complement of the string.
  // The next values are used for score-based matching only.
  double *scores;		// scores of node for width w; index 0 is scores[w-min_width]
} Subtreedata;

// Added by TLB
// The auxiliary information about each subtree is stored in a stack.
DECLAREARRAYSTRUCT(Subtreedata);
typedef Subtreedata *Subtreedataptr;
DECLAREARRAYSTRUCT(Subtreedataptr);

struct Suffixtreetype
{
// Added by TLB
  ArraySubtreedata subtreedata;	// stack that holds information on subtrees	
  Uint *subtreedataindex;	// index of subtreedata indexed by headposition
  Uint *leafdataindex;		// index of subtreedata indexed by leafnumber
//
  Uint textlen,               // the length of the input string
       *leaftab,              // stores the brother-references of the leafs
       *branchtab,            // table TBranch
       *rootchildren;         // references to successors of root
  SYMBOL *text,               // points to the input string
         *sentinel;           // points to the position of the \(\$\)-symbol

  Uint nextfreeleafnum,       // the number of the next leaf
       headnodedepth,         // the depth of the headnode
       insertnode,            // the node the split edge leads to 
       insertprev,            // the edge preceeding the split edge
       smallnotcompleted,     // the number of small nodes in the current chain
       nextfreebranchnum,     // the number of the next free branch node
       onsuccpath,            // refers to node on success path of headnode
       currentdepth,          // depth of the new branch node
       branchnodeoffset,      // number of leafs in tree
       alphasize,             // the number of different characters in t
       maxbranchdepth,        // maximal depth of branching node
       largenode,             // number of large nodes
       smallnode,             // number of small nodes
       *setlink,              // address of a nil-reference
       *nextfreeleafptr,      // points to next free entry in leaftab
       *chainstart,           // address of the node, current chains starts at
       *nextfreebranch,       // reference to next free base addr. in branchtab
       *headnode,             // left component of head location
       currentbranchtabsize,  // current number of cells in branchtab
       *firstnotallocated,    // refers to the last address, such that at
                              // least \texttt{LARGEINTS} integers are 
                              // available. So a large node can be stored in 
                              // the available amount of space.
       *nonmaximal,           // bit table: if node with headposition \(i\) is 
                              // not maximal, then \(nonmaximal[i]\) is set.
       *leafcounts;           // holds counts of the number of leafs in subtree
                              // indexed by headposition
  BOOL setatnewleaf;          // nil-reference is stored in new leaf
  SYMBOL *headstart,          // these references represent the right component
         *headend,            // of the head location \((\overline{u},v)\). 
                              // \emph{headstart} refers to the first character
                              // of \(v\), and \emph{headend} to the last
                              // character. In case, \(v=\varepsilon\),
                              // \(\emph{headend}=\emph{NULL}\).
         *tailptr;            // points to the tail

#ifdef ST_DEBUG
  char * (*showsymbolstree)(SYMBOL,Uchar *);
  Uchar *alphabet;
  Uint splitleafedge,
       splitinternaledge,
       artificial,
       insertleafcalls,
       largelinks,
       largelinkwork,
       largelinklinkwork,
       multiplications,
       nodecount,
       *maxset;
  void *generalcounter;
#endif
#if (SYMBOLBYTES == 2) || (SYMBOLBYTES == 4)
  Sint lastcharindex;
#endif

};

typedef struct Suffixtreetype Suffixtree;

DECLAREARRAYSTRUCT(Bref);

//}

/*
  A location is implemented by the type \texttt{Location}.
*/

typedef struct 
{
  Stringtype locstring; // string represented by location
  Bref previousnode;    // reference to previous node (which is branching)
  SYMBOL *firstptr;     // pointer to first character of edge label
  Uint edgelen,         // length of edge
       remain;          // number of remaining characters on edge
  Reference nextnode;   // reference to node the edge points to
} Location;             // \Typedef{Location}

/*
  If a location is a node \(\overline{u}\), we set \texttt{remain} to 0, and
  store a reference to \(\overline{u}\) in \texttt{nextnode}. Moreover, we 
  store a position where \(u\) starts and its length in \texttt{locstring}.
  If the location is of the form \((\overline{u},v,w,\overline{uvw})\),
  then the components of the location satisfies the following values:
  \begin{enumerate}
  \item
  \texttt{previousnode} is a reference to \(\overline{u}\)
  \item
  \texttt{firstptr} points to the first symbol of the edge label \(vw\).
  \item
  \(\texttt{edgelen}=\Size{vw}\)
  \item
  \(\texttt{remain}=\Size{w}\)
  \item
  \texttt{nextnode} is a reference to \(\overline{uvw}\).
  \end{enumerate}
  Since \(w\) is not empty, a location is a node location if and only if
  \texttt{remain} is 0.
*/

//\Ignore{

/*
  A simple location stores just a part of information stored in a suffix tree.
*/

typedef struct 
{
  Uint remain,
       textpos;  // these last two items are redundant and can be computed
  Reference nextnode;
} Simpleloc;     // \Typedef{Simpleloc}

DECLAREARRAYSTRUCT(Simpleloc);

/*
  A path in the suffix tree is stored as an array of \texttt{Pathinfo}-records.
*/

typedef struct
{
  Uint depth, headposition;
  Bref ref;
} Pathinfo;      // \Typedef{Pathinfo}

DECLAREARRAYSTRUCT(Pathinfo);

typedef struct
{
  BOOL secondtime;
  ArrayBref stack;
} DFSstate;      // \Typedef{DFSstate}

#endif

//}
