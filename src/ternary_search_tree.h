#ifndef TERNARY_SEARCH_TREE
#define TERNARY_SEARCH_TREE

#define TST_BUFFERS 256 // number of buffers
#define TST_BUFSIZE 512 // inital buffer size

template <class KEY> struct Tnode {
  KEY splitkey;
  int level; /* level of tree */
  unsigned long long
      count; /* count how often the key leading to this node exists */

  Tnode *lokid;
  Tnode *eqkid;
  Tnode *hikid;
};

// class ternary_search_tree stores multikey-data with fixed key length
template <class KEY> class ternary_search_tree {
protected:
  typedef Tnode<KEY> *Tptr;

  // these four variables are used to accelerate allocation of Tnode objects

  Tptr buf;           // pointer to current buffer
  long next_buf_size; // size of the buffer that will be allocated next
  long old_buf_size;  // size of the buffer that was allocated last
  long bufn;          // number of next buffer
  long freen;         // number of free nodes in current buffer

  Tptr freearr[TST_BUFFERS];

  const long len; // key length
  Tptr root;      // tree root node

  long nr_total_nodes;

public:
  ternary_search_tree(const long keylength);
  ~ternary_search_tree();
  int insert(const KEY *const key); // insert key vector
  long total_nodes() const { return nr_total_nodes; }
  template <class Evaluater> void traverse(Evaluater &eval);
};

template <class KEY>
ternary_search_tree<KEY>::ternary_search_tree(const long keylength)
    : bufn(0), freen(0), root(0), next_buf_size(TST_BUFSIZE), old_buf_size(0),
      len(keylength), nr_total_nodes(0) {}

// return 0 on SUCCESS
template <class KEY>
int ternary_search_tree<KEY>::insert(const KEY *const key) {
  KEY d;
  Tptr pp;

  Tptr *p = &root;
  long level = 0; // level goes up to len-1

  while ((pp = *p)) { // as long as we encounter already exisiting nodes, we
                      // stay inside this while loop
    if (((d = key[level] - pp->splitkey)) == 0) { // go to next tree level
      pp->count++;
      p = &(pp->eqkid);
      if ((++level) == len)
        return 0;
    } else if (d < 0) {
      p = &(pp->lokid); /* move left in the current level */
    } else {
      p = &(pp->hikid); /* move right in the current level */
    }
  }
  for (;;) { /* once we find a node that is not allocated (==0), we must create
                every next node */
    if (freen-- == 0) {
#ifdef VERBOSE
      mexPrintf("Allocating buffer %ld of size %ld  \n", bufn, next_buf_size);
#endif

      const long temp_buf_size = next_buf_size;
      if (bufn == TST_BUFFERS) {
        // mexErrMsgTxt("Ran out of available buffers for tree nodes");
        return -1; // FAILURE
      }
      buf = new Tnode<KEY>[next_buf_size];
      freearr[bufn++] = buf;
      freen = next_buf_size - 1;
      next_buf_size += old_buf_size; // fibonacci like increasing of buffer size
                                     // (this keeps overall number of
                                     // allocations small)
      old_buf_size = temp_buf_size;
    }

    *p = buf++;
    pp = *p;
    pp->splitkey = key[level];
    pp->count = 1; // this node is newly created, so count is set to one
    nr_total_nodes++;
    pp->level = level;
    pp->lokid = pp->eqkid = pp->hikid = 0;
    if ((++level) == len)
      return 0;

    p = &(pp->eqkid);
  }
}

// traverse tree in an arbitrary order, execute the given function object on
// every node that is not empty
template <class KEY>
template <class Evaluater>
void ternary_search_tree<KEY>::traverse(Evaluater &eval) {
  long last_buf_size = 0;
  long buf_size =
      TST_BUFSIZE; // the actual size of the buffer is increasing each iteration

  // traverse through all buffers that are completely filled
  for (long i = 0; i < bufn; i++) {
    const long temp_buf_size = buf_size;
    long number_nodes = buf_size;
    const Tptr b = (Tptr)freearr[i];

    if (i == (bufn - 1)) // last buffer may not be completely filled
      number_nodes = buf_size - freen;

    for (long j = 0; j < number_nodes; j++) {
      const Tptr p = b + j;
      eval(p->count, p->level);
    }

    buf_size += last_buf_size;
    last_buf_size = temp_buf_size;
  }
}

template <class KEY> ternary_search_tree<KEY>::~ternary_search_tree() {
  for (long i = 0; i < bufn; i++)
    delete[] freearr[i];
}

#endif
