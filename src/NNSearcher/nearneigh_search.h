#ifndef NEARNEIGH_SEARCH_H
#define NEARNEIGH_SEARCH_H

#include <vector>
#include <Rcpp.h>

// This file contains a class collection for nearest neighbor (NN) searching.
// It is possible to search for k (k=1..N) nearest neighbors to a given query point ("fixed-mass" approach)
// or to search all neighbors within a distance r (range) around the query point ("fixes size" approach)
// ATRIA : A class that implements an advanced triangle inequality algorithm
// During preprocessing, a search tree is constructed by dividing the set of data nearneigh_searcher<POINT_SET>::points
// in two (sub)clusters. Each cluster is than subdivided until a minimum number of nearneigh_searcher<POINT_SET>::points is reached
// During search, the triangle inequality is used to exclude cluster from further search
// ATRIA might be a good choice for unevenly distributed points in very high dimensional spaces.
// Here are some parameters that might need some tuning for your special case, but for general purposes these are OK
// Parameters for the ATRIA nearest neighbor search
// A cluster will not be further subdivided if it contains less
// than ATRIAMINPOINTS nearneigh_searcher<POINT_SET>::points
// A smaller value might accelerate actual search but increase pre-processing time.
// Memory consumption will not change very much when choosing a smaller  value of ATRIAMINPOINTS.
// Christian Merkwirth, DPI Goettingen, 1998 - 2000

#define ATRIAMINPOINTS 64

#include "nn_aux.h"
#include "utilities.h"

// Base class for nearest neighbor searchers.
template <class POINT_SET>
class nearneigh_searcher : protected My_Utilities {
protected:
  int err; // error state, == 0 means OK, every other value is a failure
  const POINT_SET points;

  long Nused; // Number of points of the data set actually used

  SortedNeighborTable table;

  // test point number #index of nearneigh_searcher<POINT_SET>::points
  template <class ForwardIterator>
  void test(const long index, ForwardIterator qp, const double thresh) {
#ifdef PARTIAL_SEARCH
    const double d = points.distance(index, qp, thresh);
#else
    const double d = points.distance(index, qp);
#endif
    if (d < thresh)
      table.insert(neighbor(index, d));
#ifdef PROFILE
    points_searched++;
#endif
  }

public:
  typedef POINT_SET point_set;

  // prepare searching for a point set nearneigh_searcher<POINT_SET>::points
  // excl gives the number of samples from the end of
  // nearneigh_searcher<POINT_SET>::points that should be omitted from searching
  nearneigh_searcher(POINT_SET &&p, const long excl = 0);
  ~nearneigh_searcher();

  long number_of_points() const { return Nused; }
  const POINT_SET& get_point_set() const { return points; }

  // Return the error state, = 0 means OK. An error will only occur when the object
  // was not initialized correctly.
  inline int geterr() const { return err; }
};

// Advanced triangle inequaltity algorithm
template <class POINT_SET> class ATRIA : public nearneigh_searcher<POINT_SET> {
protected:
  const long MINPOINTS;
  cluster* root;

  neighbor* const permutation_table;
  typedef typename POINT_SET::Metric METRIC;
  typedef searchitem SearchItem;

  priority_queue<SearchItem, vector<SearchItem>, searchitemCompare>
      search_queue;
  stack<SearchItem, vector<SearchItem> > SearchStack; // used for range searches/counts

  long total_clusters;
  long terminal_nodes;
  long total_points_in_terminal_node;
  unsigned long terminal_cluster_searched;
  unsigned long points_searched;
  unsigned long number_of_queries;

  void create_tree();
  void destroy_tree();

  pair<long, long> find_child_cluster_centers(const cluster* const c,
                                              neighbor* const Section,
                                              const long c_length);
  long assign_points_to_centers(neighbor* const Section, const long c_length,
                                pair<cluster*, cluster*> childs);

  template <class ForwardIterator>
  void search(ForwardIterator query_point, const long first, const long last,
              const double epsilon);

public:
  ATRIA(POINT_SET &&p, const long excl = 0, const long minpts = ATRIAMINPOINTS);
  ~ATRIA();

  // Search for k nearest neighbors of the point query_point, excluding
  // points with indices between first and last from the search. Returns a
  // sorted vector of neighbors (by reference).
  template <class ForwardIterator>
  long search_k_neighbors(vector<neighbor> &v, const long k,
                          ForwardIterator query_point, const long first = -1,
                          const long last = -1, const double epsilon = 0);

  // Count the number of points within distance 'radius' from the query point,
  // excluding points with indices between first and last from the search.
  template <class ForwardIterator>
  long count_range(const double radius, ForwardIterator query_point,
                   const long first = -1, const long last = -1);

  // Search points within distance 'radius' from the query point,  excluding points
  // with indices between first and last  Returns an unsorted vector v of neigbors by
  // reference.
  template <class ForwardIterator>
  long search_range(vector<neighbor> &v, const double radius,
                    ForwardIterator query_point, const long first = -1,
                    const long last = -1);

  // Returns an approximation of the data set radius such that any pairwise
  // distance in the data set is smaller than twice this radius. This bound is
  // not necessarily tight.
  inline double data_set_radius() const { return root->Rmax; };
  inline long total_tree_nodes() const { return total_clusters; };
  double search_efficiency() const {
    return ((double)points_searched) / ((double)nearneigh_searcher<POINT_SET>::number_of_points() * number_of_queries);
  }
};

template <class POINT_SET>
nearneigh_searcher<POINT_SET>::nearneigh_searcher(POINT_SET&& p, const long excl)
    : err(0), points(std::move(p)), Nused(p.size() - excl)
{
  if ((nearneigh_searcher<POINT_SET>::Nused < 1) || (excl < 0)) {
    Rcpp::Rcerr << "Wrong parameters for nearest neighbour search" <<std::endl;
    Rcpp::Rcerr << "Nused : " << nearneigh_searcher<POINT_SET>::Nused <<std::endl;
    nearneigh_searcher<POINT_SET>::err = 1;
    return;
  }
}

template <class POINT_SET>
nearneigh_searcher<POINT_SET>::~nearneigh_searcher() {}

template <class POINT_SET>
ATRIA<POINT_SET>::ATRIA(POINT_SET &&p, const long excl, const long minpts)
    : nearneigh_searcher<POINT_SET>(std::move(p), excl), MINPOINTS(minpts), root(nullptr),
      permutation_table(new neighbor[nearneigh_searcher<POINT_SET>::Nused]),
      total_clusters(1), terminal_nodes(0), total_points_in_terminal_node(0),
      terminal_cluster_searched(0), points_searched(0), number_of_queries(0) {

#ifdef VERBOSE
  Rcpp::Rcout << "ATRIA Constructor" <<std::endl;
  Rcpp::Rcout << "Size of point set : " << p.size() << "  points of dimension " << p.dimension() <<std::endl;
  Rcpp::Rcout << "Number of points used : " << nearneigh_searcher<POINT_SET>::number_of_points() <<std::endl;
  Rcpp::Rcout << "MINPOINTS : " << MINPOINTS <<std::endl;
#endif
  terminal_cluster_searched = 0;
  if (nearneigh_searcher<POINT_SET>::err) {
    Rcpp::Rcerr << "Error initializing parent object" <<std::endl;
    return;
  }

  if (permutation_table == 0) {
    Rcpp::Rcerr << "Out of memory" <<std::endl;
    nearneigh_searcher<POINT_SET>::err = 1;
    return;
  }
  create_tree();

#ifdef VERBOSE
  Rcpp::Rcout << "Created tree structure for ATRIA searcher" <<std::endl;
#endif
}

template <class POINT_SET> ATRIA<POINT_SET>::~ATRIA() {
#ifdef VERBOSE
  Rcpp::Rcout << "ATRIA Destructor" <<std::endl;
  Rcpp::Rcout << "Total_clusters : " << total_clusters <<std::endl;
  Rcpp::Rcout << "Total number of points in terminal nodes : "
        << total_points_in_terminal_node <<std::endl;
  Rcpp::Rcout << "Average number of points in a terminal node : "
       << ((double)total_points_in_terminal_node) / terminal_nodes <<std::endl;
  if (number_of_queries == 0)
    Rcpp::Rcout << "No queries were done" <<std::endl;
  else
    Rcpp::Rcout << "Average percentage of points searched "
         << (100.0 * (double)points_searched) /
    ((double)nearneigh_searcher<POINT_SET>::Nused *
    number_of_queries)
    << "% (" << ceil(((double)points_searched) / number_of_queries)
    << ")" <<std::endl;
    Rcpp::Rcout << "Average number of terminal nodes visited : "
         << ((double)terminal_cluster_searched) / number_of_queries <<std::endl;
#endif

  destroy_tree();

  delete[] permutation_table;
}

template <class POINT_SET>
pair<long, long> ATRIA<POINT_SET>::find_child_cluster_centers(
    const cluster* const c, neighbor* const Section, const long length) {
  pair<long, long> centers(-1, -1);

  if (c->Rmax == 0) { // if all data nearneigh_searcher<POINT_SET>::points seem
                      // to be identical
#ifdef VERBOSE
    Rcpp::Rcout << "ATRIA : Data seem to be singular, search may be very inefficient"
          <<std::endl;
#endif
    terminal_nodes++;
    total_points_in_terminal_node += length;

    return centers; // indicate that there's no need to further divide this data
                    // set
  }

  // Compute right center, the point that is farthest away from the c->center.
  // We assume that distances in Section are with respect to the c->center.
  long index = 0;
  long center_right = Section[index].index();
  double dist = Section[index].dist();
  // nearneigh_searcher<POINT_SET>::points.distance(c->center,// Section[index].index());
  for (long i = 1; i < length; i++) {
    const double d = Section[i].dist();
    //cout << d <<std::endl;
    if (d > dist) {
      dist = d;
      center_right = Section[i].index();
      index = i;
    }
  }
  centers.second = center_right;
  // move this center the the last (rightmost) element of this Section
  swap(Section, index, length - 1);
  // Compute left center, the point that is farthest away from the center_right.
  // We also overwrite distances in Section to distances wrt the right center.
  index = 0;
  long center_left = Section[index].index();
  dist = nearneigh_searcher<POINT_SET>::points.distance(center_right,
                                                        Section[index].index());
  //cout << dist <<std::endl;
  Section[index].dist() = dist;
  for (long i = 1; i < length - 1; i++) {
    const double d = nearneigh_searcher<POINT_SET>::points.distance(
        center_right, Section[i].index());
    Section[i].dist() = d;
    //cout << center_right << " " << i << " " <<  Section[i].index() << "  " << d << " " << Section[i].dist() <<std::endl;
    if (d > dist) {
      dist = d;
      center_left = Section[i].index();
      index = i;
    }
  }
  // move this center the the first (leftmost) element of this Section
  swap(Section, index, 0);
  centers.first = center_left;

  //cout << "Centers: " << center_left << " " << center_right <<std::endl;
  return centers;
}

// assign each point to the nearest center, using a kind of quicksort like
// sorting procedure
template <class POINT_SET>
long ATRIA<POINT_SET>::assign_points_to_centers(
    neighbor *const Section, const long c_length,
    pair<cluster *, cluster *> childs) {
  const long center_left = childs.first->center;
  long i = 0;
  long j = c_length - 1;

  // maximal distance fron one cluster's center to
  // nearneigh_searcher<POINT_SET>::points belonging to this cluster
  double Rmax_left = 0;
  double Rmax_right = 0;

  while (1) {
    short i_belongs_to_left = 1;
    short j_belongs_to_right = 1;

    while (i + 1 < j) {
      i++;
      const double dl = nearneigh_searcher<POINT_SET>::points.distance(
          center_left, Section[i].index());
      // reuse information instead of calculating dr =
      //const double dr = nearneigh_searcher<POINT_SET>::points.distance(center_right,
      //  Section[i].index());
      const double dr = Section[i].dist();

      if (dl > dr) {
        // point belongs to the right corner
        Section[i].dist() = dr;
        i_belongs_to_left = 0;

        Rmax_right = max(Rmax_right, dr);
        break;
      }
      // point belongs to the left corner
      Section[i].dist() = dl;
      Rmax_left = max(Rmax_left, dl);
    }

    while (j - 1 > i) { // either we reached the start of the array or this
                        // element was already checked by the i loop
      --j;

      const double dr =
          Section[j]
              .dist(); // nearneigh_searcher<POINT_SET>::points.distance(center_right,
                       // Section[j].index());
      const double dl = nearneigh_searcher<POINT_SET>::points.distance(
          center_left, Section[j].index());

      if (dr >= dl) {
        // point belongs to the left corner
        Section[j].dist() = dl;
        j_belongs_to_right = 0;

        Rmax_left = max(Rmax_left, dl);
        break;
      }
      // point belongs to the right corner
      Section[j].dist() = dr;
      Rmax_right = max(Rmax_right, dr);
    }

    if (i == j - 1) {
      if ((!i_belongs_to_left) && (!j_belongs_to_right)) {
        swap(Section, i, j);
      } else if (!i_belongs_to_left) {
        i--;
        j--;
      } else if (!j_belongs_to_right) {
        i++;
        j++;
      }
      break; // finished walking through the array
    } else {
      swap(Section, i, j);
    }
  }

  childs.first->Rmax = Rmax_left;
  childs.second->Rmax = Rmax_right;

  return j;
}

template <class POINT_SET> void ATRIA<POINT_SET>::create_tree() {
  long k;

  if (nearneigh_searcher<POINT_SET>::err)
    return;

  std::stack<cluster_pointer, cluster_pointer_vector> Stack; // used for tree construction

  // select random center for root cluster, move this to first position of the
  // indices array
  root = new cluster(1, nearneigh_searcher<POINT_SET>::Nused - 1);
  root->center = My_Utilities::randindex(nearneigh_searcher<POINT_SET>::Nused);
  permutation_table[0] = neighbor(root->center, 0);

  root->Rmax = 0;
  for (k = 0; k < root->center; k++) {
    const double d =
        nearneigh_searcher<POINT_SET>::points.distance(k, root->center);
    permutation_table[k + 1] = neighbor(k, d);
    if (d > root->Rmax)
      root->Rmax = d;
  }
  for (k = root->center + 1; k < nearneigh_searcher<POINT_SET>::Nused; k++) {
    const double d =
        nearneigh_searcher<POINT_SET>::points.distance(k, root->center);
    permutation_table[k] = neighbor(k, d);
    if (d > root->Rmax)
      root->Rmax = d;
  }

#ifdef VERBOSE
  Rcpp::Rcout << "Root center : " << root->center <<std::endl;
  Rcpp::Rcout << "Root starting index  : " << root->start <<std::endl;
  Rcpp::Rcout << "Root length : " << root->length <<std::endl;
  Rcpp::Rcout << "Root Rmax : " << root->Rmax <<std::endl;
#endif

  // now create the tree
  Stack.push(root); // push root cluster on Stack

  while (!Stack.empty()) {
    cluster* const c = Stack.top();
    Stack.pop();

    const long c_start = c->start;
    const long c_length = c->length;

    neighbor* const Section = permutation_table + c_start;

    if (c->length >= MINPOINTS) { // Further divide this cluster ?
      pair<long, long> new_child_centers =
          find_child_cluster_centers((const cluster*) c, Section, c_length);

      if ((new_child_centers.first == -1) || (new_child_centers.second == -1))
        continue; // cluster could not be divided further and has already been
                  // marked as terminal cluster/node

      c->left = new cluster(new_child_centers.first);
      c->right = new cluster(new_child_centers.second);

      // create two subclusters and set properties
      const long j = assign_points_to_centers(Section, c_length, pair<cluster *, cluster *>(c->left, c->right));

      c->left->start = c_start + 1; // leave centers out
      c->left->length = j - 1;

      c->right->start = c_start + j;
      c->right->length = c_length - j - 1;

      // process new subclusters (use stacks to avoid recursive call of this
      // function)
      Stack.push(c->right);
      Stack.push(c->left);

      total_clusters += 2;
#ifdef DEBUG
      Rcpp::Rcout << "Cluster start " << c_start << " len: " << c_length << " radius: " << c->Rmax << " ";
      Rcpp::Rcout << "left child ctr & len: " << (c->left)->center << " " << (c->left)->length << "  ";
      Rcpp::Rcout << "right child ctr & len: " << (c->right)->center << " "  << (c->right)->length <<std::endl;
#endif
    } else {              // this is going to be a terminal node
      c->Rmax = -c->Rmax; // a Rmax value <= 0 marks this cluster as a terminal
                          // node of the search tree
      terminal_nodes++;
      total_points_in_terminal_node += c_length;
#ifdef DEBUG
      Rcpp::Rcout << "Terminal node " << c->length << "  " << c->center << "  "
           << c->Rmax <<std::endl;
#endif
    }
  }
}

template <class POINT_SET> void ATRIA<POINT_SET>::destroy_tree() {
  stack<cluster_pointer, cluster_pointer_vector>
      Stack; // used for tree construction

  if (!root->is_terminal()) {
    Stack.push(root->left);
    Stack.push(root->right);
  }

  while (!Stack.empty()) {
    cluster *c = Stack.top();
    Stack.pop();

    if (c != 0) {
      if (!c->is_terminal()) {
        Stack.push(c->left);
        Stack.push(c->right);
      }
      delete c;
    }
  }
}

template <class POINT_SET>
template <class ForwardIterator>
long ATRIA<POINT_SET>::search_k_neighbors(vector<neighbor> &v, const long k,
                                          ForwardIterator query_point,
                                          const long first, const long last,
                                          const double epsilon) {
  number_of_queries++;
  nearneigh_searcher<POINT_SET>::table.init_search(k);

  search(query_point, first, last, epsilon);

  return nearneigh_searcher<POINT_SET>::table.finish_search(
      v); // append nearneigh_searcher<POINT_SET>::table items to v, v should be
          // empty, afterwards nearneigh_searcher<POINT_SET>::table is empty
}

template <class POINT_SET>
template <class ForwardIterator>
void ATRIA<POINT_SET>::search(ForwardIterator query_point, const long first,
                              const long last, const double epsilon) {
  points_searched++;
  const double root_dist =
      nearneigh_searcher<POINT_SET>::points.distance(root->center, query_point);

  while (!search_queue.empty())
    search_queue.pop(); // clear search queue

  // push root cluster as search item into the PR-QUEUE
  search_queue.push(SearchItem(root, root_dist));

  while (!search_queue.empty()) {
    const SearchItem si = search_queue.top();
    search_queue.pop();
    const cluster *const c = si.clusterp();

    if ((nearneigh_searcher<POINT_SET>::table.highdist() > si.dist()) &&
        ((c->center < first) || (c->center > last)))
      nearneigh_searcher<POINT_SET>::table.insert(
          neighbor(c->center, si.dist()));

    if (nearneigh_searcher<POINT_SET>::table.highdist() >=
        (si.d_min() *
         (1.0 +
          epsilon))) { // approximative (epsilon > 0) queries are supported
      if (c->is_terminal()) {
        const neighbor *const Section = permutation_table + c->start;
        terminal_cluster_searched++;

        if (c->Rmax == 0.0) {
          for (long i = 0; i < c->length; i++) {
            const long j = Section[i].index();

            if (nearneigh_searcher<POINT_SET>::table.highdist() <= si.dist())
              break;

            if ((j < first) || (j > last))
              nearneigh_searcher<POINT_SET>::table.insert(
                  neighbor(j, si.dist()));
          }
        } else {
          for (long i = 0; i < c->length; i++) {
            const long j = Section[i].index();

            if ((j < first) || (j > last)) {
              if (nearneigh_searcher<POINT_SET>::table.highdist() >
                  fabs(si.dist() - Section[i].dist()))
                nearneigh_searcher<POINT_SET>::test(j, query_point, nearneigh_searcher<POINT_SET>::table.highdist());
            }
          }
        }
      } else { // this is an internal node
        const double dl = nearneigh_searcher<POINT_SET>::points.distance(
            c->left->center, query_point);
        const double dr = nearneigh_searcher<POINT_SET>::points.distance(
            c->right->center, query_point);
        points_searched += 2;
        // create child cluster search items
        SearchItem si_left = SearchItem(c->left, dl, dr, si);
        SearchItem si_right = SearchItem(c->right, dr, dl, si);

        // priority based search
        search_queue.push(si_right);
        search_queue.push(si_left);
      }
    }
  }
}

template <class POINT_SET>
template <class ForwardIterator>
long ATRIA<POINT_SET>::search_range(vector<neighbor> &v, const double radius,
                                    ForwardIterator query_point,
                                    const long first, const long last) {
  long count = 0;

  number_of_queries++;
  points_searched++;

  while (!SearchStack.empty())
    SearchStack.pop(); // make shure stack is empty

  SearchStack.push(
      SearchItem(root, nearneigh_searcher<POINT_SET>::points.distance(
                            root->center, query_point)));

  while (!SearchStack.empty()) {
    const SearchItem si = SearchStack.top();
    SearchStack.pop();

    if (radius >= si.d_min()) {
      const cluster *const c = si.clusterp();

      if (((c->center < first) || (c->center > last)) &&
          (si.dist() <= radius)) {
        v.push_back(neighbor(c->center, si.dist()));
        count++;
      }

      if (c->is_terminal()) { // this is a terminal node
        const neighbor *const Section = permutation_table + c->start;

        if (c->Rmax == 0.0) { // cluster has zero radius, so all
                              // nearneigh_searcher<POINT_SET>::points inside
                              // will have the same distance to q
          if (radius >= si.dist()) {
            for (long i = 0; i < c->length; i++) {
              const long j = Section[i].index();

              if ((j < first) || (j > last)) {
                v.push_back(neighbor(j, si.dist()));
                count++;
              }
            }
          }
        } else {
          for (long i = 0; i < c->length; i++) {
            const long j = Section[i].index();

            if (((j < first) || (j > last)) &&
                (radius >= fabs(si.dist() - Section[i].dist()))) {
#ifdef PARTIAL_SEARCH
              const double d = nearneigh_searcher<POINT_SET>::points.distance(
                  j, query_point, radius);
#else
              const double d = nearneigh_searcher<POINT_SET>::points.distance(
                  j, query_point);
#endif
              if (d <= radius) {
                v.push_back(neighbor(j, d));
                count++;
              }
              points_searched++;
            }
          }
        }
        terminal_cluster_searched++;
      } else { // this is an internal node
        const double dl = nearneigh_searcher<POINT_SET>::points.distance(
            c->left->center, query_point);
        const double dr = nearneigh_searcher<POINT_SET>::points.distance(
            c->right->center, query_point);
        points_searched += 2;
        const SearchItem x = SearchItem(c->left, dl, dr, si);
        const SearchItem y = SearchItem(c->right, dr, dl, si);

        SearchStack.push(x);
        SearchStack.push(y);
      }
    }
  }

  return count;
}

template <class POINT_SET>
template <class ForwardIterator>
long ATRIA<POINT_SET>::count_range(const double radius,
                                   ForwardIterator query_point,
                                   const long first, const long last) {
  long count = 0;

#ifdef PROFILE
  number_of_queries++;
  points_searched++;
#endif

  while (!SearchStack.empty())
    SearchStack.pop(); // make shure stack is empty

  SearchStack.push(
      SearchItem(root, nearneigh_searcher<POINT_SET>::points.distance(
                            root->center, query_point)));

  while (!SearchStack.empty()) {
    const SearchItem si = SearchStack.top();

    SearchStack.pop();

    if (radius >= si.d_min()) {
      const cluster *const c = si.clusterp();

      if (((c->center < first) || (c->center > last)) &&
          (si.dist() <= radius)) {
        count++;
      }

      if (c->is_terminal()) { // this is a terminal terminal node
        const neighbor *const Section = permutation_table + c->start;

        if (c->Rmax == 0.0) { // cluster has zero radius, so all
                              // nearneigh_searcher<POINT_SET>::points inside
                              // will have the same distance to q
          if (radius >= si.dist()) {
            for (long i = 0; i < c->length; i++) {
              const long j = Section[i].index();

              if ((j < first) || (j > last)) {
                count++;
              }
            }
          }
        } else {
          for (long i = 0; i < c->length; i++) {
            const long j = Section[i].index(); // index of Vergleichspunkt

            if (((j < first) || (j > last)) &&
                (radius >= fabs(si.dist() - Section[i].dist()))) {
#ifdef PARTIAL_SEARCH
              if (nearneigh_searcher<POINT_SET>::points.distance(
                      j, query_point, radius) <= radius)
                count++;
#else
              if (nearneigh_searcher<POINT_SET>::points.distance(
                      j, query_point) <= radius)
                count++;
#endif
              points_searched++;
            }
          }
        }
        terminal_cluster_searched++;
      } else { // this is an internal node
        const double dl = nearneigh_searcher<POINT_SET>::points.distance(
            c->left->center, query_point);
        const double dr = nearneigh_searcher<POINT_SET>::points.distance(
            c->right->center, query_point);
        points_searched += 2;

        const SearchItem x = SearchItem(c->left, dl, dr, si);
        const SearchItem y = SearchItem(c->right, dr, dl, si);

        SearchStack.push(x);
        SearchStack.push(y);
      }
    }
  }

  return count;
}

#endif // ifdef NEARNEIGH_SEARCH_H
