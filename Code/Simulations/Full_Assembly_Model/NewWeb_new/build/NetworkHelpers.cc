// $Id: NetworkHelpers.cc 2351 2013-12-18 21:10:51Z axel $

#include <set>
#include <list>
#include <stack>
#include <algorithm>
#include <cmath>
#include "NetworkAnalysis.h"
#include "CMatrix.h"
#include "Statistics.h"
#include "random.h"

// We also use the BOOST generic functions 
#include <boost/graph/vector_as_graph.hpp> // Boost bug #2119

#ifndef BOOST_VERSION 
#include <boost/version.hpp> 
#endif 
#if BOOST_VERSION > 104100 
#warning --------- FIXME ----------
#warning YOU HAVE A BOOST VERSION LATER THAN 1.41
#warning THIS CODE WILL COMPILE BUT NETWORK HELPERS MIGHT GIVE 
#warning AN ERROR MESSEAG WHEN EXECUTED.
#warning --------- END FIXME ------
#include <boost/graph/adjacency_matrix.hpp>
#else // adjacency_matrix needs some fixing in older versions:
#if BOOST_VERSION <= 104100 
#include "my_adjacency_matrix.hpp"
#else
#error my_adjacency_matrix.hpp will not work with latest versions of boost :(
#endif
#endif

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104000 
#include <boost/property_map/vector_property_map.hpp>
#else
#include <boost/vector_property_map.hpp>
#endif
#include <boost/graph/transitive_closure.hpp>

///////// This does not seem to be necessary anymore???
// // Finally we include the very beautiful 
// #if BOOST_VERSION == 104000 
// #include </home/arossberg/include/boost/foreach.hpp>
// #else
#include <boost/foreach.hpp>
// #endif
#include <boost/version.hpp>
  
namespace boost {
    
#if BOOST_VERSION != 104900
  namespace BOOST_FOREACH = foreach;
#endif
 
} // namespace boost


// and rename it to
#define foreach BOOST_FOREACH
// and alternatively to
#define forall BOOST_FOREACH


// Here we (nearly) allways use the following implementation of a
// graph:
#include <boost/graph/adjacency_list.hpp>
typedef boost::
adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS> Graph;
typedef boost::
graph_traits<Graph>::edge_descriptor edge;
typedef boost::
graph_traits<Graph>::vertex_descriptor node;

// ... but sometimes we shall also need undirected graphs ...
typedef boost::
adjacency_list<boost::listS, boost::vecS, boost::undirectedS> ugraph;
typedef boost::
graph_traits<ugraph>::edge_descriptor uedge;
typedef boost::
graph_traits<ugraph>::vertex_descriptor unode;

// ... perhaps made from directed graphs ...
ugraph to_ugraph(const Graph& G){
  ugraph ug(num_vertices(G));
  edge e;
  forall(e,edges(G)){
    add_edge(unode(source(e,G)),unode(target(e,G)),ug);
  }
  return ug;
}


inline void NODE_REPORT(const node & centernode,const Graph & G){
#if 0
  edge e;
  REPORT(centernode);
  forall(e,in_edges(centernode,G)){
    REPORT(source(e,G));
  }
  forall(e,out_edges(centernode,G)){
    REPORT(target(e,G));
  }
#endif
}


// #ifndef REVERSE_GRAPH_DWA092300_H_
// asdf;
// #endif
// #ifdef BOOST_SUBGRAPH_HPP
// asdfads;
// #endif
// #ifdef BOOST_FILTERED_GRAPH_HPP
// asdfkdk;
// #endif
// #ifdef BOOST_VECTOR_AS_GRAPH_HPP
// kkcj;
// #endif


#if defined(ON_SX5FSV) || defined(SX)
#include <ieeefp.h>
//#define isinf(X) (((X)+1)==(X)) //quick and dirty
inline bool isinf(double f){
  int cl=fpclass(f);
  return cl==FP_NINF || cl==FP_PINF;
}
#endif

//#define forall_set_elements(E,SET) for((E)=(SET).choose();!(SET).empty();(SET).del(E),(E)=(SET).choose())


Interaction_Matrix forgiving_tsort(Interaction_Matrix &	m){
  // Implements the "greedy" algorithm to approximate maximal acyclic
  // subgraphs. There are better algorithms, but they are more
  // complex.  For now, this might be enough.
  
  Graph G(m.size());

  permutation p(m.size());

  // Add edges one after the other, but only if graph remains acyclic:
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats){
	edge e=add_edge(i,j,G).first;

	try {
	  topological_sort(G,p.begin());
	}catch(boost::not_a_dag){
	  //graph is not cyclic
	  remove_edge(e,G);
	}	  

      }
    }
  }
    
  // G is now an approximate maximal acyclic subgraph of the food web
  topological_sort(G,p.begin());

  permutation p2(m.size());
  for(int i=0;i<m.size();i++){
    p2[p[(m.size()-1)-i]]=i;
  }

  return m.permute(p2);
}

typedef Interaction_Matrix::histogram histogram;

// These are global static for speed:
static std::simple_vector<int> passed;
static Graph G; 

//histogram accumulate_chains(const node n)__attribute__((pure));

// helper for loop_free_chain_hist:
histogram accumulate_chains(const node n){
  // histogram of chain lengths of chains starting from n
  // and not passing through passed.

  int & passed_here = passed[n];

  if(passed_here) 
    return histogram(); // drop chains with loops

  histogram sum;
  edge e;
  passed_here=true;

  foreach(e,in_edges(n, G)){
    sum+=accumulate_chains(source(e,G));
  }
  passed_here=false;

  sum.prepend(1);
  return sum;
}


histogram loop_free_chain_hist(Interaction_Matrix & m){
  G=Graph(m.size());
  
  // build the graph
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }

  std::set<node> bottom_species;
  
  // get the bottom species
  node n;
  forall(n,vertices(G)){
    if(out_degree(n,G)==0)
      bottom_species.insert(n);
  }
    
  histogram sum;
  passed=std::simple_vector<int>(num_vertices(G));


  forall(n,bottom_species){
    // the main work is done in this recursive function:
    sum+=accumulate_chains(n);
  }

  // this is probably wrong:
//   // by convention, there are no chains of length 0;
//   histogram result(sum.size()+1,0);
//   copy(sum.begin(),sum.end(),&result[1]);
//   return result; 
  return sum;
}
			

	
// This referes to loops in the graph theoretical sense:
// but the Food-Web sense seems to be different.
	
////////// LEDA manual citation:
// int STRONG_COMPONENTS(graph& G, node_array<int>& compnum) 

// STRONG_COMPONENTS takes a directed graph G(V,E) as argument and
// computes for every node v V an integer compnum[v] from [0 ...c - 1]
// where c is the number of strongly connected components of G and v
// belongs to the i-th strongly connected component iff compnum[v] =
// i. STRONG_COMPONENTS returns c.  The algorithm ([36]) has running
// time O(|V| + |E|).


double real_loop_species_fraction(Interaction_Matrix & m){
  G=Graph(m.size());
  
  // build the graph
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }

  // compute strongly connected components:
  boost::vector_property_map<int> compnum(num_vertices(G));
  strong_components(G, compnum);
  
  // compute the size of each component
  node n;
  histogram comp_size;
  forall(n,vertices(G)){
    int c=compnum[n];
    comp_size[c]++;
  }

  
  // species contained in loops are in components with more than one
  // element:
  int count=0;
  forall(n,vertices(G)){
    int c=compnum[n];
    if(comp_size[c]>1)
      count++;
  }

  return double(count)/m.size();
}
	
// helper for omnivore_fraction:

//// This definition is counterintuitive, but Williams et al Nature, 2000:
// "Omnivory [5] is the fraction of species that consume two or more
// species and have food chains of different lengths (Omniv)."

static std::simple_vector<int> unique_chain_length;
enum unique_chain_length_special_value {not_unique=-1,unset=0};

int check_unique_chain_length(const node n){
  int & passed_here=passed[n];
  // here is a design option:
  if(0){
    if(passed_here) return unique_chain_length[n]=unset;
  }else{
    if(passed_here) return unique_chain_length[n]=not_unique;
  }
  if(unique_chain_length[n]!=unset) return unique_chain_length[n];

//   // since unset==0, we do not need this:
//   if(G.outdeg(n)==0)
//     return unique_chain_length[n]=1;

  edge e;
  int unique_length=unset;
  passed_here=true;
  forall(e,out_edges(n,G)){
    int length=
      check_unique_chain_length(target(e,G));
    if(length!=unset){
      if(length==not_unique || (unique_length!=unset && length!=unique_length)){
	passed_here=false;
	return unique_chain_length[n]=not_unique;
      }
      unique_length=length;
    }
  }
  passed_here=false;
  
  return unique_chain_length[n]=1+unique_length;
}

		
double omnivore_fraction(Interaction_Matrix & m){
  G=Graph(m.size());
  
  // build the graph
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }

  passed=std::simple_vector<int>(m.size());
  unique_chain_length=std::simple_vector<int>(m.size(),unset);
  node n;
  forall(n,vertices(G)){
    check_unique_chain_length(n);
  }
  
  int count=0;
  forall(n,vertices(G)){
    if(out_degree(n,G)>1 && unique_chain_length[n]==not_unique)
      count++;
  }
  
  return double(count)/m.size();
}

enum {loopy=-1, not_loopy=1};
static std::simple_vector<int> & loopiness =unique_chain_length; 
// reuse this variable

int check_loopy(const node n){
  int & passed_here=passed[n];
  if(passed_here) return loopiness[n]=loopy;
  if(loopiness[n]!=unset) return loopiness[n];

  edge e;
  passed_here=true;
  forall(e,out_edges(n,G)){
    if(target(e,G)!=n){
      int l=
	check_loopy(target(e,G));
      if(l==loopy){
	passed_here=false;
	return loopiness[n]=loopy;
      }
    }
  }
  passed_here=false;
  
  return loopiness[n]=not_loopy;
}


double loop_species_fraction(Interaction_Matrix & m){
  G=Graph(m.size());
  
  // build the graph
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }

  passed=std::simple_vector<int>(m.size(),false);
  loopiness=std::simple_vector<int>(m.size(),unset);
  node n;
  int count=0;
  forall(n,vertices(G)){
    if(check_loopy(n)==loopy) count++;
  }
  
  return double(count)/m.size();
}

bool is_connected(Interaction_Matrix & m){
  G=Graph(m.size());
  
  // build the graph
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }
  
  // compute connected components:
  boost::vector_property_map<int> compnum(num_vertices(G));
  return 1==connected_components(to_ugraph(G), compnum);
}

Interaction_Matrix the_largest_connected_subweb(const Interaction_Matrix & m){
  Graph G(m.size()); 

  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }
  
  // compute connected components:
  boost::vector_property_map<int> compnum(num_vertices(G));
  int n_components=connected_components(to_ugraph(G), compnum);

  // if only one component, return:
  if(n_components==1) return m;

  // find the largest component:
  histogram component_size(n_components,0);
  node n;
  int largest=0;
  double max_size=0;
  forall(n,vertices(G)){
    if(++component_size[compnum[n]]>max_size){
      largest=compnum[n];
      max_size=component_size[largest];
    }
  }
  
  // reduce interaction matrix to largest component:
  std::simple_vector<bool> selection(m.size(),false);
  forall(n,vertices(G)){
    if(compnum[n]==largest)
      selection[n]=true;
  }
  
  return m.select(selection);
}

// helpers for graph_analyze:

typedef std::list< std::pair<node,node> > mini_loop_list_t;
typedef Interaction_Matrix::distribution distribution;

histogram directed_food_chains_from(node n,
				    const Graph & ASG,
				    std::simple_vector<histogram> & 
				    chains_from){
  // perhaps we already computed this
  if(chains_from[n].size()) return chains_from[n]; 

  edge e;
  histogram sum;
  forall(e,in_edges(n,ASG)){
    sum+=
      directed_food_chains_from(source(e,ASG),ASG,
				chains_from);
  }
  // add one to all chain length and one chain of length 0 to the
  // histogram:
  sum.prepend(1);
  return chains_from[n]=sum;
}

histogram directed_food_chains_to(node n,
				  const Graph & ASG,
				  std::simple_vector<histogram> & 
				  chains_to){
  // perhaps we already computed this
  if(chains_to[n].size()) return chains_to[n]; 
  
  edge e;
  histogram sum;
  forall(e,out_edges(n,ASG)){
    sum+=
      directed_food_chains_to(target(e,ASG),ASG,
			      chains_to);
  }
  // add one to all chain lengths and one chain of length 0 to the
  // histogram:
  sum.prepend(1);
  return chains_to[n]=sum;
}

CMatrix<double> table1(0,0); //make this global for convenice;

void make_table1(int maxdegree){
  table1=CMatrix<double>(maxdegree+1,maxdegree+2);
  ASSERT(table1[0][0]==0);
  for(int d=1;d<=maxdegree;d++){
    for(int t=d%2;t<=d;t+=2){
      double val=
	(1.0/(d+1))*
	(t+(d+t)/2.0*table1[d-1][abs(t-1)]+
	 (d-t)/2.0*table1[d-1][t+1] );
      table1[d][t]=val;
    }
  }
  // for debugging:
//   printf("table1:\n");
//   for(int d=0;d<=maxdegree;d++){
//     for(int t=0;t<=maxdegree;t++){
//       printf("%3g ",table1[d][t]);
//     }
//     printf("\n");
//   }
}
 
 
int rho(node v,const Graph & g){
  int i=in_degree(v,g);
  int o=out_degree(v,g);
  return (i>o ? i : o);
}
int outdegree_after(node v1, const Graph & G,node v){
  edge e;
  int count=0;
  forall(e,out_edges(v1,G)){
    if(target(e,G)!=v) count++;
  }
  return count;
}
      
int indegree_after(node v1, const Graph & G,node v){
  edge e;
  int count=0;
  forall(e,in_edges(v1,G)){
    if(source(e,G)!=v) count++;
  }
  return count;
}


double table2_val(node v2, const Graph & G,node v1){
  int dout=outdegree_after(v1,G,v2);
  int din=indegree_after(v1,G,v2);
  int d=din+dout;
  int t=abs(din-dout);
  return 0.25*d + 0.5*table1[d][t];
}


void fix_tables_at(node n1,
		   const Graph & G,
		   node n2,
		   CMatrix<double> & table2,
		   std::simple_vector<double> & sum){
  if(n1!=n2){
    int i1=n1;
    int i2=n2;
    sum[i1]-=table2[i1][i2];
    sum[i1]+=
      table2[i1][i2]=
      table2_val(n1,G,n2);
  }
}

void drop_and_fix_tables(edge e,
			 Graph & G,
			 CMatrix<double> & table2,
			 std::simple_vector<double> & sum){
  node s=source(e,G);
  node t=target(e,G);
  node n; //helper
  sum[s]-=rho(s,G);
  sum[t]-=rho(t,G);
#ifdef DEBUGGING
  forall(n,vertices(G)){
    if(n!=s) ASSERT(fabs(table2[n][s]-table2_val(n,G,s))<0.001);
    if(n!=t) ASSERT(fabs(table2[n][t]-table2_val(n,G,t))<0.001);
  }
#endif
  remove_edge(e,G);
  sum[s]+=rho(s,G);
  sum[t]+=rho(t,G);

  forall(n,vertices(G)){
    fix_tables_at(n,G,s,table2,sum);
    fix_tables_at(n,G,t,table2,sum);
  }

}


void collect_and_fix_tables(edge e,
			    Graph & G,
			    Graph & ASG,
			    CMatrix<double> & table2,
			    std::simple_vector<double> & sum){
  add_edge(source(e,G),target(e,G),ASG);
  drop_and_fix_tables(e,G,table2,sum);
}
	       
int loop_count(Interaction_Matrix & m){
  int count=0;
  for(int i=0;i<m.size();i++){
    for(int j=0;j<i;j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats){
	count++;
      }
    }
  }
  return count;
}

void check_consistency(const Graph & G){
  WARNING("check_consistency(Graph & G) is disabled");
  return;
//   LEDA_list<edge> el=G.all_edges();
//   node n;
//   forall_nodes(n,G){
//     edge e;
//     forall_in_edges(e,n)
//       ASSERT(el.search(e)!=nil);
//     forall_out_edges(e,n)
//       ASSERT(el.search(e)!=nil);
//   }
}

class random_const {
  const unsigned int _value;
public:
  random_const(unsigned int dummy=0) : _value(rand()){};
  operator const unsigned int &() const{
    return _value;
  }
  //be immune to assignent:
  random_const(const random_const & other) : _value(rand()){};
  random_const& operator=(const random_const&){return *this;};
};
  

unsigned int graph_hash(const Graph & ASG){
  static sequence<random_const> id;
  id.resize(num_vertices(ASG)); //make sure id values won't change
  unsigned int h=0;
  const unsigned int top_bit=(unsigned(1)<<(8*sizeof(h)-1));
  for(unsigned int i=0;i<id.size();i++){
    node n;
    forall(n,adjacent_vertices(i,ASG)){
      h^=id[n];
    }
    h+=i;
    if(h | top_bit)
      h=(h<<1)|1;
    else
      h<<=1;
  }
  return h;
}

// void mark_followers(const node n, const graph & ASG,
// 		    follows_t & follows){
//   if(passed[n]) return;
  
//   node m;
//   forall_adj_nodes(m,n){
//     mark_followers(m,ASG,follows);
//     node i;
//     forall_nodes(i,ASG){
//       follows(n,i)|=follows(m,i);
// //       std::cout << follows(n,i);
//     }
// //     std::cout << std::endl;
//     follows(n,m)=true;
//   }
//   passed[n]=true;
// }
  

// follows_t transitive_closure(const GRAPH<int,int> & ASG){

//   follows_t follows(ASG,false);
//   passed=node_array<bool>(ASG,false);
  
//   node n;
//   forall_nodes(n,ASG){
//     if(ASG.indeg(n)==0){
//       mark_followers(n,ASG,follows);
//     }
//   }
// //   forall_nodes(n,ASG){
// //     node m;
// //     forall_nodes(m,ASG){
// //       std::cout << follows(n,m) ;
// //     }
// //     std::cout << " " << ASG.inf(n) << std::endl;
// //   }
//   return follows;
// }


typedef boost::adjacency_matrix<> follows_t;

follows_t transitive_closure(const Graph & ASG){
  follows_t follows(0);
  boost::transitive_closure(ASG,follows);
  return follows;
}

typedef std::list<int> left_multiplicities_list_t;

graphs_stats_data graph_analyze(Interaction_Matrix & m, 
				const bool go_through_multiplicities){
  //m=m.randomize(); // turn this on to get many multiplicities
  //m.Print();
  average_meter oChnLg, oChnSD, oChnNo, oLoop, oOmniv;
  graphs_stats_data stats;

  left_multiplicities_list_t multiplicities_left;
  std::set<unsigned int> hashs_of_graphs_found_so_far;
  unsigned int Ghash=0; // to verify if we allways start with the same G

  // This loops goes through all possible outcomes of the Berger-Shor
  // Algorithm in order to average over them.  In most cases, however,
  // the outcome is unique.
  do{
    unsigned int multiplicity_index=0;

    Graph G(m.size()),ASG(m.size()); 
    std::vector<bool> vertex_alife(m.size(),true);
    int vertices_alife=m.size();
  
    mini_loop_list_t mini_loops;
  
    for(int i=0;i<m.size();i++){
      for(int j=0;j<m.size();j++){
	if((int)m[i][j]==(int)NetworkAnalysis::eats){
	  if((int)m[j][i]==(int)NetworkAnalysis::eats){
	    if(i<=j){
	      // 2-loops and 1-loops:
	      mini_loops.push_back(std::pair<node,node>(i,j));
	    }
	  }else{
	    add_edge(i,j,G);
	  }
	}
      }
    }

    // check consistency:
    if(multiplicities_left.size()==0)
      Ghash=graph_hash(G);
    else
      ALWAYS_ASSERT(Ghash==graph_hash(G));


    // find maximal degree:
    int maxdegree=0;
    node n;
    forall(n,vertices(G)){
      // dead vertices have degree zero, so they will not interfere here
      int degree=out_degree(n,G)+in_degree(n,G);
      if(degree>maxdegree) maxdegree=degree;
    }

    make_table1(maxdegree);

    CMatrix<double> table2(m.size(),m.size());
    // prepare table2:
    for(int i=0;i<m.size();i++){
      for(int j=0;j<m.size();j++){
	table2[i][j]=table2_val(i,G,j);
      }
    }

    //std::simple_vector<double> sum(m.size());
    sequence<double> sum(m.size());
    // prepare array sum:
    for(int i=0;i<m.size();i++){
      double s=0;
      s+=rho(i,G);
      for(int j=0;j<m.size();j++){
	if(j!=i)
	  s+=table2[i][j];
      }
      sum[i]=s;
    }

    //main loop:
    while(vertices_alife){
      //     node n1,n2;
      //     forall_nodes(n1,G){
      //       forall_nodes(n2,G){
      // 	printf("%6g ",table2[G.inf(n1)][G.inf(n2)]);
      //       }
      //       printf("\n");
      //     }
      
      // separate strong components
      boost::vector_property_map<int> compnum(num_vertices(G));
      int number_of_components=
	strong_components(G, compnum);

      if(number_of_components > 1){
	// we go a bit humpty dumpty here, because removal of an edge
	// of G, as it happens in collect_and_fix_tables( ),
	// invalidates the edge range edges(G) in boost.  There should
	// definitely be a better way to do this!
	edge e;
	typedef std::pair<node,node> simple_edge;
	typedef std::list< simple_edge > edge_list_t;
	edge_list_t to_remove;
	forall(e,edges(G)){
	  if(compnum[source(e,G)]!=compnum[target(e,G)]){
	    to_remove.push_back(simple_edge(source(e,G),target(e,G)));
	  }
	}
	forall(simple_edge se,to_remove){
	  collect_and_fix_tables(boost::edge(se.first,se.second,G).first,
				 G,ASG,table2,sum);
	}
      }
    
      // in order to canonicalize multiplicity iteration, find strong
      // component with node with smallest index:
      node n;
      int minindex=m.size(); // index of node with smallest index
      int min_compnum=0; // component containing this node
      forall(n,vertices(G)){
	if(vertex_alife[n]){
	  if(n < minindex){
	    minindex=n;
	    min_compnum=compnum[n];
	  }
	}
      }

      //find best nodes to process next (usually only one)
      std::list<node> best_nodes;
      double maxE=-1;
      //half max difference between "equal" values:
      const double delta=1e-10;
      forall(n,vertices(G)){
	if(vertex_alife[n]){
	  if(compnum[n]==min_compnum){
	    if(sum[n]>maxE-delta){
	      if(sum[n]>maxE+delta){
		maxE=sum[n];
		best_nodes.clear();
		best_nodes.push_back(n);
	      }else{
		best_nodes.push_back(n);
	      }
	    }
	  }
	}
      }

      //REPORT(sum);
      sequence<int> is_in_min_comp(m.size(),false);
      forall(n,vertices(G)){
	is_in_min_comp[n]=(compnum[n]==min_compnum);
      }
      //REPORT(is_in_min_comp);
      
      
      
      node max_n=*(best_nodes.begin()); //nominally unique maximal node
      ASSERT(best_nodes.begin()!=best_nodes.end());
      unsigned int bns=best_nodes.size();
      if(bns>1 && go_through_multiplicities){
// 	REPORT(bns);
	//////we detected a multiplicity
	int node_to_pick; // this we have to decide now
	const unsigned int multiplicity_randomization_threshold=0;
	unsigned int mls=multiplicities_left.size();
	if(mls==multiplicity_index){
	  //we have to start iterating over this multiplicity
	  if(multiplicity_index>=multiplicity_randomization_threshold){
	    WARNING("multiplicity_index>=multiplicity_randomization_threshold");
	    multiplicities_left.push_back(0);
	    node_to_pick=random_integer(bns);
	  }else{
	    std::cout << "mult " << multiplicity_index 
		      << " start " << bns 
		      << std::endl;
	    multiplicities_left.push_back(bns-1);
	  }
	}
	////decide which node to pick:
	if(multiplicity_index>=multiplicity_randomization_threshold){
	  node_to_pick=random_integer(bns);
	}else{
	  left_multiplicities_list_t::iterator mult=
	    multiplicities_left.begin();
	  //find the multipliticy we are working on 
	  for(int c=multiplicity_index;c>0;c--) mult++;
	  node_to_pick=(*mult);
	}
	//find the node
	std::list<node>::iterator best_node=best_nodes.begin();
	for(int c=node_to_pick;c>0;c--){
	  best_node++;
	  if(best_node==best_nodes.end()){
	    WARNING("multiplicity scanning scheme is garbled!!");
	    best_node--;
	  }
	}
	max_n=*best_node;
	
	//// do we need to step?
	if(mls==multiplicity_index+1){
	  //we have to step iteration over this multiplicity
	  left_multiplicities_list_t::reverse_iterator tail=
	    multiplicities_left.rbegin();
	  (*tail)--;
	  std::cout << "mult " << multiplicity_index 
		    << " step " << *tail 
		    << std::endl;
	}
	multiplicity_index++;
      }
      

    
      edge e;
      if(in_degree(max_n,G) >= out_degree(max_n,G)){
	// since we manipulate the edge list, this is a bit tricky:
	while(in_degree(max_n,G)){
	  e=*(in_edges(max_n,G).first);
	  NODE_REPORT(max_n,G);
	  collect_and_fix_tables(e,G,ASG,table2,sum);
	}
	while(out_degree(max_n,G)){
	  e=*(out_edges(max_n,G).first);
	  NODE_REPORT(max_n,G);
	  drop_and_fix_tables(e,G,table2,sum);
	}
      }else{
	// since we manipulate the edge list, this is a bit tricky:
	while(out_degree(max_n,G)){
	  e=*(out_edges(max_n,G).first);
	  NODE_REPORT(max_n,G);
	  collect_and_fix_tables(e,G,ASG,table2,sum);
	}
	while(in_degree(max_n,G)){
	  e=*(in_edges(max_n,G).first);
	  NODE_REPORT(max_n,G);
	  drop_and_fix_tables(e,G,table2,sum);
	}
      }
      ASSERT(in_degree(max_n,G)==0);
      ASSERT(out_degree(max_n,G)==0);

      //remove_vertex(max_n,G);
      vertex_alife[max_n]=false;
      vertices_alife--;
    }// while(vertices_alife)

    unsigned int h=graph_hash(ASG);
    if(hashs_of_graphs_found_so_far.find(h)!=
       hashs_of_graphs_found_so_far.end() ){
      std::cout << "multiplicity, but not new result" << std::endl;
    }else{
      std::vector<int>  ord(num_vertices(ASG)); // use later
      // we refrain from checking if graphs are equal when hashs are equal.
      hashs_of_graphs_found_so_far.insert(h);

      // ASG is now an approximately maximal acyclic subgraph of the
      // food web still excluding two-loops

      topological_sort(ASG,ord.rbegin());
//       ASG.sort_nodes(ord);
//       ASG.print();

      // get the ordering so far determined:
      follows_t follows=transitive_closure(ASG);

      // insert decided two-loops and list undecied two-loops:
      mini_loop_list_t undecied_two_loops;
      for(mini_loop_list_t::iterator l=mini_loops.begin();
	  l!=mini_loops.end();
	  l++  ){
	if(l->first != l->second){ // only two-loops
// 	  std::cout << "<" << ASG.inf(l->first) << "," 
// 		    << ASG.inf(l->second) << "> " ;
	  if(ord[l->second]<ord[l->first]){
	    // reorder:
// 	    std::cout << "reordered to ";
	    *l=std::pair<node,node>(l->second,l->first);
// 	    std::cout << "<" << ASG.inf(l->first) << "," 
// 		      << ASG.inf(l->second) << "> " ;
	  }
	  if(boost::edge(l->first,l->second,follows).second){
	    // the two-loop is decided, insert it
// 	    std::cout << "inserted" ;
// 	      ASG.new_edge(l->first,l->second);
	  }else{
// 	    std::cout << "kept" ;
	    // the two-loop is undecided:
	    undecied_two_loops.push_back(*l);
	  }
// 	  std::cout << std::endl;
	}
      }

//       REPORT(TOPSORT(ASG,ord));


      // sample undecided two-loop ordering at random?
      int n_undecided=
	(go_through_multiplicities ? undecied_two_loops.size() : 0 );
      const int randomization_threshold=8;
      WARN_IF(n_undecided>=randomization_threshold,"randomizing");
      const bool randomizing=(n_undecided>=randomization_threshold);
      if(randomizing) ;REPORT(n_undecided);
      
      //loop over all/many orderings of two-loops:
      typedef long long int orderbits;
      int n_acyclic=0;
      for(orderbits ordering_index=0;
	  ordering_index < 
	    (orderbits(1)<<(randomizing?randomization_threshold:n_undecided))
	    ;
	  ordering_index++){
	orderbits ordering=
	  (randomizing?
	   orderbits(random_double(orderbits(1)<<n_undecided)) :
	   ordering_index);
	
	//this ordering always goes through, so we have at least one,
	//even when randomizing:
	if(ordering_index==0) ordering=0;
	
	//insert edges according to ordering:
	std::list<edge> inserted_twoloop_edges;
	int position=0;
// 	REPORT(ordering);
	for(mini_loop_list_t::iterator i=
	      undecied_two_loops.begin();
	    i!=undecied_two_loops.end();
	    i++,position++ ){
// 	  REPORT(TOPSORT(ASG,ord));
// 	  ASG.print();
// 	  std::cout << position << "<" << ASG.inf(i->first) << "," 
// 		    << ASG.inf(i->second) << "> inserted " ;
	  if( (orderbits(1)<<position) & ordering ){
// 	    std::cout << "reverse";
	    inserted_twoloop_edges.
	      push_back(add_edge(i->second,i->first,ASG).first);
	  }else{
	    inserted_twoloop_edges.
	      push_back(add_edge(i->first,i->second,ASG).first);
	  }
// 	  std::cout << std::endl;
// 	  REPORT(TOPSORT(ASG,ord));
// 	  ASG.print();
	}

	//check if result is acyclic:
	bool is_acyclic=true;
	try{
	  topological_sort(ASG,ord.rbegin());
	}catch(boost::not_a_dag){
	  is_acyclic=false;
	}
	if(is_acyclic){
	  n_acyclic++;
	  // is acyclic. evaluate:
	  //	  ASG.sort_nodes(ord);
      
	  // compute permutation of indices:
	  permutation p(m.size());
	  forall(n,vertices(ASG)){
	    p[ord[n]]=n;
	  }
    
	  Interaction_Matrix nm=m.permute(p);
	  stats.sorted=nm;
// 	  nm.Print();
// 	  FATAL_ERROR("check_this");

	  std::simple_vector<histogram> 
	    chains_from(num_vertices(ASG),histogram());
	  histogram h; // chain length distribution
	  bool found=false;
	  forall(n,vertices(ASG)){
	    if(out_degree(n,ASG)==0){
	      found=true;
	      h+=directed_food_chains_from(n,ASG,chains_from);
	    }
	  };
	  if(!found){
	    FATAL_ERROR("could not find a lowest trophic level");
	  }
	  weighted_average_meter length;
	  long double sum=0;
	  for(unsigned int i=0;i<h.size();i++){
	    if(h[i]<0) FATAL_ERROR("histogram overflow!");
	    length.sample(i,h[i]);
	    sum+=h[i];
	  }
	  oChnLg.sample(length.readout());
	  if(!finite(length.sample_std())){
	    std::cout << h << std::endl;
	    WARNING("overflow in average meter??");
	  }
	  else
	    oChnSD.sample(length.sample_std());
	  if(!finite(log10(sum)))
	    WARNING("overflow computing sum??");
	  else
	    oChnNo.sample(log10(sum));
	  oLoop.sample(loop_count(nm));
      
	  average_meter std;
	  std::simple_vector<histogram> 
	    chains_to(num_vertices(ASG),histogram());
	  forall(n,vertices(ASG)){
	    histogram h=
	      directed_food_chains_to(n,ASG,chains_to);
	    weighted_average_meter l;
	    for(unsigned int i=0;i<h.size();i++){
	      l.sample(i,h[i]);
	    }
	    double lstd=l.sample_std();
	    if(my_isnan(lstd) || my_isinf(lstd))
	      WARNING("disconnected node encountered when computing oOmniv");
	    else
	      std.sample(lstd);
	  }
	  double stdm=std.readout();
	  if(my_isnan(stdm) || my_isinf(stdm))
	    FATAL_ERROR("could not compute oOmniv");
	  else
	    oOmniv.sample(stdm);
	} // if(TOPSORT(ASG))
	
	//remove undecided two-loop branches:
	for(std::list<edge>::iterator i=
	      inserted_twoloop_edges.begin();
	    i!=inserted_twoloop_edges.end();
	    i++){
	  remove_edge(*i,ASG);
	}
      } // loop over orderings of two-loops
      if(n_acyclic > 1){
	std::cout << n_acyclic << " graphs of " << 
	  (orderbits(1)<<(randomizing?randomization_threshold:n_undecided))
		  << " acyclic " << std::endl;
      }
    }// condition that graph is unique so far
    
    // drop finished multiplicities...
    //ASSERT(multiplicities_left.size()==multiplicity_index);
    if(multiplicities_left.size()!=multiplicity_index){
      WARN_IF(multiplicities_left.size()!=multiplicity_index,
	      "HIT THIS BUG AGAIN!");
      REPORT(multiplicities_left.size());
      REPORT(multiplicity_index);
      std::cout << "multiplicities_left:";
      while(!multiplicities_left.empty()){
	std::cout << *multiplicities_left.begin() << " " ;
	multiplicities_left.pop_front();
      }
      std::cout << std::endl;
      break;
    }
    while(!multiplicities_left.empty() &&
	  *(multiplicities_left.rbegin())<=0) 
      multiplicities_left.pop_back();
    // ... and continue if some are left
  }while(!multiplicities_left.empty());
  
  stats.oChnLg=oChnLg.readout();
  stats.oChnSD=oChnSD.readout();
  stats.oChnNo=oChnNo.readout();
  stats.oLoop=oLoop.readout();
  stats.oOmniv=oOmniv.readout();

  return stats;
}

void find_lowest_level(Interaction_Matrix & m,std::simple_vector<bool> & lowest){
  FATAL_ERROR("find_lowest_level works only with debugging options set!");
  //(LEDA bug??)
  Graph G(m.size()); 
  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }
  
  // compute strongly connected components:
  boost::vector_property_map<int> compnum(num_vertices(G));
  int number_of_components=
    strong_components(G, compnum);
  
  edge e;
  std::simple_vector<bool> lowestcomp(num_vertices(G),true);
  forall(e,edges(G)){
    if(compnum[source(e,G)]!=compnum[target(e,G)]){
      lowestcomp[source(e,G)]=false;
    }
  }
  
  for(int i=0;i<m.size();i++){
    REPORT(i);
    REPORT(lowestcomp[i]);
    lowest[i]=lowestcomp[i];
  }
  return;
}

distribution trophic_height_vector(const Interaction_Matrix & m){
  //!!! we assume G to be acyclic!!
  Graph G(m.size()); 

  for(int i=0;i<m.size();i++){
    for(int j=0;j<m.size();j++){
      if((int)m[i][j]==(int)NetworkAnalysis::eats)
	add_edge(i,j,G);
    }
  }
  
  distribution heights(m.size());
  for(int i=0;i<m.size();i++){
    std::simple_vector<histogram> 
      chains_to(num_vertices(G),histogram());
    histogram h=
      directed_food_chains_to(i,G,chains_to);
    double lsum=0,nsum=0;
    for(int j=h.size();j-->0;){
      lsum+=j*h[j];
      nsum+=h[j];
    }
    heights[i]=lsum/nsum+1;
  }
  return heights;
}

int Cy4(const Interaction_Matrix & m){
  Graph G(m.size());  // niche overlap graph
  std::simple_vector<node>  ind(m.size());
  //CMatrix<bool> adjacent(m.size(),m.size());
  Interaction_Matrix adjacent(m.size());
  
  for(int i=0;i<m.size();i++){
    for(int j=i+1;j<m.size();j++){
      int k=m.size();
      while(k--){
	if((int)m[i][k]==(int)NetworkAnalysis::eats &&
	   (int)m[j][k]==(int)NetworkAnalysis::eats ){
	  // niche overlap detected;
	  add_edge(i,j,G);	  
	  adjacent[i][j]=NetworkAnalysis::eats;
	  adjacent[j][i]=NetworkAnalysis::eats;
	  break;
	}
      }
    }
  }
  // niche overlap graph constructed
  
//   cout << "DEBUG:" << std::endl;
//   m.Print();
//   adjacent.Print();
  
  // Chordless 4-cycles are only counted if:
  // 
  // 1. n1 has the smallest index
  // 2. the index of n2 is smaller than the index of n3
  //
  // This should be sufficient to avoid double counting.
  
  int count=0;
  node n1,n2,n3,n4;
  edge e2,e3,e4;
  forall(n1,vertices(G)){
    // try n1 as starting point of chordless cycle
    forall(e2,out_edges(n1,G)){
      n2=target(e2,G);
      forall(e3,out_edges(n1,G)){
	n3=target(e3,G);
	if(n2<n3 && !adjacent[n2][n3]){
	  forall(e4,out_edges(n2,G)){
	    n4=target(e4,G);
	    if(n4 > n1
	       &&
	       !adjacent[n1][n4]
	       &&
	       adjacent[n3][n4] ){
// 	      cout << "cycle: " 
// 		   << G.inf(n1) << " "
// 		   << G.inf(n2) << " "
// 		   << G.inf(n4) << " "
// 		   << G.inf(n3) << std::endl;
	      count++;
	    }
	  }
	  forall(e4,in_edges(n2,G)){
	    n4=source(e4,G);
	    if(n4 > n1
	       &&
	       (!adjacent[n1][n4])
	       &&
	       adjacent[n3][n4] ){
	      std::cout << "";//inhibit buggy optimizations by g++
//  	      std::cout << "cycle: " 
// 		   << G.inf(n1) << " "
// 		   << G.inf(n2) << " "
// 		   << G.inf(n4) << " "
// 		   << G.inf(n3) << std::endl;
	      count++;
	    }
	  }
	}
      }
    }
  }
  return count;
}

// Lex-BFS:

typedef std::vector<node> list_to_partition;

class set_of_nodes{
public:
  list_to_partition::iterator Begin;
  list_to_partition::iterator End;
  list_to_partition::iterator Begin_Split;
  set_of_nodes(list_to_partition::iterator b,list_to_partition::iterator e):
    Begin(b),End(e),Begin_Split(e){};
};

typedef std::list<set_of_nodes> set_list;


static std::simple_vector<int> Lex_BFS_Ordering(const ugraph & G){
  // implementation of Algorithm 2 of Habib et al. (2000),
  std::simple_vector<int> pi(num_vertices(G),-1);
  std::simple_vector<list_to_partition::iterator> position_of(num_vertices(G));
  std::simple_vector<set_list::iterator> set_of(num_vertices(G));
  
  list_to_partition l(num_vertices(G));
  set_list L;
  L.push_back(set_of_nodes(l.begin(),l.end()));
//   REPORT(L.begin()->Begin_Split-l.begin());
  std::stack<set_list::iterator> might_need_fixing;

  node n;
  list_to_partition::iterator i=l.begin();
  forall(n,vertices(G)){
    position_of[n]=i;
    *i=n;
    i++;
    set_of[n]=L.begin();
  }
  
  int index=num_vertices(G);
  set_list::iterator singletons_start=L.end();
  while(singletons_start!=L.begin()){
    set_list::iterator last_non_visited_class=singletons_start;
    last_non_visited_class--;
    last_non_visited_class->End--;
    last_non_visited_class->Begin_Split--;
    node x=*(last_non_visited_class->End);
    if(last_non_visited_class->End==last_non_visited_class->Begin){
      singletons_start--;
    }
    pi[x]=--index;
    node y;
//     REPORT((position_of[x]-l.begin()));
    forall(y,adjacent_vertices(x,G)){
      if(position_of[y]<position_of[x]){
	set_list::iterator chi_b=set_of[y];
	might_need_fixing.push(chi_b);
// 	REPORT((chi_b->Begin-l.begin()));
// 	REPORT((chi_b->Begin_Split-l.begin()));
// 	REPORT((position_of[y]-l.begin()));
// 	ASSERT(position_of[y]>=chi_b->Begin);
// 	ASSERT(position_of[y]<chi_b->End);
// 	ASSERT(position_of[y]<chi_b->Begin_Split);
	--(chi_b->Begin_Split);
// 	REPORT((chi_b->Begin_Split-l.begin()));
	position_of[*(chi_b->Begin_Split)]=position_of[y];
	iter_swap(position_of[y],chi_b->Begin_Split);
	// fix position_of:
	position_of[y]=(chi_b->Begin_Split);
#if xDEBUGGING
	node z;
	forall_nodes(z,G){
	std::cout << position_of[z]-l.begin() << " ";
	}
	std::cout << std::endl;
#endif 
      }
    }

    while(!might_need_fixing.empty()){
      set_list::iterator & chi_b=might_need_fixing.top();
      if(chi_b->Begin_Split!=chi_b->End){
	if(chi_b->Begin_Split==chi_b->Begin){
	  chi_b->Begin_Split=chi_b->End;
	}else{
	  set_list::iterator chi_new=chi_b;
	  chi_new++;
#if DEBUGGING
	  set_list::iterator chi_new_hold=chi_new;
#endif
	  chi_new=L.insert(chi_new,set_of_nodes(chi_b->Begin_Split,chi_b->End));
	  ASSERT(chi_new_hold!=chi_new);
	  chi_b->End=chi_b->Begin_Split;
	  // fix set_of:
	  for(list_to_partition::iterator y=chi_new->Begin;
	      y!=chi_new->End;
	      y++){
	    set_of[*y]=chi_new;
	  }
	}
      }
      might_need_fixing.pop();
    }
  }
  ASSERT(index==0);
  return pi;
}

// helper for is_chordal below

static bool 
is_chordal_check_order(const std::vector< std::list <int> > & children_of,
		       const std::vector< int > & parent,
		       const std::vector< std::list <int> > & RN,
		       int x ,
		       int &n_cliques,
		       std::vector< bool > & makes_maxial_clique
		       ){
  if(x<parent.size()){
  
    std::list<int>::const_iterator j=RN[parent[x]].begin();
    const std::list<int>::const_iterator end_j=RN[parent[x]].end();
    // for all right neighbours *i of x, excluding its parent...
    
//     REPORT(x);
//     REPORT(*(RN[x].begin()));
    
    bool sets_identical=true;
    for(std::list<int>::const_iterator i=++(RN[x].begin());
	i!=RN[x].end();
	i++){
      
//       REPORT(*i);
      // Here is the Algorithm 3 part:
      if(j==end_j) return false;
      // ... search for *i in RN[parent[x]]:
      while(*i!=*(j++)){
	if(j==end_j){
	  // not found!
	  return false;
	}else{
	  sets_identical=false;
	}
      }
    }
      
    // Here is the Algorithm 4 part:
    if(sets_identical){
      if(parent[x]==parent.size()){
	n_cliques++;
      }else{
	if(makes_maxial_clique[parent[x]]){
	  makes_maxial_clique[parent[x]]=false;
	}else{
	  n_cliques++;
	}
      }
    }else{
      n_cliques++;
    }
    makes_maxial_clique[x]=true;
  }
  
  // the recursion:
  for(std::list< int >::const_iterator i=children_of[x].begin();
      i!=children_of[x].end();
      i++){
    if(!is_chordal_check_order(children_of,parent,RN,*i,
			       n_cliques,makes_maxial_clique ))
      return false;
  }

  return true;
}


// these three types are used for by Algorithm 9:
typedef std::vector<int> cliques_to_partition;

class set_of_cliques{
public:
  cliques_to_partition::iterator Begin;
  cliques_to_partition::iterator End;
  cliques_to_partition::iterator Begin_Split;
  set_of_cliques(cliques_to_partition::iterator b,
		 cliques_to_partition::iterator e):
    Begin(b),End(e),Begin_Split(e){};
};

typedef std::list<set_of_cliques> clique_set_list;



static int is_chordal(const ugraph & G,const std::simple_vector<int> & index){

  // implementation of Algorithms 3, 4, and 9 of Habib et al. (2000):
  const int S=num_vertices(G);
  node x,y;
  
  // generate and sort right neighbour (RN) sets:
  std::vector< std::stack< int > > arn(S);
  forall(x,vertices(G)){
    forall(y,adjacent_vertices(x,G)){
      if(index[y] > index[x]){
	arn[index[y]].push(index[x]);
      }
    }
  }
  
  std::vector< std::list< int > > RN(S+1);
  for(int i=0;i<S;i++){
    std::stack<int> & s=arn[i];
    while(!s.empty()){
      RN[s.top()].push_back(i);
      s.pop();
    }
  }
  
  // generate tree structure:
  std::vector< std::list <int> > children_of(S+1);
  std::vector< int > parent(S);
  for(int i=0;i<S;i++){
    // make sure every node has a right neighbour:
    RN[i].push_back(S);
    parent[i]=RN[i].front();
    children_of[parent[i]].push_back(i);
  }
  
  std::vector< bool > makes_maxial_clique(S,false);
  int n_cliques=0;
  
  if(!is_chordal_check_order(children_of,parent,RN,S,
			     n_cliques,makes_maxial_clique)){
    return 0;
  }

  CMatrix<NetworkAnalysis::Interaction> cm(S,n_cliques);

//   REPORT(n_cliques);

  int c=0;
  for(int i=S;i-->0;){
//     REPORT(i);
//     REPORT(makes_maxial_clique[i]);
    if(makes_maxial_clique[i]){
      for(std::list<int>::iterator m=RN[i].begin();
	  m!=RN[i].end();
	  m++){
	if(*m!=S){
	  cm[*m][c]=NetworkAnalysis::eats;
	}
      }
      cm[i][c]=NetworkAnalysis::eats;
      c++;
      ASSERT(c<=n_cliques);
    }
  }
  
  std::cout  << "CLIQUES START" << std::endl;
  for(int j=0;j<n_cliques;j++){
    for(int i=0;i<S;i++){
      std::cout << (cm[i][j]==NetworkAnalysis::eats ? "X" : " ");
    }
    std::cout << std::endl;
  }
  std::cout  << "CLIQUES END" << std::endl;
  
  if(has_consecutive_ones(cm))
    return 2;
  else
    return 1;
  
//   // Here starts the Algorithm 9 part:

//   cliques_to_partition l(n_cliques);
//   clique_set_list L;
//   L.push_back(set_of_cliques(l.begin(),l.end()));
//   std::stack<clique_set_list::iterator> might_need_fixing;
//   set_list::iterator singletons_start=L.end();

//   std::stack< int > pivots;
//   std::vector< bool > pivot_processed(n_cliques,false);
  
//   while(true){
//     set_list::iterator chi_c=singletons_start;
//     chi_c--;
//     while(chi_c is singleton && chi_c!=L.begin()){
//       singletons_start=chi_c--;
//     }
//     if(chi_c is singleton)
//       break;

//     if(pivots.empty()){
//       cliques_to_partition::iterator last_clique=chi_c.Begin;
//       for(cliques_to_partition::iterator c=
// 	    chi_c.Begin;
// 	  c!=chi_c.End;
// 	  chi_c++){
// 	if(*c>*last_clique){
// 	  last_clique=c;
// 	}
//       }
      
//       ... separate c as a new partition ...;

//     }else{
//       while(pivot_processed[pivots.top()]){
// 	pivots.pop();
//       }
//       int x=pivots.top()
//       pivot_processed[pivots.top()]=true;
//       pivots.pop();

}

int niche_overlap_graph_analysis(const Interaction_Matrix & m){
  ugraph G; 
  
  // find niche overlap graph:
  for(int i=0;i<m.size();i++){
    for(int j=0;j<i;j++){
      for(int k=0;k<m.size();k++){
	if(m.eats(i,k) && m.eats(j,k)){
	  // there is niche overlap:
	  add_edge(i,j,G);
	  break;
	}
      }
    }
  }
  
  std::simple_vector<int> index=Lex_BFS_Ordering(num_vertices(G));

  return is_chordal(G,index);
}
  


#undef forall_set_elements
 
