#include "util.h"
#include <set>

class SiteExplorer: public TreeAlignmentMatch {
  typedef vector< deque<int> >::iterator pathIter;
  typedef map< int, set< pair<int, int> > >::iterator modeIter;
public:
  SiteExplorer(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqs
  );
  map< string, set<string> > getSitePath(int mode);
private:
  ListOf<CharacterVector> ancestralSeqs;
  vector< deque<int> > evolPath;
  set<int> divPoints;
  map< int, set<string> > allele;
  map< int, set< pair<int, int> > > mutNodes;
  map< set< pair<int, int> >, vector<int> > coEvol;
  void fixedMutationFilter(
      vector<int> &linked,
      map< string, set<string> > &linkages
  );
  void coSiteFilter(
      vector<int> &linked,
      map< string, set<string> > &linkages
  );
  void nonFixedFilter(
      vector<int> &linked,
      map< string, set<string> > &linkages
  );
  void addToLinkage(
      vector<int> &linked,
      map< string, set<string> > &linkages,
      const int &node1,
      const int &node2,
      const int &siteIndex
  );
  void getEvolPath();
};
