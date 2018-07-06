#include "util.h"
#include <set>

class SiteExplorer: public TreeAlignmentMatch {
  typedef std::vector< std::deque<int> >::iterator pathIter;
  typedef std::map< int, std::set< std::pair<int, int> > >::iterator modeIter;
public:
  SiteExplorer(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqs
  );
  std::map< std::string, std::set<std::string> > getSitePath(int mode);
private:
  ListOf<CharacterVector> ancestralSeqs;
  std::vector< std::deque<int> > evolPath;
  std::set<int> divPoints;
  std::map< int, std::set<std::string> > allele;
  std::map< int, std::set< std::pair<int, int> > > mutNodes;
  std::map< std::set< std::pair<int, int> >, std::vector<int> > coEvol;
  void fixedMutationFilter(
      std::vector<int> &linked,
      std::map< std::string, std::set<std::string> > &linkages
  );
  void coSiteFilter(
      std::vector<int> &linked,
      std::map< std::string, std::set<std::string> > &linkages
  );
  void nonFixedFilter(
      std::vector<int> &linked,
      std::map< std::string, std::set<std::string> > &linkages
  );
  void addToLinkage(
      std::vector<int> &linked,
      std::map< std::string, std::set<std::string> > &linkages,
      const int &node1,
      const int &node2,
      const int &siteIndex
  );
  void getEvolPath();
};
