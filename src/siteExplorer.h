#include "util.h"
#include <set>

class MutationExplorer: public TreeAlignmentMatch {
public:
  MutationExplorer(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqs
  );
  std::set<int> const getDivPoints();
  std::vector< std::set<std::string> > const getMutationList();
private:
  bool rootDivTrue;
  std::set<int> divPoints;
  ListOf<CharacterVector> ancestralSeqs;
  std::vector< std::deque<int> > evolPath;
};

class SiteExplorer: public TreeAlignmentMatch {
public:
  SiteExplorer(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqs
  );
  std::set<int> const getDivPoints();
  std::vector< std::deque<int> > const getPath();
  std::map< std::string, std::set<std::string> > const getSitePath(const int mode);
private:
  bool rootDivTrue;
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
