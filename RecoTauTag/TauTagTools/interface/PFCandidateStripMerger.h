
/** PFCandidateStripMerger

Class that creates strips from Particle Flow Candidates
And outputs a Collection of Candidate Lists
Michail Bachtis
----------------------
University of Wisconsin
bachtis@cern.ch
**/

#include "RecoTauTag/TauTagTools/interface/PFCandidateMergerBase.h"
#include "RecoTauTag/TauTagTools/interface/TauTagTools.h"

class PFCandidateStripMerger : public PFCandidateMergerBase
{
 public:
  PFCandidateStripMerger();
  PFCandidateStripMerger(const edm::ParameterSet&);
  ~PFCandidateStripMerger() override;

  std::vector<std::vector<reco::CandidatePtr>> mergeCandidates(const std::vector<reco::CandidatePtr>&) override;

 private:

  std::vector<int> inputPdgIds_; //type of candidates to clusterize
  double etaAssociationDistance_;//eta Clustering Association Distance
  double phiAssociationDistance_;//phi Clustering Association Distance
  

  //Private Methods
  bool candidateMatches(const reco::CandidatePtr&);


};
