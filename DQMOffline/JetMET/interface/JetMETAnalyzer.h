#ifndef JetMETAnalyzer_H
#define JetMETAnalyzer_H


/** \class JetMETAnalyzer
 *
 *  DQM jetMET analysis monitoring
 *
 *  $Date: 2010/02/24 19:08:53 $
 *  $Revision: 1.15 $
 *  \author F. Chlebana - Fermilab
 *          K. Hatakeyama - Rockefeller University
 */


#include <memory>
#include <fstream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DQMOffline/JetMET/interface/JetAnalyzer.h"
#include "DQMOffline/JetMET/interface/JetPtAnalyzer.h"
#include "DQMOffline/JetMET/interface/PFJetAnalyzer.h"
#include "DQMOffline/JetMET/interface/JPTJetAnalyzer.h"
#include "DQMOffline/JetMET/interface/CaloMETAnalyzer.h"
#include "DQMOffline/JetMET/interface/METAnalyzer.h"
#include "DQMOffline/JetMET/interface/PFMETAnalyzer.h"
#include "DQMOffline/JetMET/interface/HTMHTAnalyzer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class JetMETAnalyzer : public edm::EDAnalyzer {
 public:

  /// Constructor
  JetMETAnalyzer(const edm::ParameterSet&);
  
  /// Destructor
  virtual ~JetMETAnalyzer();
  
  /// Inizialize parameters for histo binning
  void beginJob(void);

  /// Get the analysis
  void analyze(const edm::Event&, const edm::EventSetup&);

  /// Save the histos
  void endJob(void);

  /// Initialize run-based parameters
  void beginRun(const edm::Run&,  const edm::EventSetup&);

  /// Finish up a run
  void endRun(const edm::Run&,  const edm::EventSetup&);

 private:
  // ----------member data ---------------------------
  
  DQMStore* dbe;
  edm::ParameterSet parameters;
  std::string metname;

  edm::InputTag theCaloJetCollectionLabel; 
  edm::InputTag theAKJetCollectionLabel;
  edm::InputTag theSCJetCollectionLabel;
  edm::InputTag theICJetCollectionLabel;
  edm::InputTag thePFJetCollectionLabel;
  edm::InputTag theJPTJetCollectionLabel;
  edm::InputTag theTriggerResultsLabel;


  //Cleaning parameters
  edm::ParameterSet theCleaningParameters;
  edm::InputTag _theVertexLabel;
  edm::InputTag _theGTLabel;
  std::string _hlt_PhysDec;

  std::vector<unsigned > _techTrigsAND;
  std::vector<unsigned > _techTrigsOR;
  std::vector<unsigned > _techTrigsNOT;

  bool _doPVCheck;
  bool _doHLTPhysicsOn;

  bool _tightBHFiltering;
  bool _tightHcalFiltering;

  int _nvtx_min;
  int _vtxndof_min;
  int _nvtxtrks_min;
  double _vtxchi2_max;
  double _vtxz_max;
  //

  int _LSBegin;
  int _LSEnd;

  HLTConfigProvider hltConfig_;
  std::string processname_;

  MonitorElement* hltpathME;
  MonitorElement* physdecME;
  MonitorElement* lumisecME;

  std::string LoJetTrigger;
  std::string HiJetTrigger;
  
  bool theJetAnalyzerFlag;  
  bool theIConeJetAnalyzerFlag;
  bool theJetCleaningFlag;

  bool theJetPtAnalyzerFlag;
  bool theJetPtCleaningFlag;

  bool thePFJetAnalyzerFlag;
  bool thePFJetCleaningFlag;

  bool theDiJetSelectionFlag;

  bool theJPTJetAnalyzerFlag;
  bool theJPTJetCleaningFlag;

  bool theCaloMETAnalyzerFlag;

  bool theTcMETAnalyzerFlag;

  bool theMuCorrMETAnalyzerFlag;

  bool thePfMETAnalyzerFlag;

  bool theHTMHTAnalyzerFlag;

  // the jet analyzer
  JetAnalyzer       * theJetAnalyzer;
  JetAnalyzer       * theAKJetAnalyzer;
  JetAnalyzer       * theSCJetAnalyzer;
  JetAnalyzer       * theICJetAnalyzer;  
  JetAnalyzer       * theCleanedAKJetAnalyzer;
  JetAnalyzer       * theCleanedSCJetAnalyzer;
  JetAnalyzer       * theCleanedICJetAnalyzer;
  JetAnalyzer       * theDiJetAnalyzer;

  JPTJetAnalyzer    * theJPTJetAnalyzer;
  JPTJetAnalyzer    * theCleanedJPTJetAnalyzer;

  PFJetAnalyzer     * thePFJetAnalyzer; 
  PFJetAnalyzer     * theCleanedPFJetAnalyzer; 

  JetPtAnalyzer     * thePtAKJetAnalyzer;
  JetPtAnalyzer     * thePtSCJetAnalyzer;
  JetPtAnalyzer     * thePtICJetAnalyzer;
  JetPtAnalyzer     * theCleanedPtAKJetAnalyzer;
  JetPtAnalyzer     * theCleanedPtSCJetAnalyzer;
  JetPtAnalyzer     * theCleanedPtICJetAnalyzer;

  CaloMETAnalyzer   * theCaloMETAnalyzer;
  CaloMETAnalyzer   * theCaloMETNoHFAnalyzer;
  CaloMETAnalyzer   * theCaloMETHOAnalyzer;
  CaloMETAnalyzer   * theCaloMETNoHFHOAnalyzer;

  METAnalyzer       * theTcMETAnalyzer;
  METAnalyzer       * theMuCorrMETAnalyzer;

  PFMETAnalyzer     * thePfMETAnalyzer;

  HTMHTAnalyzer     * theHTMHTAnalyzer;
  
  //
  int AKNjets_HB;
  int AKNjets_BE;
  int AKNjets_HE;
  int AKNjets_EF;
  int AKNjets_HF;
  int SCNjets_HB;
  int SCNjets_BE;
  int SCNjets_HE;
  int SCNjets_EF;
  int SCNjets_HF;
  int ICNjets_HB;
  int ICNjets_BE;
  int ICNjets_HE;
  int ICNjets_EF;
  int ICNjets_HF;
  int AKNCleanedjets_HB;
  int AKNCleanedjets_BE;
  int AKNCleanedjets_HE;
  int AKNCleanedjets_EF;
  int AKNCleanedjets_HF;
  int SCNCleanedjets_HB;
  int SCNCleanedjets_BE;
  int SCNCleanedjets_HE;
  int SCNCleanedjets_EF;
  int SCNCleanedjets_HF;
  int ICNCleanedjets_HB;
  int ICNCleanedjets_BE;
  int ICNCleanedjets_HE;
  int ICNCleanedjets_EF;
  int ICNCleanedjets_HF;
  int Dijets_HB;
  int Dijets_BE;
  int Dijets_HE;
  int DijetsCleaned_HB;
  int DijetsCleaned_BE;
  int DijetsCleaned_HE;

};
#endif  
