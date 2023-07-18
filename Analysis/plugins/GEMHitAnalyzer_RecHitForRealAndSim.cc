#ifndef GEMHitAnalyzer_RecHitForRealAndSim_H
#define GEMHitAnalyzer_RecHitForRealAndSim_H
// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// GEM
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
// Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
typedef tuple<int> Key1;
// typedef tuple<int, int> Key2;
// typedef tuple<int, int, int> Key3;

class GEMHitAnalyzer_RecHitForRealAndSim : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer_RecHitForRealAndSim(const edm::ParameterSet&);
  ~GEMHitAnalyzer_RecHitForRealAndSim();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  // edm::EDGetTokenT<GEMDigiSimLink> gemDigiSimLink_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  map<Key1, TH2I*> cls_vs_ie_;
  map<Key1, TH1I*> nrechit_;
  map<Key1, TH2I*> cls_vs_nrechit_;

};

GEMHitAnalyzer_RecHitForRealAndSim::GEMHitAnalyzer_RecHitForRealAndSim(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  // gemDigiSimLink_ = consumes<GEMDigiSimLink>(iConfig.getParameter<edm::InputTag>("gemDigiSimLinkLabel"));

//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
}

#endif


GEMHitAnalyzer_RecHitForRealAndSim::~GEMHitAnalyzer_RecHitForRealAndSim(){}


void
GEMHitAnalyzer_RecHitForRealAndSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
//  iSetup.getByToken(hGEMGeom_, hGEMGeom);

//  edm::ESHandle<GEMGeometry> hGEMGeom;
//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<GEMDigiCollection> gemDigis;
  iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  // edm::Handle<GEMDigiSimLink> gemDigiSimLink;
  // iEvent.getByToken(gemDigiSimLink_, gemDigiSimLink);

  int nrechit_st1 = 0;
  int nrechit_st2 = 0;
  int nrechit_st0 = 0;

  for (const GEMRegion* Region : GEMGeometry_->regions()){
    int re = Region->region();
    for (const GEMStation* Station : Region->stations()){
      int st = Station->station();
      if (st != 1) continue;
      for (const GEMRing* Ring : Station->rings()){
        int ri = Ring->ring();
        for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){
          for (const GEMChamber* Chamber : SuperChamber->chambers()){
            int la = Chamber->id().layer();
            int ch = Chamber->id().chamber();
            for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){
              GEMDetId etaPartId = etaPart->id();
              int ie = etaPartId.ieta();
            
              Key1 key1(st);

              auto RecHitRange = gemRecHits->get(etaPartId);
              for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                int firstStrip = rechit->firstClusterStrip();
                int clsSize = rechit->clusterSize();

                cout << "clsSize: " << clsSize << "(" << re << "/" << st << "/" << la << "/" << ch << "/" << ie << ")" << endl;

                if (st == 1){
                  nrechit_st1 += 1;
                }
                else if (st == 2){
                  nrechit_st2 += 1;
                }
                else if (st == 0){
                  nrechit_st0 += 1;
                }
                cls_vs_ie_[key1]->Fill(ie, clsSize);
              }
            } // eta partition loop
          } // chamber loop
        } // super chamber loop
      } // ring loop
    } // station loop
  } // region loop

  // float RecHits = gemRecHits->size();
  if (nrechit_st1 != 0){
    cout << "ge11 nrechit: " << nrechit_st1 << endl;
    Key1 key1_st1(1);
    nrechit_[key1_st1]->Fill(nrechit_st1);
  }
  if (nrechit_st2 != 0){
    cout << "ge21 nrechit: " << nrechit_st2 << endl;
    Key1 key1_st2(2);
    nrechit_[key1_st2]->Fill(nrechit_st2);
  }
  if (nrechit_st0 != 0){
    cout << "me0 nrechit: " << nrechit_st0 << endl;
    Key1 key1_st0(0);
    nrechit_[key1_st0]->Fill(nrechit_st0);
  }

  for (const GEMRegion* Region : GEMGeometry_->regions()){
    int re = Region->region();
    for (const GEMStation* Station : Region->stations()){
      int st = Station->station();
      if (st!=1) continue;
      for (const GEMRing* Ring : Station->rings()){
        int ri = Ring->ring();
        for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){
          for (const GEMChamber* Chamber : SuperChamber->chambers()){
            int la = Chamber->id().layer();
            int ch = Chamber->id().chamber();
            for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){
              GEMDetId etaPartId = etaPart->id();
              int ie = etaPartId.ieta();

              auto RecHitRange = gemRecHits->get(etaPartId);
              for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                int firstStrip = rechit->firstClusterStrip();
                int clsSize = rechit->clusterSize();

                Key1 key1(st);
                if (st == 1){
                  cls_vs_nrechit_[key1]->Fill(nrechit_st1, clsSize);
                }
                else if (st == 2){
                  cls_vs_nrechit_[key1]->Fill(nrechit_st2, clsSize);
                }
                else if (st == 0){
                  cls_vs_nrechit_[key1]->Fill(nrechit_st0, clsSize);
                }
              }
            } // eta partition loop
          } // chamber loop
        } // super chamber loop
      } // ring loop
    } // station loop
  } // region loop
}

void GEMHitAnalyzer_RecHitForRealAndSim::beginJob(){}
void GEMHitAnalyzer_RecHitForRealAndSim::endJob(){}

void GEMHitAnalyzer_RecHitForRealAndSim::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  for (int st=0; st<3; st++){
    if (st!=1) continue;
    Key1 key1(st);
    nrechit_[key1] = fs->make<TH1I>(Form("nrechit_S%d", st),
                                    Form("The number of rechit for each event : S%d", st),
                                    11, 0, 11);
    cls_vs_ie_[key1] = fs->make<TH2I>(Form("cls_vs_ie_S%d", st),
                                      Form("Cluster size for each rechit : S%d", st),
                                      17, 0, 17,
                                      385, 0, 385);
    cls_vs_nrechit_[key1] = fs->make<TH2I>(Form("cls_vs_nrechit_S%d", st),
                                            Form("Cluster size for each rechit : S%d", st),
                                            11, 0, 11,
                                            385, 0, 385);
  } // station loop
}
void GEMHitAnalyzer_RecHitForRealAndSim::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer_RecHitForRealAndSim);
