#ifndef GEMSimClsAnalyzer_H
#define GEMSimClsAnalyzer_H
// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <map>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
// GEM
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
// Muon
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "Validation/MuonHits/interface/MuonHitHelper.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"


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
using namespace edm;
typedef tuple<int> Key1;
typedef tuple<int, int> Key2;
typedef tuple<int, int, int> Key3;

class GEMSimClsAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMSimClsAnalyzer(const edm::ParameterSet&);
  ~GEMSimClsAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // void initMuonValue();
  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;
  // edm::EDGetTokenT<edm::View<reco::Muon>> muonHandle_;
  // edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttb_;

  // edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  // const reco::VertexCollection* vertexes_;
  // edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  TTree* t_Total;
  int b_ToNEv;

  TTree* t_Event;
  float b_EvNHit, b_EvAvgCls;

  TTree* t_Hit;
  int b_HitRe, b_HitSt, b_HitRi, b_HitLa, b_HitCh, b_HitIe;
  int b_HitCls;

  TTree* t_Muon;
  double b_MuP, b_MuEnergy, b_MuMass, b_MuPt, b_MuPhi, b_MuTheta, b_MuEta;
  int b_MuGlo;
  int b_MuTight;
  int b_MuGE11;

  // TTree* t_GEMMuon;
  // double b_GMuP, b_GMuEnergy, b_GMuMass, b_GMuPt, b_GMuPhi, b_GMuTheta, b_GMuEta;
  // int b_GMuGlo;
  // int b_GMuTight;

  TTree* t_MuonHit;
  int b_MuHitRe, b_MuHitSt, b_MuHitRi, b_MuHitLa, b_MuHitCh, b_MuHitIe;
  int b_MuHitCls;

  int b_MuHitGlo;
  int b_MuHitTight;
  double b_MuHitPt;




  //////////////////////////////////////////////////////////////////////////////
  // const data members initialized in the member initializer list
  // mainly retrieved from edm::ParameterSet
  //////////////////////////////////////////////////////////////////////////////
  // ES
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> kGEMGeometryTokenBeginRun_;
  // ED
  const edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  const edm::EDGetTokenT<edm::View<reco::Muon> > muonHandle_;
  
};

GEMSimClsAnalyzer::GEMSimClsAnalyzer(const edm::ParameterSet& iConfig)
  // : hGEMGeom_(esConsumes()),
  //   hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>()),
    : kGEMGeometryTokenBeginRun_(esConsumes<edm::Transition::BeginRun>()),
      gemRecHits_(consumes<GEMRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("recHitTag"))),
      muonHandle_ (consumes<edm::View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")))
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));

  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>();
  hGEMGeomBeginRun_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 


  t_Total = fs->make<TTree>("Total", "event_total_info");
  #define ToBRANCH(name, suffix) t_Total->Branch(#name, & b_##name, #name "/" #suffix);
  ToBRANCH(ToNEv, I);

  t_Event = fs->make<TTree>("Event", "gem_hits_per_event");
  #define EvBRANCH(name, suffix) t_Event->Branch(#name, & b_##name, #name "/" #suffix);
  EvBRANCH(EvNHit, F);
  EvBRANCH(EvAvgCls, F);

  t_Hit = fs->make<TTree>("Hit", "gem_hits_per_hit");
  #define HitBRANCH(name, suffix) t_Hit->Branch(#name, & b_##name, #name "/" #suffix);
  HitBRANCH(HitRe, I);
  HitBRANCH(HitSt, I);
  HitBRANCH(HitRi, I);
  HitBRANCH(HitLa, I);
  HitBRANCH(HitCh, I);
  HitBRANCH(HitIe, I);
  HitBRANCH(HitCls, I);

  t_Muon = fs->make<TTree>("Muon", "muon_obj");
  #define MuBRANCH(name, suffix) t_Muon->Branch(#name, & b_##name, #name "/" #suffix);
  MuBRANCH(MuP, D);
  MuBRANCH(MuEnergy, D);
  MuBRANCH(MuMass, D);
  MuBRANCH(MuPt, D);
  MuBRANCH(MuPhi, D);
  MuBRANCH(MuTheta, D);
  MuBRANCH(MuEta, D);
  MuBRANCH(MuGlo, I);
  MuBRANCH(MuTight, I);
  MuBRANCH(MuGE11, I);

  // t_GEMMuon = fs->make<TTree>("GEMMuon", "muon_obj_on_gem");
  // #define GMuBRANCH(name, suffix) t_GEMMuon->Branch(#name, & b_##name, #name "/" #suffix);
  // GMuBRANCH(GMuP, D);
  // GMuBRANCH(GMuEnergy, D);
  // GMuBRANCH(GMuMass, D);
  // GMuBRANCH(GMuPt, D);
  // GMuBRANCH(GMuPhi, D);
  // GMuBRANCH(GMuTheta, D);
  // GMuBRANCH(GMuEta, D);
  // GMuBRANCH(GMuGlo, I);
  // GMuBRANCH(GMuTight, I);

  t_MuonHit = fs->make<TTree>("MuonHit", "gem_hits_by_muon");
  #define MuHitBRANCH(name, suffix) t_MuonHit->Branch(#name, & b_##name, #name "/" #suffix);
  MuHitBRANCH(MuHitRe, I);
  MuHitBRANCH(MuHitSt, I);
  MuHitBRANCH(MuHitRi, I);
  MuHitBRANCH(MuHitLa, I);
  MuHitBRANCH(MuHitCh, I);
  MuHitBRANCH(MuHitIe, I);
  MuHitBRANCH(MuHitCls, I);
  MuHitBRANCH(MuHitGlo, I);
  MuHitBRANCH(MuHitTight, I);
  MuHitBRANCH(MuHitPt, D);
}

#endif


GEMSimClsAnalyzer::~GEMSimClsAnalyzer(){}


void
GEMSimClsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
  // iSetup.getByToken(hGEMGeom_, hGEMGeom);

  // edm::ESHandle<GEMGeometry> hGEMGeom;
  // iSetup.get<MuonGeometryRecord>().get(hGEMGeom);

  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  // const GEMGeometry* gem = hGEMGeom.product();

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByToken(muonHandle_, muonHandle);

  // edm::Handle<GEMDigiCollection> gemDigis;
  // iEvent.getByToken(gemDigis_, gemDigis);


  float EvNHit = 0;
  float EvSumCls = 0;
  for (const GEMRegion* Region : GEMGeometry_->regions()){ // +, -
    int re = Region->region();
    for (const GEMStation* Station : Region->stations()){ // GE1/1, GE2/1, ME0
      int st = Station->station();
      if (st != 1) continue; // see only GE1/1 for now
      for (const GEMRing* Ring : Station->rings()){
        int ri = Ring->ring();
        for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){ // GE1/1:2, GE2/1:2, ME0:6
          for (const GEMChamber* Chamber : SuperChamber->chambers()){ // GE1/1:36(long/short), GE2/1:18, ME0:18
            GEMDetId chId = Chamber->id();
            int la = chId.layer();
            int ch = chId.chamber();
            for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){ // GE1/1:8, GE2/1:16, ME0:8
              GEMDetId ieId = etaPart->id();
              int ie = ieId.ieta();
              auto RecHitRange = gemRecHits->get(ieId);
              for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                int firstStrip = rechit->firstClusterStrip();
                int clsSize = rechit->clusterSize();
                b_HitRe = re;
                b_HitSt = st;
                b_HitRi = ri;
                b_HitLa = la;
                b_HitCh = ch;
                b_HitIe = ie;
                b_HitCls = clsSize;
                t_Hit->Fill();
                EvNHit++;
                EvSumCls += clsSize;
              }
            } // eta partition loop
          } // chamber loop
        } // super chamber loop
      } // ring loop
    } // station loop
  } // region loop
  if (EvNHit != 0){
    b_EvNHit = EvNHit;
    b_EvAvgCls = EvSumCls/EvNHit;
    t_Event->Fill();
  }
  b_ToNEv = 1;
  t_Total->Fill();


  //////////////////////////////////////////////////////////////////////////////
  // get data from Event
  //////////////////////////////////////////////////////////////////////////////
  // const GEMRecHitCollection* rechit_collection = nullptr;
  // if (auto handle = iEvent.getHandle(kGEMRecHitCollectionToken_)) {
  //   rechit_collection = handle.product();
  // } else {
  //   edm::LogError("ReadoutError") << "failed to get GEMRecHitCollection";
  //   return;
  // }

  // const edm::View<reco::Muon>* muon_view = nullptr;
  // // cout << "h: " << kMuonViewToken_ << endl;
  // if (auto handle = iEvent.getHandle(kMuonViewToken_)) {
  //   cout << "nMuon: " << handle->size() << endl;
  //   muon_view = handle.product();
  // } else {
  //   edm::LogError("ReadoutError") << "failed to get View<Muon>";
  //   return;
  // }


  //////////////////////////////////////////////////////////////////////////////
  //  Main loop
  //////////////////////////////////////////////////////////////////////////////
  for (const reco::Muon& muon : *(muonHandle.product())) {
    b_MuGlo = 1;
    b_MuTight = 1;

    b_MuHitGlo = 1;
    b_MuHitTight = 1;

    if (!muon.isGlobalMuon()) {
      b_MuGlo = 0;
      b_MuHitGlo = 0;
    }
    if (!muon.passed(reco::Muon::CutBasedIdTight)) {
      b_MuTight = 0;
      b_MuHitTight = 0;
    }
    b_MuHitPt = muon.pt(); // [GeV]


    b_MuP = muon.p();
    b_MuEnergy = muon.energy();
    b_MuMass = muon.mass();
    b_MuPt = muon.pt();
    b_MuPhi = muon.phi();
    b_MuTheta = muon.theta();
    b_MuEta = muon.eta();
    b_MuGE11 = 0;

    const reco::Track* muonTrack = 0;
    if ( muon.globalTrack().isNonnull() ) muonTrack = muon.globalTrack().get();
    else if ( muon.outerTrack().isNonnull()  ) muonTrack = muon.outerTrack().get();
    if (!muonTrack) continue;
    // b_nMuons++;

    // nGEMHitInMuontrack = 0;
    for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
      const DetId id = (*hit)->geographicalId();
      if (id.det() == 2 /* Muon */ && id.subdetId() == 4 /* GEM */) {
        GEMDetId gemid = id.rawId();

        if (gemid.station() != 1) continue;
        else b_MuGE11 = 1;

        auto RecHitRange = gemRecHits->get(gemid);
        for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
          int cls = rechit->clusterSize();

          const int re = gemid.region();
          const int st = gemid.station();
          const int ri = gemid.ring();
          const int la = gemid.layer();
          const int ch = gemid.chamber();
          const int ie = gemid.ieta();

          b_MuHitRe = re;
          b_MuHitSt = st;
          b_MuHitRi = ri;
          b_MuHitLa = la;
          b_MuHitCh = ch;
          b_MuHitIe = ie;
          b_MuHitCls = cls;
          t_MuonHit->Fill();
        }
      }     // nGEMHitInMuontrack++;
    }       // Rechit
    t_Muon->Fill();
  }        // Muon
}

void GEMSimClsAnalyzer::beginJob(){}
void GEMSimClsAnalyzer::endJob(){}

void GEMSimClsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
}
// void GEMSimClsAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
//   /* GEM Geometry */
//   // edm::ESHandle<GEMGeometry> hGEMGeom;
//   // hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

// //  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
//   // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

//   // h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);
// }

void GEMSimClsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){
}

DEFINE_FWK_MODULE(GEMSimClsAnalyzer);