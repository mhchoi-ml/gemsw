#ifndef GEMHitAnalyzerForData_H
#define GEMHitAnalyzerForData_H
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
#include "DataFormats/GEMDigi/interface/GEMOHStatusCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVFATStatusCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
#include "DataFormats/TCDS/interface/TCDSRecord.h"
// #include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"
// Muon
// #include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
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


constexpr size_t max_trigger = 16;
using namespace std;
typedef tuple<int> Key1;
// typedef tuple<int, int> Key2;
// typedef tuple<int, int, int> Key3;

class GEMHitAnalyzerForData : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzerForData(const edm::ParameterSet&);
  ~GEMHitAnalyzerForData();

private:
  int maskChamberWithError(const GEMDetId& chamber_id, const edm::Handle<GEMVFATStatusCollection>, const edm::Handle<GEMOHStatusCollection>);
  int maskBigClusterEvent(const edm::Handle<GEMRecHitCollection> gemRecHits);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<GEMOHStatusCollection> oh_status_collection_;
  edm::EDGetTokenT<GEMVFATStatusCollection> vfat_status_collection_;
  edm::EDGetTokenT<TCDSRecord> tcdsRecord_;
  // edm::EDGetTokenT<OnlineLuminosityRecord> onlineLumiRecord_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;
  // edm::EDGetTokenT<edm::View<reco::Muon> >     muonHandle_;
  // edm::EDGetTokenT<GEMDigiCollection> gemDigis_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  TTree* t_Total;
  int b_ToNEv;
  int b_ToBigCls;
  int b_ToBunchId, b_ToOrbitNumber;
  long b_ToEvent, b_ToEventTime;

  TTree* t_Event;
  float b_EvNHit, b_EvAvgCls;
  int b_EvBigCls;
  int b_EvBunchId, b_EvOrbitNumber;
  long b_EvEvent, b_EvEventTime;
  // float b_EvInstLumi;

  TTree* t_Hit;
  int b_HitRe, b_HitSt, b_HitRi, b_HitLa, b_HitCh, b_HitIe;
  int b_HitCls;
  int b_HitChamErr, b_HitBigCls, b_HitBunchId;
};

GEMHitAnalyzerForData::GEMHitAnalyzerForData(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  oh_status_collection_ = consumes<GEMOHStatusCollection>(iConfig.getParameter<edm::InputTag>("OHInputLabel"));
  vfat_status_collection_ = consumes<GEMVFATStatusCollection>(iConfig.getParameter<edm::InputTag>("VFATInputLabel"));
  tcdsRecord_ = consumes<TCDSRecord>(iConfig.getParameter<edm::InputTag>("tcdsRecord"));
  // onlineLumiRecord_ = consumes<OnlineLuminosityRecord>(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"));
  // muonHandle_    = consumes<edm::View<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonLabel"));

//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>();

  t_Total = fs->make<TTree>("Toal", "event_total_info");
  #define ToBRANCH(name, suffix) t_Total->Branch(#name, & b_##name, #name "/" #suffix);
  ToBRANCH(ToNEv, I);

  ToBRANCH(ToBigCls, I);
  ToBRANCH(ToBunchId, I);
  ToBRANCH(ToOrbitNumber, I);
  ToBRANCH(ToEvent, l);
  ToBRANCH(ToEventTime, l);

  t_Event = fs->make<TTree>("Event", "gem_hits_per_event");
  #define EvBRANCH(name, suffix) t_Event->Branch(#name, & b_##name, #name "/" #suffix);
  EvBRANCH(EvNHit, F);
  EvBRANCH(EvAvgCls, F);

  EvBRANCH(EvBigCls, I);
  EvBRANCH(EvBunchId, I);
  EvBRANCH(EvOrbitNumber, I);
  EvBRANCH(EvEvent, l);
  EvBRANCH(EvEventTime, l);
  // HitBRANCH(instLumi, F);

  t_Hit = fs->make<TTree>("Hit", "gem_hits_per_hit");
  #define HitBRANCH(name, suffix) t_Hit->Branch(#name, & b_##name, #name "/" #suffix);
  HitBRANCH(HitRe, I);
  HitBRANCH(HitSt, I);
  HitBRANCH(HitRi, I);
  HitBRANCH(HitLa, I);
  HitBRANCH(HitCh, I);
  HitBRANCH(HitIe, I);
  HitBRANCH(HitCls, I);

  HitBRANCH(HitChamErr, I);
  HitBRANCH(HitBigCls, I);
  HitBRANCH(HitBunchId, I);
}

#endif


GEMHitAnalyzerForData::~GEMHitAnalyzerForData(){}

int GEMHitAnalyzerForData::maskBigClusterEvent(const edm::Handle<GEMRecHitCollection> gemRecHits) {
  // std::vector<std::vector<int>> n_hits_each_chamber(8, std::vector<int>(36, 0)); // giving max nDigis in one chamber
  vector<vector<int>> n_hits_each_etaPart(8, vector<int>(576, 0)); // giving max nDigis in one etaPartition

  for (const GEMRecHit& cluster : *gemRecHits) {
    // ++n_clusters;
    GEMDetId hit_id = cluster.gemId();
    int layer_index = (hit_id.region()+1)/2 + 2*(hit_id.station()-1) + 4*(hit_id.layer()-1);
    // n_hits = cluster.clusterSize();
    // n_hits_each_chamber[layer_index][hit_id.chamber() - 1] += n_hits;
    int eta_index = 16 * (hit_id.chamber() - 1) + hit_id.ieta() - 1;
    n_hits_each_etaPart[layer_index][eta_index] += cluster.clusterSize();
  }

  int mask = 0;
  int max_val;
  // for (const auto& row : n_hits_each_chamber) {
  for (const auto& row : n_hits_each_etaPart) {
    max_val = *std::max_element(row.begin(), row.end());
    // if (max_val > 384) return; // big cluster event filter
    if (max_val > 48) { // big cluster event filter
      mask = 1;
      return mask;
    }
  }
  return mask;
}

int GEMHitAnalyzerForData::maskChamberWithError(const GEMDetId& chamber_id,
                                                      const edm::Handle<GEMVFATStatusCollection> vfat_status_collection,
                                                      const edm::Handle<GEMOHStatusCollection> oh_status_collection) {
  int mask = 0;
  for (auto iter = oh_status_collection->begin(); iter != oh_status_collection->end(); iter++) {
    const auto [oh_id, range] = (*iter);
    if (chamber_id != oh_id) {
      continue;
    }

    for (auto oh_status = range.first; oh_status != range.second; oh_status++) {
      if (oh_status->isBad()) {
        // GEMOHStatus is bad. Mask this chamber.
        mask = 1;
        return mask;
      }  // isBad
    }  // range
  }  // collection
  for (auto iter = vfat_status_collection->begin(); iter != vfat_status_collection->end(); iter++) {
    const auto [vfat_id, range] = (*iter);
    if (chamber_id != vfat_id.chamberId()) {
      continue;
    }
    for (auto vfat_status = range.first; vfat_status != range.second; vfat_status++) {
      if (vfat_status->isBad()) {
        mask = 1;
        return mask;
      }  // isBad
    }  // range
  }  // collection
  return mask;
}

void
GEMHitAnalyzerForData::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
//  iSetup.getByToken(hGEMGeom_, hGEMGeom);

//  edm::ESHandle<GEMGeometry> hGEMGeom;
//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  // edm::Handle<GEMDigiCollection> gemDigis;
  // iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<GEMVFATStatusCollection> vfat_status_collection;
  iEvent.getByToken(vfat_status_collection_, vfat_status_collection);

  edm::Handle<GEMOHStatusCollection> oh_status_collection;
  iEvent.getByToken(oh_status_collection_, oh_status_collection);

  edm::Handle<TCDSRecord> record;
  iEvent.getByToken(tcdsRecord_, record);

  // edm::Handle<OnlineLuminosityRecord> onlineLumiRecord;
  // iEvent.getByToken(onlineLumiRecord_, onlineLumiRecord);

  // edm::Handle<edm::View<reco::Muon> > muonHandle;
  // iEvent.getByToken(muonHandle_, muonHandle);




  cout << record.isValid() << ", " << gemRecHits.isValid() << endl;
  if (!record.isValid() || !gemRecHits.isValid()) {
    cout << "Invalid Event" << endl;
    return;
  }

    /*  we don't want to use hit based cut
  if ((n_clusters > 650.) || ((n_clusters - 50.) * (2000. / 600.) > n_hits)) {
    return;
   }
  */

  /* W.Heo's flower event cut
  int max_val = *std::max_element(n_hits_each_chamber[0].begin(),n_hits_each_chamber[0].end());
  for (const auto& row : n_hits_each_chamber) {
    max_val = std::max(max_val, *std::max_element(row.begin(), row.end()));
    if (max_val > 384) return;
  }
  */

  /*  Laurant's method */
  for (size_t i = 0; i < max_trigger; ++i) {
    long l1a_diff = 3564 * (record->getOrbitNr() - record->getL1aHistoryEntry(i).getOrbitNr())
        + record->getBXID() - record->getL1aHistoryEntry(i).getBXID();

    if ((l1a_diff > 150) && (l1a_diff < 200)) {
      cout << "Flower Event" << endl;
      return;
    }
  }

  int BigCls = maskBigClusterEvent(gemRecHits);
  int BunchId = iEvent.bunchCrossing();
  int OrbitNumber = iEvent.orbitNumber();
  long EventTime = iEvent.time().unixTime();
  long Event = iEvent.id().event();

  cout << "hit size: " << gemRecHits->size() << endl;
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
            int ChamErr = maskChamberWithError(chId, vfat_status_collection, oh_status_collection);
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
                b_HitBigCls = BigCls;
                b_HitChamErr = ChamErr;
                b_HitBunchId = BunchId;
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
    b_EvBigCls = BigCls;
    // b_instLumi = onlineLumiRecord->instLumi();
    b_EvBunchId = BunchId;
    b_EvOrbitNumber = OrbitNumber;
    b_EvEventTime = EventTime;
    b_EvEvent = Event;
    t_Event->Fill();
  }
  b_ToNEv = 1;
  b_ToBigCls = BigCls;
  // b_instLumi = onlineLumiRecord->instLumi();
  b_ToBunchId = BunchId;
  b_ToOrbitNumber = OrbitNumber;
  b_ToEventTime = EventTime;
  b_ToEvent = Event;
  t_Total->Fill();
}

void GEMHitAnalyzerForData::beginJob(){}
void GEMHitAnalyzerForData::endJob(){}

void GEMHitAnalyzerForData::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  // h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);

}
void GEMHitAnalyzerForData::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzerForData);
