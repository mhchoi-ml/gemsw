#ifndef GEMHitAnalyzer_HitForReal_H
#define GEMHitAnalyzer_HitForReal_H
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
// #include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMOHStatusCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVFATStatusCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
// #include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
#include "DataFormats/TCDS/interface/TCDSRecord.h"
// #include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"
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


constexpr size_t max_trigger = 16;
using namespace std;
typedef tuple<int> Key1;
// typedef tuple<int, int> Key2;
// typedef tuple<int, int, int> Key3;

class GEMHitAnalyzer_HitForReal : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer_HitForReal(const edm::ParameterSet&);
  ~GEMHitAnalyzer_HitForReal();

private:
  bool maskChamberWithError(const GEMDetId& chamber_id, const edm::Handle<GEMVFATStatusCollection>, const edm::Handle<GEMOHStatusCollection>);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<GEMOHStatusCollection> oh_status_collection_;
  edm::EDGetTokenT<GEMVFATStatusCollection> vfat_status_collection_;
  edm::EDGetTokenT<TCDSRecord> tcdsRecord_;
  // edm::EDGetTokenT<OnlineLuminosityRecord> onlineLumiRecord_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  TH1I* h_nEvents;

  TTree* t_RecEvent;
  float b_RecEvNrechit;
  float b_RecEvAvgCls;

  TTree* t_RecIeta;
  float b_RecIeNrechit;
  float b_RecIeAvgCls;

  TTree* t_RecHit;
  int b_RecHitRe, b_RecHitSt, b_RecHitLa, b_RecHitCh, b_RecHitIe;
  float b_RecHitEvNrechit;
  int b_RecHitCls;

  // float b_instLumi;
  int b_bunchId, b_orbitNumber;
  long b_event, b_eventTime;
  int b_chamErr, b_bigClsEv;
};

GEMHitAnalyzer_HitForReal::GEMHitAnalyzer_HitForReal(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  oh_status_collection_ = consumes<GEMOHStatusCollection>(iConfig.getParameter<edm::InputTag>("OHInputLabel"));
  vfat_status_collection_ = consumes<GEMVFATStatusCollection>(iConfig.getParameter<edm::InputTag>("VFATInputLabel"));
  tcdsRecord_ = consumes<TCDSRecord>(iConfig.getParameter<edm::InputTag>("tcdsRecord"));
  // onlineLumiRecord_ = consumes<OnlineLuminosityRecord>(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"));

//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>();

  t_RecEvent = fs->make<TTree>("RecEvent", "gem_rechits_per_event");
  #define RecEvBRANCH(name, suffix) t_RecEvent->Branch(#name, & b_##name, #name "/" #suffix);
  RecEvBRANCH(RecEvNrechit, F);
  RecEvBRANCH(RecEvAvgCls, F);

  t_RecIeta = fs->make<TTree>("RecIeta", "gem_rechits_per_ieta");
  #define RecIeBRANCH(name, suffix) t_RecIeta->Branch(#name, & b_##name, #name "/" #suffix);
  RecIeBRANCH(RecIeNrechit, F);
  RecIeBRANCH(RecIeAvgCls, F);

  t_RecHit = fs->make<TTree>("RecHit", "gem_rechits_per_hit");
  #define RecHitBRANCH(name, suffix) t_RecHit->Branch(#name, & b_##name, #name "/" #suffix);
  RecHitBRANCH(RecHitRe, I);
  RecHitBRANCH(RecHitSt, I);
  RecHitBRANCH(RecHitLa, I);
  RecHitBRANCH(RecHitCh, I);
  RecHitBRANCH(RecHitIe, I);
  RecHitBRANCH(RecHitEvNrechit, F);
  RecHitBRANCH(RecHitCls, I);

  // RecHitBRANCH(instLumi, F);
  RecHitBRANCH(bunchId, I);
  RecHitBRANCH(orbitNumber, I);
  RecHitBRANCH(event, l);
  RecHitBRANCH(eventTime, l);
  RecHitBRANCH(chamErr, I);
  RecHitBRANCH(bigClsEv, I);

}

#endif


GEMHitAnalyzer_HitForReal::~GEMHitAnalyzer_HitForReal(){}

bool GEMHitAnalyzer_HitForReal::maskChamberWithError(const GEMDetId& chamber_id,
                                                      const edm::Handle<GEMVFATStatusCollection> vfat_status_collection,
                                                      const edm::Handle<GEMOHStatusCollection> oh_status_collection) {
  const bool mask = true;
  for (auto iter = oh_status_collection->begin(); iter != oh_status_collection->end(); iter++) {
    const auto [oh_id, range] = (*iter);
    if (chamber_id != oh_id) {
      continue;
    }

    for (auto oh_status = range.first; oh_status != range.second; oh_status++) {
      if (oh_status->isBad()) {
        // GEMOHStatus is bad. Mask this chamber.
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
        return mask;
      }  // isBad
    }  // range
  }  // collection
  return not mask;
}

void
GEMHitAnalyzer_HitForReal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  if (!record.isValid() || !gemRecHits.isValid()) {
    cout << "Error!" << endl;
    return;
  }

  // std::vector<std::vector<int>> n_hits_each_chamber(8, std::vector<int>(36, 0)); // giving max nDigis in one chamber
  std::vector<std::vector<int>> n_hits_each_etaPart(8, std::vector<int>(576, 0)); // giving max nDigis in one etaPartition

  for (const GEMRecHit& cluster : *gemRecHits) {
    // ++n_clusters;
    GEMDetId hit_id = cluster.gemId();
    int layer_index = (hit_id.region()+1)/2 + 2*(hit_id.station()-1) + 4*(hit_id.layer()-1);
    // n_hits = cluster.clusterSize();
    // n_hits_each_chamber[layer_index][hit_id.chamber() - 1] += n_hits;
    int eta_index = 16 * (hit_id.chamber() - 1) + hit_id.ieta() - 1;
    n_hits_each_etaPart[layer_index][eta_index] += cluster.clusterSize();
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

  b_bigClsEv = 0;
  int max_val;
  // for (const auto& row : n_hits_each_chamber) {
  for (const auto& row : n_hits_each_etaPart) {
    max_val = *std::max_element(row.begin(), row.end());
    // if (max_val > 384) return; // big cluster event filter
    if (max_val > 48) { // big cluster event filter
      b_bigClsEv = 1; 
    }
  }

  /*  Laurant's method */
  for (size_t i = 0; i < max_trigger; ++i) {
    long l1a_diff = 3564 * (record->getOrbitNr() - record->getL1aHistoryEntry(i).getOrbitNr())
        + record->getBXID() - record->getL1aHistoryEntry(i).getBXID();

    if ((l1a_diff > 150) && (l1a_diff < 200)) {
      cout << "Flower event!!!" << endl;
      return;
    }
  }



  float RecEvNrechits = 0;
  float RecEvSumCls = 0;
  for (const GEMRegion* Region : GEMGeometry_->regions()){
      for (const GEMStation* Station : Region->stations()){
        int st = Station->station();
        if (st != 1) continue;
        for (const GEMRing* Ring : Station->rings()){
          for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){
            for (const GEMChamber* Chamber : SuperChamber->chambers()){
              // GEMDetId chId = Chamber->id();
              // b_chamErr = maskChamberWithError(chId, vfat_status_collection, oh_status_collection);
              for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){
                GEMDetId ieId = etaPart->id();
                auto RecHitRange = gemRecHits->get(ieId);

                float RecIeNrechits = 0;
                float RecIeSumCls = 0;
                for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                  int firstStrip = rechit->firstClusterStrip();
                  int clsSize = rechit->clusterSize();

                  RecIeNrechits++;
                  RecIeSumCls += clsSize;

                  RecEvNrechits++;
                  RecEvSumCls += clsSize;
                }
                if (RecIeNrechits != 0){
                  b_RecIeNrechit = RecIeNrechits;
                  b_RecIeAvgCls = RecIeSumCls/RecIeNrechits;
                  t_RecIeta->Fill();
                }
              } // eta partition loop
            } // chamber loop
          } // super chamber loop
        } // ring loop
      } // station loop
  } // region loop
  if (RecEvNrechits != 0){
    b_RecEvNrechit = RecEvNrechits;
    b_RecEvAvgCls = RecEvSumCls/RecEvNrechits;
    t_RecEvent->Fill();
  }

  // b_instLumi = onlineLumiRecord->instLumi();
  b_bunchId = iEvent.bunchCrossing();
  b_orbitNumber = iEvent.orbitNumber();
  b_eventTime = iEvent.time().unixTime();
  b_event = iEvent.id().event();

  for (const GEMRegion* Region : GEMGeometry_->regions()){
      int re = Region->region();
      for (const GEMStation* Station : Region->stations()){
        int st = Station->station();
        if (st != 1) continue;
        for (const GEMRing* Ring : Station->rings()){
          int ri = Ring->ring();
          for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){
            for (const GEMChamber* Chamber : SuperChamber->chambers()){
              GEMDetId chId = Chamber->id();
              b_chamErr = maskChamberWithError(chId, vfat_status_collection, oh_status_collection);
              int la = chId.layer();
              int ch = chId.chamber();
              for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){
                GEMDetId ieId = etaPart->id();
                int ie = ieId.ieta();
              
                auto RecHitRange = gemRecHits->get(ieId);
                for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                  int firstStrip = rechit->firstClusterStrip();
                  int clsSize = rechit->clusterSize();

                  b_RecHitRe = re;
                  b_RecHitSt = st;
                  b_RecHitLa = la;
                  b_RecHitCh = ch;
                  b_RecHitIe = ie;
                  b_RecHitCls = clsSize;
                  b_RecHitEvNrechit = RecEvNrechits;
                  t_RecHit->Fill();
                }

              } // eta partition loop
            } // chamber loop
          } // super chamber loop
        } // ring loop
      } // station loop
  } // region loop

  h_nEvents->Fill(1);
}

void GEMHitAnalyzer_HitForReal::beginJob(){}
void GEMHitAnalyzer_HitForReal::endJob(){}

void GEMHitAnalyzer_HitForReal::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);

}
void GEMHitAnalyzer_HitForReal::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer_HitForReal);
