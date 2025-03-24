#ifndef GEMDataClsAnalyzer_H
#define GEMDataClsAnalyzer_H
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
#include "DataFormats/GEMDigi/interface/GEMAMCStatusCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMC13StatusCollection.h"
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


constexpr size_t max_trigger = 16;
using namespace std;
using namespace edm;
typedef tuple<int> Key1;
typedef tuple<int, int> Key2;
typedef tuple<int, int, int> Key3;

class GEMDataClsAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMDataClsAnalyzer(const edm::ParameterSet&);
  ~GEMDataClsAnalyzer();

private:
  int maskChamberWithError(const GEMDetId& chamber_id, const edm::Handle<GEMAMCStatusCollection>,
                                                        const edm::Handle<GEMAMC13StatusCollection>,
                                                        const edm::Handle<GEMVFATStatusCollection>,
                                                        const edm::Handle<GEMOHStatusCollection>);
  int maskFlowerEvent(const edm::Handle<TCDSRecord> record);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  // edm::EDGetTokenT<OnlineLuminosityRecord> onlineLumiRecord_;
  // edm::EDGetTokenT<edm::View<reco::Muon> >     muonHandle_;

  // ES
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> kGEMGeometryTokenBeginRun_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_; 
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;
  // ED
  const edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  const edm::EDGetTokenT<edm::View<reco::Muon> > muonHandle_;
  edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMAMCStatusCollection> amc_status_collection_;
  edm::EDGetTokenT<GEMAMC13StatusCollection> amc13_status_collection_;
  edm::EDGetTokenT<GEMOHStatusCollection> oh_status_collection_;
  edm::EDGetTokenT<GEMVFATStatusCollection> vfat_status_collection_;
  edm::EDGetTokenT<TCDSRecord> tcdsRecord_;


  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  TTree* t_Event;
  // int b_ToBunchId, b_ToOrbitNumber, b_ToLumi;
  // long b_ToEvent, b_ToEventTime;
  // // float b_EvInstLumi;
  int b_EvFlr;
  // int b_ToTrigBunchId[16], b_ToTrigOrbitNumber[16];
  // int b_ToNHit, b_ToMul;


  // TTree* t_Event;
  // float b_EvNHit, b_EvAvgCls;


  // TTree* t_Hit;
  // int b_HitRe, b_HitSt, b_HitRi, b_HitLa, b_HitCh, b_HitIe;

  // int b_HitCls;
  // int b_HitFlr, b_HitChamErr;

  // int b_HitBunchId, b_HitOrbitNumber, b_HitLumi;
  // long b_HitEvent, b_HitEventTime;


  TTree* t_Muon;
  int b_MuBunchId, b_MuFlr;

  int b_MuGE11;

  double b_MuP, b_MuEnergy, b_MuMass, b_MuPt, b_MuPhi, b_MuTheta, b_MuEta;
  int b_MuGlo, b_MuTight;

  TTree* t_MuonDigi;
  int b_MuDigiFlr;

  int b_MuDigiChamErr;

  int b_MuDigiGlo, b_MuDigiTight;
  double b_MuDigiPt;

  int b_MuDigiRe, b_MuDigiSt, b_MuDigiRi, b_MuDigiLa, b_MuDigiCh, b_MuDigiIe, b_MuDigiStp;

  TTree* t_MuonHit;
  int b_MuHitLumi, b_MuHitBunchId, b_MuHitFlr;

  int b_MuHitChamErr;

  int b_MuHitGlo, b_MuHitTight;
  double b_MuHitPt;

  int b_MuHitRe, b_MuHitSt, b_MuHitRi, b_MuHitLa, b_MuHitCh, b_MuHitIe, b_MuHitCls;



  //////////////////////////////////////////////////////////////////////////////
  // const data members initialized in the member initializer list
  // mainly retrieved from edm::ParameterSet
  //////////////////////////////////////////////////////////////////////////////

};

GEMDataClsAnalyzer::GEMDataClsAnalyzer(const edm::ParameterSet& iConfig)
  // : hGEMGeom_(esConsumes()),
    // hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
    : kGEMGeometryTokenBeginRun_(esConsumes<edm::Transition::BeginRun>()),
      gemRecHits_(consumes<GEMRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("recHitLabel"))),
      muonHandle_(consumes<edm::View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonLabel")))
{
  gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  amc_status_collection_ = consumes<GEMAMCStatusCollection>(iConfig.getParameter<edm::InputTag>("AMCInputLabel"));
  amc13_status_collection_ = consumes<GEMAMC13StatusCollection>(iConfig.getParameter<edm::InputTag>("AMC13InputLabel"));
  oh_status_collection_ = consumes<GEMOHStatusCollection>(iConfig.getParameter<edm::InputTag>("OHInputLabel"));
  vfat_status_collection_ = consumes<GEMVFATStatusCollection>(iConfig.getParameter<edm::InputTag>("VFATInputLabel"));
  tcdsRecord_ = consumes<TCDSRecord>(iConfig.getParameter<edm::InputTag>("tcdsRecord"));
  // onlineLumiRecord_ = consumes<OnlineLuminosityRecord>(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"));

  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>();
  hGEMGeomBeginRun_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 


  t_Event = fs->make<TTree>("Event", "info_per_event");
  #define EvBRANCH(name, suffix) t_Event->Branch(#name, & b_##name, #name "/" #suffix);
  // ToBRANCH(ToBunchId, I);
  // ToBRANCH(ToOrbitNumber, I);
  // ToBRANCH(ToEvent, l);
  // ToBRANCH(ToEventTime, l);
  // ToBRANCH(ToLumi, I);
  //   // HitBRANCH(instLumi, F);

  EvBRANCH(EvFlr, I);

  // ToBRANCH(ToTrigBunchId[16], I);
  // ToBRANCH(ToTrigOrbitNumber[16], I);
  // ToBRANCH(ToNHit, I);
  // ToBRANCH(ToMul, I);


  // t_Event = fs->make<TTree>("Event", "gem_hits_per_event");
  // #define EvBRANCH(name, suffix) t_Event->Branch(#name, & b_##name, #name "/" #suffix);
  // EvBRANCH(EvNHit, F);
  // EvBRANCH(EvAvgCls, F);

  // t_Hit = fs->make<TTree>("Hit", "gem_hits_per_hit");
  // #define HitBRANCH(name, suffix) t_Hit->Branch(#name, & b_##name, #name "/" #suffix);
  // HitBRANCH(HitRe, I);
  // HitBRANCH(HitSt, I);
  // HitBRANCH(HitRi, I);
  // HitBRANCH(HitLa, I);
  // HitBRANCH(HitCh, I);
  // HitBRANCH(HitIe, I);

  // HitBRANCH(HitCls, I);
  // HitBRANCH(HitFlr, I);
  // HitBRANCH(HitChamErr, I);

  // HitBRANCH(HitBunchId, I);
  // HitBRANCH(HitOrbitNumber, I);
  // HitBRANCH(HitEvent, l);
  // HitBRANCH(HitEventTime, l);
  // HitBRANCH(HitLumi, I);

  t_Muon = fs->make<TTree>("Muon", "muon_obj");
  #define MuBRANCH(name, suffix) t_Muon->Branch(#name, & b_##name, #name "/" #suffix);
  MuBRANCH(MuFlr, I);
  MuBRANCH(MuBunchId, I);

  MuBRANCH(MuGE11, I);

  MuBRANCH(MuP, D);
  MuBRANCH(MuEnergy, D);
  MuBRANCH(MuMass, D);
  MuBRANCH(MuPt, D);
  MuBRANCH(MuPhi, D);
  MuBRANCH(MuTheta, D);
  MuBRANCH(MuEta, D);
  MuBRANCH(MuGlo, I);
  MuBRANCH(MuTight, I);

  t_MuonDigi = fs->make<TTree>("MuonDigi", "gem_digis_by_muon");
  #define MuDigiBRANCH(name, suffix) t_MuonDigi->Branch(#name, & b_##name, #name "/" #suffix);
  MuDigiBRANCH(MuDigiFlr, I);

  MuDigiBRANCH(MuDigiChamErr, I);

  MuDigiBRANCH(MuDigiGlo, I);
  MuDigiBRANCH(MuDigiTight, I);
  MuDigiBRANCH(MuDigiPt, D);

  MuDigiBRANCH(MuDigiRe, I);
  MuDigiBRANCH(MuDigiSt, I);
  MuDigiBRANCH(MuDigiRi, I);
  MuDigiBRANCH(MuDigiLa, I);
  MuDigiBRANCH(MuDigiCh, I);
  MuDigiBRANCH(MuDigiIe, I);
  MuDigiBRANCH(MuDigiStp, I);

  t_MuonHit = fs->make<TTree>("MuonHit", "gem_hits_by_muon");
  #define MuHitBRANCH(name, suffix) t_MuonHit->Branch(#name, & b_##name, #name "/" #suffix);
  MuHitBRANCH(MuHitLumi, I);
  MuHitBRANCH(MuHitBunchId, I);
  MuHitBRANCH(MuHitFlr, I);

  MuHitBRANCH(MuHitChamErr, I);

  MuHitBRANCH(MuHitGlo, I);
  MuHitBRANCH(MuHitTight, I);
  MuHitBRANCH(MuHitPt, D);

  MuHitBRANCH(MuHitRe, I);
  MuHitBRANCH(MuHitSt, I);
  MuHitBRANCH(MuHitRi, I);
  MuHitBRANCH(MuHitLa, I);
  MuHitBRANCH(MuHitCh, I);
  MuHitBRANCH(MuHitIe, I);
  MuHitBRANCH(MuHitCls, I);
}

#endif


GEMDataClsAnalyzer::~GEMDataClsAnalyzer(){}

// int GEMDataClsAnalyzer::maskBigClusterEvent(const edm::Handle<GEMRecHitCollection> gemRecHits) {
//   // std::vector<std::vector<int>> n_hits_each_chamber(8, std::vector<int>(36, 0)); // giving max nDigis in one chamber
//   vector<vector<int>> n_hits_each_etaPart(8, vector<int>(576, 0)); // giving max nDigis in one etaPartition

//   for (const GEMRecHit& cluster : *gemRecHits) {
//     // ++n_clusters;
//     GEMDetId hit_id = cluster.gemId();
//     int layer_index = (hit_id.region()+1)/2 + 2*(hit_id.station()-1) + 4*(hit_id.layer()-1);
//     // n_hits = cluster.clusterSize();
//     // n_hits_each_chamber[layer_index][hit_id.chamber() - 1] += n_hits;
//     int eta_index = 16 * (hit_id.chamber() - 1) + hit_id.ieta() - 1;
//     n_hits_each_etaPart[layer_index][eta_index] += cluster.clusterSize();
//   }

//   int mask = 0;
//   int max_val;
//   // for (const auto& row : n_hits_each_chamber) {
//   for (const auto& row : n_hits_each_etaPart) {
//     max_val = *std::max_element(row.begin(), row.end());
//     // if (max_val > 384) return; // big cluster event filter
//     if (max_val > 48) { // big cluster event filter
//       mask = 1;
//       return mask;
//     }
//   }
//   return mask;
// }

int GEMDataClsAnalyzer::maskFlowerEvent(const edm::Handle<TCDSRecord> record) {

  int mask = 0;

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
      mask = 1;
      return mask;
    }
  }
  return mask;
}

int GEMDataClsAnalyzer::maskChamberWithError(const GEMDetId& chamber_id,
                                              const edm::Handle<GEMAMCStatusCollection> amc_status_collection,
                                              const edm::Handle<GEMAMC13StatusCollection> amc13_status_collection,
                                              const edm::Handle<GEMVFATStatusCollection> vfat_status_collection,
                                              const edm::Handle<GEMOHStatusCollection> oh_status_collection) {
  int mask = 0;
  for (auto iter = amc_status_collection->begin(); iter != amc_status_collection->end(); iter++) {
    const auto [amc_id, range] = (*iter);
    if (chamber_id != amc_id) {
      continue;
    }
    for (auto amc_status = range.first; amc_status != range.second; amc_status++) {
      if (amc_status->isBad()) {
        // GEMAMCStatus is bad. Mask this chamber.
        mask = 1;
        return mask;
      }  // isBad
    }  // range
  }  // collection
  for (auto iter = amc13_status_collection->begin(); iter != amc13_status_collection->end(); iter++) {
    const auto [amc13_id, range] = (*iter);
    if (chamber_id != amc13_id) {
      continue;
    }
    for (auto amc13_status = range.first; amc13_status != range.second; amc13_status++) {
      if (amc13_status->isBad()) {
        // GEMAMC13Status is bad. Mask this chamber.
        mask = 1;
        return mask;
      }  // isBad
    }  // range
  }  // collection
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

void GEMDataClsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<GEMAMCStatusCollection> amc_status_collection;
  iEvent.getByToken(amc_status_collection_, amc_status_collection);

  edm::Handle<GEMAMC13StatusCollection> amc13_status_collection;
  iEvent.getByToken(amc13_status_collection_, amc13_status_collection);

  edm::Handle<GEMVFATStatusCollection> vfat_status_collection;
  iEvent.getByToken(vfat_status_collection_, vfat_status_collection);

  edm::Handle<GEMOHStatusCollection> oh_status_collection;
  iEvent.getByToken(oh_status_collection_, oh_status_collection);

  edm::Handle<TCDSRecord> record;
  iEvent.getByToken(tcdsRecord_, record);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByToken(muonHandle_, muonHandle);

  // edm::Handle<OnlineLuminosityRecord> onlineLumiRecord;
  // iEvent.getByToken(onlineLumiRecord_, onlineLumiRecord);


  if (!record.isValid() || !gemRecHits.isValid() || !gemDigis.isValid()
      || !amc_status_collection.isValid() || !amc13_status_collection.isValid()
      || !oh_status_collection.isValid() || !vfat_status_collection.isValid()) {
    cout << "Invalid Event" << endl;
    return;
  }






  int BunchId = iEvent.bunchCrossing();
  int OrbitNumber = iEvent.orbitNumber();
  long EventTime = iEvent.time().unixTime();
  long Event = iEvent.id().event();
  int Lumi = iEvent.luminosityBlock();

  int Flr = maskFlowerEvent(record);

  // float EvNHit = 0;
  // float EvSumCls = 0;
  // for (const GEMRegion* Region : GEMGeometry_->regions()){ // +, -
  //   int re = Region->region();
  //   for (const GEMStation* Station : Region->stations()){ // GE1/1, GE2/1, ME0
  //     int st = Station->station();
  //     if (st != 1) continue; // see only GE1/1 for now
  //     for (const GEMRing* Ring : Station->rings()){
  //       int ri = Ring->ring();
  //       for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){ // GE1/1:2, GE2/1:2, ME0:6
  //         for (const GEMChamber* Chamber : SuperChamber->chambers()){ // GE1/1:36(long/short), GE2/1:18, ME0:18
  //           GEMDetId chId = Chamber->id();
  //           int la = chId.layer();
  //           int ch = chId.chamber();
  //           int ChamErr = maskChamberWithError(chId, vfat_status_collection, oh_status_collection);
  //           for (const GEMEtaPartition* etaPart : Chamber->etaPartitions()){ // GE1/1:8, GE2/1:16, ME0:8
  //             GEMDetId ieId = etaPart->id();
  //             int ie = ieId.ieta();
  //             auto RecHitRange = gemRecHits->get(ieId);
  //             for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
  //               int firstStrip = rechit->firstClusterStrip();
  //               int clsSize = rechit->clusterSize();
  //               b_HitRe = re;
  //               b_HitSt = st;
  //               b_HitRi = ri;
  //               b_HitLa = la;
  //               b_HitCh = ch;
  //               b_HitIe = ie;
  //               b_HitCls = clsSize;
  //               b_HitBunchId = BunchId;
  //               b_HitOrbitNumber = OrbitNumber;
  //               b_HitEventTime = EventTime;
  //               b_HitEvent = Event;
  //               b_HitLumi = Lumi;
  //               b_HitFlr = Flr;
  //               b_HitChamErr = ChamErr;
  //               t_Hit->Fill();
  //               EvNHit++;
  //               EvSumCls += clsSize;
  //             }
  //           } // eta partition loop
  //         } // chamber loop
  //       } // super chamber loop
  //     } // ring loop
  //   } // station loop
  // } // region loop
  // if (EvNHit != 0){
  //   b_EvNHit = EvNHit;
  //   b_EvAvgCls = EvSumCls/EvNHit;
  //   t_Event->Fill();
  // }
  // b_ToNEv = 1;
  // b_ToBunchId = BunchId;
  // b_ToOrbitNumber = OrbitNumber;
  // b_ToEventTime = EventTime;
  // b_ToEvent = Event;
  // b_ToLumi = Lumi;
  // // b_instLumi = onlineLumiRecord->instLumi();
  b_EvFlr = Flr;
  // // for (size_t i = 0; i < max_trigger; ++i) {
  // //   b_ToTrigBunchId[i] = record->getL1aHistoryEntry(i).getBXID();
  // //   b_ToTrigOrbitNumber[i] = record->getL1aHistoryEntry(i).getOrbitNr();
  // // }
  // // b_ToNHit = gemRecHits->size();
  // // int multiplicity = 0;
  // // for (auto rechit : *(gemRecHits.product())){
  // //   multiplicity += rechit.clusterSize();
  // // }
  // // b_ToMul = multiplicity;

  t_Event->Fill();


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
  // if (auto handle = iEvent.getHandle(kMuonViewToken_)) {
  //   cout << "nMuon: " << handle->size() << endl;
  //   muon_view = handle.product();
  // } else {
  //   edm::LogError("ReadoutError") << "failed to get View<Muon>";
  //   return;
  // }

  // const GEMOHStatusCollection* oh_status_collectionn = nullptr;
  // const GEMVFATStatusCollection* vfat_status_collectionn = nullptr;
  // if (true) {
  //   if (auto handle = iEvent.getHandle(oh_status_collection_)) {
  //     oh_status_collectionn = handle.product();
  //   } else {
  //     edm::LogError("ReadoutError") << "failed to get OHVFATStatusCollection";
  //     return;
  //   }

  //   if (auto handle = iEvent.getHandle(vfat_status_collection_)) {
  //     vfat_status_collectionn = handle.product();
  //   } else {
  //     edm::LogError("ReadoutError") << "failed to get GEMVFATStatusCollection";
  //     return;
  //   }
  // }


  //////////////////////////////////////////////////////////////////////////////
  //  Main loop
  //////////////////////////////////////////////////////////////////////////////
  for (const reco::Muon& muon : *(muonHandle.product())) {
    double MuP = muon.p();
    double MuEnergy = muon.energy();
    double MuMass = muon.mass();
    double MuPt = muon.pt(); // [GeV]
    double MuPhi = muon.phi();
    double MuTheta = muon.theta();
    double MuEta = muon.eta();
    int MuGlo = muon.isGlobalMuon();
    int MuTight = muon.passed(reco::Muon::CutBasedIdTight);

    b_MuFlr = Flr;
    b_MuBunchId = BunchId;

    b_MuGE11 = 0;

    b_MuP = MuP;
    b_MuEnergy = MuEnergy;
    b_MuMass = MuMass;
    b_MuPt = MuPt;
    b_MuPhi = MuPhi;
    b_MuTheta = MuTheta;
    b_MuEta = MuEta;
    b_MuGlo = (MuGlo) ? 1 : 0;
    b_MuTight = (MuTight) ? 1 : 0;

    const reco::Track* muonTrack = 0;
    if ( muon.globalTrack().isNonnull() ) muonTrack = muon.globalTrack().get();
    else if ( muon.outerTrack().isNonnull() ) muonTrack = muon.outerTrack().get();
    if (!muonTrack) continue;
    // b_nMuons++;
    // nGEMHitInMuontrack = 0;
    for (auto hit = muonTrack->recHitsBegin(); hit != muonTrack->recHitsEnd(); hit++) {
      const DetId id = (*hit)->geographicalId();
      if (id.det() == 2 /* Muon */ && id.subdetId() == 4 /* GEM */) {
        GEMDetId gemid = id.rawId();

        if (gemid.station() != 1) continue;
        else b_MuGE11 = 1;

        GEMDetId chId = gemid.chamberId();
        // cout << "(ex) chId: " << chId << endl;

        const int ChamErr = maskChamberWithError(chId, amc_status_collection, amc13_status_collection, vfat_status_collection, oh_status_collection);
        
        auto DigiRange = gemDigis->get(gemid);
        for (auto digi = DigiRange.first; digi != DigiRange.second; ++digi) {
          const int re = gemid.region();
          const int st = gemid.station();
          const int ri = gemid.ring();
          const int la = gemid.layer();
          const int ch = gemid.chamber();
          const int ie = gemid.ieta();

          int stp = digi->strip();
          int bx = digi->bx();

          b_MuDigiFlr = Flr;

          b_MuDigiChamErr = ChamErr;

          b_MuDigiGlo = (MuGlo) ? 1 : 0;
          b_MuDigiTight = (MuTight) ? 1 : 0;
          b_MuDigiPt = MuPt;

          b_MuDigiRe = re;
          b_MuDigiSt = st;
          b_MuDigiRi = ri;
          b_MuDigiLa = la;
          b_MuDigiCh = ch;
          b_MuDigiIe = ie;
          b_MuDigiStp = stp;

          t_MuonDigi->Fill();
        }
        
        auto RecHitRange = gemRecHits->get(gemid);
        for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
          int cls = rechit->clusterSize();

          const int re = gemid.region();
          const int st = gemid.station();
          const int ri = gemid.ring();
          const int la = gemid.layer();
          const int ch = gemid.chamber();
          const int ie = gemid.ieta();
          const int ChamErr = maskChamberWithError(chId, amc_status_collection, amc13_status_collection, vfat_status_collection, oh_status_collection);

          b_MuHitLumi = Lumi;
          b_MuHitBunchId = BunchId;
          b_MuHitFlr = Flr;

          b_MuHitChamErr = ChamErr;

          b_MuHitGlo = (MuGlo) ? 1 : 0;
          b_MuHitTight = (MuTight) ? 1 : 0;
          b_MuHitPt = MuPt;

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

void GEMDataClsAnalyzer::beginJob(){}
void GEMDataClsAnalyzer::endJob(){}

void GEMDataClsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
}
// void GEMDataClsAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
//   /* GEM Geometry */
//   edm::ESHandle<GEMGeometry> hGEMGeom;
//   hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

// //  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
//   // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

//   // h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);

// }
void GEMDataClsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMDataClsAnalyzer);
