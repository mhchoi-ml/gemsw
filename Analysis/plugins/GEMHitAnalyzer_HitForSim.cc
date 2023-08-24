#ifndef GEMHitAnalyzer_HitForSim_H
#define GEMHitAnalyzer_HitForSim_H
// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <map>
#include <algorithm>
#include <typeinfo>

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
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
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

class GEMHitAnalyzer_HitForSim : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer_HitForSim(const edm::ParameterSet&);
  ~GEMHitAnalyzer_HitForSim();

private:
  int track2vertex(const edm::Handle<edm::SimTrackContainer>, int);
  int vertex2parent(const edm::Handle<edm::SimVertexContainer>, int);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHits_;
  edm::EDGetTokenT<edm::SimTrackContainer> gemSimTrack_;
  edm::EDGetTokenT<edm::SimVertexContainer> gemSimVertex_;
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

  TTree* t_SimHit;
  float b_SimHitPitch, b_SimHitPid, b_SimHitProcess, b_SimHitEloss, b_SimHitP;
  float b_SimHitLength, b_SimHitWidth, b_SimHitHeight;
  float b_SimHitCls, b_SimHitElosscut;
};

GEMHitAnalyzer_HitForSim::GEMHitAnalyzer_HitForSim(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  gemSimHits_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("gemSimHitLabel"));
  gemSimTrack_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("gemSimTrackLabel"));
  gemSimVertex_ = consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("gemSimVertexLabel"));

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

  t_SimHit = fs->make<TTree>("SimHit", "gem_simhits_per_hit");
  #define SimHitBRANCH(name, suffix) t_SimHit->Branch(#name, & b_##name, #name "/" #suffix);
  SimHitBRANCH(SimHitPitch, F);
  SimHitBRANCH(SimHitPid, F);
  SimHitBRANCH(SimHitProcess, F);
  SimHitBRANCH(SimHitEloss, F);
  SimHitBRANCH(SimHitP, F);
  SimHitBRANCH(SimHitLength, F);
  SimHitBRANCH(SimHitWidth, F);
  SimHitBRANCH(SimHitHeight, F);
  SimHitBRANCH(SimHitCls, F);
  SimHitBRANCH(SimHitElosscut, F);
}

#endif


GEMHitAnalyzer_HitForSim::~GEMHitAnalyzer_HitForSim(){}

int GEMHitAnalyzer_HitForSim::track2vertex(const edm::Handle<edm::SimTrackContainer> gemSimTrack,
                                                int trkid) {
  for (const auto& simtrack : *gemSimTrack.product()) {

    if (trkid == (int)simtrack.trackId()) {
      cout << "<In simtrack> " << "trkid: " << simtrack.trackId() << ", pid: " << simtrack.type() << ", vertidx: " << simtrack.vertIndex() << endl;
      return simtrack.vertIndex();
    }
  }
  return 0;
}

int GEMHitAnalyzer_HitForSim::vertex2parent(const edm::Handle<edm::SimVertexContainer> gemSimVertex,
                                                int vertid) {
  for (const auto& simvertex : *gemSimVertex.product()) {
    if (vertid == (int)simvertex.vertexId()) {
        cout << "<In simvertex> " << "vertid: " << simvertex.vertexId() << ", prcs: " << simvertex.processType() << ", parentidx: " << simvertex.parentIndex() << endl;
      return simvertex.parentIndex();
    }
  }
  return 0;
}

void
GEMHitAnalyzer_HitForSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
//  iSetup.getByToken(hGEMGeom_, hGEMGeom);

//  edm::ESHandle<GEMGeometry> hGEMGeom;
//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  const GEMGeometry* gem = hGEMGeom.product();

  // edm::Handle<GEMDigiCollection> gemDigis;
  // iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<edm::PSimHitContainer> gemSimHits;
  iEvent.getByToken(gemSimHits_, gemSimHits);

  edm::Handle<edm::SimTrackContainer> gemSimTrack;
  iEvent.getByToken(gemSimTrack_, gemSimTrack);

  edm::Handle<edm::SimVertexContainer> gemSimVertex;
  iEvent.getByToken(gemSimVertex_, gemSimVertex);

  float RecEvNrechits = 0;
  float RecEvSumCls = 0;
  for (const GEMRegion* Region : GEMGeometry_->regions()){
      for (const GEMStation* Station : Region->stations()){
        int st = Station->station();
        if (st != 1) continue;
        for (const GEMRing* Ring : Station->rings()){
          for (const GEMSuperChamber* SuperChamber : Ring->superChambers()){
            for (const GEMChamber* Chamber : SuperChamber->chambers()){

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

  for (const auto& simhit : *gemSimHits.product()) {
    GEMDetId simhit_gemid(simhit.detUnitId());
    const GEMEtaPartition* etapart = gem->etaPartition(simhit_gemid);
    auto lp = simhit.localPosition();
    // auto strip = etapart->strip(lp);
    auto pitch = etapart->localPitch(lp); // cm
    int trkid = simhit.trackId();
    cout << "trkid: " << trkid << " " << typeid(trkid).name() << endl;
    auto pid = simhit.particleType();
    cout << "pid: " << pid << endl;
    auto eloss = simhit.energyLoss(); // GeV
    auto pabs = simhit.pabs(); // GeV
    auto process = simhit.processType();
    
    auto path = simhit.entryPoint()-simhit.exitPoint();
    auto length = sqrt(pow(path.x(), 2) + pow(path.y(), 2) + pow(path.z(), 2)); // cm
    auto width = sqrt(pow(path.x(), 2) + pow(path.y(), 2)); // cm
    auto height = path.z(); // cm [sensitive detector is at 0.2975 cm]
    auto exp_cls = width/pitch+1; // width/strip pitch and (+1) for histogram matching

    int v = track2vertex(gemSimTrack, trkid);
    int p = vertex2parent(gemSimVertex, v);
    while (p != -1) {
      v = track2vertex(gemSimTrack, p);
      p = vertex2parent(gemSimVertex, v);
    }


    // for (const auto& simtrack : *gemSimTrack.product()) {
    //   if (trkid == simtrack.trackId()) {
    //     cout << "addr?: " << simtrack << ", " << &simtrack << ", " << (&simtrack).vertIndex() << endl;
    //   }
      // if (trkid == simtrack.trackId()) {
      //   cout << simtrack.trackId() << ": " << *(simtrack.trackId()) << ", " << &(simtrack.trackId()) << endl;
      //   break;
        // for (const auto& simvertex : *gemSimVertex.product()) {
        //   if (simtrack.vertIndex() == simvertex.vertexId()) {
        //     break;
        //     simvertex.parentIndex()
        //   }
        // }
      // }
    // }

    // cout << "find: " << find(simtrack.trackId().begin(), simtrack.trackId().end(), trkid) << endl;
    // cout << "addr: " << *simtrack << endl;

    // for (const auto& simtrack : *gemSimTrack.product()) {
    //   if (trkid == simtrack.trackId()) {
    //     for (const auto& simvertex : *gemSimVertex.product()) {
    //       if (simtrack.vertIndex() == simvertex.vertexId()) {
    //         break;
    //         simvertex.parentIndex()
    //       }
    //     }
    //   }
    //   // auto vert = simtrack.vertIndex();
    //   // auto genpart = simtrack.genpartIndex();
    //   // auto pid = simtrack.type();
    //   // cout << "trackid: " << trackid << ", vert: " << vert << ", genpart: " << genpart << ", pid: " << pid << endl;
    // }
    // for (const auto& simvertex : *gemSimVertex.product()) {
    //   auto parent = simvertex.parentIndex();
    //   auto vid = simvertex.vertexId();
    //   auto prcs = simvertex.processType();
    //   cout << "parent: " << parent << ", vid: " << vid << ", prcs: " << prcs << endl;
    // }


    // const auto& simtrack = *gemSimTrack.product();
    // auto trackid_ = simtrack[trackid-1].trackId();
    // auto vert = simtrack[trackid-1].vertIndex();
    // auto genpart = simtrack[trackid-1].genpartIndex();
    // auto type = simtrack[trackid-1].type();
    // cout << "trackId: " << trackid_ << ", vert: " << vert << ", genpart: " << genpart << ", type: " << type << endl;


    // sort pid/process format
    if (pid == 22) pid = 0;          // photon
    else if (pid == 11) pid = 1;     // electron
    else if (pid == -11) pid = 2;    // positron
    else if (pid == 13) pid = 3;     // muon(-)
    else if (pid == -13) pid = 4;    // muon(+)
    else if (pid == 211) pid = 5;    // pion(+)
    else if (pid == -211) pid = 6;   // pion(-)
    else if (pid == 321) pid = 7;    // kaon(+)
    else if (pid == -321) pid = 8;   // kaon(-)
    else if (pid == 2112) pid = 9;   // proton
    else if (pid == 2212) pid = 10;  // neutron
    else pid = 11;                   // the other particles

    if (process == 0) process = 0;                                                                // Primary generator
    else if (process == 91) process = 1;                                                          // Transportation
    else if (process == 92) process = 2;                                                          // CoupleTrans
    else if (process == 1) process = 3;                                                           // CoulombScat
    else if (process == 2) process = 4;                                                           // Ionisation
    else if (process == 3) process = 5;                                                           // Brems
    else if (process == 4) process = 6;                                                           // PairProdCharged
    else if (process == 5) process = 7;                                                           // Annih
    else if (process == 6) process = 8;                                                           // AnnihToMuMu
    else if (process == 7) process = 9;                                                           // AnnihToHad
    else if (process == 8) process = 10;                                                          // NuclearStopp
    else if (process == 10) process = 11;                                                         // Msc
    else if (process == 11) process = 12;                                                         // Rayleigh
    else if (process == 12) process = 13;                                                         // PhotoElectric
    else if (process == 13) process = 14;                                                         // Compton
    else if (process == 14) process = 15;                                                         // Conv
    else if (process == 15) process = 16;                                                         // ConvToMuMu
    else if ((process >= 21 && process <= 24) || (process >= 31 && process <= 35)) process = 17;  // Undefined1
    else if (process == 40) process = 18;                                                         // muDBrem
    else if (process == 49) process = 19;                                                         // MuMuonPairProd
    else if (process >= 53 && process <= 57) process = 20;                                        // Undefined2
    else if (process == 111) process = 21;                                                        // HadElastic
    else if (process == 121) process = 22;                                                        // HadInelastic
    else if (process == 131) process = 23;                                                        // HadCapture
    else if (process == 141) process = 24;                                                        // HadFission
    else if (process == 151) process = 25;                                                        // HadAtRest
    else if (process == 161) process = 26;                                                        // Undefined3
    else if (process == 201) process = 27;                                                        // Decay
    else if (process == 202) process = 28;                                                        // DecayWSpin
    else if (process == 203) process = 29;                                                        // DecayPiWSpin
    else if (process == 210) process = 30;                                                        // DecayRadio
    else if (process == 211) process = 31;                                                        // DecayUnKnown
    else if (process == 231) process = 32;                                                        // DecayExt
    else if (process == 301) process = 33;                                                        // GFlash
    else if (process == 401) process = 34;                                                        // StepLimiter
    else if (process == 402) process = 35;                                                        // Undefined4
    else if (process == 403) process = 36;                                                        // NeutronKiller
    else process = 37;  

    b_SimHitPitch = pitch;
    b_SimHitPid = pid;
    b_SimHitProcess = process;
    b_SimHitEloss = eloss;
    b_SimHitP = pabs;
    b_SimHitLength = length;
    b_SimHitWidth = width;
    b_SimHitHeight = height;
    b_SimHitCls = exp_cls;

    // energy loss cut fo ionization minimum threshold 
    if (eloss*1E9 > 28.1) b_SimHitElosscut = 1;
    else b_SimHitElosscut = 0;

    t_SimHit->Fill();
  }


  // for ( auto ivtx=SimVtx->begin(); ivtx!=SimVtx->end(); ++ivtx )
  // {
  // }
  // for ( auto itk=gemSimTrack->begin(); itk!=gemSimTrack->end(); ++itk )
  // {
  // }

  // for (const auto& simvertex : *gemSimVertex.product()) {
  //   auto parent = simvertex.parentIndex();
  //   auto vid = simvertex.vertexId();
  //   auto prcs = simvertex.processType();
  //   cout << "parent: " << parent << ", vid: " << vid << ", prcs: " << prcs << endl;
  // }
  // for (const auto& simtrack : *gemSimTrack.product()) {
  //   // GEMDetId simhit_gemid(simtrack.detUnitId());
  //   // const GEMEtaPartition* etapart = gem->etaPartition(simhit_gemid);
  //   auto trackid = simtrack.trackId();
  //   auto vert = simtrack.vertIndex();
  //   auto genpart = simtrack.genpartIndex();
  //   auto pid = simtrack.type();
  //   cout << "trackid: " << trackid << ", vert: " << vert << ", genpart: " << genpart << ", pid: " << pid << endl;
  // }

  h_nEvents->Fill(1);
}

void GEMHitAnalyzer_HitForSim::beginJob(){}
void GEMHitAnalyzer_HitForSim::endJob(){}

void GEMHitAnalyzer_HitForSim::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);

}
void GEMHitAnalyzer_HitForSim::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer_HitForSim);