#ifndef GEMHitAnalyzer_SimHitForSim_H
#define GEMHitAnalyzer_SimHitForSim_H
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
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
// #include "SimDataFormats/Track/interface/SimTrackContainer.h"
// #include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
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

// #include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
typedef tuple<int> Key1;
typedef tuple<int, int, int> Key3;


class GEMHitAnalyzer_SimHitForSim : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer_SimHitForSim(const edm::ParameterSet&);
  ~GEMHitAnalyzer_SimHitForSim();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMDigiCollection> gemDigis_;
  edm::EDGetTokenT<edm::PSimHitContainer> gemSimHits_;
  // edm::EDGetTokenT<edm::SimTrackContainer> gemSimTrack_;
  // edm::EDGetTokenT<GEMDigiSimLink> gemDigiSimLink_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  TH1F* strip_; TH1F* pitch_; TH1I* pid_; TH1I* process_; TH1F* length_; TH1F* width_; TH1F* height_; 

  map<Key1, TH1F*> exp_cls_;
  map<Key1, TH2F*> exp_cls_vs_length_;

  map<Key3, TH2F*> length_vs_eloss_MeV_;
  map<Key3, TH2F*> length_vs_eloss_keV_;
  map<Key3, TH2F*> length_vs_pabs_GeV_;
  map<Key3, TH2F*> length_vs_pabs_MeV_;
  map<Key3, TH2F*> exp_cls_vs_eloss_MeV_;
  map<Key3, TH2F*> exp_cls_vs_eloss_keV_;
  map<Key3, TH2F*> exp_cls_vs_pabs_GeV_;
  map<Key3, TH2F*> exp_cls_vs_pabs_MeV_;
  map<Key3, TH2F*> pabs_GeV_vs_eloss_MeV_;
  map<Key3, TH2F*> pabs_MeV_vs_eloss_keV_;
};

GEMHitAnalyzer_SimHitForSim::GEMHitAnalyzer_SimHitForSim(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemSimHits_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("gemSimHitLabel"));
  // gemSimTrack_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("gemSimTrackLabel"));
  // gemDigiSimLink_ = consumes<GEMDigiSimLink>(iConfig.getParameter<edm::InputTag>("gemDigiSimLinkLabel"));

//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
}

#endif


GEMHitAnalyzer_SimHitForSim::~GEMHitAnalyzer_SimHitForSim(){}

void
GEMHitAnalyzer_SimHitForSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  edm::Handle<edm::PSimHitContainer> gemSimHits;
  iEvent.getByToken(gemSimHits_, gemSimHits);

  // edm::Handle<edm::SimTrackContainer> gemSimTrack;
  // iEvent.getByToken(gemSimTrack_, gemSimTrack);

  // edm::Handle<GEMDigiSimLink> gemDigiSimLink;
  // iEvent.getByToken(gemDigiSimLink_, gemDigiSimLink);

  for (const auto& simhit : *gemSimHits.product()) {
    GEMDetId simhit_gemid(simhit.detUnitId());
    const GEMEtaPartition* etapart = gem->etaPartition(simhit_gemid);
    auto lp = simhit.localPosition();
    auto strip = etapart->strip(lp);
    auto pitch = etapart->localPitch(lp); // cm
    //  auto nstrips = etapart->nstrips(); // 384
    std::cout << pitch << ", " << strip << std::endl;

    auto pid = simhit.particleType();
    auto eloss = simhit.energyLoss()*1E9; // GeV -> eV
    auto pabs = simhit.pabs()*1E9; // GeV -> eV
    auto process = simhit.processType();
    std::cout << pid << ", " << eloss << ", " << pabs << ", " << process << std::endl;
    
    auto path = simhit.entryPoint()-simhit.exitPoint();
    auto length = sqrt(pow(path.x(), 2) + pow(path.y(), 2) + pow(path.z(), 2)); // cm
    auto width = sqrt(pow(path.x(), 2) + pow(path.y(), 2)); // cm
    auto height = path.z(); // cm [sensitive detector is at 0.2975 cm]
    auto exp_cls = width/pitch+1; // width/strip pitch and (+1) for histogram matching
    std::cout << length << ", " << width << ", " << height << ", " << exp_cls << "\n" << std::endl;

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
    else process = 37;                                                                            // Unknown

    strip_->Fill(strip); pitch_->Fill(pitch); pid_->Fill(pid); process_->Fill(process); length_->Fill(length); width_->Fill(width); height_->Fill(height);

    if (eloss > 28.1){
      Key1 key1(1);
      Key3 key3(1,pid,process);

      exp_cls_[key1]->Fill(exp_cls);
      exp_cls_vs_length_[key1]->Fill(exp_cls,length);

      length_vs_eloss_MeV_[key3]->Fill(length, eloss);
      length_vs_eloss_keV_[key3]->Fill(length, eloss);
      length_vs_pabs_GeV_[key3]->Fill(length, pabs);
      length_vs_pabs_MeV_[key3]->Fill(length, pabs);
      exp_cls_vs_eloss_MeV_[key3]->Fill(exp_cls, eloss);
      exp_cls_vs_eloss_keV_[key3]->Fill(exp_cls, eloss/1E3);
      exp_cls_vs_pabs_GeV_[key3]->Fill(exp_cls, pabs);
      exp_cls_vs_pabs_MeV_[key3]->Fill(exp_cls, pabs);
      pabs_GeV_vs_eloss_MeV_[key3]->Fill(pabs, eloss);
      pabs_MeV_vs_eloss_keV_[key3]->Fill(pabs, eloss);
    }
    else{
      Key1 key1(0);
      Key3 key3(0,pid,process);

      exp_cls_[key1]->Fill(exp_cls);
      exp_cls_vs_length_[key1]->Fill(exp_cls,length);

      length_vs_eloss_MeV_[key3]->Fill(length, eloss);
      length_vs_eloss_keV_[key3]->Fill(length, eloss);
      length_vs_pabs_GeV_[key3]->Fill(length, pabs);
      length_vs_pabs_MeV_[key3]->Fill(length, pabs);
      exp_cls_vs_eloss_MeV_[key3]->Fill(exp_cls, eloss);
      exp_cls_vs_eloss_keV_[key3]->Fill(exp_cls, eloss);
      exp_cls_vs_pabs_GeV_[key3]->Fill(exp_cls, pabs);
      exp_cls_vs_pabs_MeV_[key3]->Fill(exp_cls, pabs);
      pabs_GeV_vs_eloss_MeV_[key3]->Fill(pabs, eloss);
      pabs_MeV_vs_eloss_keV_[key3]->Fill(pabs, eloss);
    }
  }
}

void GEMHitAnalyzer_SimHitForSim::beginJob(){}
void GEMHitAnalyzer_SimHitForSim::endJob(){}

void GEMHitAnalyzer_SimHitForSim::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  strip_ = fs->make<TH1F>(Form("strip_"), Form("Strip"), 385, 0, 385);
  pitch_ = fs->make<TH1F>(Form("pitch_"), Form("Pitch"), 10, 0, 1); // pitch: [~ 1 cm (1 mm)]
  pid_ = fs->make<TH1I>(Form("pid_"), Form("Particle ID"), 12, 0, 12);
  process_ = fs->make<TH1I>(Form("process_"), Form("Process"), 38, 0, 38);
  length_ = fs->make<TH1F>(Form("length_"), Form("XYZ"), 1000, 0, 100); // length: [~ 100 cm (1 mm)]
  width_ = fs->make<TH1F>(Form("width_"), Form("XY"), 1000, 0, 100); // width: [~ 100 cm (1 mm)]
  height_ = fs->make<TH1F>(Form("height_"), Form("Z"), 10, 0, 1); // height: [~ 1 cm (1 mm)]

  for (int elcut=0; elcut<2; elcut++){
    Key1 key1(elcut);

    exp_cls_[key1] = fs->make<TH1F>(Form("exp_cls_ELcut(%d)_", elcut),
                                    Form("Expected cluster size [ELcut(%d)]", elcut),
                                    385, 0, 385);
    exp_cls_vs_length_[key1] = fs->make<TH2F>(Form("exp_cls_vs_length_ELcut(%d)_", elcut),
					   Form("Expected cluster size [ELcut(%d)]", elcut),
					   385, 0, 385, 1000, 0, 100);

    for (int pid=0; pid<12; pid++){
      for (int prcs=0; prcs<38; prcs++){
        Key3 key3(elcut, pid, prcs);

        length_vs_eloss_MeV_[key3] = fs->make<TH2F>(Form("length_vs_eloss_MeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                 Form("XYZ vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                 1000, 0, 100, 100, 0, 1000000); // length: [~100 cm (1 mm) ], eloss: [~ 1 MeV (10 keV)]
        length_vs_eloss_keV_[key3] = fs->make<TH2F>(Form("length_vs_eloss_keV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                 Form("XYZ vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                 1000, 0, 100, 100, 0, 20000); // length: [~100 cm (1 mm) ], eloss: [~ 20 keV (200 eV)]
        length_vs_pabs_GeV_[key3] = fs->make<TH2F>(Form("length_vs_pabs_GeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                Form("XYZ vs Pabs [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                1000, 0, 100, 100, 0, 1000000000); // length: [~100 cm (1 mm) ], pbas: [~ 1 GeV (10 MeV)]
        length_vs_pabs_MeV_[key3] = fs->make<TH2F>(Form("length_vs_pabs_MeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                Form("XYZ vs Pabs [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                1000, 0, 100, 100, 0, 30000000); // length: [~100 cm (1 mm) ], pbas: [~ 30 MeV (300 keV)]
        exp_cls_vs_eloss_MeV_[key3] = fs->make<TH2F>(Form("exp_cls_vs_eloss_MeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                     Form("Exp. cls. vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                     385, 0, 385, 100, 0, 1000000); // exp_cls: [~385 (1)], eloss: [~ 1 MeV (10 keV)]
        exp_cls_vs_eloss_keV_[key3] = fs->make<TH2F>(Form("exp_cls_vs_eloss_keV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                     Form("Exp. cls. vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                     385, 0, 385, 100, 0, 20); // exp_cls: [~385 (1)], eloss: [~ 20 keV (200 eV)]
        exp_cls_vs_pabs_GeV_[key3] = fs->make<TH2F>(Form("exp_cls_vs_pabs_GeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                    Form("Exp. cls. vs Pabs [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                    385, 0, 385, 100, 0, 1000000000); // exp_cls: [~385 (1)], pbas: [~ 1 GeV (10 MeV)]
        exp_cls_vs_pabs_MeV_[key3] = fs->make<TH2F>(Form("exp_cls_vs_pabs_MeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                    Form("Exp. cls. vs Pabs [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                    385, 0, 385, 100, 0, 30000000); // exp_cls: [~385 (1)], pbas: [~ 30 MeV (300 keV)]
        pabs_GeV_vs_eloss_MeV_[key3] = fs->make<TH2F>(Form("pabs_GeV_vs_eloss_MeV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                      Form("Pabs vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                      100, 0, 1000000, 100, 0, 1000000000); // eloss: [~ 1 MeV (10 keV)], pbas: [~ 1 GeV (10 MeV)]
        pabs_MeV_vs_eloss_keV_[key3] = fs->make<TH2F>(Form("pabs_MeV_vs_eloss_keV_ELcut(%d)_PID(%d)_Prcs(%d)_", elcut, pid, prcs),
                                                      Form("Pabs vs Energy loss [ELcut(%d), PID(%d), Prcs(%d)]", elcut, pid, prcs),
                                                      100, 0, 20000, 100, 0, 30000000); // eloss: [~ 20 keV (200 eV)], pbas: [~ 30 MeV (300 keV)]
      }
    }
  }
}

void GEMHitAnalyzer_SimHitForSim::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer_SimHitForSim);
