#ifndef GEMHitAnalyzer_RecHitSimHitForSim_H
#define GEMHitAnalyzer_RecHitSimHitForSim_H
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
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
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
typedef tuple<int, int> Key2;
typedef tuple<int, int, int> Key3;


class GEMHitAnalyzer_RecHitSimHitForSim : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMHitAnalyzer_RecHitSimHitForSim(const edm::ParameterSet&);
  ~GEMHitAnalyzer_RecHitSimHitForSim();

private:
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
  // edm::EDGetTokenT<edm::SimTrackContainer> gemSimTrack_;
  // edm::EDGetTokenT<GEMDigiSimLink> gemDigiSimLink_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeom_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> hGEMGeomBeginRun_;

  // map<Key4, TH2D*> rechit_occ_;
  // map<Key4, TH2D*> digi_occ_;

  // TH1F* strip_; TH1F* pitch_; TH1I* pid_; TH1I* process_; TH1F* xyz_; TH1F* xy_; TH1F* z_; 

  // map<Key1, TH1F*> exp_cls_;
  // map<Key1, TH2F*> exp_cls_vs_xyz_;

  //  map<Key3, TH2F*> xyz_vs_eloss_MeV_;
  //  map<Key3, TH2F*> xyz_vs_eloss_keV_;
  //  map<Key3, TH2F*> xyz_vs_pabs_GeV_;
  //  map<Key3, TH2F*> xyz_vs_pabs_MeV_;
  //  map<Key3, TH2F*> exp_cls_vs_eloss_MeV_;
  //  map<Key3, TH2F*> exp_cls_vs_eloss_keV_;
  //  map<Key3, TH2F*> exp_cls_vs_pabs_GeV_;
  //  map<Key3, TH2F*> exp_cls_vs_pabs_MeV_;
  //  map<Key3, TH2F*> eloss_MeV_vs_pabs_GeV_;
  //  map<Key3, TH2F*> eloss_keV_vs_pabs_MeV_;

  map<Key2, TH2F*> avg_exp_cls_vs_nsimhit_;
  TH2I* eta_vs_nsimhit_;
  TH2I* eta_vs_nrechit_;
  TH2F* eta_vs_simrec_ratio_;
  TH1F* cls_;
  TH1F* exp_cls_;
};

GEMHitAnalyzer_RecHitSimHitForSim::GEMHitAnalyzer_RecHitSimHitForSim(const edm::ParameterSet& iConfig)
  : hGEMGeom_(esConsumes()),
    hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
{
  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  gemSimHits_ = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("gemSimHitLabel"));
  // gemSimTrack_ = consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("gemSimTrackLabel"));
  // gemDigiSimLink_ = consumes<GEMDigiSimLink>(iConfig.getParameter<edm::InputTag>("gemDigiSimLinkLabel"));

//  hGEMGeomBegin_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
//  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 
}

#endif


GEMHitAnalyzer_RecHitSimHitForSim::~GEMHitAnalyzer_RecHitSimHitForSim(){}

void
GEMHitAnalyzer_RecHitSimHitForSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeom_);
  // iSetup.getByToken(hGEMGeom_, hGEMGeom);

  // edm::ESHandle<GEMGeometry> hGEMGeom;
  // iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;
  const GEMGeometry* gem = hGEMGeom.product();

  // edm::Handle<GEMDigiCollection> gemDigis;
  // iEvent.getByToken(gemDigis_, gemDigis);

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<edm::PSimHitContainer> gemSimHits;
  iEvent.getByToken(gemSimHits_, gemSimHits);

  // edm::Handle<edm::SimTrackContainer> gemSimTrack;
  // iEvent.getByToken(gemSimTrack_, gemSimTrack);

  // edm::Handle<GEMDigiSimLink> gemDigiSimLink;
  // iEvent.getByToken(gemDigiSimLink_, gemDigiSimLink);

  // int nrechit_st1 = 0;
  // int nrechit_st2 = 0;
  // int nrechit_st0 = 0;

  int nrechit[8] = {0,0,0,0,0,0,0,0};
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

              // auto SimHitRange = gemSimHits->get(etaPartId);
              // for (auto simhit = gemSimHits->begin(); simhit != gemSimHits->end(); ++simhit) {
              //   auto pid = simhit.particleType();
              //   cout << "pid: " << pid << endl;
              // }
              int n = 0;
              auto RecHitRange = gemRecHits->get(etaPartId);
              for (auto rechit = RecHitRange.first; rechit != RecHitRange.second; ++rechit) {
                int firstStrip = rechit->firstClusterStrip();
                int clsSize = rechit->clusterSize();
                cls_->Fill(clsSize);
                n++;
                // cout << "clsSize: " << clsSize << "(" << re << "/" << st << "/" << la << "/" << ch << "/" << ie << ")" << endl;
              }
              nrechit[ie-1] += n;
            } // eta partition loop
          } // chamber loop
        } // super chamber loop
      } // ring loop
    } // station loop
  } // region loop


  float exp_cls_sum[8] = {0,0,0,0,0,0,0,0};
  float nsimhit[8] = {0,0,0,0,0,0,0,0};
  float avg_exp_cls[8] = {0,0,0,0,0,0,0,0};
  for (const auto& simhit : *gemSimHits.product()) {
    GEMDetId simhit_gemid(simhit.detUnitId());
    const GEMEtaPartition* etapart = gem->etaPartition(simhit_gemid);
    auto lp = simhit.localPosition();
    auto strip = etapart->strip(lp);
    auto pitch = etapart->localPitch(lp); // cm
    //  auto nstrips = etapart->nstrips(); // 384
    cout << "pitch: " << pitch << ", strip: " << strip << endl;

    auto pid = simhit.particleType();
    auto eloss = simhit.energyLoss()*1E9; // GeV -> eV
    auto pabs = simhit.pabs()*1E9; // GeV -> eV
    auto process = simhit.processType();
    cout << "pid: " << pid << ", eloss: " << eloss << ", pabs: " << pabs << ", process: " << process << endl;
    
    auto path = simhit.entryPoint()-simhit.exitPoint();
    auto xyz = sqrt(pow(path.x(), 2) + pow(path.y(), 2) + pow(path.z(), 2)); // cm
    auto xy = sqrt(pow(path.x(), 2) + pow(path.y(), 2)); // cm
    auto z = path.z(); // cm [sensitive detector is at 0.2975 cm]
    auto exp_cls = xy/pitch+1; // xy/strip pitch and (+1) for histogram matching
    cout << "xyz: " << xyz << ", xy: " << xy << ", z: " << z << ", exp cls: " << exp_cls << endl;

    // if (pid != 11) continue;

    exp_cls_->Fill(exp_cls);

    GEMDetId etaPartId = etapart->id();

    auto ie = etaPartId.ieta();
    cout << "ieta: " << ie << endl;

    exp_cls_sum[ie-1] += exp_cls;
    nsimhit[ie-1] += 1;  
  }

  for (int ie=0; ie<8; ie++){
    if (nsimhit[ie]==0) continue;
    avg_exp_cls[ie] = exp_cls_sum[ie]/nsimhit[ie];
    cout << "avg exp. cls.: " << avg_exp_cls[ie] << ", nsimhit: " << nsimhit[ie] << ", eta: " << ie << endl;
    
    eta_vs_nsimhit_->Fill(ie,nsimhit[ie]);
    eta_vs_nrechit_->Fill(ie,nrechit[ie]);
    eta_vs_simrec_ratio_->Fill(ie,nrechit[ie]/nsimhit[ie]);
    // Key1 key1(ie);
    // avg_exp_cls_vs_nsimhit_[key1]->Fill(avg_exp_cls[ie],nsimhit[ie]);
    if (nrechit[ie] > 3) nrechit[ie] = 3;
    Key2 key2(ie,nrechit[ie]);
    avg_exp_cls_vs_nsimhit_[key2]->Fill(avg_exp_cls[ie],nsimhit[ie]);
  }
}

void GEMHitAnalyzer_RecHitSimHitForSim::beginJob(){}
void GEMHitAnalyzer_RecHitSimHitForSim::endJob(){}

void GEMHitAnalyzer_RecHitSimHitForSim::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

//  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  cls_ = fs->make<TH1F>(Form("cls_"), Form("cls."), 385, 0, 385);
  exp_cls_ = fs->make<TH1F>(Form("exp_cls_"), Form("exp. cls."), 385, 0, 385);

  eta_vs_nsimhit_ = fs->make<TH2I>(Form("eta_vs_nsimhit_"), Form("Eta vs nsimhit"), 8, 0, 8, 100, 0, 100);
  eta_vs_nrechit_ = fs->make<TH2I>(Form("eta_vs_nrechit_"), Form("Eta vs nrechit"), 8, 0, 8, 100, 0, 100);
  eta_vs_simrec_ratio_ = fs->make<TH2F>(Form("eta_vs_simrec_ratio_"), Form("Eta vs rec/sim"), 8, 0, 8, 10, 0, 1);

  // for (int ie=0; ie<8; ie++){
  //   Key1 key1(ie);
  //   avg_exp_cls_vs_nsimhit_[key1] = fs->make<TH2F>(Form("avg_exp_cls_vs_nsimhit_ie(%d)_", ie),
  //                                                   Form("average exp. cls. vs nsimhit [ie(%d)]", ie),
  //                                                   385, 0, 385, 100, 0, 100);
  // }

  for (int nrechit=0; nrechit<4; nrechit++){
    for (int ie=0; ie<8; ie++){
      Key2 key2(ie, nrechit);
      avg_exp_cls_vs_nsimhit_[key2] = fs->make<TH2F>(Form("avg_exp_cls_vs_nsimhit_ie(%d)_nrechit(%d)", ie, nrechit),
                                                      Form("average exp. cls. vs nsimhit [ie(%d), nrechit(%d)]", ie, nrechit),
                                                      385, 0, 385, 100, 0, 100);
    }
  }
}

void GEMHitAnalyzer_RecHitSimHitForSim::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMHitAnalyzer_RecHitSimHitForSim);
