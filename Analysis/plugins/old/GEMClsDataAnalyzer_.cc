#ifndef GEMClsDataAnalyzer_H
#define GEMClsDataAnalyzer_H
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

class GEMClsDataAnalyzer : public edm::one::EDAnalyzer<edm::one::WatchRuns> {  
public:
  explicit GEMClsDataAnalyzer(const edm::ParameterSet&);
  ~GEMClsDataAnalyzer();

  // currently only for STA muons
  enum class StartingStateType {
    kOutermostMeasurementState = 0,
    kInnermostMeasurementState,
    kStateOnSurfaceWithCSCSegment,
    kAlignmentStyle,
  };

  // Define the metric as the smaller the absolute value, the better the matching.
  enum class MatchingMetric {
    kDeltaPhi = 0,  // computeDeltaPhi
    kRdPhi,         // computeRdPhi
  };

  // https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre3/Configuration/Applications/python/ConfigBuilder.py#L35
  enum class ScenarioOption {
    kPP = 0,
    kCosmics,
    kHeavyIons,
  };

  struct GEMLayer {
    GEMLayer(Disk::DiskPointer disk, std::vector<const GEMChamber *> chambers, GEMDetId id)
        : disk(disk), chambers(chambers), id(id) {}
    Disk::DiskPointer disk;
    std::vector<const GEMChamber *> chambers;
    GEMDetId id;
  };

  using StartingState = std::tuple<bool, TrajectoryStateOnSurface, DetId>;

private:
  int maskChamberWithError(const GEMDetId& chamber_id, const edm::Handle<GEMVFATStatusCollection>, const edm::Handle<GEMOHStatusCollection>);
  int maskFlowerEvent(const edm::Handle<TCDSRecord> record);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  // edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  // edm::EDGetTokenT<GEMOHStatusCollection> oh_status_collection_;
  // edm::EDGetTokenT<GEMVFATStatusCollection> vfat_status_collection_;
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
  int b_ToBunchId, b_ToOrbitNumber;
  long b_ToEvent, b_ToEventTime;
  int b_ToLumi;
  // float b_EvInstLumi;
  int b_ToFlr;

  TTree* t_Event;
  float b_EvNHit, b_EvAvgCls;

  TTree* t_Hit;
  int b_HitRe, b_HitSt, b_HitRi, b_HitLa, b_HitCh, b_HitIe;
  int b_HitCls;
  int b_HitBunchId, b_HitOrbitNumber;
  long b_HitEvent, b_HitEventTime;
  int b_HitLumi;
  int b_HitFlr;
  int b_HitChamErr;

  TTree* t_Muon;
  int b_MuRe, b_MuSt, b_MuRi, b_MuLa, b_MuCh, b_MuIe;
  int b_MuCls;
  int b_MuLumi;
  int b_MuFlr;
  int b_MuChamErr;


  StartingStateType getStartingStateType(const std::string);
  MatchingMetric getMatchingMetric(const std::string);
  reco::Muon::MuonTrackType getMuonTrackType(const std::string);
  ScenarioOption getScenarioOption(const std::string);

  void buildGEMLayers(const GEMGeometry *);
  bool skipGEMStation(const int);

  bool checkPropagationDirection(const reco::Track *, const GEMLayer &);

  StartingState buildStartingState(const reco::Muon &, const reco::TransientTrack &, const GEMLayer &);
  StartingState getInnermostMeasurementState(const reco::TransientTrack &);
  StartingState getOutermostMeasurementState(const reco::TransientTrack &);
  StartingState buildStateOnSurfaceWithCSCSegment(const reco::Muon &, const reco::TransientTrack &, const GEMLayer &);
  StartingState buildStartingStateAlignmentStyle(const reco::Muon &, const reco::TransientTrack &, const GEMLayer &);

  // for kStateOnSurfaceWithCSCSegment and AlignmentStyle
  const CSCSegment *findCSCSegment(const reco::Muon &, const reco::TransientTrack &, const GEMLayer &);
  const CSCSegment *findCSCSegmentBeam(const reco::TransientTrack &, const GEMLayer &);
  const CSCSegment *findCSCSegmentCosmics(const reco::Muon &, const GEMLayer &);
  bool isMuonSubdetAllowed(const DetId &, const int);
  bool isCSCAllowed(const CSCDetId &, const int);

  bool checkBounds(const Plane &, const GlobalPoint &);
  bool checkBounds(const Plane &, const GlobalPoint &, const GlobalError &, float);
  const GEMEtaPartition *findEtaPartition(const GlobalPoint &,
                                          const GlobalError &,
                                          const std::vector<const GEMChamber *> &);

  float computeRdPhi(const GlobalPoint &, const LocalPoint &, const GEMEtaPartition *);
  float computeDeltaPhi(const GlobalPoint &, const LocalPoint &, const GEMEtaPartition *);
  float computeMatchingMetric(const GlobalPoint &, const LocalPoint &, const GEMEtaPartition *);

  // std::pair<const GEMRecHit *, float> findClosestHit(const GlobalPoint &,
  //                                                    const GEMRecHitCollection::range &,
  //                                                    const GEMEtaPartition *);
  std::tuple<const GEMRecHit *, float, int> findClosestHit(const GlobalPoint &,
                                                     const GEMRecHitCollection::range &,
                                                     const GEMEtaPartition *);

  // some helpers
  inline bool isInsideOut(const reco::Track &);

  //////////////////////////////////////////////////////////////////////////////
  // const data members initialized in the member initializer list
  // mainly retrieved from edm::ParameterSet
  //////////////////////////////////////////////////////////////////////////////
  // ES
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> kGEMGeometryTokenBeginRun_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> kTransientTrackBuilderToken_;
  // ED
  const edm::EDGetTokenT<GEMRecHitCollection> kGEMRecHitCollectionToken_;
  const edm::EDGetTokenT<edm::View<reco::Muon> > kMuonViewToken_;
  
  const std::string kMuonTrackTypeName_;
  const reco::Muon::MuonTrackType kMuonTrackType_;
  // const TString kMuonName_;
  // const std::string kFolder_;
  const ScenarioOption kScenario_;
  // cuts
  const StartingStateType kStartingStateType_;
  const std::vector<std::vector<int> > kMuonSubdetForGEM_;
  const std::vector<std::vector<int> > kCSCForGEM_;  // when using StartingStateType::kStateOnSurfaceWithCSCSegment
  const float kMuonSegmentMatchDRCut_;               // for cosmics

  const std::vector<double> kMuonPtMinCuts_;   // station as index
  const std::vector<double> kMuonEtaMinCuts_;  // station as index
  const std::vector<double> kMuonEtaMaxCuts_;  // station as index
  const float kPropagationErrorRCut_;          // cm
  const float kPropagationErrorPhiCut_;        // degree
  const float kBoundsErrorScale_;              // TODO docc
  // matching
  const MatchingMetric kMatchingMetric_;
  const float kMatchingCut_;
  // // for MinotorElement
  // const std::vector<double> kMuonPtBins_;  // station as index
  // const std::vector<int> kMuonEtaNbins_;   // station as index
  // const std::vector<double> kMuonEtaLow_;  // station as index
  // const std::vector<double> kMuonEtaUp_;   // station as index

  const bool kMaskChamberWithError_;

  // const
  const bool kModeDev_;

  const bool kMonitorGE11_;
  const bool kMonitorGE21_;
  const bool kMonitorGE0_;

  //////////////////////////////////////////////////////////////////////////////
  // const data members
  // FIXME static?
  //////////////////////////////////////////////////////////////////////////////
  // https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre3/DataFormats/CSCRecHit/interface/CSCSegment.h#L60
  const int kCSCSegmentDimension_ = 4;

  //////////////////////////////////////////////////////////////////////////////
  // non-const data members
  //////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<MuonServiceProxy> muon_service_;
  std::vector<GEMLayer> gem_layers_;
};

GEMClsDataAnalyzer::GEMClsDataAnalyzer(const edm::ParameterSet& iConfig)
  // : hGEMGeom_(esConsumes()),
    // hGEMGeomBeginRun_(esConsumes<edm::Transition::BeginRun>())
    : kGEMGeometryTokenBeginRun_(esConsumes<edm::Transition::BeginRun>()),
    kTransientTrackBuilderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
    kGEMRecHitCollectionToken_(consumes<GEMRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("recHitTag"))),
    kMuonViewToken_(consumes<edm::View<reco::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))),
    kMuonTrackTypeName_(iConfig.getUntrackedParameter<std::string>("muonTrackType")),
    kMuonTrackType_(getMuonTrackType(kMuonTrackTypeName_)),
    kScenario_(getScenarioOption(iConfig.getUntrackedParameter<std::string>("scenario"))),
    kStartingStateType_(getStartingStateType(iConfig.getUntrackedParameter<std::string>("startingStateType"))),
    kMuonSubdetForGEM_({
        iConfig.getUntrackedParameter<std::vector<int> >("muonSubdetForGE0"),
        iConfig.getUntrackedParameter<std::vector<int> >("muonSubdetForGE11"),
        iConfig.getUntrackedParameter<std::vector<int> >("muonSubdetForGE21"),
    }),
    kCSCForGEM_({
        iConfig.getUntrackedParameter<std::vector<int> >("cscForGE0"),
        iConfig.getUntrackedParameter<std::vector<int> >("cscForGE11"),
        iConfig.getUntrackedParameter<std::vector<int> >("cscForGE21"),
    }),
    kMuonSegmentMatchDRCut_(static_cast<float>(iConfig.getUntrackedParameter<double>("muonSegmentMatchDRCut"))),
    kMuonPtMinCuts_({
      iConfig.getUntrackedParameter<double>("muonPtMinCutGE0"),
      iConfig.getUntrackedParameter<double>("muonPtMinCutGE11"),
      iConfig.getUntrackedParameter<double>("muonPtMinCutGE21"),
    }),
    kMuonEtaMinCuts_({
      iConfig.getUntrackedParameter<double>("muonEtaMinCutGE0"),
      iConfig.getUntrackedParameter<double>("muonEtaMinCutGE11"),
      iConfig.getUntrackedParameter<double>("muonEtaMinCutGE21"),
    }),
    kMuonEtaMaxCuts_({
      iConfig.getUntrackedParameter<double>("muonEtaMaxCutGE0"),
      iConfig.getUntrackedParameter<double>("muonEtaMaxCutGE11"),
      iConfig.getUntrackedParameter<double>("muonEtaMaxCutGE21"),
    }),
    kPropagationErrorRCut_(static_cast<float>(iConfig.getUntrackedParameter<double>("propagationErrorRCut"))),
    kPropagationErrorPhiCut_(static_cast<float>(iConfig.getUntrackedParameter<double>("propagationErrorPhiCut"))),    
    kBoundsErrorScale_(static_cast<float>(iConfig.getUntrackedParameter<double>("boundsErrorScale"))),
    kMatchingMetric_(getMatchingMetric(iConfig.getUntrackedParameter<std::string>("matchingMetric"))),
    kMatchingCut_(static_cast<float>(iConfig.getUntrackedParameter<double>("matchingCut"))),
    kMaskChamberWithError_(iConfig.getUntrackedParameter<bool>("maskChamberWithError")),
    kModeDev_(iConfig.getUntrackedParameter<bool>("modeDev")),
    kMonitorGE11_(iConfig.getUntrackedParameter<bool>("monitorGE11")),
    kMonitorGE21_(iConfig.getUntrackedParameter<bool>("monitorGE21")),
    kMonitorGE0_(iConfig.getUntrackedParameter<bool>("monitorGE0"))
{
  muon_service_ = std::make_unique<MuonServiceProxy>(iConfig.getParameter<edm::ParameterSet>("ServiceParameters"), consumesCollector());


  // gemDigis_ = consumes<GEMDigiCollection>(iConfig.getParameter<edm::InputTag>("gemDigiLabel"));
  // gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  // oh_status_collection_ = consumes<GEMOHStatusCollection>(iConfig.getParameter<edm::InputTag>("OHInputLabel"));
  // vfat_status_collection_ = consumes<GEMVFATStatusCollection>(iConfig.getParameter<edm::InputTag>("VFATInputLabel"));
  tcdsRecord_ = consumes<TCDSRecord>(iConfig.getParameter<edm::InputTag>("tcdsRecord"));
  // onlineLumiRecord_ = consumes<OnlineLuminosityRecord>(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"));
  // muonHandle_    = consumes<edm::View<reco::Muon> >(iConfig.getParameter<edm::InputTag>("muonLabel"));

  hGEMGeom_ = esConsumes<GEMGeometry, MuonGeometryRecord>();
  hGEMGeomBeginRun_ = esConsumes<GEMGeometry, MuonGeometryRecord>(); 

  t_Total = fs->make<TTree>("Total", "event_total_info");
  #define ToBRANCH(name, suffix) t_Total->Branch(#name, & b_##name, #name "/" #suffix);
  ToBRANCH(ToNEv, I);

  ToBRANCH(ToBunchId, I);
  ToBRANCH(ToOrbitNumber, I);
  ToBRANCH(ToEvent, l);
  ToBRANCH(ToEventTime, l);
  ToBRANCH(ToLumi, I);
    // HitBRANCH(instLumi, F);

  ToBRANCH(ToFlr, I);

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

  HitBRANCH(HitBunchId, I);
  HitBRANCH(HitOrbitNumber, I);
  HitBRANCH(HitEvent, l);
  HitBRANCH(HitEventTime, l);
  HitBRANCH(HitLumi, I);

  HitBRANCH(HitChamErr, I);
  HitBRANCH(HitFlr, I);

  t_Muon = fs->make<TTree>("Muon", "gem_hits");
  #define MuBRANCH(name, suffix) t_Muon->Branch(#name, & b_##name, #name "/" #suffix);
  MuBRANCH(MuRe, I);
  MuBRANCH(MuSt, I);
  MuBRANCH(MuRi, I);
  MuBRANCH(MuLa, I);
  MuBRANCH(MuCh, I);
  MuBRANCH(MuIe, I);
  MuBRANCH(MuCls, I);

  MuBRANCH(MuLumi, I);

  MuBRANCH(MuChamErr, I);
  MuBRANCH(MuFlr, I);
}

#endif


GEMClsDataAnalyzer::~GEMClsDataAnalyzer(){}


// convert a string to enum
GEMClsDataAnalyzer::MatchingMetric GEMClsDataAnalyzer::getMatchingMetric(const std::string name) {
  MatchingMetric method;

  if (name == "DeltaPhi") {
    method = MatchingMetric::kDeltaPhi;

  } else if (name == "RdPhi") {
    method = MatchingMetric::kRdPhi;

  } else {
    edm::LogError("ReadoutError") << "received an unexpected MatchingMetric: " << name
                                 << " -> MatchingMetric::kDeltaPhi will be used instead.";
    method = MatchingMetric::kDeltaPhi;
  }

  return method;
}

// convert a string to enum
GEMClsDataAnalyzer::StartingStateType GEMClsDataAnalyzer::getStartingStateType(const std::string name) {
  StartingStateType type;

  if (name == "InnermostMeasurementState") {
    type = StartingStateType::kInnermostMeasurementState;

  } else if (name == "OutermostMeasurementState") {
    type = StartingStateType::kOutermostMeasurementState;

  } else if (name == "StateOnSurfaceWithCSCSegment") {
    type = StartingStateType::kStateOnSurfaceWithCSCSegment;

  } else if (name == "AlignmentStyle") {
    type = StartingStateType::kAlignmentStyle;

  } else {
    edm::LogError("ReadoutError") << "received an unexpected StartingStateType: " << name
                                 << " -> StartingStateType::kOutermostMeasurementState will be used instead.";
    type = StartingStateType::kOutermostMeasurementState;
  }

  return type;
}

// convert a string to enum
reco::Muon::MuonTrackType GEMClsDataAnalyzer::getMuonTrackType(const std::string name) {
  reco::Muon::MuonTrackType muon_track_type;

  // DO NOT ALLOW TYPO
  if (name == "InnerTrack") {
    muon_track_type = reco::Muon::MuonTrackType::InnerTrack;

  } else if (name == "OuterTrack") {
    muon_track_type = reco::Muon::MuonTrackType::OuterTrack;

  } else if (name == "CombinedTrack") {
    muon_track_type = reco::Muon::MuonTrackType::CombinedTrack;

  } else {
    edm::LogError("ReadoutError") << "received an unexpected reco::Muon::MuonTrackType: " << name
                                 << " --> OuterTrack will be used instead.";

    muon_track_type = reco::Muon::MuonTrackType::OuterTrack;
  }

  return muon_track_type;
}

GEMClsDataAnalyzer::ScenarioOption GEMClsDataAnalyzer::getScenarioOption(const std::string name) {
  ScenarioOption scenario;
  if (name == "pp") {
    scenario = ScenarioOption::kPP;

  } else if (name == "cosmics") {
    scenario = ScenarioOption::kCosmics;

  } else if (name == "HeavyIons") {
    scenario = ScenarioOption::kHeavyIons;

    edm::LogInfo("ReadoutError") << "The scenario is set to \"HeavyIons\""
                                << " but there is no strategy dedicated to"
                                << "\"HeavyIons\" scenario. The strategy for "
                                << "the \"pp\" scenario will be used insteqad.";

  } else {
    scenario = ScenarioOption::kPP;

    edm::LogError("ReadoutError") << "received an unexpected ScenarioOption: " << name
                                 << ". Choose from (\"pp\", \"cosmics\", \"HeavyIons\")"
                                 << " --> pp will be used instead.";
  }

  return scenario;
}

// In the `cosmics` scenario, TODO doc
bool GEMClsDataAnalyzer::isInsideOut(const reco::Track& track) {
  return track.innerPosition().mag2() > track.outerPosition().mag2();
}

void GEMClsDataAnalyzer::buildGEMLayers(const GEMGeometry* gem) {
  std::map<GEMDetId, std::vector<const GEMChamber*> > chambers_per_layer;

  for (const GEMStation* station : gem->stations()) {
    const int region_id = station->region();
    const int station_id = station->station();
    const bool is_ge11 = station_id == 1;

    if (skipGEMStation(station_id)) {
      continue;
    }

    for (const GEMSuperChamber* superchamber : station->superChambers()) {
      // GE11: chamber == 0 for even chambers, chamber == 1 for odd chambers
      // GE21 and GE0: chamber == 0 for all chambers
      const int chamber_id = is_ge11 ? superchamber->id().chamber() % 2 : 0;

      for (const GEMChamber* chamber : superchamber->chambers()) {
        const int layer_id = chamber->id().layer();

        const GEMDetId key{region_id, 1, station_id, layer_id, chamber_id, 0};

        if (chambers_per_layer.find(key) == chambers_per_layer.end()) {
          chambers_per_layer.insert({key, std::vector<const GEMChamber*>()});
        }
        chambers_per_layer.at(key).push_back(chamber);
      }  // GEMChamber => iterate over layer ids
    }    // GEMSuperChamber => iterate over chamber ids
  }      // GEMStation

  gem_layers_.reserve(chambers_per_layer.size());
  for (auto [gem_id, chambers] : chambers_per_layer) {
    // layer position and rotation
    const float z_origin = chambers.front()->position().z();
    Surface::PositionType position{0.f, 0.f, z_origin};
    Surface::RotationType rotation;

    // eta partitions should have same R and Z spans.
    // XXX is it true?
    auto [r_min, r_max] = chambers.front()->surface().rSpan();
    auto [z_min, z_max] = chambers.front()->surface().zSpan();

    z_min -= z_origin;
    z_max -= z_origin;

    // the bounds from min and max R and Z in the local coordinates.
    SimpleDiskBounds* bounds = new SimpleDiskBounds(r_min, r_max, z_min, z_max);
    const Disk::DiskPointer layer = Disk::build(position, rotation, bounds);

    gem_layers_.emplace_back(layer, chambers, gem_id);

    LogDebug("ReadoutError") << gem_id
                            << Form(" ==> (z_origin, z_min, z_max) = (%.2f, %.2f, %.2f)", z_origin, z_min, z_max);
  }  // ring
}

// TODO doc
// See https://twiki.cern.ch/twiki/pub/CMS/GEMPPDOfflineDQM/check-muon-direction.pdf
bool GEMClsDataAnalyzer::checkPropagationDirection(const reco::Track* track, const GEMLayer& layer) {
  const bool is_same_region = track->eta() * layer.id.region() > 0;

  bool skip = false;
  if (kScenario_ == ScenarioOption::kCosmics) {
    float p2_in = track->innerMomentum().mag2();
    float p2_out = track->outerMomentum().mag2();
    if (isInsideOut(*track))
      std::swap(p2_in, p2_out);
    const bool is_outgoing = p2_in > p2_out;

    skip = (is_outgoing xor is_same_region);

  } else {
    // beam scenario
    skip = not is_same_region;
  }

  return skip;
}

GEMClsDataAnalyzer::StartingState GEMClsDataAnalyzer::buildStartingState(
    const reco::Muon& muon, const reco::TransientTrack& transient_track, const GEMLayer& gem_layer) {
  bool found = false;
  TrajectoryStateOnSurface state;
  DetId det_id;

  switch (kStartingStateType_) {
    case StartingStateType::kOutermostMeasurementState: {
      std::tie(found, state, det_id) = getOutermostMeasurementState(transient_track);
      break;
    }
    case StartingStateType::kInnermostMeasurementState: {
      std::tie(found, state, det_id) = getInnermostMeasurementState(transient_track);
      break;
    }
    case StartingStateType::kStateOnSurfaceWithCSCSegment: {
      std::tie(found, state, det_id) = buildStateOnSurfaceWithCSCSegment(muon, transient_track, gem_layer);
      break;
    }
    case StartingStateType::kAlignmentStyle: {
      std::tie(found, state, det_id) = buildStartingStateAlignmentStyle(muon, transient_track, gem_layer);
      break;
    }
    default: {
      edm::LogError("ReadoutError") << "got an unexpected StartingStateType";
      break;
    }
  }

  found &= state.isValid();

  if (found and (det_id.det() == DetId::Detector::Muon)) {
    found &= isMuonSubdetAllowed(det_id, gem_layer.id.station());
  }

  if (found) {
    if (MuonHitHelper::isGEM(det_id)) {
      const GEMDetId start_id{det_id};

      const bool are_same_region = gem_layer.id.region() == start_id.region();
      const bool are_same_station = gem_layer.id.station() == start_id.station();
      const bool are_same_layer = gem_layer.id.layer() == start_id.layer();
      if (are_same_region and are_same_station and are_same_layer) {
        LogDebug("ReadoutError")
            << "The starting detector of the muon propagation is same with the destination. Skip this propagation.";
        found = false;
      }
    }  // isGEM
  }    // found

  return std::make_tuple(found, state, det_id);
}

// Use the innermost measurement state as an initial state for the muon propagation.
// NOTE If the analyzer uses global or standalone muons and GEM hits are used in the
// muon reconstruction, the result should be biased.
// In 12_4_0_pre3, GEM hits are used in the pp scenario, but not in the cosmics scenario.
// https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre3/RecoMuon/StandAloneMuonProducer/python/standAloneMuons_cfi.py#L111-L127
// https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre3/RecoMuon/CosmicMuonProducer/python/cosmicMuons_cfi.py
GEMClsDataAnalyzer::StartingState GEMClsDataAnalyzer::getInnermostMeasurementState(
    const reco::TransientTrack& transient_track) {
  TrajectoryStateOnSurface state;
  DetId det_id;

  const reco::Track& track = transient_track.track();
  // get real innermost state
  if (isInsideOut(track)) {
    state = transient_track.outermostMeasurementState();
    det_id = track.outerDetId();

  } else {
    state = transient_track.innermostMeasurementState();
    det_id = track.innerDetId();
  }

  return std::make_tuple(true, state, det_id);
}

// Use the outermost measurement state as an initial state for the muon propagation.
GEMClsDataAnalyzer::StartingState GEMClsDataAnalyzer::getOutermostMeasurementState(
    const reco::TransientTrack& transient_track) {
  const reco::Track& track = transient_track.track();

  TrajectoryStateOnSurface state;
  DetId det_id;

  // get real innermost state
  if (isInsideOut(track)) {
    state = transient_track.innermostMeasurementState();
    det_id = track.innerDetId();

  } else {
    state = transient_track.outermostMeasurementState();
    det_id = track.outerDetId();
  }

  return std::make_tuple(true, state, det_id);
}

// Find the nearest CSC segment to the given GEM layer and then use a trajectory
// state on the surface with the segment as an initial state.
// XXX This method results in the residual phi distribution with two peaks
// because the muon and antimuon make different peaks.
GEMClsDataAnalyzer::StartingState GEMClsDataAnalyzer::buildStateOnSurfaceWithCSCSegment(
    const reco::Muon& muon, const reco::TransientTrack& transient_track, const GEMLayer& gem_layer) {
  bool found = false;
  TrajectoryStateOnSurface state;
  DetId det_id;

  if (const CSCSegment* csc_segment = findCSCSegment(muon, transient_track, gem_layer)) {
    const GeomDet* det = muon_service_->trackingGeometry()->idToDet(csc_segment->cscDetId());
    const GlobalPoint global_position = det->toGlobal(csc_segment->localPosition());

    found = true;
    state = transient_track.stateOnSurface(global_position);
    det_id = csc_segment->geographicalId();
  }

  return std::make_tuple(found, state, det_id);
}

// Find an ME11 segment and the build an initial state using the location and
// direction of the ME11 segment. If the muon has an inner track, the outerP of
// the inner track is used as the momentum magnitude. If not, the momentum
// magnitude is set to 1 GeV.
// https://github.com/gem-sw/alignment/blob/713e8fa/GEMCSCBendingAnalyzer/MuonAnalyser/plugins/analyser.cc#L435-L446
GEMClsDataAnalyzer::StartingState GEMClsDataAnalyzer::buildStartingStateAlignmentStyle(
    const reco::Muon& muon, const reco::TransientTrack& transient_track, const GEMLayer& gem_layer) {
  bool found = false;
  TrajectoryStateOnSurface state;
  DetId det_id;

  if (const CSCSegment* csc_segment = findCSCSegment(muon, transient_track, gem_layer)) {
    found = true;
    det_id = csc_segment->geographicalId();

    // position
    const LocalPoint position = csc_segment->localPosition();
    // momentum
    const reco::TrackRef inner_track = muon.innerTrack();
    const float momentum_magnitude = inner_track.isNonnull() ? inner_track.get()->outerP() : 1.0f;
    const LocalVector momentum = momentum_magnitude * csc_segment->localDirection();

    // trajectory parameter
    const LocalTrajectoryParameters trajectory_parameters{position, momentum, muon.charge()};

    // trajectory error
    const LocalTrajectoryError trajectory_error =
        asSMatrix<5>(csc_segment->parametersError().similarityT(csc_segment->projectionMatrix()));

    // surface
    const Plane& surface = muon_service_->trackingGeometry()->idToDet(det_id)->surface();

    state =
        TrajectoryStateOnSurface{trajectory_parameters, trajectory_error, surface, &*muon_service_->magneticField()};
  }

  return std::make_tuple(found, state, det_id);
}

// for beam scenario
const CSCSegment* GEMClsDataAnalyzer::findCSCSegmentBeam(const reco::TransientTrack& transient_track,
                                                            const GEMLayer& gem_layer) {
  const CSCSegment* best_csc_segment = nullptr;
  double min_z_distance = std::numeric_limits<double>::infinity();  // in cm

  for (trackingRecHit_iterator tracking_rechit_iter = transient_track.recHitsBegin();
       tracking_rechit_iter != transient_track.recHitsEnd();
       tracking_rechit_iter++) {
    const TrackingRecHit* tracking_rechit = *tracking_rechit_iter;
    if (not tracking_rechit->isValid()) {
      LogDebug("ReadoutError") << "got an invalid trackingRecHit_iterator from transient_track. skip it.";
      continue;
    }

    const DetId det_id = tracking_rechit->geographicalId();
    if (not MuonHitHelper::isCSC(det_id)) {
      continue;
    }

    if (tracking_rechit->dimension() != kCSCSegmentDimension_) {
      continue;
    }

    const CSCDetId csc_id{det_id};
    if (not isCSCAllowed(csc_id, gem_layer.id.station())) {
      continue;
    }

    if (auto csc_segment = dynamic_cast<const CSCSegment*>(tracking_rechit)) {
      const GeomDet* det = muon_service_->trackingGeometry()->idToDet(csc_id);
      if (det == nullptr) {
        edm::LogError("ReadoutError") << "GlobalTrackingGeometry::idToDet returns nullptr; CSCDetId=" << csc_id;
        continue;
      }
      const GlobalPoint global_position = det->toGlobal(csc_segment->localPosition());
      const float z_distance = std::abs(gem_layer.disk->localZclamped(global_position));

      if (z_distance < min_z_distance) {
        best_csc_segment = csc_segment;
        min_z_distance = z_distance;
      }

    } else {
      edm::LogError("ReadoutError")
          << "failed to perform the conversion from `const TrackingRechit*` to `const CSCSegment*`";
    }
  }  // trackingRecHit_iterator

  return best_csc_segment;
}

const CSCSegment* GEMClsDataAnalyzer::findCSCSegmentCosmics(const reco::Muon& muon, const GEMLayer& gem_layer) {
  const CSCSegment* best_csc_segment = nullptr;

  for (const reco::MuonChamberMatch& chamber_match : muon.matches()) {
    if (not MuonHitHelper::isCSC(chamber_match.id)) {
      continue;
    }

    const CSCDetId csc_id{chamber_match.id};
    if (not isCSCAllowed(csc_id, gem_layer.id.station())) {
      continue;
    }

    const float x_track = chamber_match.x;
    const float y_track = chamber_match.y;

    for (const reco::MuonSegmentMatch& segment_match : chamber_match.segmentMatches) {
      if (not segment_match.isMask(reco::MuonSegmentMatch::BestInStationByDR)) {
        continue;
      }

      const float dr = std::hypot(x_track - segment_match.x, y_track - segment_match.y);
      std::cout << "ReadoutError" << ": dr=" << dr << std::endl;

      if (dr > kMuonSegmentMatchDRCut_) {
        LogDebug("ReadoutError") << "too large dR(muon, segment)";
        break;
      }

      if (segment_match.cscSegmentRef.isNonnull()) {
        best_csc_segment = segment_match.cscSegmentRef.get();
      }
    }  // MuonSegmentMatch
  }    // MuonChamberMatch

  return best_csc_segment;
}

// just thin wrapper
const CSCSegment* GEMClsDataAnalyzer::findCSCSegment(const reco::Muon& muon,
                                                        const reco::TransientTrack& transient_track,
                                                        const GEMLayer& gem_layer) {
  if (kScenario_ == ScenarioOption::kCosmics) {
    return findCSCSegmentCosmics(muon, gem_layer);
  } else {
    // pp or HI
    return findCSCSegmentBeam(transient_track, gem_layer);
  }
}

bool GEMClsDataAnalyzer::isMuonSubdetAllowed(const DetId& det_id, const int gem_station) {
  if ((gem_station < 0) or (gem_station > 2)) {
    edm::LogError("ReadoutError") << "got unexpected gem station " << gem_station;
    return false;
  }

  if (det_id.det() != DetId::Detector::Muon) {
    edm::LogError("ReadoutError") << Form(
        "(Detector, Subdetector) = (%d, %d)", static_cast<int>(det_id.det()), det_id.subdetId());
    return false;
  }

  const std::vector<int> allowed = kMuonSubdetForGEM_.at(gem_station);
  return allowed.empty() or (std::find(allowed.begin(), allowed.end(), det_id.subdetId()) != allowed.end());
}

// Returns a bool value indicating whether or not the CSC detector can be used
// as a start detector for a given GEM station.
// See https://github.com/cms-sw/cmssw/blob/CMSSW_12_4_0_pre3/DataFormats/MuonDetId/interface/CSCDetId.h#L187-L193
// This method is used when using `buildStateOnSurfaceWithCSCSegment` or
// `buildStartingStateAlignmentStyle`
bool GEMClsDataAnalyzer::isCSCAllowed(const CSCDetId& csc_id, const int gem_station) {
  if ((gem_station < 0) or (gem_station > 2)) {
    edm::LogError("ReadoutError") << "got unexpected gem station " << gem_station;
    return false;
  }

  // unsigned short to int
  const int csc_chamber_type = static_cast<int>(csc_id.iChamberType());

  const std::vector<int> allowed = kCSCForGEM_.at(gem_station);
  return allowed.empty() or (std::find(allowed.begin(), allowed.end(), csc_chamber_type) != allowed.end());
}

bool GEMClsDataAnalyzer::checkBounds(const Plane& plane, const GlobalPoint& global_point) {
  const LocalPoint local_point = plane.toLocal(global_point);
  const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
  return plane.bounds().inside(local_point_2d);
}

// TODO comment on the scale
// https://github.com/cms-sw/cmssw/blob/CMSSW_12_0_0_pre3/DataFormats/GeometrySurface/src/SimpleDiskBounds.cc#L20-L35
bool GEMClsDataAnalyzer::checkBounds(const Plane& plane,
                                        const GlobalPoint& global_point,
                                        const GlobalError& global_error,
                                        const float scale) {
  const LocalPoint local_point = plane.toLocal(global_point);
  const LocalError local_error = ErrorFrameTransformer::transform(global_error, plane);

  const LocalPoint local_point_2d{local_point.x(), local_point.y(), 0.0f};
  return plane.bounds().inside(local_point_2d, local_error, scale);
}

const GEMEtaPartition* GEMClsDataAnalyzer::findEtaPartition(const GlobalPoint& global_point,
                                                               const GlobalError& global_error,
                                                               const std::vector<const GEMChamber*>& chamber_vector) {
  const GEMEtaPartition* bound = nullptr;
  for (const GEMChamber* chamber : chamber_vector) {
    if (not checkBounds(chamber->surface(), global_point, global_error, kBoundsErrorScale_)) {
      continue;
    }

    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      if (checkBounds(eta_partition->surface(), global_point, global_error, kBoundsErrorScale_)) {
        bound = eta_partition;
        break;
      }
    }  // GEMEtaPartition
  }    // GEMChamber

  return bound;
}

// Borrowed from https://github.com/gem-sw/alignment/blob/713e8fa/GEMCSCBendingAnalyzer/MuonAnalyser/plugins/analyser.cc#L321-L327
float GEMClsDataAnalyzer::computeRdPhi(const GlobalPoint& prop_global_pos,
                                          const LocalPoint& hit_local_pos,
                                          const GEMEtaPartition* eta_partition) {
  const StripTopology& topology = eta_partition->specificTopology();
  const LocalPoint prop_local_pos = eta_partition->toLocal(prop_global_pos);

  const float dx = prop_local_pos.x() - hit_local_pos.x();
  const float dy = prop_local_pos.y() - hit_local_pos.y();
  const float hit_strip = eta_partition->strip(hit_local_pos);
  const float hit_phi = topology.stripAngle(hit_strip);
  const float rdphi = std::cos(hit_phi) * dx - std::sin(hit_phi) * dy;
  return rdphi;
}

// Returns a global delta phi between a propagated muon and a reconstructed hit.
float GEMClsDataAnalyzer::computeDeltaPhi(const GlobalPoint& prop_global_pos,
                                             const LocalPoint& hit_local_pos,
                                             const GEMEtaPartition* eta_partition) {
  const GlobalPoint hit_global_pos = eta_partition->toGlobal(hit_local_pos);
  const float dphi = Geom::convertRadToDeg(prop_global_pos.phi() - hit_global_pos.phi());
  return dphi;
}

// a thin wrapper to hide a messy conditional statement
float GEMClsDataAnalyzer::computeMatchingMetric(const GlobalPoint& prop_global_pos,
                                                   const LocalPoint& hit_local_pos,
                                                   const GEMEtaPartition* eta_partition) {
  float metric;
  switch (kMatchingMetric_) {
    case MatchingMetric::kDeltaPhi: {
      metric = computeDeltaPhi(prop_global_pos, hit_local_pos, eta_partition);
      break;
    }
    case MatchingMetric::kRdPhi: {
      metric = computeRdPhi(prop_global_pos, hit_local_pos, eta_partition);
      break;
    }
    default: {
      edm::LogError("ReadoutError") << "unknown MatchingMetric.";  // TODO
      metric = std::numeric_limits<float>::infinity();
    }
  }

  return metric;
}

// This method finds the closest hit to a propagated muon in the eta partition
// with that propagated muon. Adjacent eta partitions are excluded from the area
// of interst to avoid ambiguity in defining the detection efficiency of each
// eta partition.
// std::pair<const GEMRecHit*, float> GEMClsDataAnalyzer::findClosestHit(const GlobalPoint& prop_global_pos,
//                                                                          const GEMRecHitCollection::range& rechit_range,
//                                                                          const GEMEtaPartition* eta_partition) {
std::tuple<const GEMRecHit*, float, int> GEMClsDataAnalyzer::findClosestHit(const GlobalPoint& prop_global_pos,
                                                                         const GEMRecHitCollection::range& rechit_range,
                                                                         const GEMEtaPartition* eta_partition) {
  const GEMRecHit* closest_hit = nullptr;
  float min_metric = std::numeric_limits<float>::infinity();
  int cls = -1;

  for (auto hit = rechit_range.first; hit != rechit_range.second; ++hit) {
    const LocalPoint hit_local_pos = hit->localPosition();

    const float metric = computeMatchingMetric(prop_global_pos, hit_local_pos, eta_partition);

    int clusterSize = hit->clusterSize();
    if (std::abs(metric) < std::abs(min_metric)) {
      min_metric = metric;
      closest_hit = &(*hit);
      cls = clusterSize;
    }
  }

  // return std::make_pair(closest_hit, min_metric);
  return std::make_tuple(closest_hit, min_metric, cls);
}

bool GEMClsDataAnalyzer::skipGEMStation(const int station) {
  bool skip = false;

  if (station == 0) {
    skip = not kMonitorGE0_;

  } else if (station == 1) {
    skip = not kMonitorGE11_;

  } else if (station == 2) {
    skip = not kMonitorGE21_;

  } else {
    edm::LogError("ReadoutError") << "got an unexpected GEM station " << station << ". skip this station.";
    skip = true;
  }

  return skip;
}

// int GEMClsDataAnalyzer::maskBigClusterEvent(const edm::Handle<GEMRecHitCollection> gemRecHits) {
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

int GEMClsDataAnalyzer::maskFlowerEvent(const edm::Handle<TCDSRecord> record) {

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

int GEMClsDataAnalyzer::maskChamberWithError(const GEMDetId& chamber_id,
                                                      const edm::Handle<GEMVFATStatusCollection> vfat_status_collection,
                                                      const edm::Handle<GEMOHStatusCollection> oh_status_collection) {
  int mask = 0;
  cout << "a" << endl;
  for (auto iter = oh_status_collection->begin(); iter != oh_status_collection->end(); iter++) {
    const auto [oh_id, range] = (*iter);
    if (chamber_id != oh_id) {
      cout << "aa" << endl;
      continue;
    }
    cout << "b" << endl;
    for (auto oh_status = range.first; oh_status != range.second; oh_status++) {
      if (oh_status->isBad()) {
        // GEMOHStatus is bad. Mask this chamber.
        mask = 1;
        return mask;
      }  // isBad
    }  // range
    cout << "c" << endl;
  }  // collection
  for (auto iter = vfat_status_collection->begin(); iter != vfat_status_collection->end(); iter++) {
    const auto [vfat_id, range] = (*iter);
    if (chamber_id != vfat_id.chamberId()) {
      continue;
    }
    cout << "d" << endl;
    for (auto vfat_status = range.first; vfat_status != range.second; vfat_status++) {
      if (vfat_status->isBad()) {
        mask = 1;
        return mask;
      }  // isBad
    }  // range
    cout << "e" << endl;
  }  // collection
  return mask;
}

void
GEMClsDataAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  iEvent.getByToken(kGEMRecHitCollectionToken_, gemRecHits);

  // edm::Handle<GEMVFATStatusCollection> vfat_status_collection;
  // iEvent.getByToken(vfat_status_collection_, vfat_status_collection);

  // edm::Handle<GEMOHStatusCollection> oh_status_collection;
  // iEvent.getByToken(oh_status_collection_, oh_status_collection);

  edm::Handle<TCDSRecord> record;
  iEvent.getByToken(tcdsRecord_, record);

  // edm::Handle<OnlineLuminosityRecord> onlineLumiRecord;
  // iEvent.getByToken(onlineLumiRecord_, onlineLumiRecord);

  // edm::Handle<edm::View<reco::Muon> > muonHandle;
  // iEvent.getByToken(muonHandle_, muonHandle);




  if (!record.isValid() || !gemRecHits.isValid()) {
    cout << "Invalid Event" << endl;
    return;
  }


  int BunchId = iEvent.bunchCrossing();
  int OrbitNumber = iEvent.orbitNumber();
  long EventTime = iEvent.time().unixTime();
  long Event = iEvent.id().event();
  int Lumi = iEvent.luminosityBlock();
  int Flr = maskFlowerEvent(record);

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
            // int ChamErr = maskChamberWithError(chId, vfat_status_collection, oh_status_collection);
            int ChamErr = 0;
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
                b_HitBunchId = BunchId;
                b_HitOrbitNumber = OrbitNumber;
                b_HitEventTime = EventTime;
                b_HitEvent = Event;
                b_HitLumi = Lumi;
                b_HitFlr = Flr;
                b_HitChamErr = ChamErr;
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
  b_ToBunchId = BunchId;
  b_ToOrbitNumber = OrbitNumber;
  b_ToEventTime = EventTime;
  b_ToEvent = Event;
  b_ToLumi = Lumi;
  // b_instLumi = onlineLumiRecord->instLumi();
  b_ToFlr = Flr;
  t_Total->Fill();


  //////////////////////////////////////////////////////////////////////////////
  // get data from Event
  //////////////////////////////////////////////////////////////////////////////
  const GEMRecHitCollection* rechit_collection = nullptr;
  if (auto handle = iEvent.getHandle(kGEMRecHitCollectionToken_)) {
    rechit_collection = handle.product();
  } else {
    edm::LogError("ReadoutError") << "failed to get GEMRecHitCollection";
    return;
  }

  const edm::View<reco::Muon>* muon_view = nullptr;
  if (auto handle = iEvent.getHandle(kMuonViewToken_)) {
    cout << "nMuon: " << handle->size() << endl;
    muon_view = handle.product();
  } else {
    edm::LogError("ReadoutError") << "failed to get View<Muon>";
    return;
  }

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
  // get data from EventSetup
  //////////////////////////////////////////////////////////////////////////////
  const TransientTrackBuilder* transient_track_builder = nullptr;
  if (auto handle = iSetup.getHandle(kTransientTrackBuilderToken_)) {
    transient_track_builder = handle.product();
  } else {
    edm::LogError("ReadoutError") << "failed to get TransientTrackBuilder";
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  // get more data from EventSetup using MuonServiceProxy
  //////////////////////////////////////////////////////////////////////////////
  muon_service_->update(iSetup);

  // TODO StraightLinePropagator if B < epsilon else SteppingHelixPropagatorAny
  const Propagator* propagator = nullptr;
  if (auto handle = muon_service_->propagator("SteppingHelixPropagatorAny")) {
    propagator = handle.product();
  } else {
    edm::LogError("ReadoutError") << "failed to get Propagator";
    return;
  }

  //////////////////////////////////////////////////////////////////////////////
  //  Main loop
  //////////////////////////////////////////////////////////////////////////////
  // cout << "main loop" << endl;
  for (const reco::Muon& muon : *muon_view) {

    // cout << "main loop in" << endl;
    // cout << "vaild (track,match): (" << muon.isAValidMuonTrack(kMuonTrackType_) << "," << muon.isMatchesValid() << ")" << endl;
    // cout << "valid (energy,quality,time): (" << muon.isEnergyValid() << "," << muon.isQualityValid() << "," << muon.isTimeValid() << ")" << endl;
    // cout << "valid (isol,pfisol): (" << muon.isIsolationValid() << "," << muon.isPFIsolationValid() << ")" << endl;
    // cout << "#chs: " << muon.numberOfChambers() << endl;
    // cout << "what mu: (" << muon.isMuon() << "," << muon.isGlobalMuon() << "," << muon.isTrackerMuon() << "," << muon.isStandAloneMuon() << "," << muon.isCaloMuon() << "," << muon.isPFMuon() << "," << muon.isRPCMuon() << "," << muon.isGEMMuon() << "," << muon.isME0Muon() << ")" << endl;
    // cout << "type: " << kMuonTrackType_ << endl;
    // auto tt = muon.muonTrack(kMuonTrackType_);
    // cout << "trackref: " << tt.isNull() << "," << tt.isNonnull() << "," << tt.isAvailable() << "," << tt.get() << "!" << endl;
    const reco::Track* track = muon.muonTrack(kMuonTrackType_).get();
    // cout << "track: " << track << "," << &track << endl;
    
    if (track == nullptr) {
      // cout << "trackended" << endl;
      LogDebug("ReadoutError") << "failed to get a " << kMuonTrackTypeName_;
      continue;
    }

    const reco::TransientTrack transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      // cout << "transient_track ended" << endl;
      edm::LogError("ReadoutError") << "failed to build TransientTrack";
      continue;
    }
    // cout << "layer loop" << endl;
    for (const GEMLayer& layer : gem_layers_) {
      // cout << "layer loop in" << endl;
      if (checkPropagationDirection(track, layer)) {
        LogDebug("ReadoutError") << "bad flight path. skip this propagation.";
        continue;
      }
      // cout << "build startingstate" << endl;
      const auto [found_start_state, start_state, start_id] = buildStartingState(muon, transient_track, layer);
      if (not found_start_state) {
        LogDebug("ReadoutError") << "propagation starting state not found";
        continue;
      }
      // cout << "propagator" << endl;
      // the trajectory state on the destination surface
      const auto [propagated_state, prop_path_length] = propagator->propagateWithPath(start_state, *(layer.disk));
      if (not propagated_state.isValid()) {
        LogDebug("ReadoutError") << "failed to propagate a muon from "
                                << Form("(Detector, Subdetector) = (%d, %d)",
                                        static_cast<int>(start_id.det()),
                                        start_id.subdetId())
                                << " to " << layer.id << ". The path length is " << prop_path_length;
        continue;
      }

      const GlobalPoint prop_global_pos = propagated_state.globalPosition();
      const GlobalError& prop_global_err = ErrorFrameTransformer::transform(propagated_state.localError().positionError(), *layer.disk);
      // cout << "checkbound" << endl;
      if (not checkBounds(*layer.disk, prop_global_pos, prop_global_err, kBoundsErrorScale_)) {
        LogDebug("ReadoutError") << "failed to pass checkBounds";
        continue;
      }
      // cout << "find eta partition" << endl;
      const GEMEtaPartition* eta_partition = findEtaPartition(prop_global_pos, prop_global_err, layer.chambers);
      if (eta_partition == nullptr) {
        LogDebug("ReadoutError") << "failed to find an eta partition";
        continue;
      }
      // cout << "aa" << endl;
      const GEMDetId gem_id = eta_partition->id();

      // if (true) {
      //   const bool has_error = maskChamberWithError(gem_id.chamberId(), vfat_status_collection, oh_status_collection);
      //   b_MuChamErr = has_error;
      //   if (has_error) {
      //     LogDebug("ReadoutError") << gem_id.chamberId() << " has an erorr. Skip this propagation.";
      //     continue;
      //   }
      // }

      //////////////////////////////////////////////////////////////////////////
      //
      //////////////////////////////////////////////////////////////////////////
      // const GEMDetId rs_key = Key2(gem_id);     // region-station
      // const GEMDetId rsl_key = Key3(gem_id);  // region-station-layer
      // const GEMDetId rse_key = Key3(gem_id);  // region-station-ieta
      const int re = gem_id.region();
      const int st = gem_id.station();
      const int ri = gem_id.ring();
      const int la = gem_id.layer();
      const int ch = gem_id.chamber();
      const int ie = gem_id.ieta();

      const double muon_pt = muon.pt();
      const double muon_eta = std::fabs(muon.eta());
      const double muon_phi = Geom::convertRadToDeg(muon.phi());

      const double prop_global_err_r = std::sqrt(prop_global_err.rerr(prop_global_pos));
      const double prop_global_err_phi = Geom::convertRadToDeg(std::sqrt(prop_global_err.phierr(prop_global_pos)));

      // cuts
      const bool passed_prop_err_r_cut = (prop_global_err_r < kPropagationErrorRCut_);
      const bool passed_prop_err_phi_cut = (prop_global_err_phi < kPropagationErrorPhiCut_);
      const bool passed_pt_cut = muon_pt > kMuonPtMinCuts_.at(st);
      const bool passed_eta_cut =
          (muon_eta > kMuonEtaMinCuts_.at(st)) and (muon_eta < kMuonEtaMaxCuts_.at(st));

      const bool passed_prop_err_cuts = passed_prop_err_r_cut and passed_prop_err_phi_cut;
      const bool passed_all_cuts = passed_prop_err_cuts and passed_pt_cut and passed_eta_cut;

      const int cutflow_last = not kModeDev_                 ? 0
                               : not passed_prop_err_r_cut   ? 1
                               : not passed_prop_err_phi_cut ? 2
                               : not passed_pt_cut           ? 3
                               : not passed_eta_cut          ? 4
                                                             : 5;

      //////////////////////////////////////////////////////////////////////////
      // Find a closet hit
      //////////////////////////////////////////////////////////////////////////
      // cout << "closest hit" << endl;
      // const auto [hit, matching_metric] =
      //     findClosestHit(prop_global_pos, rechit_collection->get(gem_id), eta_partition);
      const auto [hit, matching_metric, cls] =
          findClosestHit(prop_global_pos, rechit_collection->get(gem_id), eta_partition);

      if (hit == nullptr) {
        LogDebug("ReadoutError") << "hit not found";
        continue;
      }

      // if (kModeDev_) {
      //   fillMEWithinLimits(me_matching_metric_all_, rse_key, matching_metric);
      // }

      if (std::abs(matching_metric) > kMatchingCut_) {
        LogDebug("ReadoutError") << "failed to pass the residual rphi cut";
        continue;
      }


      //////////////////////////////////////////////////////////////////////////
      // Fill resolutions
      //////////////////////////////////////////////////////////////////////////
      // cout << "cut filtered" << endl;
      if (passed_all_cuts) {
        // cout << "passed" << endl;
        const LocalPoint hit_local_pos = hit->localPosition();
        const GlobalPoint& hit_global_pos = eta_partition->toGlobal(hit_local_pos);
        const float residual_phi = Geom::convertRadToDeg(prop_global_pos.phi() - hit_global_pos.phi());

        // fillMEWithinLimits(me_residual_phi_, rse_key, residual_phi);
        // cout << "modeDev: " << kModeDev_ << endl;
        if (kModeDev_) {
          const LocalPoint prop_local_pos = eta_partition->toLocal(prop_global_pos);
          const StripTopology& topology = eta_partition->specificTopology();

          const float residual_x = prop_local_pos.x() - hit_local_pos.x();
          const float residual_y = prop_local_pos.y() - hit_local_pos.y();
          const float residual_strip = topology.strip(prop_local_pos) - topology.strip(hit_local_pos);


          if (muon.charge() < 0) {
            cout << "negative muon" << endl;
            // fillMEWithinLimits(me_residual_phi_muon_, rse_key, residual_phi);
          } else {
            cout << "positive muon" << endl;
            // fillMEWithinLimits(me_residual_phi_antimuon_, rse_key, residual_phi);
          }
          // cout << "hit: " << *hit << endl;
          cout << "cls: " << cls << endl;
          cout << re << ", " << st << ", " << ri << ", " << la << ", " << ch << ", " << ie << endl;

          b_MuRe = re;
          b_MuSt = st;
          b_MuRi = ri;
          b_MuLa = la;
          b_MuCh = ch;
          b_MuIe = ie;
          b_MuCls = cls;
          b_MuLumi = Lumi;
          b_MuFlr = Flr;
          t_Muon->Fill();

          // cout << "end!" << endl;
        }  // kModeDev_
      }    // passed_all_cuts
    }      // destination
  }        // Muon
  // cout << "\n" << endl;

}

void GEMClsDataAnalyzer::beginJob(){}
void GEMClsDataAnalyzer::endJob(){}

void GEMClsDataAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup) {
  const GEMGeometry* gem = nullptr;
  if (auto handle = iSetup.getHandle(kGEMGeometryTokenBeginRun_)) {
    gem = handle.product();
  } else {
    edm::LogError("ReadoutError") << "failed to get GEMGeometry";
    return;
  }

  buildGEMLayers(gem);
}
// void GEMClsDataAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) { 
//   /* GEM Geometry */
//   edm::ESHandle<GEMGeometry> hGEMGeom;
//   hGEMGeom = iSetup.getHandle(hGEMGeomBeginRun_);

// //  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
//   // const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

//   // h_nEvents = fs->make<TH1I>("nEvents", "The number of events", 2, 0, 2);

// }
void GEMClsDataAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){
}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(GEMClsDataAnalyzer);
