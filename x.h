/* Auto-generated header file */
#ifndef x_h
#define x_h
#include <math.h>
#include <set>

namespace mini {

  unsigned const NMAX(512);

  float
  deltaR(float _eta1, float _phi1, float _eta2, float _phi2)
  {
    float dEta(_eta1 - _eta2);
    float dPhi(TVector2::Phi_mpi_pi(_phi1 - _phi2));
    return sqrt(dEta * dEta + dPhi * dPhi);
  }

  float 
  mt(float _met,float _metPhi,float _lPt,float _lPhi)
  {
  return sqrt(2.*_met*_lPt*(1-cos(_metPhi-_lPhi)));
  }

  class photon {
  public:
    photon() : size(0) {}
    ~photon() {}
    void setAddress(TTree&);
    void clear() { size = 0; }

    unsigned size;
    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float hOverE[NMAX];
    float sigmaIetaIeta[NMAX];
    float sigmaIphiIphi[NMAX];
    float etaWidth[NMAX];
    float phiWidth[NMAX];
    float r9[NMAX];
    float r5[NMAX];
    float trackerIso[NMAX];
    float ecalIso[NMAX];
    float hcalIso[NMAX];
    float chargedHadronIso[NMAX];
    float neutralHadronIso[NMAX];
    float photonIso[NMAX];
    float caloX[NMAX];
    float caloY[NMAX];
    float caloZ[NMAX];
    short iSubdet[NMAX];
    short superClusterIndex[NMAX];
    unsigned char nPixelSeeds[NMAX];
    unsigned char nClusters[NMAX];
    bool hasMatchedElectron[NMAX];
    bool electronVetoBit[NMAX];
    bool looseElectronVetoBit[NMAX];
    bool isLoose[NMAX];
    bool isMedium[NMAX];
    bool isTight[NMAX];
    bool isLoosePix[NMAX];
    bool isMediumPix[NMAX];
    bool isTightPix[NMAX];
    bool isLooseLV[NMAX];
    bool isMediumLV[NMAX];
    bool isTightLV[NMAX];
  };

  class electron {
  public:
    electron() : size(0) {}
    ~electron() {}
    void setAddress(TTree&);
    void clear() { size = 0; }

    unsigned size;
    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float combRelSubdetIso[NMAX];
    float combRelIso[NMAX];
    float deltaEta[NMAX];
    float deltaPhi[NMAX];
    float sigmaIetaIeta[NMAX];
    float sigmaIphiIphi[NMAX];
    float r9[NMAX];
    float r5[NMAX];
    float etaWidth[NMAX];
    float phiWidth[NMAX];
    float hOverE[NMAX];
    float d0[NMAX];
    float dz[NMAX];
    float epDiff[NMAX];
    float vtxFitProb[NMAX];
    float dCot[NMAX];
    float dist[NMAX];
    float caloX[NMAX];
    float caloY[NMAX];
    float caloZ[NMAX];
    short iSubdet[NMAX];
    short superClusterIndex[NMAX];
    unsigned char nClusters[NMAX];
    unsigned char nPixelHits[NMAX];
    unsigned char nMissingHits[NMAX];
    bool passConversionVeto[NMAX];
    bool isVeto[NMAX];
    bool isLoose[NMAX];
    bool isMedium[NMAX];
    bool isTight[NMAX];
  };

  class muon {
  public:
    muon() : size(0) {}
    ~muon() {}
    void setAddress(TTree&);
    void clear() { size = 0; }

    unsigned size;
    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float normChi2[NMAX];
    float dxy[NMAX];
    float dz[NMAX];
    float combRelSubdetIso[NMAX];
    float combRelIso[NMAX];
    short iSubdet[NMAX];
    unsigned char nMatchedStations[NMAX];
    unsigned char nLayersWithMmt[NMAX];
    unsigned char nValidMuonHits[NMAX];
    unsigned char nValidPixelHits[NMAX];
    bool isGlobalMuon[NMAX];
    bool isPFMuon[NMAX];
    bool hasInnerTrack[NMAX];
    bool hasGlobalTrack[NMAX];
    bool hasBestTrack[NMAX];
    bool isLoose[NMAX];
    bool isTight[NMAX];
  };

  class jet {
  public:
    jet() : size(0) {}
    ~jet() {}
    void setAddress(TTree&);
    void clear() { size = 0; }

    unsigned size;
    float pt[NMAX];
    float eta[NMAX];
    float phi[NMAX];
    float px[NMAX];
    float py[NMAX];
    float pz[NMAX];
    float energy[NMAX];
    float jecScale[NMAX];
    float chFraction[NMAX];
    float nhFraction[NMAX];
    float ceFraction[NMAX];
    float neFraction[NMAX];
    short iSubdet[NMAX];
    unsigned char nConstituents[NMAX];
    unsigned char nCharged[NMAX];
    bool passPUJetIdLoose[NMAX];
    bool isLoose[NMAX];
  };

  class vertex {
  public:
    vertex() : size(0) {}
    ~vertex() {}
    void setAddress(TTree&);
    void clear() { size = 0; }

    unsigned size;
    float x[NMAX];
    float y[NMAX];
    float z[NMAX];
    float rho[NMAX];
    float sumPt2[NMAX];
    float chi2[NMAX];
    float ndof[NMAX];
    unsigned short nTracks[NMAX];
    bool isGood[NMAX];
  };

  void
  photon::setAddress(TTree& _tree)
  {
    _tree.SetBranchAddress("photon.size", &size);
    if(_tree.GetBranch("photon.pt")) _tree.SetBranchAddress("photon.pt", pt);
    if(_tree.GetBranch("photon.eta")) _tree.SetBranchAddress("photon.eta", eta);
    if(_tree.GetBranch("photon.phi")) _tree.SetBranchAddress("photon.phi", phi);
    if(_tree.GetBranch("photon.px")) _tree.SetBranchAddress("photon.px", px);
    if(_tree.GetBranch("photon.py")) _tree.SetBranchAddress("photon.py", py);
    if(_tree.GetBranch("photon.pz")) _tree.SetBranchAddress("photon.pz", pz);
    if(_tree.GetBranch("photon.energy")) _tree.SetBranchAddress("photon.energy", energy);
    if(_tree.GetBranch("photon.hOverE")) _tree.SetBranchAddress("photon.hOverE", hOverE);
    if(_tree.GetBranch("photon.sigmaIetaIeta")) _tree.SetBranchAddress("photon.sigmaIetaIeta", sigmaIetaIeta);
    if(_tree.GetBranch("photon.sigmaIphiIphi")) _tree.SetBranchAddress("photon.sigmaIphiIphi", sigmaIphiIphi);
    if(_tree.GetBranch("photon.etaWidth")) _tree.SetBranchAddress("photon.etaWidth", etaWidth);
    if(_tree.GetBranch("photon.phiWidth")) _tree.SetBranchAddress("photon.phiWidth", phiWidth);
    if(_tree.GetBranch("photon.r9")) _tree.SetBranchAddress("photon.r9", r9);
    if(_tree.GetBranch("photon.r5")) _tree.SetBranchAddress("photon.r5", r5);
    if(_tree.GetBranch("photon.trackerIso")) _tree.SetBranchAddress("photon.trackerIso", trackerIso);
    if(_tree.GetBranch("photon.ecalIso")) _tree.SetBranchAddress("photon.ecalIso", ecalIso);
    if(_tree.GetBranch("photon.hcalIso")) _tree.SetBranchAddress("photon.hcalIso", hcalIso);
    if(_tree.GetBranch("photon.chargedHadronIso")) _tree.SetBranchAddress("photon.chargedHadronIso", chargedHadronIso);
    if(_tree.GetBranch("photon.neutralHadronIso")) _tree.SetBranchAddress("photon.neutralHadronIso", neutralHadronIso);
    if(_tree.GetBranch("photon.photonIso")) _tree.SetBranchAddress("photon.photonIso", photonIso);
    if(_tree.GetBranch("photon.caloX")) _tree.SetBranchAddress("photon.caloX", caloX);
    if(_tree.GetBranch("photon.caloY")) _tree.SetBranchAddress("photon.caloY", caloY);
    if(_tree.GetBranch("photon.caloZ")) _tree.SetBranchAddress("photon.caloZ", caloZ);
    if(_tree.GetBranch("photon.iSubdet")) _tree.SetBranchAddress("photon.iSubdet", iSubdet);
    if(_tree.GetBranch("photon.superClusterIndex")) _tree.SetBranchAddress("photon.superClusterIndex", superClusterIndex);
    if(_tree.GetBranch("photon.nPixelSeeds")) _tree.SetBranchAddress("photon.nPixelSeeds", nPixelSeeds);
    if(_tree.GetBranch("photon.nClusters")) _tree.SetBranchAddress("photon.nClusters", nClusters);
    if(_tree.GetBranch("photon.hasMatchedElectron")) _tree.SetBranchAddress("photon.hasMatchedElectron", hasMatchedElectron);
    if(_tree.GetBranch("photon.electronVetoBit")) _tree.SetBranchAddress("photon.electronVetoBit", electronVetoBit);
    if(_tree.GetBranch("photon.looseElectronVetoBit")) _tree.SetBranchAddress("photon.looseElectronVetoBit", looseElectronVetoBit);
    if(_tree.GetBranch("photon.isLoose")) _tree.SetBranchAddress("photon.isLoose", isLoose);
    if(_tree.GetBranch("photon.isMedium")) _tree.SetBranchAddress("photon.isMedium", isMedium);
    if(_tree.GetBranch("photon.isTight")) _tree.SetBranchAddress("photon.isTight", isTight);
    if(_tree.GetBranch("photon.isLoosePix")) _tree.SetBranchAddress("photon.isLoosePix", isLoosePix);
    if(_tree.GetBranch("photon.isMediumPix")) _tree.SetBranchAddress("photon.isMediumPix", isMediumPix);
    if(_tree.GetBranch("photon.isTightPix")) _tree.SetBranchAddress("photon.isTightPix", isTightPix);
    if(_tree.GetBranch("photon.isLooseLV")) _tree.SetBranchAddress("photon.isLooseLV", isLooseLV);
    if(_tree.GetBranch("photon.isMediumLV")) _tree.SetBranchAddress("photon.isMediumLV", isMediumLV);
    if(_tree.GetBranch("photon.isTightLV")) _tree.SetBranchAddress("photon.isTightLV", isTightLV);
  }

   void
   electron::setAddress(TTree& _tree)
   {
     _tree.SetBranchAddress("electron.size", &size);
     if(_tree.GetBranch("electron.pt")) _tree.SetBranchAddress("electron.pt", pt);
     if(_tree.GetBranch("electron.eta")) _tree.SetBranchAddress("electron.eta", eta);
     if(_tree.GetBranch("electron.phi")) _tree.SetBranchAddress("electron.phi", phi);
     if(_tree.GetBranch("electron.px")) _tree.SetBranchAddress("electron.px", px);
     if(_tree.GetBranch("electron.py")) _tree.SetBranchAddress("electron.py", py);
     if(_tree.GetBranch("electron.pz")) _tree.SetBranchAddress("electron.pz", pz);
     if(_tree.GetBranch("electron.energy")) _tree.SetBranchAddress("electron.energy", energy);
     if(_tree.GetBranch("electron.combRelSubdetIso")) _tree.SetBranchAddress("electron.combRelSubdetIso", combRelSubdetIso);
     if(_tree.GetBranch("electron.combRelIso")) _tree.SetBranchAddress("electron.combRelIso", combRelIso);
     if(_tree.GetBranch("electron.deltaEta")) _tree.SetBranchAddress("electron.deltaEta", deltaEta);
     if(_tree.GetBranch("electron.deltaPhi")) _tree.SetBranchAddress("electron.deltaPhi", deltaPhi);
     if(_tree.GetBranch("electron.sigmaIetaIeta")) _tree.SetBranchAddress("electron.sigmaIetaIeta", sigmaIetaIeta);
     if(_tree.GetBranch("electron.sigmaIphiIphi")) _tree.SetBranchAddress("electron.sigmaIphiIphi", sigmaIphiIphi);
     if(_tree.GetBranch("electron.r9")) _tree.SetBranchAddress("electron.r9", r9);
     if(_tree.GetBranch("electron.r5")) _tree.SetBranchAddress("electron.r5", r5);
     if(_tree.GetBranch("electron.etaWidth")) _tree.SetBranchAddress("electron.etaWidth", etaWidth);
     if(_tree.GetBranch("electron.phiWidth")) _tree.SetBranchAddress("electron.phiWidth", phiWidth);
     if(_tree.GetBranch("electron.hOverE")) _tree.SetBranchAddress("electron.hOverE", hOverE);
     if(_tree.GetBranch("electron.d0")) _tree.SetBranchAddress("electron.d0", d0);
     if(_tree.GetBranch("electron.dz")) _tree.SetBranchAddress("electron.dz", dz);
     if(_tree.GetBranch("electron.epDiff")) _tree.SetBranchAddress("electron.epDiff", epDiff);
     if(_tree.GetBranch("electron.vtxFitProb")) _tree.SetBranchAddress("electron.vtxFitProb", vtxFitProb);
     if(_tree.GetBranch("electron.dCot")) _tree.SetBranchAddress("electron.dCot", dCot);
     if(_tree.GetBranch("electron.dist")) _tree.SetBranchAddress("electron.dist", dist);
     if(_tree.GetBranch("electron.caloX")) _tree.SetBranchAddress("electron.caloX", caloX);
     if(_tree.GetBranch("electron.caloY")) _tree.SetBranchAddress("electron.caloY", caloY);
     if(_tree.GetBranch("electron.caloZ")) _tree.SetBranchAddress("electron.caloZ", caloZ);
     if(_tree.GetBranch("electron.iSubdet")) _tree.SetBranchAddress("electron.iSubdet", iSubdet);
     if(_tree.GetBranch("electron.superClusterIndex")) _tree.SetBranchAddress("electron.superClusterIndex", superClusterIndex);
     if(_tree.GetBranch("electron.nClusters")) _tree.SetBranchAddress("electron.nClusters", nClusters);
     if(_tree.GetBranch("electron.nPixelHits")) _tree.SetBranchAddress("electron.nPixelHits", nPixelHits);
     if(_tree.GetBranch("electron.nMissingHits")) _tree.SetBranchAddress("electron.nMissingHits", nMissingHits);
     if(_tree.GetBranch("electron.passConversionVeto")) _tree.SetBranchAddress("electron.passConversionVeto", passConversionVeto);
     if(_tree.GetBranch("electron.isVeto")) _tree.SetBranchAddress("electron.isVeto", isVeto);
     if(_tree.GetBranch("electron.isLoose")) _tree.SetBranchAddress("electron.isLoose", isLoose);
     if(_tree.GetBranch("electron.isMedium")) _tree.SetBranchAddress("electron.isMedium", isMedium);
     if(_tree.GetBranch("electron.isTight")) _tree.SetBranchAddress("electron.isTight", isTight);
   }
   
  void
  muon::setAddress(TTree& _tree)
  {
    _tree.SetBranchAddress("muon.size", &size);
    if(_tree.GetBranch("muon.pt")) _tree.SetBranchAddress("muon.pt", pt);
    if(_tree.GetBranch("muon.eta")) _tree.SetBranchAddress("muon.eta", eta);
    if(_tree.GetBranch("muon.phi")) _tree.SetBranchAddress("muon.phi", phi);
    if(_tree.GetBranch("muon.px")) _tree.SetBranchAddress("muon.px", px);
    if(_tree.GetBranch("muon.py")) _tree.SetBranchAddress("muon.py", py);
    if(_tree.GetBranch("muon.pz")) _tree.SetBranchAddress("muon.pz", pz);
    if(_tree.GetBranch("muon.energy")) _tree.SetBranchAddress("muon.energy", energy);
    if(_tree.GetBranch("muon.normChi2")) _tree.SetBranchAddress("muon.normChi2", normChi2);
    if(_tree.GetBranch("muon.dxy")) _tree.SetBranchAddress("muon.dxy", dxy);
    if(_tree.GetBranch("muon.dz")) _tree.SetBranchAddress("muon.dz", dz);
    if(_tree.GetBranch("muon.combRelSubdetIso")) _tree.SetBranchAddress("muon.combRelSubdetIso", combRelSubdetIso);
    if(_tree.GetBranch("muon.combRelIso")) _tree.SetBranchAddress("muon.combRelIso", combRelIso);
    if(_tree.GetBranch("muon.iSubdet")) _tree.SetBranchAddress("muon.iSubdet", iSubdet);
    if(_tree.GetBranch("muon.nMatchedStations")) _tree.SetBranchAddress("muon.nMatchedStations", nMatchedStations);
    if(_tree.GetBranch("muon.nLayersWithMmt")) _tree.SetBranchAddress("muon.nLayersWithMmt", nLayersWithMmt);
    if(_tree.GetBranch("muon.nValidMuonHits")) _tree.SetBranchAddress("muon.nValidMuonHits", nValidMuonHits);
    if(_tree.GetBranch("muon.nValidPixelHits")) _tree.SetBranchAddress("muon.nValidPixelHits", nValidPixelHits);
    if(_tree.GetBranch("muon.isGlobalMuon")) _tree.SetBranchAddress("muon.isGlobalMuon", isGlobalMuon);
    if(_tree.GetBranch("muon.isPFMuon")) _tree.SetBranchAddress("muon.isPFMuon", isPFMuon);
    if(_tree.GetBranch("muon.hasInnerTrack")) _tree.SetBranchAddress("muon.hasInnerTrack", hasInnerTrack);
    if(_tree.GetBranch("muon.hasGlobalTrack")) _tree.SetBranchAddress("muon.hasGlobalTrack", hasGlobalTrack);
    if(_tree.GetBranch("muon.hasBestTrack")) _tree.SetBranchAddress("muon.hasBestTrack", hasBestTrack);
    if(_tree.GetBranch("muon.isLoose")) _tree.SetBranchAddress("muon.isLoose", isLoose);
    if(_tree.GetBranch("muon.isTight")) _tree.SetBranchAddress("muon.isTight", isTight);
  }
  
  void
  jet::setAddress(TTree& _tree)
  {
    _tree.SetBranchAddress("jet.size", &size);
    if(_tree.GetBranch("jet.pt")) _tree.SetBranchAddress("jet.pt", pt);
    if(_tree.GetBranch("jet.eta")) _tree.SetBranchAddress("jet.eta", eta);
    if(_tree.GetBranch("jet.phi")) _tree.SetBranchAddress("jet.phi", phi);
    if(_tree.GetBranch("jet.px")) _tree.SetBranchAddress("jet.px", px);
    if(_tree.GetBranch("jet.py")) _tree.SetBranchAddress("jet.py", py);
    if(_tree.GetBranch("jet.pz")) _tree.SetBranchAddress("jet.pz", pz);
    if(_tree.GetBranch("jet.energy")) _tree.SetBranchAddress("jet.energy", energy);
    if(_tree.GetBranch("jet.jecScale")) _tree.SetBranchAddress("jet.jecScale", jecScale);
    if(_tree.GetBranch("jet.chFraction")) _tree.SetBranchAddress("jet.chFraction", chFraction);
    if(_tree.GetBranch("jet.nhFraction")) _tree.SetBranchAddress("jet.nhFraction", nhFraction);
    if(_tree.GetBranch("jet.ceFraction")) _tree.SetBranchAddress("jet.ceFraction", ceFraction);
    if(_tree.GetBranch("jet.neFraction")) _tree.SetBranchAddress("jet.neFraction", neFraction);
    if(_tree.GetBranch("jet.iSubdet")) _tree.SetBranchAddress("jet.iSubdet", iSubdet);
    if(_tree.GetBranch("jet.nConstituents")) _tree.SetBranchAddress("jet.nConstituents", nConstituents);
    if(_tree.GetBranch("jet.nCharged")) _tree.SetBranchAddress("jet.nCharged", nCharged);
    if(_tree.GetBranch("jet.passPUJetIdLoose")) _tree.SetBranchAddress("jet.passPUJetIdLoose", passPUJetIdLoose);
    if(_tree.GetBranch("jet.isLoose")) _tree.SetBranchAddress("jet.isLoose", isLoose);
  }  
  
  void
  vertex::setAddress(TTree& _tree)
  {
    _tree.SetBranchAddress("vertex.size", &size);
    if(_tree.GetBranch("vertex.x")) _tree.SetBranchAddress("vertex.x", x);
    if(_tree.GetBranch("vertex.y")) _tree.SetBranchAddress("vertex.y", y);
    if(_tree.GetBranch("vertex.z")) _tree.SetBranchAddress("vertex.z", z);
    if(_tree.GetBranch("vertex.rho")) _tree.SetBranchAddress("vertex.rho", rho);
    if(_tree.GetBranch("vertex.sumPt2")) _tree.SetBranchAddress("vertex.sumPt2", sumPt2);
    if(_tree.GetBranch("vertex.chi2")) _tree.SetBranchAddress("vertex.chi2", chi2);
    if(_tree.GetBranch("vertex.ndof")) _tree.SetBranchAddress("vertex.ndof", ndof);
    if(_tree.GetBranch("vertex.nTracks")) _tree.SetBranchAddress("vertex.nTracks", nTracks);
    if(_tree.GetBranch("vertex.isGood")) _tree.SetBranchAddress("vertex.isGood", isGood);
  }  

  std::set<unsigned>
  LooseJet(jet& _p) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if((_p.eta[i]<2.6) && (_p.eta[i]>-2.6)) {
        if(_p.nhFraction[i] >= 0.99) continue;
        if(_p.neFraction[i] >= 0.99) continue;
        if(_p.nConstituents[i] <= 1) continue;
        if(_p.chFraction[i] <= 0) continue;
        if(_p.ceFraction[i] >= 0.99) continue;
        if(_p.nCharged[i] <= 0) continue;
        loose.insert(i);
      } 
    }
    return loose;
  }

  std::set<unsigned> 
  LooseJetNew(jet& _p) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if((_p.eta[i]<2.4) && (_p.eta[i]>-2.4)) {
        if(_p.nhFraction[i] >= 0.99) continue;
        if(_p.neFraction[i] >= 0.99) continue;
        if(_p.nConstituents[i] <= 1) continue;
        if(_p.chFraction[i] <= 0) continue;
        if(_p.ceFraction[i] >= 0.99) continue;
        if(_p.nCharged[i] <= 0) continue;
        loose.insert(i);
      } else  {
        if(_p.nhFraction[i] >= 0.99) continue;
        if(_p.neFraction[i] >= 0.99) continue;
        if(_p.nConstituents[i] <= 1) continue;
        loose.insert(i);
      }
    }
    return loose;
  }

  std::set<unsigned> 
  LoosePhoton(photon& _p) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.iSubdet[i]==0) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if(_p.chargedHadronIso[i] > 2.6) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        loose.insert(i);
      } else if(_p.iSubdet[i]==1) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if(_p.chargedHadronIso[i] > 2.3) continue;
        if(_p.neutralHadronIso[i] > 2.9) continue;
        loose.insert(i);
      }
    }
    return loose;
  }

  std::set<unsigned> 
  LoosePhoton(photon& _p, float _pt) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.pt[i]<_pt) continue;
      if(_p.iSubdet[i]==0) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if(_p.chargedHadronIso[i] > 2.6) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        loose.insert(i);
      } else if(_p.iSubdet[i]==1) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if(_p.chargedHadronIso[i] > 2.3) continue;
        if(_p.neutralHadronIso[i] > 2.9) continue;
        loose.insert(i);
      }
    }
    return loose;
  }

  std::set<unsigned> 
  LoosePhoton(photon& _p, float _pt, bool _NoEndCap) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.pt[i]<_pt) continue;
      if(_NoEndCap && _p.iSubdet[i]==1) continue;
      if(_p.iSubdet[i]==0) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if(_p.chargedHadronIso[i] > 2.6) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        loose.insert(i);
      } else if(_p.iSubdet[i]==1) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if(_p.chargedHadronIso[i] > 2.3) continue;
        if(_p.neutralHadronIso[i] > 2.9) continue;
        loose.insert(i);
      }
    }
    return loose;
  }

  std::set<unsigned> 
  LoosePhoton(photon& _p, std::set<unsigned> _loose, float _pt=70.0) {
    for (std::set<unsigned>::iterator it=_loose.begin(); it!=_loose.end(); ++it) {
      if(_p.pt[*it]<_pt) _loose.erase(it);
    }
    return _loose;
  }


  std::set<unsigned>
  EMObject(photon& _p,float _pt) {
    std::set<unsigned> obj;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.pt[i] < _pt) continue;
      if(_p.iSubdet[i]==0) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if(_p.chargedHadronIso[i] > 2.6) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        obj.insert(i);
      } else if(_p.iSubdet[i]==1) {
        //if(_p.nPixelSeeds[i] > 0) continue;
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if(_p.chargedHadronIso[i] > 2.3) continue;
        if(_p.neutralHadronIso[i] > 2.9) continue;
        obj.insert(i);
      }
    }
    return obj;
  }

  std::set<unsigned>
  YLoosePhoton(photon& _p) {
    std::set<unsigned> loose;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.isLoose[i]) {
        loose.insert(i);
      }
    }
    return loose;
  }


  std::set<unsigned>
  YTightMuon(muon& _p) {
    std::set<unsigned> tight;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.isTight[i]) {
        tight.insert(i);
      } 
    }
    return tight;
  }


  std::set<unsigned> 
  JetFakePhoton(photon& _p) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      //if(_p.pt[i]<10.) continue;
      if(_p.iSubdet[i]==0) {
        if(_p.nPixelSeeds[i] > 0) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] < 0.012) continue;
        if(_p.chargedHadronIso[i] < 2.6) continue;  
        if(_p.neutralHadronIso[i] < 3.5) continue;
        if(_p.photonIso[i] < 1.3) continue;
        fake.insert(i);
      } else if(_p.iSubdet[i]==1) {
        if(_p.nPixelSeeds[i] > 0) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] < 0.034) continue;
        if(_p.chargedHadronIso[i] < 2.3) continue;  
        if(_p.neutralHadronIso[i] < 2.9) continue;
        fake.insert(i);
      }
    } 
    return fake;
  }

  std::set<unsigned> 
  JetFake(photon& _p) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.iSubdet[i]==0) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if((_p.chargedHadronIso[i] < 2.6)&&(_p.chargedHadronIso[i] > 15.0)) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        fake.insert(i);
      } else if(_p.iSubdet[i]==1) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if((_p.chargedHadronIso[i] < 2.3)&&(_p.chargedHadronIso[i] > 15.0)) continue;        
        if(_p.neutralHadronIso[i] > 2.9) continue;
        fake.insert(i);
      }
    }
    return fake;
  }

  std::set<unsigned> 
  JetFakeObject(photon& _p) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.iSubdet[i]==0) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.012) continue;
        //if((_p.chargedHadronIso[i] < 2.6)&&(_p.chargedHadronIso[i] > 15.0)) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        fake.insert(i);
      } else if(_p.iSubdet[i]==1) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.034) continue;
        //if((_p.chargedHadronIso[i] < 2.3)&&(_p.chargedHadronIso[i] > 15.0)) continue;        
        if(_p.neutralHadronIso[i] > 2.9) continue;
        fake.insert(i);
      }
    }
    return fake;
  }


  std::set<unsigned> 
  JetFake1(photon& _p, float _pt=0.0, bool _NoEndCap=false) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.pt[i]<_pt) continue;
      if(_NoEndCap && _p.iSubdet[i]==1) continue;
      if(_p.iSubdet[i]==0) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if((_p.chargedHadronIso[i] < 2.6)&&(_p.chargedHadronIso[i] > 15.0)) continue;
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        fake.insert(i);
      } else if(_p.iSubdet[i]==1) {
        if(!(_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        //if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if((_p.chargedHadronIso[i] < 2.3)&&(_p.chargedHadronIso[i] > 15.0)) continue;        
        if(_p.neutralHadronIso[i] > 2.9) continue;
        fake.insert(i);
      }
    }
    return fake;
  }

  std::set<unsigned> 
  JetPhoton(photon& _p,jet& _j) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.iSubdet[i]==1) continue; 
      for (unsigned k=0;k<_j.size;k++) {
        if(deltaR(_p.eta[i],_p.phi[i],_j.eta[k],_j.phi[k])<0.05) {
          fake.insert(i);
          break;
        }
      }
    } 
    return fake;
  }

  std::set<unsigned> 
  ElectronFakePhoton(photon& _p) {
    std::set<unsigned> fake;
    for(unsigned i=0;i<_p.size;i++) {
      if(_p.iSubdet[i]==0) {
        //if(_p.nPixelSeeds[i] == 0) continue;
        if((_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.012) continue;
        if(_p.chargedHadronIso[i] > 2.6) continue;  
        if(_p.neutralHadronIso[i] > 3.5) continue;
        if(_p.photonIso[i] > 1.3) continue;
        fake.insert(i);
      } else if(_p.iSubdet[i]==1) {
        //if(_p.nPixelSeeds[i] == 0) continue;
        if((_p.electronVetoBit[i])) continue;
        if(_p.hOverE[i] > 0.05) continue;
        if(_p.sigmaIetaIeta[i] > 0.034) continue;
        if(_p.chargedHadronIso[i] > 2.3) continue;  
        if(_p.neutralHadronIso[i] > 2.9) continue;
        fake.insert(i);
      }
    }
    return fake; 
  }


  //fake photon from jets
  //if(photon.sigmaIetaIeta < 0.012 && chIso < 2.6) goodPhotons.push_back(&photon);
  //else if(photon.sigmaIetaIeta < 0.014 && chIso < 15.) fakePhotons.push_back(&photon);
  std::set<unsigned>
  TightMuon(muon & _m) {
    std::set<unsigned> tight;
    for (unsigned i=0;i<_m.size;i++) {
        if(_m.pt[i] <= 200.) {
          if(_m.iSubdet[i] == -1) continue;
          if(!_m.isGlobalMuon[i]) continue;
          if(!_m.isPFMuon[i]) continue; 
          if(_m.normChi2[i] >= 10.) continue;
          if(_m.nValidMuonHits[i] <= 0) continue;
          if(_m.nMatchedStations[i] <= 1) continue;
          if(_m.dxy[i] >= 0.2) continue;
          if(_m.dz[i] >= 0.5) continue;
          if(_m.nValidPixelHits[i] <= 0) continue;
          if(_m.nLayersWithMmt[i] <=5) continue;
          if(_m.combRelIso[i] >= 0.12) continue;
          tight.insert(i);
        } else {
          if(_m.iSubdet[i] == -1) continue;
          if(!_m.isGlobalMuon[i]) continue;
          if(_m.nValidMuonHits[i] <= 0) continue;
          if(_m.nMatchedStations[i] <= 1) continue;
          if(_m.dxy[i] >= 0.2) continue;
          if(_m.dz[i] >= 0.5) continue;
          if(_m.nValidPixelHits[i] <= 0) continue;
          if(_m.nLayersWithMmt[i] <=8) continue;
          if(_m.combRelIso[i] >= 0.12) continue;
          tight.insert(i);
        }
      
    }
    return tight;
  }


  std::set<unsigned>
  FakeMuon(muon & _m) {
    std::set<unsigned> tight;
    for (unsigned i=0;i<_m.size;i++) {
        if(_m.pt[i] <= 200.) {
          if(_m.iSubdet[i] == -1) continue;
          if(!_m.isGlobalMuon[i]) continue;
          if(!_m.isPFMuon[i]) continue;
          if(_m.normChi2[i] >= 10.) continue;
          if(_m.nValidMuonHits[i] <= 0) continue;
          if(_m.nMatchedStations[i] <= 1) continue;
          if(_m.dxy[i] >= 0.2) continue;
          if(_m.dz[i] >= 0.5) continue;
          if(_m.nValidPixelHits[i] <= 0) continue;
          if(_m.nLayersWithMmt[i] <=5) continue;
          if(_m.combRelIso[i] < 0.12 && _m.combRelIso[i] > 0.25) continue;
          tight.insert(i);
        } else {
          if(_m.iSubdet[i] == -1) continue;
          if(!_m.isGlobalMuon[i]) continue;
          if(_m.nValidMuonHits[i] <= 0) continue;
          if(_m.nMatchedStations[i] <= 1) continue;
          if(_m.dxy[i] >= 0.2) continue;
          if(_m.dz[i] >= 0.5) continue;
          if(_m.nValidPixelHits[i] <= 0) continue;
          if(_m.nLayersWithMmt[i] <=8) continue;
          if(_m.combRelIso[i] < 0.12 && _m.combRelIso[i] > 0.25) continue;
          tight.insert(i);
        }

    }
    return tight;
  }


  std::set<unsigned>
  MediumElectron(electron& _e) {
    std::set<unsigned> medium;
    for (unsigned i=0;i<_e.size;i++) {
      if(_e.iSubdet[i]==0) {
        if(_e.deltaEta[i] > 0.004) continue;
        if(_e.deltaPhi[i] > 0.06) continue;
        if(_e.sigmaIetaIeta[i] > 0.01) continue;
        if(_e.hOverE[i] > 0.12) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] > 0.15) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        medium.insert(i);
      } else if(_e.iSubdet[i]==1){
        if(_e.deltaEta[i] > 0.007) continue;
        if(_e.deltaPhi[i] > 0.03) continue;
        if(_e.sigmaIetaIeta[i] > 0.03) continue;
        if(_e.hOverE[i] > 0.1) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] > 0.15) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        medium.insert(i);
      }
    }
    return medium;
  }

  std::set<unsigned>
  MediumElectron(electron& _e,photon& _p) {
    std::set<unsigned> medium;
    std::set<unsigned>::iterator it;

    for (unsigned i=0;i<_e.size;i++) {
      if(_e.iSubdet[i]==0) {
        if(_e.deltaEta[i] > 0.004) continue;
        if(_e.deltaPhi[i] > 0.06) continue;
        if(_e.sigmaIetaIeta[i] > 0.01) continue;
        if(_e.hOverE[i] > 0.12) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] > 0.15) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        medium.insert(i);
      } else if(_e.iSubdet[i]==1){
        if(_e.deltaEta[i] > 0.007) continue;
        if(_e.deltaPhi[i] > 0.03) continue;
        if(_e.sigmaIetaIeta[i] > 0.03) continue;
        if(_e.hOverE[i] > 0.1) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] > 0.15) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        medium.insert(i);
      }
    }

    for (it=medium.begin(); it!=medium.end(); ++it) {
      for (unsigned k=0;k<_p.size;k++) {
        if(deltaR(_p.eta[k],_p.phi[k],_e.eta[*it],_e.phi[*it])<0.1) {
          medium.erase(it);
          break;
        }
      }
    }

    return medium;
  }



  std::set<unsigned>
  FakeElectron(electron& _e) {
    std::set<unsigned> fake;
    for (unsigned i=0;i<_e.size;i++) {
      if(_e.iSubdet[i]==0) {
        if(_e.deltaEta[i] > 0.004) continue;
        if(_e.deltaPhi[i] > 0.06) continue;
        if(_e.sigmaIetaIeta[i] > 0.01) continue;
        if(_e.hOverE[i] > 0.12) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] < 0.15) continue;
        if(_e.combRelIso[i] > 0.30) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        fake.insert(i);
      } else if(_e.iSubdet[i]==1){
        if(_e.deltaEta[i] > 0.007) continue;
        if(_e.deltaPhi[i] > 0.03) continue;
        if(_e.sigmaIetaIeta[i] > 0.03) continue;
        if(_e.hOverE[i] > 0.1) continue;
        if(_e.d0[i] > 0.02) continue;
        if(_e.dz[i] > 0.1) continue;
        if(_e.epDiff[i] > 0.05) continue;
        if(_e.combRelIso[i] < 0.15) continue;
        if(_e.combRelIso[i] > 0.30) continue;
        if(!_e.passConversionVeto[i]) continue;
        if(_e.vtxFitProb[i] > 0.000001) continue;
        if(_e.nMissingHits[i] > 1) continue;
        fake.insert(i);
      }
    }
    return fake;
  }


}

#endif
