#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <iostream>
#include "x.h"
#include "TCanvas.h"
#include <math.h>       
#include <set>       
using namespace mini;
const double PI = 4.0*atan(1.0); 

void phoHadInfo()
{
  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("../*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/D*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("../*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/D*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("../*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/D*.root/allObjects");

  allObjects.AddFriend("selectedObjects");
  allObjects.AddFriend("eventVars");

  //define variables
  bool hlt0;
  bool hlt1;
  float met;
  float metPhi;
  photon   p, rp;
  muon     m, rm;
  electron e, re;
  jet      j, rj; 

  eventVars.SetBranchAddress("met", &met);
  eventVars.SetBranchAddress("metPhi", &metPhi);
  eventVars.SetBranchAddress("HLT_Photon70_CaloIdXL_PFHT400", &hlt0);
  eventVars.SetBranchAddress("HLT_Photon70_CaloIdXL_PFNoPUHT400", &hlt1);

  p.setAddress(selectedObjects);
  m.setAddress(selectedObjects);
  e.setAddress(selectedObjects);
  j.setAddress(selectedObjects);

  rp.setAddress(allObjects);
  re.setAddress(allObjects);
  rm.setAddress(allObjects);
  rj.setAddress(allObjects);
 
  //selectedObjects.Print();
  //define histograms
  TH1F* h_st = new TH1F("h_st", "St", 600, 0., 6000.);

  unsigned int nCnt[90] = {0};
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    unsigned n(0);
    nCnt[n]++;

    if(!hlt0 && !hlt1) continue;
    nCnt[++n]++;

    std::set<unsigned> ljB;
    std::set<unsigned> ljE;
        
    std::set<unsigned> lm; 
    std::set<unsigned> tm;    

    std::set<unsigned> leB;
    std::set<unsigned> leE;
    std::set<unsigned> meB;          
    std::set<unsigned> meE;    
    std::set<unsigned> teB;
    std::set<unsigned> teE;

    std::set<unsigned> lcB;
    std::set<unsigned> lcE;
    std::set<unsigned> mcB;//medium clean electron barrel
    std::set<unsigned> mcE;
    std::set<unsigned> tcB;//tight clean electron barrel
    std::set<unsigned> tcE;              

    std::set<unsigned> lpB;
    std::set<unsigned> lpE;
    std::set<unsigned> mpB;
    std::set<unsigned> mpE;
    std::set<unsigned> tpB; 
    std::set<unsigned> tpE;

   
    //jet
    for(unsigned i(0);i<j.size;i++) {
      if(j.iSubdet[i]==0) {//barrel
        if(j.isLoose[i]) ljB.insert(i);
      } else if (j.iSubdet[i]==1) {//endcap
        if(j.isLoose[i]) ljE.insert(i);
      }
    }

    //muon 
    for(unsigned i(0);i<m.size;i++) {
      if(m.isLoose[i]) lm.insert(i);      
      if(m.isTight[i]) tm.insert(i);      
    }

    //electron
    for(unsigned i(0);i<e.size;i++) {
      int dirty(0);
      for (unsigned k=0;k<p.size;k++) {
        if(mini::deltaR(p.eta[k],p.phi[k],e.eta[i],e.phi[i])<0.1) {
          dirty++; 
          break;
        }
      }

      if(e.iSubdet[i]==0) {//barrel
        if(e.isLoose[i])  leB.insert(i);
        if(e.isMedium[i]) meB.insert(i);
        if(e.isTight[i])  teB.insert(i);
        if(!dirty) {
          if(e.isLoose[i])  lcB.insert(i);
          if(e.isMedium[i]) mcB.insert(i);
          if(e.isTight[i])  tcB.insert(i);
        }
      } else if (e.iSubdet[i]==1) {//endcap
        if(e.isLoose[i])  leE.insert(i);
        if(e.isMedium[i]) meE.insert(i);
        if(e.isTight[i])  teE.insert(i);
        if(!dirty) {
          if(e.isLoose[i])  lcE.insert(i);
          if(e.isMedium[i]) mcE.insert(i);
          if(e.isTight[i])  tcE.insert(i);
        }

      }
    }

    //photon
    for(unsigned i(0);i<p.size;i++) {
      if(p.iSubdet[i]==0) {//barrel
        if(p.isLoose[i])  lpB.insert(i);
        if(p.isMedium[i]) mpB.insert(i);
        if(p.isTight[i])  tpB.insert(i);
      } else if (p.iSubdet[i]==1) {//endcap
        if(p.isLoose[i])  lpE.insert(i);
        if(p.isMedium[i]) mpE.insert(i);
        if(p.isTight[i])  tpE.insert(i);
      }
    }

    if(ljB.size()>0) nCnt[++n]++;
    if(ljE.size()>0) nCnt[++n]++; 
    
    if(lm.size()>0) nCnt[++n]++;
    if(tm.size()>0) nCnt[++n]++;
    
    if(leB.size()>0) nCnt[++n]++;
    if(leE.size()>0) nCnt[++n]++;
    if(meB.size()>0) nCnt[++n]++;      
    if(meE.size()>0) nCnt[++n]++;
    if(teB.size()>0) nCnt[++n]++;      
    if(teE.size()>0) nCnt[++n]++;
    
    if(lcB.size()>0) nCnt[++n]++;
    if(lcE.size()>0) nCnt[++n]++;
    if(mcB.size()>0) nCnt[++n]++;      
    if(mcE.size()>0) nCnt[++n]++;
    if(tcB.size()>0) nCnt[++n]++;      
    if(tcE.size()>0) nCnt[++n]++;    
    
    if(lpB.size()>0) nCnt[++n]++;
    if(lpE.size()>0) nCnt[++n]++;
    if(mpB.size()>0) nCnt[++n]++;      
    if(mpE.size()>0) nCnt[++n]++;
    if(tpB.size()>0) nCnt[++n]++;      
    if(tpE.size()>0) nCnt[++n]++;   




  }//while

int q(0);
std::cout<< "Total events                                 : "<< nCnt[q] <<std::endl;
std::cout<< "HLT passed                                   : "<< nCnt[++q] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose jet Barrel                          : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose jet EndCap                          : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose muon                                : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight muon                                : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose electron Barrel                     : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose electron EndCap                     : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron Barrel                    : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron EndCap                    : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron Barrel                     : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron EndCap                     : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose 'clean' electron Barrel             : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose 'clean' electron EndCap             : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron Barrel            : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron EndCap            : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron Barrel             : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron EndCap             : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose photon Barrel                       : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose photon EndCap                       : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium photon Barrel                      : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium photon EndCap                      : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight photon Barrel                       : "<< nCnt[++q]<<" , "<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight photon EndCap                       : "<< nCnt[++q]<<" , "<<q<<nCnt[q]/float(nCnt[1]) <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
//hlist.Add( h_mt );

TFile fout("info.root", "recreate");
hlist.Write();
fout.Close();



}

