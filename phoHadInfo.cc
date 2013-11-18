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
  eventVars.Add("./A1.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/D*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("./A1.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/D*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("./A1.root/allObjects");
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
    nCnt[0]++;

    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;

    std::set<unsigned int> tmB;
    std::set<unsigned int> lmB;
    std::set<unsigned int> teB;
    std::set<unsigned int> meB;
    std::set<unsigned int> leB;
    std::set<unsigned int> tpB;
    std::set<unsigned int> mpB;
    std::set<unsigned int> lpB;
    std::set<unsigned int> ljB;

    std::set<unsigned int> tmE;
    std::set<unsigned int> lmE;
    std::set<unsigned int> teE;
    std::set<unsigned int> meE;
    std::set<unsigned int> leE;
    std::set<unsigned int> tpE;
    std::set<unsigned int> mpE;
    std::set<unsigned int> lpE;
    std::set<unsigned int> ljE;

    std::set<unsigned int> tcB;//tight clean electron barrel
    std::set<unsigned int> mcB;//medium clean electron barrel
    std::set<unsigned int> lcB;
    std::set<unsigned int> tcE;
    std::set<unsigned int> mcE;
    std::set<unsigned int> lcE;
   
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
      if(m.iSubdet[i]==0) {//barrel
        if(m.isLoose[i]) lmB.insert(i);      
        if(m.isTight[i]) tmB.insert(i);      
      } else if (m.iSubdet[i]==1) {//endcap
        if(m.isLoose[i]) lmE.insert(i);
        if(m.isTight[i]) tmE.insert(i);
      }
    }

    //electron
    for(unsigned i(0);i<e.size;i++) {
      for (unsigned k=0;k<p.size;k++) {
        int clean(0);
        if(mini::deltaR(p.eta[k],p.phi[k],e.eta[i],e.phi[i])<0.1) {
          clean++; 
          break;
        }
      }

      if(e.iSubdet[i]==0) {//barrel
        if(e.isLoose[i])  leB.insert(i);
        if(e.isMedium[i]) meB.insert(i);
        if(e.isTight[i])  teB.insert(i);
      } else if (e.iSubdet[i]==1) {//endcap
        if(e.isLoose[i])  leE.insert(i);
        if(e.isMedium[i]) meE.insert(i);
        if(e.isTight[i])  teE.insert(i);
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




  }//while


std::cout<< "Total events                                 : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                   : "<< nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose jet Barrel                          : "<< nCnt[2]<<" , "<<nCnt[2]/nCnt[1] <<std::endl;
std::cout<< "w/ loose jet EndCap                          : "<< nCnt[3]<<" , "<<nCnt[3]/nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose muon Barrel                         : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ loose muon EndCap                         : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ tight muon Barrel                         : "<< nCnt[6]<<" , "<<nCnt[6]/nCnt[1] <<std::endl;
std::cout<< "w/ tight muon EndCap                         : "<< nCnt[7]<<" , "<<nCnt[7]/nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose electron Barrel                     : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ loose electron EndCap                     : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ medium electron Barrel                    : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ mudium electron EndCap                    : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ tight electron Barrel                     : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ tight electron EndCap                     : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose 'clean' electron Barrel             : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ loose 'clean' electron EndCap             : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ medium 'clean' electron Barrel            : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ mudium 'clean' electron EndCap            : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ tight 'clean' electron Barrel             : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ tight 'clean' electron EndCap             : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ loose photon Barrel                     : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ loose photon EndCap                     : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ medium photon Barrel                    : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ medium photon EndCap                    : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "w/ tight photon Barrel                     : "<< nCnt[4]<<" , "<<nCnt[4]/nCnt[1] <<std::endl;
std::cout<< "w/ tight photon EndCap                     : "<< nCnt[5]<<" , "<<nCnt[5]/nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;



std::cout<< "w/ selected electron , medium electron       : "<< nCnt[3]<<" , "<<nCnt[7] <<std::endl;
std::cout<< "w/ selected muon, tight muon                 : "<< nCnt[4]<<" , "<<nCnt[8] <<std::endl;
std::cout<< "w/ selected photon, loose photon             : "<< nCnt[5]<<" , "<<nCnt[9] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ 'clean' medium electron                   : "<< nCnt[10] <<std::endl;
std::cout<< "Loose photon w/ pt > 70                      : "<< nCnt[30] <<std::endl;
std::cout<< "Medium electron & loose photon 70            : "<< nCnt[31] <<std::endl;
std::cout<< "'Clean' medium electron & loose photon 70    : "<< nCnt[32] <<std::endl;
std::cout<< "Tight muon & loose photon 70                 : "<< nCnt[40] <<std::endl;
std::cout<< "Fake electron & EM Object 70                 : "<< nCnt[50] <<std::endl;
std::cout<< "Fake muon     & EM Object 70                 : "<< nCnt[60] <<std::endl;


TFile fout("st.root", "recreate");
hlist.Write();
fout.Close();


}

