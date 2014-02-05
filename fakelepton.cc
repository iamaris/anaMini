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
#include "TH1.h"

void fakelepton()
{
  TChain eventVars("eventVars");//define eventVars tree
  //eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");
  eventVars.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/QCD/*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");
  selectedObjects.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/QCD/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  //allObjects.Add("/export/cmss/acalamba/photonhad/*.root/allObjects");
  allObjects.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/QCD/*.root/allObjects");

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
  //link the branches to the variables defined above
  //this is where you add variables relevant to your analysis
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
  
//  TH1F* h_st1 = new TH1F("h_st1", "St j=1", 50, 0., 2500.);
//  TH1F* h_st2 = new TH1F("h_st2", "St j=2", 50, 0., 2500.);  
//  h_st1->Sumw2();  
//  h_st2->Sumw2();  
//  TH1F* h_njet = new TH1F("h_njet","# of jets (corrected)",10,0.,10.);
//  TH1F* h_njet_raw = new TH1F("h_njet_raw","# of jets (uncorrected)",10,0.,10.);

  unsigned nCnt[20] = {0};
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;


    unsigned nLooseE = 0.0;
    unsigned nLooseM = 0.0;
    unsigned nTightE = 0.0;
    unsigned nTightM = 0.0;

    for(unsigned k(0);k<re.size;k++) {
      if(re.isLoose[k]) {
        nLooseE++;
      }
      if(re.isTight[k]) {
        nTightE++;
      }
    }

    for(unsigned k(0);k<rm.size;k++) {
      if(rm.isLoose[k]) {
        nLooseM++;
      }
      if(rm.isTight[k]) {
        nTightM++;
      }
    }


    if(nLooseE!=1 && nLooseM!=1) continue;

    if(nLooseE==1) {
       nCnt[2]++;
       if(nTightE==1) nCnt[3]++;
    }

    if(nLooseM==1) {
       nCnt[4]++;
       if(nTightM==1) nCnt[5]++;
    }


  }//while

std::cout<< "Total events                                 : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                   : "<< nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "event w/ 1 loose electron                    : "<< nCnt[2]<<" , "<<nCnt[2] <<std::endl;
std::cout<< "event w/ 1 loose&tight electron              : "<< nCnt[3]<<" , "<<nCnt[3] <<std::endl;
std::cout<< "event w/ 1 loose muon                        : "<< nCnt[4]<<" , "<<nCnt[4] <<std::endl;
std::cout<< "event w/ 1 loose&tight muon                  : "<< nCnt[5]<<" , "<<nCnt[5] <<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
//hlist.Add(h_st1);

TFile fout("fakelepton.root", "recreate");
hlist.Write();
fout.Close();


}

