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

void met()
{

  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");
  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");
  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/photonhad/*.root/allObjects");

  allObjects.AddFriend("selectedObjects");
  allObjects.AddFriend("eventVars");

  TChain eventVarsWg("eventVarsWg");//define eventVars tree
  eventVarsWg.Add("/export/cmss/acalamba/Wg/*.root/eventVars");
  TChain selectedObjectsWg("selectedObjectsWg");  //define eventVars tree
  selectedObjectsWg.Add("/export/cmss/acalamba/Wg/*.root/selectedObjects");
  TChain allObjectsWg("allObjectsWg");  //define eventVars tree
  allObjectsWg.Add("/export/cmss/acalamba/Wg/*.root/allObjects");

  allObjectsWg.AddFriend("selectedObjectsWg");
  allObjectsWg.AddFriend("eventVarsWg");

  //define variables
  bool hlt0;
  bool hlt1;
  float met;
  float metPhi;
  photon   p, rp;
  muon     m, rm;
  electron e, re;
  jet      j, rj; 


  bool hlt0Wg;
  bool hlt1Wg;
  float metWg;
  float metPhiWg;
  photon   pWg, rpWg;
  muon     mWg, rmWg;
  electron eWg, reWg;
  jet      jWg, rjWg;

  //link the branches to the variables defined above
  //this is where you add variables relevant to your analysis
  eventVars.SetBranchAddress("met", &met);
  eventVars.SetBranchAddress("metPhi", &metPhi);
  eventVars.SetBranchAddress("HLT_Photon70_CaloIdXL_PFHT400", &hlt0);
  eventVars.SetBranchAddress("HLT_Photon70_CaloIdXL_PFNoPUHT400", &hlt1);

  eventVarsWg.SetBranchAddress("met", &metWg);
  eventVarsWg.SetBranchAddress("metPhi", &metPhiWg);
  eventVarsWg.SetBranchAddress("HLT_Photon70_CaloIdXL_PFHT400", &hlt0Wg);
  eventVarsWg.SetBranchAddress("HLT_Photon70_CaloIdXL_PFNoPUHT400", &hlt1Wg);

  p.setAddress(selectedObjects);
  m.setAddress(selectedObjects);
  e.setAddress(selectedObjects);
  j.setAddress(selectedObjects);

  rp.setAddress(allObjects);
  re.setAddress(allObjects);
  rm.setAddress(allObjects);
  rj.setAddress(allObjects);

  pWg.setAddress(selectedObjectsWg);
  mWg.setAddress(selectedObjectsWg);
  eWg.setAddress(selectedObjectsWg);
  jWg.setAddress(selectedObjectsWg);

  rpWg.setAddress(allObjectsWg);
  reWg.setAddress(allObjectsWg);
  rmWg.setAddress(allObjectsWg);
  rjWg.setAddress(allObjectsWg);


  float bin[19] = {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,500};
  TH1F* h_met_d          = new TH1F("h_met_d", "MET Data", 18, bin);
  TH1F* h_met_b          = new TH1F("h_met_b", "MET Wg MC", 18, bin);
  TH1F* h_bin          = new TH1F("h_bin", "Divider", 18, bin);
 
  TH1F* h_mt_de = new TH1F("h_mt_de", "e Mt Data", 50, 0., 500.);
  TH1F* h_mt_be = new TH1F("h_mt_be", "e Mt Wg MC", 50, 0., 500.);

  TH1F* h_mt_dm = new TH1F("h_mt_dm", "m Mt Data", 50, 0., 500.);
  TH1F* h_mt_bm = new TH1F("h_mt_bm", "m Mt Wg MC", 50, 0., 500.);
  
  THStack* stack_met = new THStack("stack_met","MET");

  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    std::set<unsigned int> ce = MediumElectron(re,p);

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0 && !hlt1) continue;

    //Loose photon requirement
    if(p.size<1) continue;

    //Signal
    if(p.pt[0]<70.) continue;

    float mt_m = 0.0;
    for(unsigned int k(0);k<m.size;k++) {
      mt_m = mt(met,metPhi,m.pt[k],m.phi[k]);
      h_mt_dm->Fill(mt_m);
    }

    float mt_e = 0.0;
    for(std::set<unsigned int>::iterator it=ce.begin();it!=ce.end();++it) {
      mt_e = mt(met,metPhi,re.pt[*it],re.phi[*it]);
      h_mt_de->Fill(mt_e);
    }

    if(m.size>0 || ce.size()>0) {
      if(mt_e>100.0 || mt_m>100.0) {
        h_met_d->Fill(met);
      }
    }

  }//while

  long iEntryWg = 0;
  while(allObjectsWg.GetEntry(iEntryWg++) != 0){
    std::set<unsigned int> ceWg = MediumElectron(reWg,pWg);

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0Wg && !hlt1Wg) continue;

    //Loose photon requirement
    if(pWg.size<1) continue;

    //Signal
    if(pWg.pt[0]<70.) continue;
    float mt_m = 0.0;
    for(unsigned int k(0);k<mWg.size;k++) {
      mt_m = mt(metWg,metPhiWg,mWg.pt[k],mWg.phi[k]);
      h_mt_bm->Fill(mt_m,1.8);
    }

    float mt_e = 0.0;
    for(std::set<unsigned int>::iterator it=ceWg.begin();it!=ceWg.end();++it) {
      mt_e = mt(metWg,metPhiWg,reWg.pt[*it],reWg.phi[*it]);
      h_mt_be->Fill(mt_e,1.8);
    }

    if(mWg.size>0 || ceWg.size()>0) {
      if(mt_e>100.0 || mt_m>100.0) {
        h_met_b->Fill(metWg,1.8);
      }
    }

  }//while


TCanvas *hmet = new TCanvas("hmet", "MET plots",800, 800);
hmet->Divide(2,2);
hmet->cd(1);
int len = (sizeof(bin)/sizeof(*bin));
for (int i=1;i<len;i++) {
  int dx = bin[i]-bin[i-1];
  h_bin->Fill(bin[i-1],dx);
}
h_met_d->Divide(h_bin);
h_met_b->Divide(h_bin);
h_met_d->Draw();
h_met_b->Draw("same");
hmet->cd(2);
h_mt_dm->Draw();
h_mt_bm->Draw("same");
h_mt_de->Draw("same");
h_mt_be->Draw("same");
hmet->cd(3);
TH1F *h1 = (TH1F*)h_met_d->Clone("h1");
h1->Divide(h_met_b);
h1->Draw("e");
hmet->cd(4);
TH1F *h2 = (TH1F*)h_mt_dm->Clone("h2");
h2->Divide(h_mt_bm);
h2->Draw("e");
TH1F *h3 = (TH1F*)h_mt_de->Clone("h3");
h3->Divide(h_mt_be);
h3->Draw("samee");

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
hlist.Add( h_met_d );
hlist.Add( h_met_b );
hlist.Add( h_mt_dm );
hlist.Add( h_mt_de );
hlist.Add( h_mt_bm );
hlist.Add( h_mt_be);

TFile fout("met.root", "recreate");
hlist.Write();
fout.Close();


}

