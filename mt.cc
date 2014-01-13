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

void mt()
{

  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/Wg/Wg9.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/Wg/Wg9.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/Wg/Wg9.root/allObjects");

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
 
  TH1F* h_p_pt0      = new TH1F("h_p_pt0", "Raw Photon", 50, 0., 500.);
  TH1F* h_p_pt1      = new TH1F("h_p_pt1", "Loose Photon", 50, 0., 500.);
 
  TH1F* h_e_pt0      = new TH1F("h_e_pt0", "Raw Electron", 50, 0., 500.);
  TH1F* h_e_pt1      = new TH1F("h_e_pt1", "Medium Electron", 50, 0., 500.);

  TH1F* h_m_pt0      = new TH1F("h_m_pt0", "Raw Muon", 50, 0., 500.);
  TH1F* h_m_pt1      = new TH1F("h_m_pt1", "Tight Muon", 50, 0., 500.);

  TH1F* h_j_pt0      = new TH1F("h_j_pt0", "Raw Jet", 50, 0., 500.);
  TH1F* h_j_pt1      = new TH1F("h_j_pt1", "Loose Jet", 50, 0., 500.);
  TH1F* h_n_jets      = new TH1F("h_n_jets", "N Loose Jets", 8, 0., 8.);
  TH1F* h_n_jets2      = new TH1F("h_n_jets2", "N Raw Jets", 8, 0., 8.);
  float xx[8] = {2.,2.,2.,2.,4.,4.,4.,4.};
  //loop over events
  long iEntry = 0;//change to 0
  while(allObjects.GetEntry(iEntry++) != 0){

    h_n_jets->Fill(j.size);
    h_n_jets2->Fill(rj.size);
    for(unsigned int i(0);i<rp.size;i++) h_p_pt0->Fill(rp.pt[i]);
    for(unsigned int i(0);i<p.size;i++)  h_p_pt1->Fill(p.pt[i]);
    
    for(unsigned int i(0);i<re.size;i++) h_e_pt0->Fill(re.pt[i]);
    for(unsigned int i(0);i<e.size;i++)  h_e_pt1->Fill(e.pt[i]);

    for(unsigned int i(0);i<rm.size;i++) h_m_pt0->Fill(rm.pt[i]);
    for(unsigned int i(0);i<m.size;i++)  h_m_pt1->Fill(m.pt[i]);

    for(unsigned int i(0);i<rj.size;i++) h_j_pt0->Fill(rj.pt[i]);
    for(unsigned int i(0);i<j.size;i++)  h_j_pt1->Fill(j.pt[i]);

  }//while
TCanvas *x = new TCanvas("x","plots",800,800);
x->Divide(2,2);
x->cd(1);
h_n_jets->Scale(xx);
h_n_jets->Draw();
x->cd(2);
h_n_jets2->Draw();
x->cd(3);
TH1F *h1 = (TH1F*)h_n_jets->Clone("h1");
TH1F *h2 = (TH1F*)h_n_jets2->Clone("h2");
h1->Divide(h_n_jets2);
h1->Draw();
x->cd(4);
h2->Divide(h_n_jets);
h2->Draw();
/*
TCanvas *x = new TCanvas("x", "plots",800, 800);
x->Divide(2,2);
x->cd(1);
h_p_pt0->Draw();
h_p_pt1->Draw("same");
x->cd(2);
h_e_pt0->Draw();
h_e_pt1->Draw("same");
x->cd(3);
h_m_pt0->Draw();
h_m_pt1->Draw("same");
x->cd(4);
h_j_pt0->Draw();
h_j_pt1->Draw("same");
*/
}

