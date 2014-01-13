#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <iostream>
#include "x.h"
#include "TCanvas.h"
#include <math.h>       
#include <set>       
#include <THStack.h>
using namespace mini;
const double PI = 4.0*atan(1.0); 
#include "TH1.h"

void signal()
{
  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/photonhad/D*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/photonhad/D*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/photonhad/D*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/ttbarjetgamma/*.root/allObjects");

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

  TH1F* h_met   = new TH1F("h_met", "met signal", 50, 0., 500.);
  TH1F* h_met_jet2gamma   = new TH1F("h_met_jet2gamma", "met j->#gamma", 50, 0., 500.);
  TH1F* h_met_e2gamma   = new TH1F("h_met_e2gamma", "met e->#gamma", 50, 0., 500.);
  TH1F* h_met_QCD   = new TH1F("h_met_QCD", "met QCD", 50, 0., 500.);
  THStack* stack_met = new THStack("stack_met","MET");

  
  TH1F* h_eg_pt_e = new TH1F("h_eg_pt_e", "e#gamma: electron pt", 50, 0., 500.);
  TH1F* h_eg_pt_fe = new TH1F("h_eg_pt_fe", "e#gamma: fake electron pt", 50, 0., 500.);
  TH1F* h_feg_pt_e = new TH1F("h_feg_pt_e", "e#gamma: electron pt", 50, 0., 500.);
  TH1F* h_feg_pt_fe = new TH1F("h_feg_pt_fe", "e#gamma: fake electron pt", 50, 0., 500.);
  TH1F* h_ug_pt_u = new TH1F("h_ug_pt_u", "#mu#gamma: muon pt", 50, 0., 500.);
  TH1F* h_ug_pt_fu = new TH1F("h_ug_pt_fu", "#mu#gamma: fake muon pt", 50, 0., 500.);
  TH1F* h_fug_pt_u = new TH1F("h_fug_pt_u", "#mu#gamma: muon pt", 50, 0., 500.);
  TH1F* h_fug_pt_fu = new TH1F("h_fug_pt_fu", "#mu#gamma: fake muon pt", 50, 0., 500.);


  unsigned nCnt[90] = {0};
  unsigned Nfg = 0;
  unsigned Ng = 0;
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;
    std::set<unsigned> fe = FakeElectron(re,5.0); //fake e
    std::set<unsigned> fu = FakeMuon(rm,5.0); //fake m
    std::set<unsigned> fg = ElectronFakePhoton(rp,70.0); //electron identified as photon?
    std::set<unsigned> fo = FakeableObject(rp,70.0); //gamma-e 
    std::set<unsigned> em = EMObject(rp,70.0); //EMObject 

    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;

    //l+FO (l+j->g) sample
    if(fo.size()>0) {
        //e+FO
        if(e.size>0 && e.pt[0]>5.0) {
          //TODO:Add deltaR requirement?
          h_met_jet2gamma->Fill(met);
        }

        //u+FO
       if(m.size>0 && m.pt[0]>5.0) {
          //TODO:Add deltaR requirement?
          h_met_jet2gamma->Fill(met);
       }      

    }

    //(l+e->g) sample    
    if(fg.size()>0) {
        //e+e->g
        if(e.size>0 && e.pt[0]>5.0) {
          //TODO:Add deltaR requirement?
           h_met_e2gamma->Fill(met);       
        }

        //u+e->g
       if(m.size>0 && m.pt[0]>5.0) {
          //TODO:Add deltaR requirement?
          h_met_e2gamma->Fill(met);
       }   
    }

    //QCD(fl + EM)
    if(em.size()>0) {
        //fe-g sample
        if(fe.size()>0) {
          if(e.size>0) h_feg_pt_e->Fill(e.pt[0]);
          h_feg_pt_fe->Fill(re.pt[*(fe.begin())]);
          h_met_QCD->Fill(met);
        }
        //fu-g sample
        if(fu.size()>0) {
          if(m.size>0) h_fug_pt_u->Fill(m.pt[0]);
          h_fug_pt_fu->Fill(rm.pt[*(fu.begin())]);
          h_met_QCD->Fill(met);
        }
    }

    
    if(p.size==0 && p.pt[0]>70.) continue;
   
    
    //signal
    if(e.size==0 || m.size==0) continue;
    Nfg = Nfg + fg.size();
    Ng  = Ng + p.size;

    //e-g sample
    if(e.size>0 && e.pt[0]>5.0) {
       if(deltaR(p.eta[0],p.phi[0],e.eta[0],e.phi[0])> 0.4) {
         h_eg_pt_e->Fill(e.pt[0]);
          h_met->Fill(met);
         if(fe.size()>0) h_eg_pt_fe->Fill(re.pt[*(fe.begin())]);
       }
    }

    //u-g sample
    if(m.size>0 && m.pt[0]>5.0) {
       if(deltaR(p.eta[0],p.phi[0],m.eta[0],m.phi[0])> 0.4) {
         h_ug_pt_u->Fill(m.pt[0]);
         h_met->Fill(met);
         if(fu.size()>0)h_ug_pt_fu->Fill(rm.pt[*(fu.begin())]);
       }
    }

  }//while

TCanvas *he = new TCanvas("he", "plots",800, 800);
he->Divide(2,2);
he->cd(1);
h_eg_pt_e->Draw();
h_feg_pt_e->Draw("same");
he->cd(2);
h_eg_pt_fe->Draw();
h_feg_pt_fe->Draw("same");
he->cd(3);
TH1F *e1 = (TH1F*)h_eg_pt_e->Clone("e1");
e1->Divide(h_feg_pt_e);
e1->Draw("e");
he->cd(4);
TH1F *e2 = (TH1F*)h_eg_pt_fe->Clone("e2");
e2->Divide(h_feg_pt_fe);
e2->Draw("e");


TCanvas *hu = new TCanvas("hu", "plots",800, 800);
hu->Divide(2,2);
hu->cd(1);
h_ug_pt_u->Draw();
h_fug_pt_u->Draw("same");
hu->cd(2);
h_ug_pt_fu->Draw();
h_fug_pt_fu->Draw("same");
hu->cd(3);
TH1F *u1 = (TH1F*)h_ug_pt_u->Clone("u1");
u1->Divide(h_fug_pt_u);
u1->Draw("e");
hu->cd(4);
TH1F *u2 = (TH1F*)h_ug_pt_fu->Clone("u2");
u2->Divide(h_fug_pt_fu);
u2->Draw("e");

TCanvas *hmet = new TCanvas("hmet","MET Stacked plots",800,800);
hmet->cd(1);
h_met_jet2gamma->Scale(0.03);
h_met_e2gamma->Scale(0.016);
h_met_QCD->Scale(0.000001);
stack_met->Add(h_met_jet2gamma);
stack_met->Add(h_met_e2gamma);
stack_met->Add(h_met_QCD);
h_met_jet2gamma->SetFillColor(kBlue);
h_met_e2gamma->SetFillColor(kRed);
h_met_QCD->SetFillColor(kGreen);
stack_met->Draw();
h_met->Draw("same");

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);

TFile fout("signal.root", "recreate");
hlist.Write();
fout.Close();


}

