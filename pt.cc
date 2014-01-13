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

void pt()
{
  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/photonhad/A*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/photonhad/A*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/photonhad/A*.root/allObjects");

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
  TH1F* h_ec_pt = new TH1F("h_ec_pt", "control: e Pt", 25, 0., 500.);
  TH1F* h_es_pt = new TH1F("h_es_pt", "signal: e Pt", 25, 0., 500.);
  TH1F* h_fec_pt = new TH1F("h_fec_pt", "control: fake e Pt", 25, 0., 500.);
  TH1F* h_fes_pt = new TH1F("h_fes_pt", "signal: fake e Pt", 25, 0., 500.);

  TH1F* h_mc_pt = new TH1F("h_mc_pt", "control: m Pt", 25, 0., 500.);
  TH1F* h_ms_pt = new TH1F("h_ms_pt", "signal: m Pt", 25, 0., 500.);
  TH1F* h_fmc_pt = new TH1F("h_fmc_pt", "control: fake m Pt", 25, 0., 500.);
  TH1F* h_fms_pt = new TH1F("h_fms_pt", "signal: fake m Pt", 25, 0., 500.);

  unsigned int nCnt[90] = {0};
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;
    std::set<unsigned int> fe = FakeElectron(re);
    std::set<unsigned int> fm = FakeMuon(rm);
    std::set<unsigned int> em = EMObject(rp,70.);
    std::set<unsigned int> me = MediumElectron(re);
    std::set<unsigned int> tm = TightMuon(rm);
    std::set<unsigned int> lp = LoosePhoton(rp);
    std::set<unsigned int> ce = MediumElectron(re,p);

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;
 
    //ControlSample
    if(em.size()>0) {
      if(fe.size()>0 || fm.size()>0) {
        for (unsigned int i(0); i<m.size; i++) {
          h_mc_pt->Fill(m.pt[i]);
        }
        if(ce.size()>0) {
          for (std::set<unsigned int>::iterator it=ce.begin();it!=ce.end();++it) {
            h_ec_pt->Fill(re.pt[*it]);
          }
        }
        if(fe.size()>0) {
          for (std::set<unsigned int>::iterator it=fe.begin();it!=fe.end();++it) {
            h_fec_pt->Fill(re.pt[*it]);
          }
        }
        if(fm.size()>0) {
          for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
            h_fmc_pt->Fill(rm.pt[*it]);
          }
        }
      }
    }

    //Loose photon requirement
    if(p.size<1) continue;

    //Signal
    if(p.pt[0]<70.) continue;
    if(m.size<1 && ce.size()<1) continue;

    for (unsigned int i(0); i<m.size; i++) {
      h_ms_pt->Fill(m.pt[i]);
    }
    if(ce.size()>0) {
      for (std::set<unsigned int>::iterator it=ce.begin();it!=ce.end();++it) {
        h_es_pt->Fill(re.pt[*it]);
       }
    }
    if(fe.size()>0) {
      for (std::set<unsigned int>::iterator it=fe.begin(); it!=fe.end(); ++it) {
        h_fes_pt->Fill(re.pt[*it]);
      }
    }
    if(fm.size()>0) {
      for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
        h_fms_pt->Fill(rm.pt[*it]);
      }
    }

  }//while

TCanvas *elec = new TCanvas("elec", "pt",800, 800);
elec->Divide(2,2);
elec->cd(1);
h_ec_pt->Draw("e");
h_es_pt->Draw("samee");
elec->cd(2);
h_fec_pt->Draw("e");
h_fes_pt->Draw("samee");
elec->cd(3);
TH1F *e1 = (TH1F*)h_es_pt->Clone("e1");
e1->Divide(h_ec_pt);
e1->Draw("e");
elec->cd(4);
TH1F *e2 = (TH1F*)h_fes_pt->Clone("e2");
e2->Divide(h_fec_pt);
e2->Draw("e");

TCanvas *muon = new TCanvas("muon", "pt",800, 800);
muon->Divide(2,2);
muon->cd(1);
h_mc_pt->Draw("e");
h_ms_pt->Draw("samee");
muon->cd(2);
h_fmc_pt->Draw("e");
h_fms_pt->Draw("samee");
muon->cd(3);
TH1F *m1 = (TH1F*)h_ms_pt->Clone("m1");
m1->Divide(h_mc_pt);
m1->Draw("e");
muon->cd(4);
TH1F *m2 = (TH1F*)h_fms_pt->Clone("m2");
m2->Divide(h_fmc_pt);
m2->Draw("e");

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);

TFile fout("delete.root", "recreate");
hlist.Write();
fout.Close();


}

