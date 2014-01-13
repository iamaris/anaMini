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

void st_phohad()
{
  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/photonhad/*.root/allObjects");

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
  
  TH1F* h_st1 = new TH1F("h_st1", "St j=1", 50, 0., 2500.);
  TH1F* h_st2 = new TH1F("h_st2", "St j=2", 50, 0., 2500.);  
  TH1F* h_st3 = new TH1F("h_st3", "St j=3", 50, 0., 2500.);
  TH1F* h_st4 = new TH1F("h_st4", "St j=4", 50, 0., 2500.);
  TH1F* h_st5 = new TH1F("h_st5", "St j=5", 50, 0., 2500.);
  TH1F* h_st6 = new TH1F("h_st6", "St j=6", 50, 0., 2500.);
  TH1F* h_st7 = new TH1F("h_st7", "St j=7", 50, 0., 2500.);
  TH1F* h_st8 = new TH1F("h_st8", "St j>=8", 50, 0., 2500.);
  TH1F* h_st12 = new TH1F("h_st12", "St 1-2 jets", 50, 0., 2500.);
  TH1F* h_st34 = new TH1F("h_st34", "St 3-4 jets", 50, 0., 2500.);
  TH1F* h_st55 = new TH1F("h_st55", "St >4 jets", 50, 0., 2500.);
  
  TH1F* h_njet = new TH1F("h_njet","# of jets",10,0.,10.);
  TH1F* h_njet_raw = new TH1F("h_njet_raw","# of jets (uncorrected)",10,0.,10.);


  unsigned nCnt[90] = {0};
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;

    std::set<unsigned> tm;
    std::set<unsigned> te;
    std::set<unsigned> mp;

    float st(0);

    if (p.size>0) {
      for (unsigned k(0);k<rp.size;k++) {
        if(rp.isMedium[k]) {
          st = st + rp.pt[k];
          mp.insert(k); 
        }
      }
    }

    if (e.size>0) {
      for(unsigned k(0);k<re.size;k++) {
        if(re.isTight[k]) {
          st = st + re.pt[k];
          te.insert(k); 
        }
      }
    }

    if (m.size>0) {
      for(unsigned k(0);k<rm.size;k++) {
        if(rm.isTight[k]) {
          st = st + rm.pt[k];
          tm.insert(k);
        }
      }
    }


    int jSize(0);

    if (j.size>0) {
      for(unsigned k(0);k<j.size;k++) {
        if(j.pt[k]<30.) continue;
        bool same = false;
        for (std::set<unsigned>::iterator it=tm.begin(); it!=tm.end(); ++it) {
          if(deltaR(j.eta[k],j.phi[k],rm.eta[*it],rm.phi[*it])<0.5) {
            same = true; 
            break;
          }
        }
        if(same) continue;
        if (same==false) {
          for (std::set<unsigned>::iterator it=te.begin(); it!=te.end(); ++it) {
            if(deltaR(j.eta[k],j.phi[k],re.eta[*it],re.phi[*it])<0.5) {
              same = true;
              break;
            }
          }
        }
        if(same) continue;
        if (same==false) {
          for (std::set<unsigned>::iterator it=mp.begin(); it!=mp.end(); ++it) {
            if(deltaR(j.eta[k],j.phi[k],rp.eta[*it],rp.phi[*it])<0.5) {
              same = true;
              break;
            }
          }
        }
        if(same) continue;
        jSize++;
        st = st + j.pt[k];

      }//k
    }//j.size
    
    if(met>20) st = st + met;
    if(p.size==0) continue;
    if(p.pt[0]<70.) continue;

    if(e.size==0 && m.size==0) continue;
    h_njet->Fill(jSize);
    h_njet_raw->Fill(j.size);
    if(jSize==1) h_st1->Fill(st);
    if(jSize==2) h_st2->Fill(st);
    if(jSize==3) h_st3->Fill(st);
    if(jSize==4) h_st4->Fill(st);
    if(jSize==5) h_st5->Fill(st);
    if(jSize==6) h_st6->Fill(st);
    if(jSize==7) h_st7->Fill(st);
    if(jSize>=8) h_st8->Fill(st);

    if(jSize==1||jSize==2) h_st12->Fill(st);
    if(jSize==3||jSize==4) h_st34->Fill(st);
    if(jSize>=5) h_st55->Fill(st);
  }//while

TCanvas *hh = new TCanvas("hh", "St plots",800, 800);
hh->Divide(2,2);
hh->cd(1);
h_st3->Draw();
h_st1->Draw("same");
h_st2->Draw("same");
h_st4->Draw("same");
h_st5->Draw("same");
h_st6->Draw("same");
h_st7->Draw("same");
h_st8->Draw("same");
hh->cd(2);
h_njet_raw->Draw();
h_njet->Draw("same");


TCanvas *h = new TCanvas("h", "St plots",800, 800);
h->cd(1);
h_st1->Draw();
h_st2->Draw("same");
h->cd(2);
h_st2->Draw();
h_st3->Draw("same");
h->cd(3);
h_st3->Draw();
h_st4->Draw("same");
h->cd(4);
TH1F *h1 = (TH1F*)h_st2->Clone("h1");
h1->Divide(h_st1);
h1->Draw("e");
h->cd(5);
TH1F *h2 = (TH1F*)h_st3->Clone("h2");
h2->Divide(h_st2);
h2->Draw("e");
h->cd(6);
TH1F *h3 = (TH1F*)h_st4->Clone("h3");
h3->Divide(h_st3);
h3->Draw("e");

TCanvas *c = new TCanvas("c", "St plots",800, 800);
c->Divide(2,2);
c->cd(1);
h_st12->Draw();
h_st34->Draw("same");
c->cd(3);
TH1F *h4 = (TH1F*)h_st34->Clone("h4");
h4->Divide(h_st12);
h4->Draw("e");
c->cd(2);
h_st34->Draw();
h_st55->Draw("same");
c->cd(4);
TH1F *h5 = (TH1F*)h_st55->Clone("h5");
h5->Divide(h_st34);
h5->Draw("e");


/*
h_st1->Scale(1/Norm(h_st1,700,2000));
h_st2->Scale(1/Norm(h_st2,700,2000));
h_st3->Scale(1/Norm(h_st3,700,2000));
h_st4->Scale(1/Norm(h_st4,700,2000));
h_st5->Scale(1/Norm(h_st5,700,2000));
h_st6->Scale(1/Norm(h_st6,700,2000));
h_st7->Scale(1/Norm(h_st7,700,2000));
*/

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
hlist.Add(h_st1);
hlist.Add(h_st2);
hlist.Add(h_st3);
hlist.Add(h_st4);
hlist.Add(h_st5);
hlist.Add(h_st6);
hlist.Add(h_st7);
hlist.Add(h_st8);
hlist.Add(h_st12);
hlist.Add(h_st34);
hlist.Add(h_st55);

TFile fout("st_plots_photonhad.root", "recreate");
hlist.Write();
fout.Close();


}

