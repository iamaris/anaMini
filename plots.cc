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

void plots()
{
  TChain eventVars("eventVars");//define eventVars tree
  //eventVars.Add("./A1.root/eventVars");
  //eventVars.Add("./*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/A*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/B*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/C*.root/eventVars");
  eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");
  eventVars.Add("/home/acalamba/Desktop/d/*.root/eventVars");


  TChain selectedObjects("selectedObjects");  //define eventVars tree
  //selectedObjects.Add("./A1.root/selectedObjects");
  //selectedObjects.Add("./*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/A*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/B*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/C*.root/selectedObjects");
  selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");
  selectedObjects.Add("/home/acalamba/Desktop/d/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  //allObjects.Add("./A1.root/allObjects");
  //allObjects.Add("./*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/A*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/B*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/C*.root/allObjects");
  allObjects.Add("/export/cmss/acalamba/photonhad/*.root/allObjects");
  allObjects.Add("/home/acalamba/Desktop/d/*.root/allObjects");


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
  TH1F* h_p0_pt = new TH1F("h_p0_pt", "Raw Photon Pt", 50, 0., 500.);
  TH1F* h_p1_pt = new TH1F("h_p1_pt", "HLT + Raw Photon Pt", 50, 0., 500.);
  TH1F* h_p2_pt = new TH1F("h_p2_pt", "Loose Photon Pt", 50, 0., 500.);
  TH1F* h_p3_pt = new TH1F("h_p3_pt", "Loose Photon Pt 70", 50, 0., 500.);
  TH1F* h_p4_pt = new TH1F("h_p4_pt", "Loose Photon Pt + Clean e", 50, 0., 500.);
  TH1F* h_p5_pt = new TH1F("h_p5_pt", "Loose Photon Pt + Tight Muon", 50, 0., 500.);
  TH1F* h_p6_pt = new TH1F("h_p6_pt", "Photon Pt", 50, 0., 500.);
  TH1F* h_p7_pt = new TH1F("h_p7_pt", "Photon Pt", 50, 0., 500.);

  TH1F* h_e0_pt = new TH1F("h_e0_pt", "Raw electron Pt", 50, 0., 500.);
  TH1F* h_e1_pt = new TH1F("h_e1_pt", "HLT + raw electron Pt", 50, 0., 500.);
  TH1F* h_e2_pt = new TH1F("h_e2_pt", "medium electron Pt", 50, 0., 500.);
  TH1F* h_e3_pt = new TH1F("h_e3_pt", "medium electron Pt w/ loose photon", 50, 0., 500.);
  TH1F* h_e4_pt = new TH1F("h_e4_pt", "medium electron Pt w/ loose p 70", 50, 0., 500.);
  TH1F* h_e5_pt = new TH1F("h_e5_pt", "clean electron Pt w/ loose p 70", 50, 0., 500.);

  TH1F* h_m0_pt = new TH1F("h_m0_pt", "Raw muon Pt", 50, 0., 500.);
  TH1F* h_m1_pt = new TH1F("h_m1_pt", "HLT + raw muon Pt", 50, 0., 500.);
  TH1F* h_m2_pt = new TH1F("h_m2_pt", "tight muon Pt", 50, 0., 500.);
  TH1F* h_m3_pt = new TH1F("h_m3_pt", "tight muon Pt w/ loose p", 50, 0., 500.);
  TH1F* h_m4_pt = new TH1F("h_m4_pt", "tight muon Pt w/ loose p 70", 50, 0., 500.);

  TH1F* h_p0_met = new TH1F("h_p0_met", "Raw photon met", 50, 0., 500.);
  TH1F* h_p1_met = new TH1F("h_p1_met", "HLT + Raw photon met", 50, 0., 500.);
  TH1F* h_p2_met = new TH1F("h_p2_met", "Loose photon met", 50, 0., 500.);
  TH1F* h_p3_met = new TH1F("h_p3_met", "Loose photon 70 met", 50, 0., 500.);

  TH1F* h_e0_met = new TH1F("h_e0_met", "Raw Electron met", 50, 0., 500.);
  TH1F* h_e1_met = new TH1F("h_e1_met", "HLT + Raw Electron met", 50, 0., 500.);
  TH1F* h_e2_met = new TH1F("h_e2_met", "Medium Electron met", 50, 0., 500.);
  TH1F* h_e3_met = new TH1F("h_e3_met", "Medium Electron + Loose Photon met", 50, 0., 500.);
  TH1F* h_e4_met = new TH1F("h_e4_met", "Medium e + Loose Photon 70 met", 50, 0., 500.);
  TH1F* h_e5_met = new TH1F("h_e5_met", "Clean electron + loose g 70 met", 50, 0., 500.);

  TH1F* h_m0_met = new TH1F("h_m0_met", "Raw muon met", 50, 0., 500.);
  TH1F* h_m1_met = new TH1F("h_m1_met", "HLT + Raw muon met", 50, 0., 500.);
  TH1F* h_m2_met = new TH1F("h_m2_met", "Tight muon met", 50, 0., 500.);
  TH1F* h_m3_met = new TH1F("h_m3_met", "Tight muon + loose photon met", 50, 0., 500.);
  TH1F* h_m4_met = new TH1F("h_m4_met", "Tight muon + loose photon 70 met", 50, 0., 500.);

  TH1F* h_e0_mt = new TH1F("h_e0_mt", "Raw Electron Mt", 50, 0., 500.);
  TH1F* h_e1_mt = new TH1F("h_e1_mt", "HLT + Raw Electron Mt", 50, 0., 500.);
  TH1F* h_e2_mt = new TH1F("h_e2_mt", "Medium Electron Mt", 50, 0., 500.);
  TH1F* h_e3_mt = new TH1F("h_e3_mt", "Medium Electron Mt + Loose Photon", 50, 0., 500.);
  TH1F* h_e4_mt = new TH1F("h_e4_mt", "Medium Electron Mt + Loose Photon 70", 50, 0., 500.);
  TH1F* h_e5_mt = new TH1F("h_e5_mt", "Clean Electron Mt + Loose Photon", 50, 0., 500.);

  TH1F* h_m0_mt = new TH1F("h_m0_mt", "Raw Muon Mt", 50, 0., 500.);
  TH1F* h_m1_mt = new TH1F("h_m1_mt", "Raw Muon Mt + HLT", 50, 0., 500.);
  TH1F* h_m2_mt = new TH1F("h_m2_mt", "Tight Muon Mt", 50, 0., 500.);
  TH1F* h_m3_mt = new TH1F("h_m3_mt", "Tight Muon Mt + Loose Photon", 50, 0., 500.);
  TH1F* h_m4_mt = new TH1F("h_m4_mt", "Tight Muon Mt + Loose Photon > 70GeV", 50, 0., 500.);
 

  TH1F* h_m_pt  = new TH1F("h_m_pt", "Tight muon pt", 50, 0., 100.);
  TH1F* h_e_pt  = new TH1F("h_e_pt", "Loose electron pt", 50, 0., 100.);
  TH1F* h_ce_pt = new TH1F("h_ce_pt", "Clean electron pt", 50, 0., 100.);


  TH1F* h_em0_pt  = new TH1F("h_em0_pt", "EM Object pt", 50, 0., 500.);
  TH1F* h_em1_pt  = new TH1F("h_em1_pt", "EM Object pt", 50, 0., 500.);
  TH1F* h_em2_pt  = new TH1F("h_em2_pt", "EM Object pt", 50, 0., 500.);

  TH1F* h_met   = new TH1F("h_met", "met", 50, 0., 500.);
  TH1F* h_met0   = new TH1F("h_met0", "fe-g met", 50, 0., 500.);
  TH1F* h_met1   = new TH1F("h_met1", "fm-g met", 50, 0., 500.);
  TH1F* h_met2   = new TH1F("h_met2", "fe-em met", 50, 0., 500.);
  TH1F* h_met3   = new TH1F("h_met3", "fm-em met", 50, 0., 500.);

  TH1F* h_mt0 = new TH1F("h_mt0", "fe mt", 50, 0., 500.);
  TH1F* h_mt1 = new TH1F("h_mt1", "fm mt", 50, 0., 500.);
  TH1F* h_mt2 = new TH1F("h_mt2", "fe+em mt", 50, 0., 500.);
  TH1F* h_mt3 = new TH1F("h_mt3", "fm+em mt", 50, 0., 500.);
  TH1F* h_mt4 = new TH1F("h_mt4", "fe+g mt", 50, 0., 500.);
  TH1F* h_mt5 = new TH1F("h_mt5", "fm+g mt", 50, 0., 500.);

  unsigned int nCnt[90] = {0};
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;
    h_met->Fill(met);

    std::set<unsigned int> fe = FakeElectron(re);
    std::set<unsigned int> fm = FakeMuon(rm);
    std::set<unsigned int> em = EMObject(rp,70.);
    std::set<unsigned int> me = MediumElectron(re);
    std::set<unsigned int> tm = TightMuon(rm);
    std::set<unsigned int> lp = LoosePhoton(rp);
    std::set<unsigned int> ce = MediumElectron(re,p);
    std::set<unsigned int> lp70 = LoosePhoton(rp,lp,70.0);

    if (rp.size>0) {
      h_p0_met->Fill(met);
      for (unsigned int k(0);k<rp.size;k++) h_p0_pt->Fill(rp.pt[k]);
    }
    if (re.size>0) {
      h_e0_met->Fill(met);
      for(unsigned int k(0);k<re.size;k++) {
        h_e0_pt->Fill(re.pt[k]);
        h_e0_mt->Fill(mt(met,metPhi,re.pt[k],re.phi[k]));
      }  
    }
    if (rm.size>0) {
      h_m0_met->Fill(met);
      for(unsigned int k(0);k<rm.size;k++) {
        h_m0_pt->Fill(rm.pt[k]);
        h_m0_mt->Fill(mt(met,metPhi,rm.pt[k],rm.phi[k]));
      }  
    }

    ///HLT/////std::cout << hlt << std::endl;///
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;
 
    if (rp.size>0) {
      for (unsigned int k(0);k<rp.size;k++) h_p1_pt->Fill(rp.pt[k]);   
      h_p1_met->Fill(met);
    }

    if (re.size>0) {
      h_e1_met->Fill(met);
      for (unsigned int k(0);k<re.size;k++) {
        h_e1_pt->Fill(re.pt[k]);
        h_e1_mt->Fill(mt(met,metPhi,re.pt[k],re.phi[k]));
      }
    }

    if (rm.size>0) {
      h_m1_met->Fill(met);
      for(unsigned int k(0);k<rm.size;k++) {
        h_m1_pt->Fill(rm.pt[k]);
        h_m1_mt->Fill(mt(met,metPhi,rm.pt[k],rm.phi[k]));
      }  
    }

    ///Fake-Electron-Muon///
    if (fe.size()>0) {
      for (std::set<unsigned int>::iterator it=fe.begin(); it!=fe.end(); ++it) {
        h_mt0->Fill(mt(met,metPhi,re.pt[*it],re.phi[*it]));
      }
    }

    if (fm.size()>0) {
      for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
        h_mt1->Fill(mt(met,metPhi,rm.pt[*it],rm.phi[*it]));
      }
    }
    
    ///Loose-Medium-Tight///
    if (e.size>0) {
       h_e2_met->Fill(met);
       for(unsigned int k(0);k<e.size;k++) {
         h_e2_pt->Fill(e.pt[k]);
         h_e2_mt->Fill(mt(met,metPhi,e.pt[k],e.phi[k]));
       }  
    }

    if (m.size>0) {
      h_m2_met->Fill(met);
      for(unsigned int k(0);k<m.size;k++) {
        h_m2_pt->Fill(m.pt[k]);
        h_m2_mt->Fill(mt(met,metPhi,m.pt[k],m.phi[k]));
      }  
    }


    //EMObject//Fake-Signal
    if(em.size()>0) {
      for (std::set<unsigned int>::iterator it=em.begin(); it!=em.end(); ++it) {
        h_em0_pt->Fill(rp.pt[*it]);
      }
      if(m.size>0) {
        for (std::set<unsigned int>::iterator it=em.begin(); it!=em.end(); ++it) {
          h_em2_pt->Fill(rp.pt[*it]);
        }
      }
      if (ce.size()>0) {
        for (std::set<unsigned int>::iterator it=em.begin(); it!=em.end(); ++it) {
          h_em1_pt->Fill(rp.pt[*it]);
        }
      }
      if(fe.size()>0) {
        h_met2->Fill(met);
        for (std::set<unsigned int>::iterator it=fe.begin(); it!=fe.end(); ++it) {
          h_mt2->Fill(mt(met,metPhi,re.pt[*it],re.phi[*it]));
        }
      }
      if(fm.size()>0) {
        h_met3->Fill(met);  
        for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
          h_mt3->Fill(mt(met,metPhi,rm.pt[*it],rm.phi[*it]));
        }
      }
    }


    //Loose photon requirement
    if(p.size<1) continue;
    h_p2_pt->Fill(p.pt[0]); 
    h_p2_met->Fill(met);  

    if (e.size>0) {
      h_e3_met->Fill(met);
      for(unsigned int k(0);k<e.size;k++) {
        h_e3_pt->Fill(e.pt[k]);
        h_e3_mt->Fill(mt(met,metPhi,e.pt[k],e.phi[k]));
      }
    }

    if (m.size>0) {
      h_m3_met->Fill(met);
      for(unsigned int k(0);k<m.size;k++) {
        h_m3_pt->Fill(m.pt[k]);
        h_m3_mt->Fill(mt(met,metPhi,m.pt[k],m.phi[k]));
      }
    }


    //Signal
    //if(p.pt[0]<70.) continue;
    if(lp70.size()==0) continue;
    h_p3_pt->Fill(p.pt[0]);
    h_p3_met->Fill(met);   

    if(m.size>0) {
      for (unsigned int k(0);k<m.size;k++) {
        h_m4_pt->Fill(m.pt[k]);
        h_m_pt->Fill(m.pt[k]); 
        h_m4_mt->Fill(mt(met,metPhi,m.pt[k],m.phi[k]));
      }
      h_m4_met->Fill(met);  
      for(unsigned int k(0);k<p.size;k++) if(p.pt[k]>=70.)h_p5_pt->Fill(p.pt[k]);
    }

    if(e.size>0) {
      for(unsigned int k(0);k<e.size;k++) {
        h_e4_pt->Fill(e.pt[k]);
        h_e_pt->Fill(e.pt[k]);
        h_e4_mt->Fill(mt(met,metPhi,e.pt[k],e.phi[k]));
      }
      h_e4_met->Fill(met);
    }

    if(ce.size()>0) {
      for (std::set<unsigned int>::iterator it=ce.begin(); it!=ce.end(); ++it) {
        h_e5_pt->Fill(re.pt[*it]);
        h_ce_pt->Fill(re.pt[*it]);
        h_e5_mt->Fill(mt(met,metPhi,re.pt[*it],re.phi[*it]));
      }
      h_e5_met->Fill(met);
      for(unsigned int k(0);k<p.size;k++) if(p.pt[k]>=70.) h_p4_pt->Fill(p.pt[k]);
    }

    //FakeSignal
    if(fe.size()>0) {
      h_met0->Fill(met);
      for (std::set<unsigned int>::iterator it=fe.begin(); it!=fe.end(); ++it) {
        h_mt4->Fill(mt(met,metPhi,re.pt[*it],re.phi[*it]));
      }
      for(unsigned int k(0);k<p.size;k++) if(p.pt[k]>=70.) h_p6_pt->Fill(p.pt[k]);
    }
    
   
    if(fm.size()>0) {
      h_met1->Fill(met);
      for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
        h_mt5->Fill(mt(met,metPhi,re.pt[*it],re.phi[*it]));
      }
      for(unsigned int k(0);k<p.size;k++) if(p.pt[k]>=70.) h_p7_pt->Fill(p.pt[k]);
    }

    //if(e.size <1 && m.size <1) continue;

  }//while

TCanvas *hmet = new TCanvas("hmet", "MET plots",800, 800);
hmet->Divide(2,2);
hmet->cd(1);
h_p0_met->Draw();
h_p1_met->Draw("same");
h_p2_met->Draw("same");
h_p3_met->Draw("same");
hmet->cd(2);
h_e0_met->Draw();
h_e1_met->Draw("same");
h_e2_met->Draw("same");
h_e3_met->Draw("same");
h_e4_met->Draw("same");
h_e5_met->Draw("same");
hmet->cd(3);
h_m0_met->Draw();
h_m1_met->Draw("same");
h_m2_met->Draw("same");
h_m3_met->Draw("same");
h_m4_met->Draw("same");

TCanvas *hpt = new TCanvas("hpt", "Pt plots",800, 800);
hpt->Divide(2,2);
hpt->cd(1);
h_p0_pt->Draw();
h_p1_pt->Draw("same");
h_p2_pt->Draw("same");
h_p3_pt->Draw("same");
h_p4_pt->Draw("same");
h_p5_pt->Draw("same");
hpt->cd(2);
h_e0_pt->Draw();
h_e1_pt->Draw("same");
h_e2_pt->Draw("same");
h_e3_pt->Draw("same");
h_e4_pt->Draw("same");
h_e5_pt->Draw("same");
hpt->cd(3);
h_m0_pt->Draw();
h_m1_pt->Draw("same");
h_m2_pt->Draw("same");
h_m3_pt->Draw("same");
h_m4_pt->Draw("same");
hpt->cd(4);
h_e_pt->Draw();
h_ce_pt->Draw("same");
h_m_pt->Draw("same");

TCanvas *hmt = new TCanvas("hmt", "Mt plots",800, 800);
hmt->Divide(2,2);
hmt->cd(1);
h_e0_mt->Draw();
h_e1_mt->Draw("same");
h_e2_mt->Draw("same");
h_e3_mt->Draw("same");
h_e4_mt->Draw("same");
h_e5_mt->Draw("same");
hmt->cd(2);
h_m0_mt->Draw();
h_m1_mt->Draw("same");
h_m2_mt->Draw("same");
h_m3_mt->Draw("same");
h_m4_mt->Draw("same");

TCanvas *hf = new TCanvas("hf", "Fake-Object plots",800, 800);
hf->Divide(2,2);
hf->cd(1);
h_em0_pt->Draw();
h_em1_pt->Draw("same");
h_em2_pt->Draw("same");
hf->cd(2);
h_met0->Draw();
h_met1->Draw("same");
h_met2->Draw("same");
h_met3->Draw("same");
h_m4_met->Draw("same");
h_e5_met->Draw("same");
hf->cd(3);
h_mt0->Draw();
h_mt1->Draw("same");
h_mt2->Draw("same");
h_mt3->Draw("same");
h_mt4->Draw("same");
h_mt5->Draw("same");
h_e5_mt->Draw("same");
h_m4_mt->Draw("same");
hf->cd(4);


std::cout<< "Total events                               : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                 : "<< nCnt[1] <<std::endl;
std::cout<< "-------------------------------------------: "<<std::endl;
std::cout<< "Photon                                     : "<< nCnt[10] <<std::endl;
std::cout<< "HLT + photon                               : "<< nCnt[11] <<std::endl;
std::cout<< "HLT + Loose photon                         : "<< nCnt[12] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV                 : "<< nCnt[13] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV+ e              : "<< nCnt[14] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV+ m              : "<< nCnt[15] <<std::endl;
std::cout<< "-------------------------------------------: "<<std::endl;
std::cout<< "Electron                                   : "<< nCnt[20] <<std::endl;
std::cout<< "HLT + electron                             : "<< nCnt[21] <<std::endl;
std::cout<< "HLT + medium electron                      : "<< nCnt[22] <<std::endl;
std::cout<< "HLT + 'clean' medium e                     : "<< nCnt[23] <<std::endl;
std::cout<< "HLT + 'clean' medium e + g                 : "<< nCnt[24] <<std::endl;
std::cout<< "HLT + 'clean' medium e + g > 70GeV         : "<< nCnt[25] <<std::endl;
std::cout<< "HLT + medium e + g > 70GeV                 : "<< nCnt[26] <<std::endl;
std::cout<< "-------------------------------------------: "<<std::endl;
std::cout<< "Muon                                       : "<< nCnt[30] <<std::endl;
std::cout<< "HLT + muon                                 : "<< nCnt[31] <<std::endl;
std::cout<< "HLT + tight muon                           : "<< nCnt[32] <<std::endl;
std::cout<< "HLT + tight muon + g                       : "<< nCnt[33] <<std::endl;
std::cout<< "HLT + tight muon + g > 70GeV               : "<< nCnt[34] <<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
hlist.Add(h_p0_met);
hlist.Add(h_p1_met);
hlist.Add(h_p2_met);
hlist.Add(h_p3_met);

hlist.Add(h_e0_met);
hlist.Add(h_e1_met);
hlist.Add(h_e2_met);
hlist.Add(h_e3_met);
hlist.Add(h_e4_met);
hlist.Add(h_e5_met);

hlist.Add(h_m0_met);
hlist.Add(h_m1_met);
hlist.Add(h_m2_met);
hlist.Add(h_m3_met);
hlist.Add(h_m4_met);

hlist.Add(h_p0_pt);
hlist.Add(h_p1_pt);
hlist.Add(h_p2_pt);
hlist.Add(h_p3_pt);
hlist.Add(h_p4_pt);
hlist.Add(h_p5_pt);

hlist.Add(h_e0_pt);
hlist.Add(h_e1_pt);
hlist.Add(h_e2_pt);
hlist.Add(h_e3_pt);
hlist.Add(h_e4_pt);
hlist.Add(h_e5_pt);

hlist.Add(h_m0_pt);
hlist.Add(h_m1_pt);
hlist.Add(h_m2_pt);
hlist.Add(h_m3_pt);
hlist.Add(h_m4_pt);

hlist.Add(h_e_pt);
hlist.Add(h_ce_pt);
hlist.Add(h_m_pt);

hlist.Add(h_e0_mt);
hlist.Add(h_e1_mt);
hlist.Add(h_e2_mt);
hlist.Add(h_e3_mt);
hlist.Add(h_e4_mt);
hlist.Add(h_e5_mt);

hlist.Add(h_m0_mt);
hlist.Add(h_m1_mt);
hlist.Add(h_m2_mt);
hlist.Add(h_m3_mt);
hlist.Add(h_m4_mt);

hlist.Add(h_em0_pt);
hlist.Add(h_em1_pt);
hlist.Add(h_em2_pt);

hlist.Add(h_met0);
hlist.Add(h_met1);
hlist.Add(h_met2);
hlist.Add(h_met3);

hlist.Add(h_mt0);
hlist.Add(h_mt1);
hlist.Add(h_mt2);
hlist.Add(h_mt3);
hlist.Add(h_mt4);
hlist.Add(h_mt5);

TFile fout("PtMtMet_plots.root", "recreate");
hlist.Write();
fout.Close();


}

