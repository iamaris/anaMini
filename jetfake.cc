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

void jetfake()
{

  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("/export/cmss/acalamba/photonhad/A*1.root/eventVars");
  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("/export/cmss/acalamba/photonhad/A*1.root/selectedObjects");
  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("/export/cmss/acalamba/photonhad/A*1.root/allObjects");

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


  TH1F* h_lp_pt = new TH1F("h_lp_pt", "Loose photon pt", 25, 0., 500.);
  TH1F* h_jf1_pt = new TH1F("h_jf_pt1", "charged hadron iso", 25, 0., 500.);
  TH1F* h_jf2_pt = new TH1F("h_jf_pt2", "neutral hadron iso", 25, 0., 500.);
  TH1F* h_jf3_pt = new TH1F("h_jf_pt3", "photon iso", 25, 0., 500.);
  TH1F* h_lpbar_pt = new TH1F("h_lpbar_pt", "(Raw photon - loose photon) pt", 25, 0., 500.);

  TH1F* h_p0 = new TH1F("h_p0", "photon object pt ", 25, 0., 500.);
  TH1F* h_p1 = new TH1F("h_p1", "pixelseed <= 0 ", 25, 0., 500.);
  TH1F* h_p2 = new TH1F("h_p2", "hOverE <= 0.05", 25, 0., 500.);
  TH1F* h_p3 = new TH1F("h_p3", "SigmaInIn <= 0.012", 25, 0., 500.);
  TH1F* h_p4 = new TH1F("h_p4", "chaHadIso <= 2.6", 25, 0., 500.);
  TH1F* h_p5 = new TH1F("h_p5", "neuHadIso <= 3.5", 25, 0., 500.);
  TH1F* h_p6 = new TH1F("h_p6", "photonIso <= 1.3", 25, 0., 500.);
  TH1F* h_p11 = new TH1F("h_p11", "pixelseed = 0 ", 25, 0., 500.);
  TH1F* h_p22 = new TH1F("h_p22", "hOverE > 0.05", 25, 0., 500.);
  TH1F* h_p33 = new TH1F("h_p33", "SigmaInIn > 0.012", 25, 0., 500.);
  TH1F* h_p44 = new TH1F("h_p44", "chaHadIso > 2.6", 25, 0., 500.);
  TH1F* h_p55 = new TH1F("h_p55", "neuHadIso > 3.5", 25, 0., 500.);
  TH1F* h_p66 = new TH1F("h_p66", "photonIso > 1.3", 25, 0., 500.);


  unsigned int nCnt[30] = {0};  
  unsigned int np[20] = {0};  
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;
    std::set<unsigned int> ce = MediumElectron(re,p);
    std::set<unsigned int> jf = JetPhoton(rp,j);
    std::set<unsigned int> lp = LoosePhoton(rp,0,true);
    std::set<unsigned int> lpbar;

    ///HLT
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;

    for (std::set<unsigned int>::iterator it=jf.begin();it!=jf.end();++it) {
      if(rp.iSubdet[*it]==0) {
         h_p0->Fill(rp.pt[*it]);
        np[19]++;
        if(rp.nPixelSeeds[*it] > 0) {
          h_p11->Fill(rp.pt[*it]);
          np[0]++;
        } else {
          h_p1->Fill(rp.pt[*it]);
          np[1]++;
        }
        if(rp.hOverE[*it] > 0.05) {
          h_p22->Fill(rp.pt[*it]);
          np[2]++;
        } else {
          h_p2->Fill(rp.pt[*it]);
          np[3]++;
        }
        if(rp.sigmaIetaIeta[*it] > 0.012) {
          h_p33->Fill(rp.pt[*it]);
          np[4]++;
        } else {
          h_p3->Fill(rp.pt[*it]);
          np[5]++;
        }
        if(rp.chargedHadronIso[*it] > 2.6) {
          h_p44->Fill(rp.pt[*it]);
          np[6]++;
        } else {
          h_p4->Fill(rp.pt[*it]);
          np[7]++;
        }
        if(rp.neutralHadronIso[*it] > 3.5) {
          h_p55->Fill(rp.pt[*it]);
          np[8]++;
        } else {
          h_p5->Fill(rp.pt[*it]);
          np[9]++;
        }
        if(rp.photonIso[*it] > 1.3) {
          h_p66->Fill(rp.pt[*it]);
          np[10]++;
        } else {
          h_p6->Fill(rp.pt[*it]);
          np[11]++;
        }
      } 
    }

    //if ((lpbar.size()+lp.size())!=rp.size) std::cout<<"OH NO!";

    if(lpbar.size()>0) {
      for (std::set<unsigned int>::iterator it=lpbar.begin();it!=lpbar.end();++it) {
        h_lpbar_pt->Fill(rp.pt[*it]);
      }
    }

    if(lp.size()>0) {
      for (std::set<unsigned int>::iterator it=lp.begin();it!=lp.end();++it) {
        h_lp_pt->Fill(rp.pt[*it]);
      }
    }

    if(m.size<1 && ce.size()<1) continue;

    if(ce.size()>0) nCnt[4]++;
    if(m.size>0) nCnt[5]++;

  }//while

TCanvas *c0 = new TCanvas("c0", "pt",800, 800);
c0->Divide(2,2);
c0->cd(1);
h_p0->Draw();
h_lp_pt->Draw("same");
h_lpbar_pt->Draw("same");
h_p1->Draw("same");
h_p2->Draw("same");
h_p3->Draw("same");
h_p4->Draw("same");
h_p5->Draw("same");
h_p6->Draw("same");
c0->cd(2);
h_p0->Draw();
h_lp_pt->Draw("same");
h_lpbar_pt->Draw("same");
h_p11->Draw("same");
h_p22->Draw("same");
h_p33->Draw("same");
h_p44->Draw("same");
h_p55->Draw("same");
h_p66->Draw("same");

TCanvas *canvas = new TCanvas("canvas", "pt",800, 800);
canvas->Divide(3,3);
canvas->cd(1);
h_lp_pt->Draw();
h_lpbar_pt->Draw("same");
h_jf1_pt->Draw("same");
canvas->cd(2);
h_lp_pt->Draw();
h_lpbar_pt->Draw("same");
h_jf2_pt->Draw("same");
canvas->cd(3);
h_lp_pt->Draw();
h_lpbar_pt->Draw("same");
h_jf3_pt->Draw("same");
canvas->cd(4);
TH1F *e1 = (TH1F*)h_jf1_pt->Clone("e1");
e1->Divide(h_lp_pt);
e1->Draw("e");
canvas->cd(5);
TH1F *e2 = (TH1F*)h_jf2_pt->Clone("e2");
e2->Divide(h_lp_pt);
e2->Draw("e");
canvas->cd(6);
TH1F *e3 = (TH1F*)h_jf3_pt->Clone("e3");
e3->Divide(h_lp_pt);
e3->Draw("e");


std::cout<< "# raw photon                      : "<< np[19] <<std::endl;
std::cout<< "# raw photon / pixelseed > 0      : "<< np[0] << "  -  "<< (float)np[0]/np[19] <<std::endl;
std::cout<< "# raw photon / pixelseed <= 0     : "<< np[1] << "  -  "<< (float)np[1]/np[19] <<std::endl;
std::cout<< "# raw photon / hOverE > 0.05      : "<< np[2] << "  -  "<< (float)np[2]/np[19] <<std::endl;
std::cout<< "# raw photon / hOverE <= 0.05     : "<< np[3] << "  -  "<< (float)np[3]/np[19] <<std::endl;
std::cout<< "# raw photon / SigmaInIn > 0.012  : "<< np[4] << "  -  "<< (float)np[4]/np[19] <<std::endl;
std::cout<< "# raw photon / SigmaInIn <= 0.012 : "<< np[5] << "  -  "<< (float)np[5]/np[19] <<std::endl;
std::cout<< "# raw photon / chaHadIso >  2.6   : "<< np[6] << "  -  "<< (float)np[6]/np[19] <<std::endl;
std::cout<< "# raw photon / chaHadIso <= 2.6   : "<< np[7] << "  -  "<< (float)np[7]/np[19] <<std::endl;
std::cout<< "# raw photon / neuHadIso  > 3.5   : "<< np[8] << "  -  "<< (float)np[8]/np[19] <<std::endl;
std::cout<< "# raw photon / neuHadIso <= 3.5   : "<< np[9] << "  -  "<< (float)np[9]/np[19] <<std::endl;
std::cout<< "# raw photon / photonIso  > 1.3   : "<< np[10] << "  -  "<< (float)np[10]/np[19] <<std::endl;
std::cout<< "# raw photon / photonIso <= 1.3   : "<< np[11] << "  -  "<< (float)np[11]/np[19] <<std::endl;

std::cout<< "Total events                               : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                 : "<< nCnt[1] <<std::endl;
std::cout<< "-------------------------------------------: "<<std::endl;
std::cout<< "HLT + Loose photon                         : "<< nCnt[2] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV                 : "<< nCnt[3] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV+ e              : "<< nCnt[4] <<std::endl;
std::cout<< "HLT + loose photon > 70GeV+ m              : "<< nCnt[5] <<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);


TFile fout("ttbar.root", "recreate");
hlist.Write();
fout.Close();


}

