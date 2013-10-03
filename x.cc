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

void x()
{

  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("./A1.root/eventVars");
  //eventVars.Add("./*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/A*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/B*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/C*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");
  //eventVars.Add("/home/acalamba/Desktop/d/*.root/eventVars");


  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("./A1.root/selectedObjects");
  //selectedObjects.Add("./*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/A*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/B*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/C*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");
  //selectedObjects.Add("/home/acalamba/Desktop/d/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("./A1.root/allObjects");
  //allObjects.Add("./*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/A*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/B*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/C*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/*.root/allObjects");
  //allObjects.Add("/home/acalamba/Desktop/d/*.root/allObjects");

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
  //float bin[31] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,255,270,300,330,360,390};
  //TH1F* h_met          = new TH1F("h_met", "MET", 30, bin);
  TCanvas *c2h = new TCanvas("c2h", "plots",800, 800);
  TCanvas *x = new TCanvas("x", "plots",800, 800);
  //TH1F* h_mass= new TH1F("h_mass", "M_{l#gamma}", 120, 0., 120.);
  //TH1F* h_mass_fake= new TH1F("h_mass_fake", "M_{l#gamma_{fake}}", 120, 0., 120.);
  //TH1F* h_mass_fake_met= new TH1F("h_mass_fake_met", "M_{#mu#mu}", 120, 0., 120.);
  //TH1F* h_pt_rawmu = new TH1F("h_pt_rawmu", "#mu Pt", 120, 0., 120.);
  TH1F* h_met      = new TH1F("h_met", "MET", 50, 0., 500.);
  TH1F* h_met_s    = new TH1F("h_met_s", "MET", 50, 0., 500.);
  TH1F* h_mt       = new TH1F("h_mt", "Mt", 50, 0., 500.);
  TH1F* h_mt_s     = new TH1F("h_mt_s", "Mt", 50, 0., 500.);
  TH1F* h_jp_pt    = new TH1F("h_jp_pt", "Photon Pt", 50, 0., 500.);
  TH1F* h_j_pt     = new TH1F("h_j_pt", "Jet Pt", 50, 0., 500.);
  TH1F* h_e_pt     = new TH1F("h_e_pt", "Electron Pt", 50, 0., 500.);
  TH1F* h_fe_pt    = new TH1F("h_fe_pt", "Fake Electron Pt", 50, 0., 500.);
  TH1F* h_ce_pt    = new TH1F("h_ce_pt", "Clean Electron Pt", 50, 0., 500.);
  TH1F* h_m_pt     = new TH1F("h_m_pt", "Muon Pt", 50, 0., 500.);
  TH1F* h_lp       = new TH1F("h_lp","Looose Photon Index", 6,0.0,6.0);
  TH1F* h_dre      = new TH1F("h_dre","e#gamma/#mu#gamma #DeltaR", 100,0.,6.);
  TH1F* h_drm      = new TH1F("h_drm","#mu#gamma #DeltaR", 100,0.,6.);
 
  TH1F* h_p_pt_raw = new TH1F("h_p_pt_raw", "Photon Pt", 50, 0., 500.);
  TH1F* h_p_pt_hlt = new TH1F("h_p_pt_hlt", "Photon Pt", 50, 0., 500.);
  TH1F* h_p_pt_lph = new TH1F("h_p_pt_lph", "Photon Pt", 50, 0., 500.);
  TH1F* h_p_pt_l70 = new TH1F("h_p_pt_l70", "Photon Pt", 50, 0., 500.);
  TH1F* h_p_pt_lep = new TH1F("h_p_pt_lep", "Photon Pt", 50, 0., 500.);
  TH1F* h_p_pt_jet = new TH1F("h_p_pt_jet", "Photon Pt", 50, 0., 500.);

  TH1F* h_e_pt_raw = new TH1F("h_e_pt_raw", "Electron Pt", 50, 0., 500.);
  TH1F* h_e_pt_hlt = new TH1F("h_e_pt_hlt", "Electron Pt", 50, 0., 500.);
  TH1F* h_e_pt_lph = new TH1F("h_e_pt_lph", "Electron Pt", 50, 0., 500.);
  TH1F* h_e_pt_l70 = new TH1F("h_e_pt_l70", "Electron Pt", 50, 0., 500.);
  TH1F* h_e_pt_lep = new TH1F("h_e_pt_lep", "Electron Pt", 50, 0., 500.);
  TH1F* h_e_pt_jet = new TH1F("h_e_pt_jet", "Electron Pt", 50, 0., 500.);

  TH1F* h_m_pt_raw = new TH1F("h_m_pt_raw", "Muon Pt", 50, 0., 500.);
  TH1F* h_m_pt_hlt = new TH1F("h_m_pt_hlt", "Muon Pt", 50, 0., 500.);
  TH1F* h_m_pt_lph = new TH1F("h_m_pt_lph", "Muon Pt", 50, 0., 500.);
  TH1F* h_m_pt_l70 = new TH1F("h_m_pt_l70", "Muon Pt", 50, 0., 500.);
  TH1F* h_m_pt_lep = new TH1F("h_m_pt_lep", "Muon Pt", 50, 0., 500.);
  TH1F* h_m_pt_jet = new TH1F("h_m_pt_jet", "Muon Pt", 50, 0., 500.);
 
  unsigned int nCnt[20] = {0};
  //loop over events
  long iEntry = 0;//change to 0
  while(allObjects.GetEntry(iEntry++) != 0){
    h_met->Fill(met);

    if(p.size>0 && e.size>0) {
       for (int k(0);k<p.size;k++) {
          for (int j(0);j<e.size;j++) {
              float dn = p.eta[k] - e.eta[j];
              //float dp = p.phi[k] - e.phi[j];
              float dp = fabs(fabs(fabs(p.phi[k] - e.phi[j]) - PI) - PI);
              h_dre->Fill(sqrt(dp*dp+dn*dn));
          }
       }
    }

    if(p.size>0 && m.size>0) {
       for (int k(0);k<p.size;k++) {
          for (int j(0);j<m.size;j++) {
              float dn = p.eta[k] - m.eta[j];
              float dp = fabs(fabs(fabs(p.phi[k] - m.phi[j]) - PI) - PI);
              h_drm->Fill(sqrt(dp*dp+dn*dn));
          }
       }
    }

    for (int k = 0;k<rp.size;k++) {
        h_p_pt_raw->Fill(rp.pt[k]);
    }
    for (int k = 0;k<re.size;k++) {
        h_e_pt_raw->Fill(re.pt[k]);
    }
    for (int k = 0;k<rm.size;k++) {
        h_m_pt_raw->Fill(rm.pt[k]);
    }


    nCnt[0]++;
    //std::cout << hlt << std::endl;
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;
    for (int k =0;k<rp.size;k++) {
        h_p_pt_hlt->Fill(rp.pt[k]);    
    }
    for (int k = 0;k<re.size;k++) {
        h_e_pt_hlt->Fill(re.pt[k]);
    }
    for (int k = 0;k<rm.size;k++) {
        h_m_pt_hlt->Fill(rm.pt[k]);
    }


    std::set<int> yp = YLoosePhoton(rp);
    std::set<int> lp = LoosePhoton(rp);
    std::set<int> ep = ElectronFakePhoton(rp);
    std::set<int> jp = JetFakePhoton(rp);
    std::set<int> lj = LooseJet(rj);
    std::set<int> me = MediumElectron(re);

    if(j.size > 0) nCnt[2]++;
    if(e.size > 0) nCnt[3]++;
    if(m.size > 0) nCnt[4]++;
    if(lp.size()>0)   nCnt[5]++;
    if(yp.size()>0)   nCnt[6]++;
    if(p.size > 0) nCnt[7]++;
    
    if(p.size<1) continue;
    h_p_pt_lph->Fill(p.pt[0]);    
    for (int k = 0;k<re.size;k++) {
        h_e_pt_lph->Fill(re.pt[k]);
    }
    for (int k = 0;k<rm.size;k++) {
        h_m_pt_lph->Fill(rm.pt[k]);
    }

    if(p.pt[0]<70.) continue;
    h_p_pt_l70->Fill(p.pt[0]);    
    nCnt[8]++;
    set<int> fakeE;
    fakeE = FakeElectron(re);
    for (std::set<int>::iterator it=fakeE.begin(); it!=fakeE.end(); ++it)
    h_fe_pt->Fill(re.pt[*it]);

    if(e.size>0) {
      for(int k(0);k<e.size;k++) h_e_pt->Fill(e.pt[k]);
    }

    for (int k = 0;k<re.size;k++) {
        h_e_pt_l70->Fill(re.pt[k]);
    }
    for (int k = 0;k<rm.size;k++) {
        h_m_pt_l70->Fill(rm.pt[k]);
    }


    if(e.size <1 && m.size <1) continue;
    if (e.size>0) {
       float mt = sqrt(2.*met*e.pt[0]*(1-cos(metPhi-e.phi[0])));
       h_mt->Fill(mt);
       for(int k(0);k<e.size;k++) {
          h_e_pt_lep->Fill(e.pt[k]);  
       }  
    }

    h_met->Fill(met);
    nCnt[9]++;
    h_p_pt_lep->Fill(p.pt[0]);    
    h_m_pt_lep->Fill(m.pt[0]);    
    h_met_s->Fill(met);

    if(j.size < 1) continue;
    nCnt[10]++;
    h_p_pt_jet->Fill(p.pt[0]);    
    h_e_pt_lep->Fill(e.pt[0]);    
    h_m_pt_lep->Fill(m.pt[0]);    
    //cout<<met<<endl;
  }//while
x->Divide(2,1);
x->cd(1);
h_dre->Draw();
h_drm->Draw("same");
x->cd(2);
h_e_pt->Draw();
h_fe_pt->Draw("same");
/*
//h_met->Scale(1/1000.);
//h_met_fake->Draw("same");
//h_mass->Draw();
//h_mass_fake->Draw();
//h_met->Draw();
c2h->Divide(2,2);
c2h->cd(1);
h_p_pt_raw->Draw(); 
h_p_pt_hlt->Draw("same");
h_p_pt_lph->Draw("same");
h_p_pt_l70->Draw("same");
h_p_pt_lep->Draw("same");
h_p_pt_jet->Draw("same");
//h_ep_pt->SetFillColor(kBlue); // Fill fill color to yellow
//h_ep_pt->Draw("same");
//h_lp_pt->SetMarkerColor(kYellow);
//h_lp_pt->SetFillColor(kYellow); // Fill fill color to yellow
//h_jp_pt->SetFillColor(kGreen); // Fill fill color to yellow
//h_jp_pt->SetMarkerColor(kGreen);
//h_jp_pt->Draw("same");
c2h->cd(2);
//float n1 = h_met->Integral();
//float n2 = h_met_s->Integral();
//h_met->Scale(1/n1);
//h_met_s->Scale(1/n2);
h_met->Draw();
h_met_s->Draw("same");
h_mt->Draw("same");
c2h->cd(3);
h_e_pt_raw->Draw();
h_e_pt_hlt->Draw("same");
h_e_pt_lph->Draw("same");
h_e_pt_l70->Draw("same");
h_e_pt_lep->Draw("same");
h_e_pt_jet->Draw("same");
c2h->cd(4);
h_m_pt_raw->Draw();
h_m_pt_hlt->Draw("same");
h_m_pt_lph->Draw("same");
h_m_pt_l70->Draw("same");
h_m_pt_lep->Draw("same");
h_m_pt_jet->Draw("same");
*/
std::cout<< "Total events                               : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                 : "<< nCnt[1] <<std::endl;
std::cout<< "Loose Jet                                  : "<< nCnt[2] <<std::endl;
std::cout<< "Medium electron                            : "<< nCnt[3] <<std::endl;
std::cout<< "Tight muon                                 : "<< nCnt[4] <<std::endl;
std::cout<< "Loose photon                               : "<< nCnt[5] <<std::endl;
std::cout<< "Y Loose photon                             : "<< nCnt[6] <<std::endl;
std::cout<< "Selected photon                            : "<< nCnt[7] <<std::endl;
std::cout<< "Loose photon w/ pt > 70                    : "<< nCnt[8] <<std::endl;
std::cout<< "Medium electron/Tight mu & Loose Photon 70 : "<< nCnt[9] <<std::endl;
std::cout<< "Medium e/Tight mu & Loose Photon 70 & jet  : "<< nCnt[10] <<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
hlist.Add( h_met );

TFile fout("outputA_sampleAnalysis.root", "recreate");
hlist.Write();
fout.Close();


}

