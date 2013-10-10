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
  //eventVars.Add("./A1.root/eventVars");
  //eventVars.Add("./*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/A*.root/eventVars");
  eventVars.Add("/export/cmss/acalamba/photonhad/B*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/C*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/*.root/eventVars");
  //eventVars.Add("/home/acalamba/Desktop/d/*.root/eventVars");


  TChain selectedObjects("selectedObjects");  //define eventVars tree
  //selectedObjects.Add("./A1.root/selectedObjects");
  //selectedObjects.Add("./*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/A*.root/selectedObjects");
  selectedObjects.Add("/export/cmss/acalamba/photonhad/B*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/C*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/*.root/selectedObjects");
  //selectedObjects.Add("/home/acalamba/Desktop/d/*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  //allObjects.Add("./A1.root/allObjects");
  //allObjects.Add("./*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/A*.root/allObjects");
  allObjects.Add("/export/cmss/acalamba/photonhad/B*.root/allObjects");
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
  //TCanvas *c2h = new TCanvas("c2h", "plots",800, 800);
  //TCanvas *x = new TCanvas("x", "plots",800, 800);
  //TH1F* h_mass= new TH1F("h_mass", "M_{l#gamma}", 120, 0., 120.);
  //TH1F* h_mass_fake= new TH1F("h_mass_fake", "M_{l#gamma_{fake}}", 120, 0., 120.);
  //TH1F* h_mass_fake_met= new TH1F("h_mass_fake_met", "M_{#mu#mu}", 120, 0., 120.);
  //TH1F* h_pt_rawmu = new TH1F("h_pt_rawmu", "#mu Pt", 120, 0., 120.);
  TH1F* h_met      = new TH1F("h_met", "MET", 50, 0., 500.);
  TH1F* h_met_s    = new TH1F("h_met_s", "MET", 50, 0., 500.);
  TH1F* h_met_e    = new TH1F("h_met_e", "MET", 50, 0., 500.);
  TH1F* h_met_m    = new TH1F("h_met_m", "MET", 50, 0., 500.);
  TH1F* h_met_fe    = new TH1F("h_met_fe", "MET", 50, 0., 500.);
  TH1F* h_met_fm    = new TH1F("h_met_fm", "MET", 50, 0., 500.);
  TH1F* h_mt       = new TH1F("h_mt", "Mt", 50, 0., 500.);
  TH1F* h_mt_e       = new TH1F("h_mt_e", "Mt", 50, 0., 500.);
  TH1F* h_mt_m       = new TH1F("h_mt_m", "Mt", 50, 0., 500.);
  TH1F* h_mt_fe       = new TH1F("h_mt_fe", "Mt", 50, 0., 500.);
  TH1F* h_mt_fm       = new TH1F("h_mt_fm", "Mt", 50, 0., 500.);
  TH1F* h_mt_s     = new TH1F("h_mt_s", "Mt", 50, 0., 500.);
  TH1F* h_jp_pt    = new TH1F("h_jp_pt", "Photon Pt", 50, 0., 500.);
  TH1F* h_j_pt     = new TH1F("h_j_pt", "Jet Pt", 50, 0., 500.);
  TH1F* h_e_pt     = new TH1F("h_e_pt", "Electron Pt", 50, 0., 500.);
  TH1F* h_fe_pt    = new TH1F("h_fe_pt", "Fake Electron Pt", 50, 0., 500.);
  TH1F* h_ce_pt    = new TH1F("h_ce_pt", "Clean Electron Pt", 50, 0., 500.);
  TH1F* h_m_pt     = new TH1F("h_m_pt", "Muon Pt", 50, 0., 500.);
  TH1F* h_fm_pt    = new TH1F("h_fm_pt", "Fake Muon Pt", 50, 0., 500.);
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
 
  unsigned int nCnt[90] = {0};
  //loop over events
  long iEntry = 0;//change to 0
  while(allObjects.GetEntry(iEntry++) != 0){
    h_met->Fill(met);

    if(p.size>0 && e.size>0) {
       for (unsigned int k(0);k<p.size;k++) {
          for (unsigned int i(0);i<e.size;i++) {
              float dn = p.eta[k] - e.eta[i];
              //float dp = p.phi[k] - e.phi[j];
              float dp = fabs(fabs(fabs(p.phi[k] - e.phi[i]) - PI) - PI);
              h_dre->Fill(sqrt(dp*dp+dn*dn));
          }
       }
    }

    if(p.size>0 && m.size>0) {
       for (unsigned int k(0);k<p.size;k++) {
          for (unsigned int i(0);i<m.size;i++) {
              float dn = p.eta[k] - m.eta[i];
              float dp = fabs(fabs(fabs(p.phi[k] - m.phi[i]) - PI) - PI);
              h_drm->Fill(sqrt(dp*dp+dn*dn));
          }
       }
    }

    for (unsigned int k = 0;k<rp.size;k++) {
        h_p_pt_raw->Fill(rp.pt[k]);
    }
    for (unsigned int k = 0;k<re.size;k++) {
        h_e_pt_raw->Fill(re.pt[k]);
    }
    for (unsigned int k = 0;k<rm.size;k++) {
        h_m_pt_raw->Fill(rm.pt[k]);
    }


    nCnt[0]++;
    //std::cout << hlt << std::endl;
    if(!hlt0 && !hlt1) continue;
    nCnt[1]++;
    for (unsigned int k =0;k<rp.size;k++) {
        h_p_pt_hlt->Fill(rp.pt[k]);    
    }
    for (unsigned int k = 0;k<re.size;k++) {
        h_e_pt_hlt->Fill(re.pt[k]);
    }
    for (unsigned int k = 0;k<rm.size;k++) {
        h_m_pt_hlt->Fill(rm.pt[k]);
    }


    std::set<unsigned int> yp = YLoosePhoton(rp);
    std::set<unsigned int> ep = ElectronFakePhoton(rp);
    std::set<unsigned int> jp = JetFakePhoton(rp);
    std::set<unsigned int> lj = LooseJet(rj);
    std::set<unsigned int> me = MediumElectron(re);
    std::set<unsigned int> tm = TightMuon(rm);
    std::set<unsigned int> lp = LoosePhoton(rp);
    std::set<unsigned int> ce = MediumElectron(re,p);
    std::set<unsigned int> fe = FakeElectron(re);
    std::set<unsigned int> fm = FakeMuon(rm);
    std::set<unsigned int> em = EMObject(rp,70.);


    if(j.size > 0) nCnt[2]++;
    if(e.size > 0) nCnt[3]++;
    if(m.size > 0) nCnt[4]++;
    if(p.size > 0) nCnt[5]++;

    if(lj.size()>0) nCnt[6]++;
    if(me.size()>0) nCnt[7]++;
    if(tm.size()>0) nCnt[8]++;
    if(lp.size()>0) nCnt[9]++;
    if(ce.size()>0) nCnt[10]++;

    if(em.size()>0 && fe.size()>0) {
      nCnt[50]++;
      h_met_fe->Fill(met);
      for (std::set<unsigned int>::iterator it=fe.begin(); it!=fe.end(); ++it) {
        h_fe_pt->Fill(re.pt[*it]);
        float mt_fe = sqrt(2.*met*re.pt[*it]*(1-cos(metPhi-re.phi[*it])));
        h_mt_fe->Fill(mt_fe);
      }
    }

    if(em.size()>0 && fm.size()>0) {
      nCnt[60]++;
      h_met_fm->Fill(met);
      for (std::set<unsigned int>::iterator it=fm.begin(); it!=fm.end(); ++it) {
        h_fm_pt->Fill(rm.pt[*it]);
        float mt_fm = sqrt(2.*met*rm.pt[*it]*(1-cos(metPhi-rm.phi[*it])));
        h_mt_fm->Fill(mt_fm);
      }
    }
    
    if(p.size<1) continue;

    if(p.pt[0]<70.) continue;
    h_p_pt_l70->Fill(p.pt[0]);    
    nCnt[30]++;


    if(e.size>0) {
      nCnt[31]++;
      for(unsigned int k(0);k<e.size;k++) h_e_pt->Fill(e.pt[k]);
    }

    if(ce.size()>0) {
      nCnt[32]++;
      h_met_e->Fill(met);
      for (std::set<unsigned int>::iterator it=ce.begin(); it!=ce.end(); ++it) {
        h_ce_pt->Fill(re.pt[*it]);
        float mt_e = sqrt(2.*met*re.pt[*it]*(1-cos(metPhi-re.phi[*it])));
        h_mt_e->Fill(mt_e);
      }
    }

    if(m.size>0) {
      nCnt[40]++;
      h_met_m->Fill(met);
      for(unsigned int k(0);k<m.size;k++) {
        h_m_pt->Fill(m.pt[k]);
        float mt_m = sqrt(2.*met*m.pt[k]*(1-cos(metPhi-m.phi[k])));
        h_mt_m->Fill(mt_m);
      }
    }

  }//while



TCanvas *x = new TCanvas("x", "plots",800, 800);
x->Divide(2,2);
x->cd(1);
h_ce_pt->Draw();
h_fe_pt->Draw("same");
x->cd(2);
h_fm_pt->Draw();
h_m_pt->Draw("same");
x->cd(3);
h_met_e->Draw();
h_met_fe->Draw("same");
x->cd(4);
h_met_fm->Draw();
h_met_m->Draw("same");


TCanvas *y = new TCanvas("y", "plots",800, 800);
y->Divide(2,2);
y->cd(1);
h_mt_e->Draw();
h_mt_fe->Draw("same");
y->cd(2);
h_mt_fm->Draw();
h_mt_m->Draw("same");


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
std::cout<< "Total events                                 : "<< nCnt[0] <<std::endl;
std::cout<< "HLT passed                                   : "<< nCnt[1] <<std::endl;
std::cout<< "-----------------------------------------------"<<std::endl;
std::cout<< "w/ selected jet , loose jet                  : "<< nCnt[2]<<" , "<<nCnt[6] <<std::endl;
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


//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
hlist.Add( h_mt );

TFile fout("outputA_sampleAnalysis.root", "recreate");
hlist.Write();
fout.Close();


}

