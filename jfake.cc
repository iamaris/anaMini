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

void jfake()
{

  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("./A*2.root/eventVars");
  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("./A*2.root/selectedObjects");
  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("./A*2.root/allObjects");

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

//  TH1F* h_p0 = new TH1F("h_p0", "photon object pt ", 25, 0., 500.);
  TH1F* h_sieta0 = new TH1F("h_sieta0", "BARREL: #gamma object #sigma_{i#etai#eta}", 70, 0., 0.035);
  TH1F* h_sieta1 = new TH1F("h_sieta1", "BARREL: #gamma->j #sigma_{i#etai#eta}", 70, 0., 0.035);  
  TH1F* h_sieta2 = new TH1F("h_sieta2", "BARREL: #gamma  #sigma_{i#etai#eta}", 70, 0., 0.035);  
  TH1F* h_sieta3 = new TH1F("h_sieta3", "ENDCAP: #gamma object #sigma_{i#etai#eta}", 70, 0.02, 0.055);
  TH1F* h_sieta4 = new TH1F("h_sieta4", "ENDCAP: #gamma->j #sigma_{i#etai#eta}", 70, 0.02, 0.055);  
  TH1F* h_sieta5 = new TH1F("h_sieta5", "ENDCAP: #gamma  #sigma_{i#etai#eta}", 70, 0.02, 0.055);    
  
  
  unsigned nCnt[30] = {0};  
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    nCnt[0]++;
    //std::set<unsigned> me = MediumElectron(re,p);
    std::set<unsigned> jf = JetFake(rp);     
    ///HLT
    if(!hlt0 && !hlt1) continue;
    for(unsigned i(0);i<rp.size;++i) {
      if(rp.iSubdet[i]==0)
      h_sieta0->Fill(rp.sigmaIetaIeta[i]);
      if(rp.iSubdet[i]==1)
      h_sieta3->Fill(rp.sigmaIetaIeta[i]);      
    }
    

    for (std::set<unsigned>::iterator it=jf.begin();it!=jf.end();++it) {
      if(rp.iSubdet[*it]==0)
      h_sieta1->Fill(rp.sigmaIetaIeta[*it]);
      if(rp.iSubdet[*it]==1)
      h_sieta4->Fill(rp.sigmaIetaIeta[*it]);      
      
    }
    
    if(!p.size) continue;

    for(unsigned i(0);i<p.size;++i) {
      if(p.iSubdet[i]==0)
      h_sieta2->Fill(p.sigmaIetaIeta[i]);
      if(p.iSubdet[i]==1)
      h_sieta5->Fill(p.sigmaIetaIeta[i]);      
    }    

    
  }//while
  
TCanvas *h = new TCanvas("h", "#sigma_{i#etai#eta}",800, 800);
h->Divide(2,2);
h->cd(1);   
h_sieta0->Draw(); 
h_sieta1->Draw("same");
h_sieta2->Draw("same");
h->cd(2);   
h_sieta3->Draw(); 
h_sieta4->Draw("same");
h_sieta5->Draw("same");

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);


TFile fout("test.root", "recreate");
hlist.Write();
fout.Close();


}

