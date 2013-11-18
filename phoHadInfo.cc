#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <iostream>
#include "x.h"
#include "TCanvas.h"
#include <math.h>       
#include <set>
#include <iomanip>  
using namespace mini;
const double PI = 4.0*atan(1.0); 

void phoHadInfo()
{
  TChain eventVars("eventVars");//define eventVars tree
  eventVars.Add("./*.root/eventVars");
  //eventVars.Add("/export/cmss/acalamba/photonhad/D*.root/eventVars");

  TChain selectedObjects("selectedObjects");  //define eventVars tree
  selectedObjects.Add("./*.root/selectedObjects");
  //selectedObjects.Add("/export/cmss/acalamba/photonhad/D*.root/selectedObjects");

  TChain allObjects("allObjects");  //define eventVars tree
  allObjects.Add("./*.root/allObjects");
  //allObjects.Add("/export/cmss/acalamba/photonhad/D*.root/allObjects");

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
 
  //selectedObjects.Print();
  //define histograms
  //TH1F* h_st = new TH1F("h_st", "St", 600, 0., 6000.);

  unsigned nCnt[90] = {0};
  unsigned ee(0),mm(0),jj(0),pp(0);
  //loop over events
  long iEntry = 0;
  while(allObjects.GetEntry(iEntry++) != 0){
    unsigned n(0);
    nCnt[n++]++;

    if(!hlt0 && !hlt1) continue;
    nCnt[n++]++;

    std::set<unsigned> ljB;
    std::set<unsigned> ljE;
        
    std::set<unsigned> lmB; 
    std::set<unsigned> lmE; 
    std::set<unsigned> tmB;    
    std::set<unsigned> tmE; 

    std::set<unsigned> leB;
    std::set<unsigned> leE;
    std::set<unsigned> meB;          
    std::set<unsigned> meE;    
    std::set<unsigned> teB;
    std::set<unsigned> teE;

    std::set<unsigned> lcB;
    std::set<unsigned> lcE;
    std::set<unsigned> mcB;//medium clean electron barrel
    std::set<unsigned> mcE;
    std::set<unsigned> tcB;//tight clean electron barrel
    std::set<unsigned> tcE;              

    std::set<unsigned> lpB;
    std::set<unsigned> lpE;
    std::set<unsigned> mpB;
    std::set<unsigned> mpE;
    std::set<unsigned> tpB; 
    std::set<unsigned> tpE;

   
    //jet
    for(unsigned i(0);i<j.size;i++) {
      if(j.iSubdet[i]==0) {//barrel
        if(j.isLoose[i]) ljB.insert(i);
      } else if (j.iSubdet[i]==1) {//endcap
        if(j.isLoose[i]) ljE.insert(i);
      }
    }
    if(j.size>0) jj++;

    //muon 
    if (m.size>0)
    for(unsigned i(0);i<m.size;i++) {
      if(m.iSubdet[i]==0) {//barrel
        if(m.isLoose[i]) lmB.insert(i);
        if(m.isTight[i]) tmB.insert(i);  
      } else if (m.iSubdet[i]==1) {//endcap
        if(m.isLoose[i]) lmE.insert(i);
        if(m.isTight[i]) tmE.insert(i);  
      }
    }
    if(m.size>0) mm++;

    //electron
    for(unsigned i(0);i<e.size;i++) {
      int dirty(0);
      for (unsigned k=0;k<p.size;k++) {
        if(mini::deltaR(p.eta[k],p.phi[k],e.eta[i],e.phi[i])<0.3) {
          dirty++; 
          break;
        }
      }

      if(e.iSubdet[i]==0) {//barrel
        if(e.isLoose[i])  leB.insert(i);
        if(e.isMedium[i]) meB.insert(i);
        if(e.isTight[i])  teB.insert(i);
        if(!dirty) {
          if(e.isLoose[i])  lcB.insert(i);
          if(e.isMedium[i]) mcB.insert(i);
          if(e.isTight[i])  tcB.insert(i);
        }
      } else if (e.iSubdet[i]==1) {//endcap
        if(e.isLoose[i])  leE.insert(i);
        if(e.isMedium[i]) meE.insert(i);
        if(e.isTight[i])  teE.insert(i);
        if(!dirty) {
          if(e.isLoose[i])  lcE.insert(i);
          if(e.isMedium[i]) mcE.insert(i);
          if(e.isTight[i])  tcE.insert(i);
        }
      }
    }
    if(e.size>0) ee++;

    //photon
    for(unsigned i(0);i<p.size;i++) {
      if(p.iSubdet[i]==0) {//barrel
        if(p.isLoose[i])  lpB.insert(i);
        if(p.isMedium[i]) mpB.insert(i);
        if(p.isTight[i])  tpB.insert(i);
      } else if (p.iSubdet[i]==1) {//endcap
        if(p.isLoose[i])  lpE.insert(i);
        if(p.isMedium[i]) mpE.insert(i);
        if(p.isTight[i])  tpE.insert(i);
      }
    }
    if(p.size>0) pp++;


    if(ljB.size()>0) nCnt[n]++;
    n++;
    if(ljE.size()>0) nCnt[n]++; 
    n++;
    if(lmB.size()>0) nCnt[n]++;
    n++;
    if(tmB.size()>0) nCnt[n]++;
    n++;
    if(lmE.size()>0) nCnt[n]++;
    n++;
    if(tmE.size()>0) nCnt[n]++;
    n++;    
    if(leB.size()>0) nCnt[n]++;
    n++;
    if(leE.size()>0) nCnt[n]++;
    n++;
    if(meB.size()>0) nCnt[n]++;      
    n++;
    if(meE.size()>0) nCnt[n]++;
    n++;
    if(teB.size()>0) nCnt[n]++;      
    n++;
    if(teE.size()>0) nCnt[n]++;
    n++;    
    if(lcB.size()>0) nCnt[n]++;
    n++;
    if(lcE.size()>0) nCnt[n]++;
    n++;
    if(mcB.size()>0) nCnt[n]++;      
    n++;
    if(mcE.size()>0) nCnt[n]++;
    n++;
    if(tcB.size()>0) nCnt[n]++;      
    n++;
    if(tcE.size()>0) nCnt[n]++;    
    n++;    
    if(lpB.size()>0) nCnt[n]++;
    n++;
    if(lpE.size()>0) nCnt[n]++;
    n++;
    if(mpB.size()>0) nCnt[n]++;      
    n++;
    if(mpE.size()>0) nCnt[n]++;
    n++;
    if(tpB.size()>0) nCnt[n]++;      
    n++;
    if(tpE.size()>0) nCnt[n]++;   
    n++;

    if(lpB.size()!=1&&lpE.size()!=0) continue;
    if(ljB.size()>0) nCnt[n]++;
    n++;
    if(ljE.size()>0) nCnt[n]++; 
    n++;
    if(lmB.size()>0) nCnt[n]++;
    n++;
    if(tmB.size()>0) nCnt[n]++;
    n++;
    if(lmE.size()>0) nCnt[n]++;
    n++;
    if(tmE.size()>0) nCnt[n]++;
    n++;    
    if(leB.size()>0) nCnt[n]++;
    n++;
    if(leE.size()>0) nCnt[n]++;
    n++;
    if(meB.size()>0) nCnt[n]++;      
    n++;
    if(meE.size()>0) nCnt[n]++;
    n++;
    if(teB.size()>0) nCnt[n]++;      
    n++;
    if(teE.size()>0) nCnt[n]++;
    n++;    
    if(lcB.size()>0) nCnt[n]++;
    n++;
    if(lcE.size()>0) nCnt[n]++;
    n++;
    if(mcB.size()>0) nCnt[n]++;      
    n++;
    if(mcE.size()>0) nCnt[n]++;
    n++;
    if(tcB.size()>0) nCnt[n]++;      
    n++;
    if(tcE.size()>0) nCnt[n]++;    
    n++;   

  }//while

int q(0);
std::cout<< "Total events                                 : "<<setw(7)<< nCnt[q] <<std::endl;
std::cout<< "HLT passed                                   : "<<setw(7)<< nCnt[++q] <<std::endl;
std::cout<< "Event with selected jet                      : "<<setw(7)<< jj <<std::endl;
std::cout<< "Event with selected muon                     : "<<setw(7)<< mm <<std::endl;
std::cout<< "Event with selected electron                 : "<<setw(7)<< ee <<std::endl;
std::cout<< "Event with selected photon                   : "<<setw(7)<< pp <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose jet Barrel                          : "<<setw(7)<< nCnt[++q]<<" , "<<setprecision(3)<<nCnt[2]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose jet EndCap                          : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[3]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose muon Barrel                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[4]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight muon Barrel                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[5]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose muon EndCap                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[6]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight muon EndCap                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[7]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose electron Barrel                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[8]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose electron EndCap                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[9]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron Barrel                    : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[10]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron EndCap                    : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[11]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron Barrel                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[12]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron EndCap                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[13]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose 'clean' electron Barrel             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[14]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose 'clean' electron EndCap             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[15]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron Barrel            : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[16]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron EndCap            : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[17]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron Barrel             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[18]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron EndCap             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[19]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose photon Barrel                       : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[20]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose photon EndCap                       : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[21]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium photon Barrel                      : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[22]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium photon EndCap                      : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[23]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight photon Barrel                       : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[24]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight photon EndCap                       : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[25]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "---Single loose barrel photon w/ Pt > 70 required-------------"<<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose jet Barrel                          : "<<setw(7)<< nCnt[++q]<<" , "<<setprecision(3)<<nCnt[26]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose jet EndCap                          : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[27]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose muon Barrel                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[28]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight muon Barrel                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[29]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose muon EndCap                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[30]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight muon EndCap                         : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[31]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose electron Barrel                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[32]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose electron EndCap                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[33]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron Barrel                    : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[34]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium electron EndCap                    : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[35]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron Barrel                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[36]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight electron EndCap                     : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[37]/float(nCnt[1]) <<std::endl;
std::cout<< "--------------------------------------------------------------"<<std::endl;
std::cout<< "w/ loose 'clean' electron Barrel             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[38]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ loose 'clean' electron EndCap             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[39]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron Barrel            : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[40]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ medium 'clean' electron EndCap            : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[41]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron Barrel             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[42]/float(nCnt[1]) <<std::endl;
std::cout<< "w/ tight 'clean' electron EndCap             : "<<setw(7)<< nCnt[++q]<<" , "<<nCnt[43]/float(nCnt[1]) <<std::endl;

//save histograms inside sampleAnalysis.root
TObjArray hlist(0);
//hlist.Add( h_mt );

TFile fout("info.root", "recreate");
hlist.Write();
fout.Close();



}

