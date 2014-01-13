final() { 

//TFile f("st_plots_photonhad.root","READ");
 TFile *f = new TFile("st_plots_photonhad.root");
 //f->ls();
 //TH1F *h1 = (TH1F*)f->Get("h_st1"); 

 //TCanvas *hh = new TCanvas("hh", "St plots",800, 800);
 //hh->Divide(2,1);
 //hh->cd(1);
 //h_st3->Draw();
 h_st1->Draw();
 //h_st2->Draw("same");
 //h_st4->Draw("same");
 //h_st5->Draw("same");
 //h_st6->Draw("same");
 //h_st7->Draw("same");
 //h_st8->Draw("same");


 f->Close();

}
