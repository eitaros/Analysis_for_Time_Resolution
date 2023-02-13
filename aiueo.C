#define aiueo_cxx
#include "aiueo.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <ctime>
#include <stdio.h>

void aiueo::Loop(int id)
{
   TH1D *h0[16], *h1[16], *ht0[16], *ht1[16], *ht2[16];
   bool frag[16];
   for(int i=0; i<16; i++){
      frag[i] = false;
   }
   string hname;
   double fitmin[16] = {200, 500, 1000, 600, 1000, 1000, 600, 600, 500, 400, 600, 600, 300, 700, 1000, 0};
   double fitmax[16];
   int bin = 128;
   int bint = 256;
   for(int i=0; i<16; i++){
      hname = "h0" + to_string(i);
      h0[i] = new TH1D(hname.c_str(),hname.c_str(),bin,0,4095);
      hname = "h1" + to_string(i);
      h1[i] = new TH1D(hname.c_str(),hname.c_str(),bin,0,4095);
      hname = "ht0" + to_string(i);
      ht0[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);
      hname = "ht1" + to_string(i);
      ht1[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);
      hname = "ht2" + to_string(i);
      ht2[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);

      fitmax[i] = 3500;
   }

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      for(int i=0; i<16; i++){
         h0[i]->Fill(adc[i]);
         ht0[i]->Fill(tdc[i]);
         if(Hit(4)){
            frag[i] = true;
            h1[i]->Fill(adc[i]);
            ht1[i]->Fill(tdc[i]);
         }
         if(Hit(5)) ht2[i]->Fill(tdc[i]);
      }
   }
   if(frag[12]!=true){
      cout << "no finger scinti signal" << endl;
      return;
   }
   double chi2[16], ndof[16];
   for(int i=0; i<16; i++){
      chi2[i] = 0.; ndof[i] = 0.;
   }


   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(4,4);
   TF1 *func;
   for(int i=0; i<16; i++){
      c1->cd(i+1);
      gPad->SetLogy();
      h0[i]->SetMinimum(0.7);
      h0[i]->SetMaximum(100000);
      h0[i]->Draw();
      if(frag[i]){
         h1[i]->SetLineColor(6);
         h1[i]->Draw("same");
         h1[i]->Fit("landau","","",fitmin[i],fitmax[i]);
         func = h1[i]->GetFunction("landau");
         chi2[i] = func->GetChisquare();
         ndof[i] = func->GetNDF();
      } else {
         chi2[i] = -1;
         ndof[i] = -1;
      }
   }
   char fname[256];
   snprintf(fname,sizeof(fname),"./ana/cut/ana%03d_adc.pdf",ananum);
   c1->Print(fname);

   /* TCanvas *c2 = new TCanvas("c2","c2");
   c2->Divide(4,4);
   TF1 *func2;
   for(int i=0; i<16; i++){
      c2->cd(i+1);
      gPad->SetLogy();
      ht0[i]->SetMinimum(0.7);
      ht0[i]->SetMaximum(10000);
      ht0[i]->Draw();
      if(frag[i]){
         ht1[i]->SetLineColor(6);
         ht1[i]->Draw("same");
         ht1[i]->Fit("gaus");
         func2 = ht1[i]->GetFunction("gaus");
         mean_tdc[i] = func2->GetParameter(1);
         sigma_tdc[i] = func2->GetParameter(2);
      } else {
         mean_tdc[i] = -1;
         sigma_tdc[i] = -1;
      }
   }
   snprintf(fname,sizeof(fname),"./ana/cut/ana%03d_tdc.pdf",ananum);
   c2->Print(fname); */

   TCanvas *c3 = new TCanvas("c3","c3");
   c3->Divide(2,2);
   for(int i=8; i<12; i++){
      c3->cd(i-7);
      gPad->SetLogy();
      ht0[i]->SetMinimum(0.5);
      ht0[i]->SetMaximum(10000);
      ht0[i]->Draw();
      if(frag[i]){
         ht2[i]->SetLineColor(3);
         ht2[i]->Draw("same");
         ht1[i]->SetLineColor(6);
         ht1[i]->Draw("same");
      }
   }
   snprintf(fname,sizeof(fname),"./ana/cut/anancasOn_tdc.pdf");
   c3->Print(fname);

   cout << fname << endl;
   cout << h1[0]->GetEntries() << endl;

   //save
   string ofname;
   ofname = "./ana/cut/chi2.txt";
   ofstream of1(ofname.c_str(),std::ios::app);
   if(of1){
      of1 << ananum << " ";
      for(int i=0; i<16; i++){
         of1 << chi2[i] << " ";
      }
      of1 << endl;
      cout << "chi2 saved" << endl;
   }
   of1.close();
   ofname = "./ana/cut/ndof.txt";
   ofstream of2(ofname.c_str(),std::ios::app);
   if(of2){
      of2 << ananum << " ";
      for(int i=0; i<16; i++){
         of2 << ndof[i] << " ";
      }
      of2 << endl;
      cout << "ndof saved" << endl;
   }
   of2.close();
   ofname = "./ana/cut/tdc_mean.txt";
   ofstream of3(ofname.c_str(),std::ios::app);
   if(of3){
      of3 << ananum << " ";
      for(int i=0; i<16; i++){
         of3 << mean_tdc[i] << " ";
      }
      of3 << endl;
      cout << "tdc mean saved" << endl;
   }
   of3.close();
   ofname = "./ana/cut/tdc_sigma.txt";
   ofstream of4(ofname.c_str(),std::ios::app);
   if(of4){
      of4 << ananum << " ";
      for(int i=0; i<16; i++){
         of4 << sigma_tdc[i] << " ";
      }
      of4 << endl;
      cout << "tdc sigma saved" << endl;
   }
   of4.close();

}

void aiueo::Search_tdc_garbage(){
   //cout << "hoge0" << endl;

   TH1D *hstg0_qdc[16], *hstg1_qdc[16], *hstg0_tdc[16], *hstg1_tdc[16], *htmax[16], *htmin[16], *hdelta[16];
   bool frag[16];
   string hname;
   double fitmin[16] = {200, 500, 1000, 600, 1000, 1000, 600, 600, 500, 400, 600, 600, 300, 700, 1000, 0};
   double fitmax[16];
   int bin = 128;
   int bint = 256;
   for(int i=0; i<16; i++){
      hname = "hstg0" + to_string(i);
      hstg0_qdc[i] = new TH1D(hname.c_str(),hname.c_str(),bin,0,4095);
      hname = "hstg1" + to_string(i);
      hstg1_qdc[i] = new TH1D(hname.c_str(),hname.c_str(),bin,0,4095);
      hname = "hstgt0" + to_string(i);
      hstg0_tdc[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);
      hname = "hstgt1" + to_string(i);
      hstg1_tdc[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);
      hname = "hdelta" + to_string(i);
      hdelta[i] = new TH1D(hname.c_str(),hname.c_str(),bint,0,4095);

      fitmax[i] = 3500;
      frag[i] = false;
      tmax[i] = mean_tdc[i] + 5*sigma_tdc[i];
      tmin[i] = mean_tdc[i] - 5*sigma_tdc[i];
      cout << tmax[i] << " " << tmin[i] << endl;
   }

   //cout << "hoge1" << endl;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //cout << "hoge2" << endl;
      for(int i=0; i<16; i++){
         hstg0_qdc[i]->Fill(adc[i]);
         hstg0_tdc[i]->Fill(tdc[i]);
         //cout << "hoge3" << endl;
         if(Hit2(i)&&Hit(4)){
            frag[i] = true;
            hstg1_qdc[i]->Fill(adc[i]);
            hstg1_tdc[i]->Fill(tdc[i]);
         }
      }
   }
   //cout << "hoge" << endl;


   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(4,4);
   for(int i=0; i<16; i++){
      c1->cd(i+1);
      gPad->SetLogy();
      hstg0_qdc[i]->SetMinimum(0.7);
      hstg0_qdc[i]->SetMaximum(100000);
      hstg0_qdc[i]->Draw();
      if(frag[i] && mean_tdc[i]!=-1 && sigma_tdc[i]!=-1){
         hstg1_qdc[i]->SetLineColor(6);
         hstg1_qdc[i]->Draw("same");
      } 
   }
   char fname[256];
   snprintf(fname,sizeof(fname),"./ana/cut/ana%03d_stg_q.pdf",ananum);
   c1->Print(fname);

   TCanvas *c2 = new TCanvas("c2","c2");
   c2->Divide(4,4);
   TF1 *func2;
   for(int i=0; i<16; i++){
      c2->cd(i+1);
      gPad->SetLogy();
      hstg0_tdc[i]->SetMinimum(0.1);
      hstg0_tdc[i]->SetMaximum(10000);
      hstg0_tdc[i]->Draw();
      if(frag[i] && mean_tdc[i]!=-1 && sigma_tdc[i]!=-1){
         hstg1_tdc[i]->SetLineColor(6);
         hstg1_tdc[i]->Draw("same");
      }
   }
   snprintf(fname,sizeof(fname),"./ana/cut/ana%03d_stg_t.pdf",ananum);
   c2->Print(fname); 

   cout << hstg1_qdc[0]->GetEntries() << endl;

}

bool aiueo::Hit(int id){
   bool frag = false;
   if(id == 0){frag = true;}
   if(id == 1){
      if(tdc[12]!=0){
         frag = true;
      }
   }
   if(id == 2){
      for(int i=0; i<16; i++){
         if(detch[i]!=-1){
            if(i==0 && tdc[i]!=0) frag = true;
            else if(i>0 && frag && tdc[i]!=0) frag = true;
            else frag = false;
         }
      }
   }
   if(id == 3){
      for(int i=0; i<16; i++){
         if(detch[i]!=-1){
            if(i==0 && tdc[i]!=0) frag = true;
            else if(i>0 && frag && tdc[i]!=0) frag = true;
            else frag = false;
         }
      }
      if(frag && adc[1]<1200) frag = true;
      else frag = false;
   }
   if(id == 4){
      //frag = true;
      //if(tdc[12]!=0) frag = true;
      /* for(int i=0; i<11; i++){
         if(detch[i]!=-1){
            if(i==0 && tdc[i]!=0) frag = true;
            else if(i>0 && frag && tdc[i]!=0) frag = true;
            else frag = false;
         }
      }  */
      for(int i=0; i<16; i++){
         if(detch[i]!=-1){
            if(i==0 && tdc[i]!=0) frag = true;
            else if(i>0 && frag && tdc[i]!=0) frag = true;
            else frag = false;
         }
      }  
      if(frag && adc[13]<1500 && adc[14]<2000) frag = true;
      else frag = false;
      if(frag && adc[6]<1500 && adc[7]<1400) frag = true;
      else frag = false;
      if(frag && adc[1]<1200) frag = true;
      else frag = false;
      /* if(frag && adc[8]<1200 && adc[9]<1200) frag = true;
      else frag = false;
      if(frag && adc[10]<1350 && adc[11]<1350) frag = true;
      else frag = false; */
   }
   if(id == 5){
      for(int i=0; i<16; i++){
         if(detch[i]!=-1){
            if(i==0 && tdc[i]!=0) frag = true;
            else if(i>0 && frag && tdc[i]!=0) frag = true;
            else frag = false;
         }
      }
   }

   return frag;
}

bool aiueo::Hit2(int id){
   bool frag = false ;

   if(tdc[id]>tmax[id] || tdc[id]<tmin[id]) frag = true;
   else frag = false;
   
   return frag;
}

void aiueo::Set_exist_det(int id){
   //return;
   char dummy[32];
   string filename;
   filename = "./ana/det" + to_string(id) + ".txt";
   ifstream file(filename.c_str());
   for (int i = 0; i < 16; i++)
   {
      file >> dummy >>detch[i];
      cout << dummy << " " << detch[i] << endl;
   } 
}

void aiueo::Check(){
   TGraph *gr1;
}

void aiueo::MipTest(){
   TH1D *h0[16], *h1[16];
   char hname[16];
   int bin = 128;
   int hmin = 0;
   int hmax = 4095;
   bool detect[16];
   for(int i=0; i<16; i++){
      snprintf(hname,sizeof(hname),"hm0%02d",i);
      h0[i] = new TH1D(hname,hname,bin,hmin,hmax);
      snprintf(hname,sizeof(hname),"hm1%02d",i);
      h1[i] = new TH1D(hname,hname,bin,hmin,hmax);
      detect[i] = false;
   }
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for(int i=0; i<16; i++){
         h0[i]->Fill(adc[i]);
         if(tdc[i]!=0){ 
            h1[i]->Fill(adc[i]);
            detect[i] = true;
         }
      }
   }
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(4,4);
   TF1 *func[16];
   int mean = 1000;
   int range = 600;
   for(int i=0; i<16; i++){
      func[i] = new TF1("func","gaus");
      func[i]->SetParameter(1,mean);
      c1->cd(i+1);
      gPad->SetLogy();
      h0[i]->Draw();
      if(detect[i]){
         h1[i]->SetLineColor(3);
         h1[i]->Draw("same");
         h1[i]->Fit(func[i],"","",mean-range,mean+range);
         mip_mean[i] = func[i]->GetParameter(1);
         mip_sigma[i] = func[i]->GetParameter(2);
      } else {
         mip_mean[i] = 1000;
         mip_sigma[i] = 100;
      }
   }
   c1->Delete();
}

void aiueo::MipCheck(){
   TH1D *h0[16], *h1[16], *h2[16];
   char hname[16];
   int bin = 128;
   int hmin = 0;
   int hmax = 4095;
   bool detect[16];
   double max1[16], max2[16];
   for(int i=0; i<16; i++){
      snprintf(hname,sizeof(hname),"hmc0%02d",i);
      h0[i] = new TH1D(hname,hname,bin,hmin,hmax);
      snprintf(hname,sizeof(hname),"hmc1%02d",i);
      h1[i] = new TH1D(hname,hname,bin,hmin,hmax);
      snprintf(hname,sizeof(hname),"hmc2%02d",i);
      h2[i] = new TH1D(hname,hname,bin,hmin,hmax);
      detect[i] = false;

      max1[i] = mip_mean[i] + 3*mip_sigma[i];
      max2[i] = mip_mean[i] + 4*mip_sigma[i];
   }
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for(int i=0; i<16; i++){
         h0[i]->Fill(adc[i]);
         if(adc[i]<max1[i]) h1[i]->Fill(adc[i]);
         if(adc[i]<max2[i]) h2[i]->Fill(adc[i]);
      }
   }
   TCanvas *cmc = new TCanvas("cmc","cmc");
   cmc->Divide(4,4);
   TF1 *func[16];
   int mean = 1000;
   int range = 600;
   for(int i=0; i<16; i++){
      func[i] = new TF1("func","gaus");
      func[i]->SetParameter(1,mean);
      cmc->cd(i+1);
      gPad->SetLogy();
      h0[i]->Draw();
      h1[i]->SetLineColor(6);
      h1[i]->Draw("same");
      h2[i]->SetLineColor(3);
      h2[i]->Draw("same");
   }
}