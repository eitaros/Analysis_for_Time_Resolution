#define ana_eito_cxx
#include "ana_eito.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define noofdet 15

bool ana_eito::hit(int t)
{
   bool flag = false;
   //ignore
   if (t == 0)
   {
      flag = true;
   }
   
   //hit
   if (t == 1){
      for (int i =0 ; i<noofdet; i++)
      {
         if ( i==0 ) 
         { 
            if (leading[detno[i]][0]!=-1) flag = true; 
         } else if (i<13) {        
            if ( flag && leading[detno[i]][0] !=-1 ) flag = true; 
            else flag = false;
         } else {}
      }
   }

   //good adc
   if (t == 2){
      for (int i = 3; i<noofdet; i++)
      {
         if ( i==3 ) { 
            if (charge[i]>0.6&&charge[i]<4.0) flag = true; 
         } else if (i<7) {        
            if ( flag && charge[i]>0.6&charge[i]<4.0) flag = true; 
            else flag = false;
         } else if (i<13) {
            if ( flag && charge[i]>0.7&charge[i]<3.0) flag = true; 
            else flag = false;
         } else {}
      }
   }

   return flag;
}

void ana_eito::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ana_eito.C
//      root> ana_eito t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

void ana_eito::ped_mip()
{
   //initiarise
   for (int i=0; i<noofdet; i++)
   {
      pedestal[i]=0.;
      ped_sigma[i]=0.;
      mip[i]=0.;
      t_offset[i]=0.;
      t_width[i]=0.;
      if(h_adcmip[i]) delete h_adcmip[i];
      if(h_adcped[i]) delete h_adcped[i];
      if(h_rawtdc[i]) delete h_rawtdc[i];
   }

   

   //make hist
   char name[32];
   for (int i = 0; i < noofdet; i++)
   {
      sprintf(name,"h_adcmip_%d",detno[i]);
      h_adcmip[i] = new TH1F(name,name,4095,0.,4095.); 
      sprintf(name,"h_adcped_%d",detno[i]);
      h_adcped[i] = new TH1F(name,name,500,0.,500.);
      sprintf(name,"h_rawtdc_%d",detno[i]);
      h_rawtdc[i] = new TH1F(name,name,10000,800000.,1000000.);
   }

   if (fChain == 0) return;

   //loop
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //Convert();
      for (int i=0; i<noofdet; i++) 
      {
        if (leading[detno[i]][0]!=-1) h_adcmip[i]->Fill(adc[detno[i]]); 
        if (leading[detno[i]][0]!=-1) h_rawtdc[i]->Fill(leading[detno[i]][0]); 
        if (leading[detno[i]][0]==-1) h_adcped[i]->Fill(adc[detno[i]]);
      }
   }
   
   //ped
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(2,2);
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->Divide(2,2);
   TCanvas *c3 = new TCanvas("c3","c3");
   if (noofdet == 11) c3->Divide(2,2);
   else c3->Divide(2,3);
   TCanvas *c10 = new TCanvas("c10","c10");
   c10->Divide(1,2);
   if (noofdet < 15) delete c10;

   TFitResultPtr results1;
   for (int i = 0; i < noofdet; i++)
   {
      if (i < 3)
      {
         c1->cd(i+1);
         results1 = h_adcped[i]->Fit("gaus","S","",0,500);
         h_adcped[i]->SetXTitle("ADCchannel(ch)");
         h_adcped[i]->SetYTitle("Counts");
         pedestal[i] = results1->Parameter(1);
         ped_sigma[i] = results1->Parameter(2);
      } else if (i>2&&i<7)
      {
         c2->cd(i-2);
         results1 = h_adcped[i]->Fit("gaus","S","",0,500);
         h_adcped[i]->SetXTitle("ADCchannel(ch)");
         h_adcped[i]->SetYTitle("Counts");
         pedestal[i] = results1->Parameter(1);
         ped_sigma[i] = results1->Parameter(2);
      } else if (i>6&&i<13)
      {
         c3->cd(i-6);
         results1 = h_adcped[i]->Fit("gaus","S","",0,500);
         h_adcped[i]->SetXTitle("ADCchannel(ch)");
         h_adcped[i]->SetYTitle("Counts");
         pedestal[i] = results1->Parameter(1);
         ped_sigma[i] = results1->Parameter(2);
      } else
      {
         c10->cd(i-12);
         results1 = h_adcped[i]->Fit("gaus","S","",300,400);
         h_adcped[i]->SetXTitle("ADCchannel(ch)");
         h_adcped[i]->SetYTitle("Counts");
         pedestal[i] = results1->Parameter(1);
         ped_sigma[i] = results1->Parameter(2);
      }
   }
   char ped_fig[64];
   sprintf(ped_fig,"pdf/run%4.4d_ped_trig.pdf",runno);
   c1->Print(ped_fig);
   sprintf(ped_fig,"pdf/run%4.4d_ped_pa.pdf",runno);
   c2->Print(ped_fig);
   sprintf(ped_fig,"pdf/run%4.4d_ped_cnc.pdf",runno);
   c3->Print(ped_fig);
   if (noofdet == 15) {
      sprintf(ped_fig,"pdf/run%4.4d_ped_ref.pdf",runno);
      c10->Print(ped_fig);
   }

   //mip
   TCanvas *c4 = new TCanvas("c4","c4");
   c4->Divide(2,2);
   TCanvas *c5 = new TCanvas("c5","c5");
   c5->Divide(2,2);
   TCanvas *c6 = new TCanvas("c6","c6");
   if (noofdet == 11) c6->Divide(2,2);
   else c6->Divide(2,3);
   TCanvas *c11 = new TCanvas("c11","c11");
   c11->Divide(1,2);
   if (noofdet < 15) delete c11;
   
   TFitResultPtr results2;
   for (int i = 0; i < noofdet; i++)
   {
      if (i < 3)
      {
         c4->cd(i+1);
         results2 = h_adcmip[i]->Fit("landau","S","",300,4095);
         h_adcmip[i]->SetXTitle("ADCchannel(ch)");
         h_adcmip[i]->SetYTitle("Counts");
         mip[i] = results2->Parameter(1);
      } else if (i>2&&i<7)
      {
         c5->cd(i-2);
         results2 = h_adcmip[i]->Fit("landau","S","",300,4095);
         h_adcmip[i]->SetXTitle("ADCchannel(ch)");
         h_adcmip[i]->SetYTitle("Counts");
         mip[i] = results2->Parameter(1);
      } else if (i>6&&i<13)
      {
         c6->cd(i-6);
         results2 = h_adcmip[i]->Fit("landau","S","",300,4095);
         h_adcmip[i]->SetXTitle("ADCchannel(ch)");
         h_adcmip[i]->SetYTitle("Counts");
         mip[i] = results2->Parameter(1);
      } else
      {
         c11->cd(i-12);
         results1 = h_adcmip[i]->Fit("landau","S","",1000,2500);
         h_adcmip[i]->SetXTitle("ADCchannel(ch)");
         h_adcmip[i]->SetYTitle("Counts");
         mip[i] = results2->Parameter(1);
      }
   }
   char mip_fig[64];
   sprintf(mip_fig,"pdf/run%4.4d_mip_trig.pdf",runno);
   c4->Print(mip_fig);
   sprintf(mip_fig,"pdf/run%4.4d_mip_pa.pdf",runno);
   c5->Print(mip_fig);
   sprintf(mip_fig,"pdf/run%4.4d_mip_cnc.pdf",runno);
   c6->Print(mip_fig);
    if (noofdet == 15) {
      sprintf(ped_fig,"pdf/run%4.4d_mip_ref.pdf",runno);
      c11->Print(ped_fig);
   }

   //tdc offset
   TCanvas *c7 = new TCanvas("c7","c7");
   c7->Divide(2,2);   
   TCanvas *c8 = new TCanvas("c8","c8");
   c8->Divide(2,2); 
   TCanvas *c9 = new TCanvas("c9","c9");
   if (noofdet == 11) c9->Divide(2,2);
   else c9->Divide(2,3);
   TCanvas *c12 = new TCanvas("c12","c12");
   c12->Divide(1,2);
   if (noofdet < 15) delete c11;
   TFitResultPtr results3;
   for (int i=0; i<noofdet; i++)
   {
      if (i < 3)
      {
         c7->cd(i+1);
         results3  = h_rawtdc[i]  ->Fit("gaus","S");
         h_rawtdc[i]->SetXTitle("TDCchannel(ch)");
         h_rawtdc[i]->SetYTitle("Counts");
         t_offset[i]   = results3->Parameter(1);
         t_width[i] = results3->Parameter(2);
      } else if (i>2&&i<7)
      {
         c8->cd(i-2);
         results3  = h_rawtdc[i]  ->Fit("gaus","S");
         h_rawtdc[i]->SetXTitle("TDCchannel(ch)");
         h_rawtdc[i]->SetYTitle("Counts");
         t_offset[i]   = results3->Parameter(1);
         t_width[i] = results3->Parameter(2);
      } else if (i>6&&i<13)
      {
         c9->cd(i-6);
         results3  = h_rawtdc[i]  ->Fit("gaus","S");
         h_rawtdc[i]->SetXTitle("TDCchannel(ch)");
         h_rawtdc[i]->SetYTitle("Counts");
         t_offset[i]   = results3->Parameter(1);
         t_width[i] = results3->Parameter(2);
      } else
      {
         c12->cd(i-12);
         results3 = h_rawtdc[i]->Fit("gaus","S");
         h_rawtdc[i]->SetXTitle("TDCchannel(ch)");
         h_rawtdc[i]->SetYTitle("Counts");
         t_offset[i]   = results3->Parameter(1);
         t_width[i] = results3->Parameter(2);
      }
   }
   char tdc_fig[64]; 
   sprintf(tdc_fig,"pdf/run%4.4d_tdc_trig.pdf",runno);
   c7->Print(tdc_fig);
   sprintf(tdc_fig,"pdf/run%4.4d_tdc_pa.pdf",runno);
   c8->Print(tdc_fig);
   sprintf(tdc_fig,"pdf/run%4.4d_tdc_cnc.pdf",runno);
   c9->Print(tdc_fig);
    if (noofdet == 15) {
      sprintf(ped_fig,"pdf/run%4.4d_tdc_ref.pdf",runno);
      c12->Print(ped_fig);
   }

   //write file
   char ofile1[64]; 
   sprintf(ofile1,"run%4.4d.ped_mip",runno);
   ofstream output_file1(ofile1);
   output_file1 << fixed << setprecision(2);
   for (int i=0; i<noofdet; i++) 
   {
     output_file1 << pedestal[i]  << " " << ped_sigma[i] << " "<< mip[i] << endl;
   }
   output_file1.close();

   char ofile2[64]; 
   sprintf(ofile2,"run%4.4d.t_offset",runno);
   ofstream output_file2(ofile2);
   output_file2 << fixed << setprecision(2);
   for (int i=0; i<noofdet; i++) 
   {
     output_file2 << t_offset[i] << " " << t_width[i] << endl;
   }
   output_file2.close(); 
}

void ana_eito::phtw(int cor1,int cor2,int ref1,int ref2)
{
   //reset
   //if (c1) delete c1;
   if (func) delete func;
   for (int i = 0; i < 2; i++)
   {
      if (dt[i]) delete dt[i];
   }
   
   //return;
   Bool_t flag = kFALSE;
   for (int i =0 ; i<noofdet; i++)
   {
      if ( i==0 ) 
      { 
         if ((mip[i]-pedestal[i])!=0) flag = kTRUE; 
      }
      else 
      {        
         if ( flag && (mip[i]-pedestal[i])!=0 ) {flag = kTRUE;} 
         else {printf("%d %f %f\n",i,mip[i],pedestal[i]); flag = kFALSE;}
      }
   }
   if(!flag)
   {
      cout << "do ped_mip() before phtw()" << endl;
      return;
   }

   if (fChain == 0) return;
   
   //set id
   // input raw adc,tdc cahnnnel to cor1,2/ref1,2. 
   // cor_r,l/ref_r,l return list's id (without ch3)
   std::vector<int> vec;
   cor_r = cor1-1;
   cor_l = cor2-1;
   ref_r = ref1-1;
   ref_l = ref2-1;
   vec.push_back(ref_r); vec.push_back(ref_l); 
   vec.push_back(cor_r); vec.push_back(cor_l); 

   //parameter intiarise
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> slew_no[i];
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
         slew_no[i] = 0;
      }
   }
   slew_pass_file.close();

   //swel_no check
   if (slew_no[cor_r]!=slew_no[cor_l]||slew_no[ref_r]!=slew_no[ref_l])
   {
      cout << "input id is bat. retry." << endl;
      cout << "cor1 : " << cor1 << " slew_no : " <<slew_no[cor_r] << endl;
      cout << "cor2 : " << cor2 << " slew_no : " <<slew_no[cor_l] << endl;
      cout << "ref1 : " << ref1 << " slew_no : " <<slew_no[ref_r] << endl;
      cout << "ref2 : " << ref2 << " slew_no : " <<slew_no[ref_l] << endl;
      return;
   }
   

   //max min define
   double hist_min[2];
   double hist_max[2];
   for (size_t i = 0; i < 2; i++)
   {
      hist_min[i] = 0.;
      hist_max[i] = 5.0;
   }
   //hist_max[0] = (4095. - pedestal[cor_r])/(mip[cor_r] - pedestal[cor_r]);
   //hist_max[1] = (4095. - pedestal[cor_l])/(mip[cor_l] - pedestal[cor_l]);
   double fit_min[2];
   double fit_max[2];
   for (size_t i = 0; i < 2; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.5;
      //fit_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      fit_max[i] = 4.0;
   }

   //make hist
   char histname[64];
   sprintf(histname,"dt:adc[%d]",detno[cor_r]);
   dt[0] = new TH2F(histname,histname,64,hist_min[0],hist_max[0],10000,-50000,50000);
   sprintf(histname,"dt:adc[%d]",detno[cor_l]);
   dt[1] = new TH2F(histname,histname,64,hist_min[1],hist_max[1],1000,-50000,50000);
   
   //define function
   func = new TF1("func","[0] + [1]/sqrt(x) + [2]*sqrt(x)",0.,4.);

   if (fChain == 0) return;

   //printf("########################\n");

   //slew
   double t_ref_r, t_ref_l, t_cor_r, t_cor_l, t_cor, t_ref, delta_time;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit2(vec))
      {
         t_cor_r = time[cor_r] - (slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r]));
         t_cor_l = time[cor_l] - (slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]));
         t_ref_r = time[ref_r] + (slew_par0[ref_r] + slew_par1[ref_r]/sqrt(charge[ref_r]) + slew_par2[ref_r]*sqrt(charge[ref_r]));
         t_ref_l = time[ref_l] + (slew_par0[ref_l] + slew_par1[ref_l]/sqrt(charge[ref_l]) + slew_par2[ref_l]*sqrt(charge[ref_l]));
         t_cor = (t_cor_r + t_cor_l)/2.0;
         t_ref = (t_ref_r + t_ref_l)/2.0;
         delta_time = t_cor - t_ref;
         if (slew_no[ref_r]==0||slew_no[ref_l]==0)
         {
            if (charge[ref_r]>1.2&&charge[ref_l]>1.2) {
               dt[0]->Fill(charge[cor_r],delta_time);
               dt[1]->Fill(charge[cor_l],delta_time);
            }
         } else
         {
            dt[0]->Fill(charge[cor_r],delta_time);
            dt[1]->Fill(charge[cor_l],delta_time);
         }
      }
   }

   //fit
   TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,1);
   char prof_name[32];
   TFitResultPtr result;
   TProfile *prof[2];
   //right
   c1->cd(1);
   sprintf(prof_name,"slew_right");
   dt[0]->ProfileX(prof_name);
   prof[0] = (TProfile*)gROOT->FindObject(prof_name);
   result = prof[0]->Fit("func","S","",fit_min[0],fit_max[0]);
   slew_par0[cor_r] += result->Parameter(0);
   slew_par1[cor_r] += result->Parameter(1);
   slew_par2[cor_r] += result->Parameter(2);
   slew_no[cor_r] += 1;
   //left
   c1->cd(2);
   sprintf(prof_name,"slew_left");
   dt[1]->ProfileX(prof_name);
   prof[1] = (TProfile*)gROOT->FindObject(prof_name);
   result = prof[1]->Fit("func","S","",fit_min[1],fit_max[1]);
   slew_par0[cor_l] += result->Parameter(0);
   slew_par1[cor_l] += result->Parameter(1);
   slew_par2[cor_l] += result->Parameter(2);
   slew_no[cor_l] += 1;
   
   //save
   //pdf
   char pdfname[64];
   if (slew_no[cor_r]==0||slew_no[cor_l]==0)
   {
      sprintf(pdfname,"pdf/run%4.4d_1stslew_%d%dvs%d%d.pdf",runno,cor1,cor2,ref1,ref2);
   }else
   {
      sprintf(pdfname,"pdf/run%4.4d_slew_%d%dvs%d%d.pdf",runno,cor1,cor2,ref1,ref2);
   }
   c1->Print(pdfname);
   //par
   char fname[64];
   sprintf(fname,"run%4.4d.slew_pass",runno);
   ofstream outputfile(fname);
   outputfile << fixed << setprecision(2);
   for (int i = 0; i < noofdet; i++)
   {
      outputfile << slew_par0[i] << " " << slew_par1[i] << " " << slew_par2[i] << " " << slew_no[i] << endl;
   }
   outputfile.close();
}

void ana_eito::dt_fit(int test)
{
   //return;
   /*
      test = 0 cdh
           = 1 cnc1
           = 2 cnc2
           = 3 ref
   */

  //define fit channel name
   int test_r,test_l,ref1_r,ref1_l,ref2_r,ref2_l;
   std::vector<int> vec;
   ref1_r = 3;
   ref1_l = 4;
   ref2_r = 5;
   ref2_l = 6;
   if (test == 0) {test_r = 7; test_l = 8;}
   if (test == 1) {test_r = 9; test_l = 10;}
   if (test == 2) {test_r = 11; test_l = 12;}
   if (test == 3) {test_r = 13; test_l = 14;}
   vec.push_back(ref1_r); vec.push_back(ref1_l);
   vec.push_back(ref2_r); vec.push_back(ref2_l);
   vec.push_back(test_r); vec.push_back(test_l);
   char test_name[8];
   if (test == 0) sprintf(test_name,"cdh");
   if (test == 1) sprintf(test_name,"cnc1");
   if (test == 2) sprintf(test_name,"cnc2"); 
   if (test == 3) sprintf(test_name,"ref"); 

   //check ped_mip
   /* bool x_conv_frag = (pedestal[0]-mip[0])==0.||(pedestal[1]-mip[1])==0.
                        ||(pedestal[2]-mip[2])==0.||(pedestal[3]-mip[3])==0.;
   if(x_conv_frag)
   {
      cout << "do ped_mip() before phtw()" << endl;
      return;
   } */

   //slew parameter intiarise
   int dummy;
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> dummy;
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
      }
   }
   slew_pass_file.close();

   //make hist
   char histname[32];
   for (size_t i = 0; i < 3; i++)
   {
      if (i == 0) sprintf(histname,"dt_pa1_vs_pa2");
      if (i == 1) sprintf(histname,"dt_pa1_vs_%s",test_name);
      if (i == 2) sprintf(histname,"dt_pa2_vs_%s",test_name);
      h_dt[i] = new TH1F(histname,histname,20000,-30000,30000);
   }

   if (fChain == 0) return;

   //loop fill
   double raw_ref1, raw_ref2, raw_test, slew_pa1, slew_pa2, slew_test;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (1)
         //hit2(vec))
         //(time[test_r]!=-1 && time[test_l]!=-1 && time[ref1_r]!=-1 && time[ref1_l]!=-1 && time[ref2_r]!=-1 && time[ref2_l]!=-1)
      {
         raw_ref1 = (time[ref1_r] + time[ref1_l])/2;
         raw_ref2 = (time[ref2_r] + time[ref2_l])/2;
         raw_test = (time[test_r] + time[test_l])/2;
         slew_pa1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]) 
                        + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2;
         slew_pa2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]) 
                        + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2;
         slew_test = (slew_par0[test_r] + slew_par1[test_r]/sqrt(charge[test_r]) + slew_par2[test_r]*sqrt(charge[test_r]) 
                        + slew_par0[test_l] + slew_par1[test_l]/sqrt(charge[test_l]) + slew_par2[test_l]*sqrt(charge[test_l]))/2;
         if (hit(0))
         {
            h_dt[0]->Fill(raw_ref1 - raw_ref2 - slew_pa1 + slew_pa2);
            h_dt[1]->Fill(raw_test - raw_ref1 - slew_test + slew_pa1);
            h_dt[2]->Fill(raw_test - raw_ref2 - slew_test + slew_pa2);
         }
      }
   }
   
   //draw fit
   TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,2);
   TFitResultPtr result;
   double mean[3];
   double sigma[3];
   for (int i = 0; i < 3; i++)
   {
      sigma[i] = 0.;
   }
   for (int i = 0; i < 3; i++)
   {
      c1->cd(i+1);
      mean[i] = h_dt[i]->GetMean(1);
      h_dt[i]->GetXaxis()->SetRangeUser(mean[i] - 2500.,mean[i] + 2500.);
      h_dt[i]->Draw();
      result = h_dt[i]->Fit("gaus","S","",mean[i] - 2500.,mean[i] + 2500.);
      sigma[i] = result->Parameter(2);
   }
   
   double pa1_res, pa2_res, test_res;
   pa1_res = 0.;
   pa2_res = 0.;
   test_res = 0.;
   pa1_res = sqrt((sigma[0]*sigma[0] + sigma[1]*sigma[1] - sigma[2]*sigma[2])/2);
   pa2_res = sqrt((sigma[0]*sigma[0] + sigma[2]*sigma[2] - sigma[1]*sigma[1])/2);
   test_res = sqrt((sigma[1]*sigma[1] + sigma[2]*sigma[2] - sigma[0]*sigma[0])/2);

   char comment[64];
   c1->cd(4);
   auto frame = c1->DrawFrame(0,0,1,1);
   frame->GetXaxis()->SetAxisColor(kWhite);
   frame->GetYaxis()->SetAxisColor(0);
   TLatex latex;
   latex.SetTextSize(0.08);
   latex.SetTextAlign(13);  //align at top
   latex.DrawLatex(0.05,.8,"############################");
   latex.DrawLatex(0.05,.7,"          time resolution           ");
   sprintf(comment,"       pa1    = %3.3f ps ",pa1_res);
   latex.DrawLatex(0.05,.6,comment);
   sprintf(comment,"       pa2    = %3.3f ps ",pa2_res);
   latex.DrawLatex(0.05,.5,comment);
   sprintf(comment,"       %s     = %3.3f ps ",test_name,test_res);
   latex.DrawLatex(0.05,.4,comment);
   latex.DrawLatex(0.05,.3,"############################");

   char time_resolution[64];
   sprintf(time_resolution,"run%4.4d_resolution%d.pdf",runno,test);
   c1->Print(time_resolution);
   
   //print result
   cout << endl;
   cout << "####################################" << endl;
   cout << "          time resolution           " << endl;
   cout << "       pa1    = " << pa1_res << " ps" << endl;
   cout << "       pa2    = " << pa2_res << " ps" << endl;
   cout << "       " << test_name << "   = " << test_res << " ps" << endl;
   //cout << "       test   = " << test_res << " ps" << endl;
   cout << "####################################" << endl; 
}

void ana_eito::changed_data()
{
   //return;
   /* //check ped_mip
   bool x_conv_frag = (pedestal[0]-mip[0])==0.||(pedestal[1]-mip[1])==0.
                        ||(pedestal[2]-mip[2])==0.||(pedestal[3]-mip[3])==0.;
   if(x_conv_frag)
   {
      cout << "do ped_mip() before changed_data()" << endl;
      return;
   } */

   //define adc hist max
   double hist_min[noofdet];     //adc
   double hist_max[noofdet];
   double tdc_min[noofdet];
   double tdc_max[noofdet];
   for (size_t i = 0; i < noofdet; i++)
   {
      hist_min[i] = 0.;
      hist_max[i] = 5.0;
      //(4095. - pedestal[i])/(mip[i] - pedestal[i]);
      tdc_min[i] = -10*t_width[i]*timeConv[i];
      tdc_max[i] = 10*t_width[i]*timeConv[i];
   }
   double fit_min[noofdet];
   double fit_max[noofdet];
   for (size_t i = 0; i < 4; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.5;
      fit_max[i] = 5.0;
      //(4095. - pedestal[i])/(mip[i] - pedestal[i]);
      //fit_max[i] = 3.0;
   }

   //make hist
   char histname[64];
   for (int i = 0; i < noofdet; i++)
   {
      sprintf(histname,"changed_adc%d",detno[i]);
      changed_adc[i] = new TH1F(histname,histname,400,hist_min[i],hist_max[i]);
      sprintf(histname,"changed_tdc%d",detno[i]);
      changed_tdc[i] = new TH1F(histname,histname,100,tdc_min[i],tdc_max[i]);
   }
   

   //fill adc
   //TProfile *prof[4];
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit(1))
      {
         for (int i = 0; i < noofdet; i++)
         {
            changed_adc[i]->Fill(charge[i]);
            changed_tdc[i]->Fill(time[i]);
         }
      }
   }
   
   //draw save hist
   /* TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,2);
   for (int i = 0; i < 3; i++)
   {
      c1->cd(i+1);
      changed_adc[i]->Draw();
      changed_adc[i]->Fit("landau","","",fit_min[i],fit_max[i]);
   }
   char pdf_name[64];
   sprintf(pdf_name,"run%4.4d_adc_trig.pdf",runno);
   c1->Print(pdf_name); */
   
   TCanvas *c2 = new TCanvas("c2","");
   c2->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c2->cd(i+1);
      changed_adc[i+3]->SetXTitle("energy(E/E_mip)");
      changed_adc[i+3]->SetYTitle("counts");
      changed_adc[i+3]->Draw();
      changed_adc[i+3]->Fit("landau","","",fit_min[i],fit_max[i]);
   }
   char pdf_name[64];
   sprintf(pdf_name,"run%4.4d_adc_pa.pdf",runno);
   c2->Print(pdf_name);

   TCanvas *c3 = new TCanvas("c3","");
   if (noofdet == 11) c3->Divide(2,2);
   if (noofdet == 13) c3->Divide(2,3);
   for (int ii = 0; ii < 6; ii++)
   {
      if (noofdet == 11&& 4<=ii) break;

      c3->cd(ii+1);
      changed_adc[ii+7]->SetXTitle("energy(E/E_mip)");
      changed_adc[ii+7]->SetYTitle("counts");
      changed_adc[ii+7]->Draw();
      changed_adc[ii+7]->Fit("landau","","",fit_min[ii],fit_max[ii]);
   }
   //char pdf_name[64];
   sprintf(pdf_name,"run%4.4d_adc_cnc.pdf",runno);
   c3->Print(pdf_name);
   
   if (noofdet == 15)
   {
      TCanvas* c4 = new TCanvas("c4","c4");
      c4->Divide(1,2);
      for (int ii = 0; ii < 2; ii++)
      {
         c4->cd(ii+1);
         changed_adc[ii+13]->SetXTitle("energy(E/E_mip)");
         changed_adc[ii+13]->SetYTitle("counts");
         changed_adc[ii+13]->Draw();
         changed_adc[ii+13]->Fit("landau","","",fit_min[ii+13],fit_max[ii+13]);
      }
      sprintf(pdf_name,"run%4.4d_adc_ref.pdf",runno);
      c4->Print(pdf_name);
   }
   
   
   
   /* TCanvas *c2 = new TCanvas("c2","");
   c2->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c2->cd(i+1);
      changed_tdc[i]->Draw();
   }
   sprintf(pdf_name,"run%4.4d_changed_data.pdf",runno);
   c2->Print(pdf_name);  */
}

void ana_eito::check()
{
   return;
   /* if (fChain == 0) return;

   //make hist
   char histname[32];
   for (int i = 0; i < noofdet; i++)
   {
      sprintf(histname,"raw_adc_%d",i);
      raw_adc[i] = new TH1F(histname,histname,4095,0.,4095.);
      sprintf(histname,"raw_tdc_%d",i);
      raw_tdc[i] = new TH1F(histname,histname,4095,0.,4095.);
      sprintf(histname,"hit_adc_%d",i);
      hit_adc[i] = new TH1F(histname,histname,4095,0.,4095.);
   }

   //loop
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //Convert();
      for (int i=0; i<noofdet; i++) 
      {
         raw_adc[i]->Fill(adc[detno[i]])
         if (tdc[noofdet]!=4095) 
         {
            hit_adc[i]->Fill(adc[detno[i]]); 
            raw_tdc[i]->Fill(tdc[detno[i]]);
         }
      }
   }

   //Draw adc
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(3,3);
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->Divide(2,2);
   for (int i = 0; i < noofdet; i++)
   {
      if (i < 7)
      {
         c1->cd(i+1);
         gPad->SetLogy(1);
         raw_adc[i]->GetXaxis()->SetTitle("ADC_channel");
         raw_adc[i]->GetYaxis()->SetTitle("Counts");
         raw_adc[i]->Draw();
         
         hit_adc[i]->GetXaxis()->SetTitle("ADC_channel");
         hit_adc[i]->GetYaxis()->SetTitle("Counts");
         hit_adc[i]->SetLineColor(2);
         hit_adc[i]->Draw("same");
      }else if (i=>7&&i<noofdet)
      {
         c2->cd(i-6);
         gPad->SetLogy(1);
         raw_adc[i]->Draw();
         hit_adc[i]->SetLineColor(2);
         hit_adc[i]->Draw("same");
      }
   }
   char name[64];
   sprintf(name,"run%4.4d_rawadc_trig_pa.pdf",runno);
   c1->Print(name);
   sprintf(name,"run%4.4d_rawadc_cnc.pdf",runno);
   c1->Print(name);

   //Draw tdc
   TCanvas *c1 = new TCanvas("c3","c3");
   c3->Divide(3,3);
   TCanvas *c2 = new TCanvas("c4","c4");
   c4->Divide(2,2);
   for (int i = 0; i < noofdet; i++)
   {
      if (i < 7)
      {
         c3->cd(i+1);
         raw_tdc[i]->GetXaxis()->SetTitle("TDC_channel");
         raw_tdc[i]->GetYaxis()->SetTitle("Counts");
         raw_tdc[i]->Draw();
      }else if (i=>7&&i<noofdet)
      {
         c4->cd(i-6);
         gPad->GetXaxis()->SetTitle("TDC_channel");
         gPad->GetYaxis()->SetTitle("Counts");
         raw_tdc[i]->Draw();
      }
   }
   char name[64];
   sprintf(name,"run%4.4d_rawtdc_trig_pa.pdf",runno);
   c3->Print(name);
   sprintf(name,"run%4.4d_rawtdc_cnc.pdf",runno);
   c4->Print(name); */
}

void ana_eito::check_slew() 
{
   if (fChain == 0) return;

   //parameter intiarise
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> slew_no[i];
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
         slew_no[i] = 0;
      }
   }
   slew_pass_file.close();

   //make hist
   char hname[64];
   double hbin = 10000.;
   double hmin = -200000.;
   double hmax = 200000.;
   for (int i = 0; i < 4; i++)
   {
      sprintf(hname,"before_pa1,pa2_vs_adc[%d]",detno[i+3]);
      raw_slew_pa[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa1,pa2_vs_adc[%d]",detno[i+3]);
      check_slew_pa[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
   }
   for (int i = 0; i < 2; i++)
   {
      sprintf(hname,"before_pa1,cdh_vs_adc[%d]",detno[i+7]);
      raw_slew_cdh[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa1,cdh_vs_adc[%d]",detno[i+7]);
      check_slew_cdh[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"before_pa1,cnc1_vs_adc[%d]",detno[i+9]);
      raw_slew_cnc1[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa1,cnc1_vs_adc[%d]",detno[i+9]);
      check_slew_cnc1[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"before_pa1,cnc2_vs_adc[%d]",detno[i+11]);
      raw_slew_cnc2[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa1,cnc2_vs_adc[%d]",detno[i+11]);
      check_slew_cnc2[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax); 
   }
   for (int i = 0; i < 2; i++)
   {
      sprintf(hname,"before_pa2,cdh_vs_adc[%d]",detno[i+7]);
      raw_slew_cdh[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa2,cdh_vs_adc[%d]",detno[i+7]);
      check_slew_cdh[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"before_pa2,cnc1_vs_adc[%d]",detno[i+9]);
      raw_slew_cnc1[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa2,cnc1_vs_adc[%d]",detno[i+9]);
      check_slew_cnc1[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"before_pa2,cnc2_vs_adc[%d]",detno[i+11]);
      raw_slew_cnc2[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa2,cnc2_vs_adc[%d]",detno[i+11]);
      check_slew_cnc2[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
   }

   //loop
   double rawpa1, rawpa2, rawcdh, rawcnc1, rawcnc2;
   double slewpa1, slewpa2, slewcdh, slewcnc1, slewcnc2;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit(2))
      {
         rawpa1 = (time[3]+time[4])/2;
         rawpa2 = (time[5]+time[6])/2;
         rawcdh = (time[7]+time[8])/2;
         rawcnc1 = (time[9]+time[10])/2;
         rawcnc2 = (time[11]+time[12])/2;
         slewpa1 = (slew_par0[3] + slew_par1[3]/sqrt(charge[3]) + slew_par2[3]*sqrt(charge[3]) + slew_par0[4] + slew_par1[4]/sqrt(charge[4]) + slew_par2[4]*sqrt(charge[4]))/2;
         slewpa2 = (slew_par0[5] + slew_par1[5]/sqrt(charge[5]) + slew_par2[5]*sqrt(charge[5]) + slew_par0[6] + slew_par1[6]/sqrt(charge[6]) + slew_par2[6]*sqrt(charge[6]))/2;
         slewcdh = (slew_par0[7] + slew_par1[7]/sqrt(charge[7]) + slew_par2[7]*sqrt(charge[7]) + slew_par0[8] + slew_par1[8]/sqrt(charge[8]) + slew_par2[8]*sqrt(charge[8]))/2;
         slewcnc1 = (slew_par0[9] + slew_par1[9]/sqrt(charge[9]) + slew_par2[9]*sqrt(charge[9]) + slew_par0[10] + slew_par1[10]/sqrt(charge[10]) + slew_par2[10]*sqrt(charge[10]))/2;
         slewcnc2 = (slew_par0[11] + slew_par1[11]/sqrt(charge[11]) + slew_par2[11]*sqrt(charge[11]) + slew_par0[12] + slew_par1[12]/sqrt(charge[12]) + slew_par2[12]*sqrt(charge[12]))/2;
         
         //pa
         raw_time = rawpa1 - rawpa2;
         slewed_time = raw_time - slewpa1 + slewpa2;
         for (int i = 0; i < 4; i++)
         {
            raw_slew_pa[i]->Fill(charge[i+3],raw_time);
            check_slew_pa[i]->Fill(charge[i+3],slewed_time );
         }
         //cdh pa
         raw_time = rawcdh - rawpa1;
         slewed_time = raw_time + slewpa1 - slewcdh;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_cdh[i]->Fill(charge[i+7],raw_time);
            check_slew_cdh[i]->Fill(charge[i+7],slewed_time);
         }
         raw_time = rawcdh - rawpa2;
         slewed_time = raw_time + slewpa2 - slewcdh;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_cdh[i+2]->Fill(charge[i+7],raw_time);
            check_slew_cdh[i+2]->Fill(charge[i+7],slewed_time);
         }
         //cnc1 pa
         raw_time = rawcnc1 - rawpa1;
         slewed_time = raw_time + slewpa1 - slewcnc1;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_cnc1[i]->Fill(charge[i+9],raw_time);
            check_slew_cnc1[i]->Fill(charge[i+9],slewed_time);
         }
         raw_time = rawcnc1 - rawpa2;
         slewed_time = raw_time + slewpa2 - slewcnc1;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_cnc1[i+2]->Fill(charge[i+9],raw_time);
            check_slew_cnc1[i+2]->Fill(charge[i+9],slewed_time);
         }
         //cnc2 pa
         if (noofdet == 13)
         {
            raw_time = rawcnc2 - rawpa1;
            slewed_time = raw_time + slewpa1 - slewcnc2;
            for (int i = 0; i < 2; i++)
            {
               raw_slew_cnc2[i]->Fill(charge[i+11],raw_time);
               check_slew_cnc2[i]->Fill(charge[i+11],slewed_time);
            }
            raw_time = rawcnc2 - rawpa2;
            slewed_time = raw_time + slewpa2 - slewcnc2;
            for (int i = 0; i < 2; i++)
            {
               raw_slew_cnc2[i+2]->Fill(charge[i+9],raw_time);
               check_slew_cnc2[i+2]->Fill(charge[i+9],slewed_time);
            }
         }
      }
   }
   //pa1 pa2
   TCanvas *c1 = new TCanvas("c1","");
   TCanvas *c2 = new TCanvas("c2","");
   c1->Divide(2,2);
   c2->Divide(2,2);
   char fname[64];
   for (int i = 0; i < 4; i++)
   {
      c1->cd(i+1);
      mean[i] = raw_slew_pa[i]->GetMean(2);
      //gPad->DrawFrame(0.,mean[i]-5000.,5.,mean[i]+5000.);
      raw_slew_pa[i]->GetYaxis()->SetRangeUser(mean[i]-10000.,mean[i]+10000.);
      raw_slew_pa[i]->SetXTitle("Charge(E/E_mip)");
      raw_slew_pa[i]->SetYTitle("Time(ps)");
      raw_slew_pa[i]->Draw();
      //gPad->RedrawAxis();
      c2->cd(i+1);
      smean[i] = check_slew_pa[i]->GetMean(2);
      //gPad->DrawFrame(0.,smean[i]-5000.,5.,smean[i]+5000.);
      check_slew_pa[i]->GetYaxis()->SetRangeUser(smean[i]-10000.,smean[i]+10000.);
      check_slew_pa[i]->SetXTitle("Charge(E/E_mip)");
      check_slew_pa[i]->SetYTitle("Time(ps)");
      check_slew_pa[i]->Draw();
      //gPad->RedrawAxis();

   }
   sprintf(fname,"slew%5.5d_pa1_vs_pa2_before.pdf",runno);
   c1->Print(fname);
   sprintf(fname,"slew%5.5d_pa1_vs_pa2_after.pdf",runno);
   c2->Print(fname);
   
   //pa cdh
   TCanvas *c3 = new TCanvas("c3","");
   TCanvas *c4 = new TCanvas("c4","");
   c3->Divide(2,2);
   c4->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c3->cd(i+1);
      mean[i+4] = raw_slew_cdh[i]->GetMean(2);
      //gPad->DrawFrame(0.,mean[i+4]-5000.,5.,mean[i+4]+5000.);
      raw_slew_cdh[i]->GetYaxis()->SetRangeUser(mean[i+4]-10000.,mean[i+4]+10000.);
      raw_slew_cdh[i]->SetXTitle("Charge(E/E_mip)");
      raw_slew_cdh[i]->SetYTitle("Time(ps)");
      raw_slew_cdh[i]->Draw();
      //gPad->RedrawAxis();
      c4->cd(i+1);
      smean[i+4] = check_slew_cdh[i]->GetMean(2);
      //gPad->DrawFrame(0.,smean[i+4]-5000.,5.,smean[i+4]+5000.);
      check_slew_cdh[i]->GetYaxis()->SetRangeUser(smean[i+4]-10000.,smean[i+4]+10000.);
      check_slew_cdh[i]->SetXTitle("Charge(E/E_mip)");
      check_slew_cdh[i]->SetYTitle("Time(ps)");
      check_slew_cdh[i]->Draw();
      //gPad->RedrawAxis();
   }
   sprintf(fname,"slew%5.5d_pa_vs_cdh_before.pdf",runno);
   c3->Print(fname);
   sprintf(fname,"slew%5.5d_pa_vs_cdh_after.pdf",runno);
   c4->Print(fname);
   
   //pa cnc1
   TCanvas *c5 = new TCanvas("c5","");
   TCanvas *c6 = new TCanvas("c6","");
   c5->Divide(2,2);
   c6->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c5->cd(i+1);
      mean[i+8] = raw_slew_cnc1[i]->GetMean(2);
      //gPad->DrawFrame(0.,mean[i+8]-5000.,5.,mean[i+8]+5000.);
      raw_slew_cnc1[i]->SetXTitle("Charge(E/E_mip)");
      raw_slew_cnc1[i]->SetYTitle("Time(ps)");
      raw_slew_cnc1[i]->GetYaxis()->SetRangeUser(mean[i+8]-10000.,mean[i+8]+10000.);
      raw_slew_cnc1[i]->Draw();
      //gPad->RedrawAxis();
      c6->cd(i+1);
      smean[i+8] = check_slew_cnc1[i]->GetMean(2);
      //gPad->DrawFrame(0.,smean[i+8]-5000.,5.,smean[i+8]+5000.);
      check_slew_cnc1[i]->GetYaxis()->SetRangeUser(smean[i+8]-10000.,smean[i+8]+10000.);
      check_slew_cnc1[i]->SetXTitle("Charge(E/E_mip)");
      check_slew_cnc1[i]->SetYTitle("Time(ps)");
      check_slew_cnc1[i]->Draw("same");
      //gPad->RedrawAxis();
   }
   sprintf(fname,"slew%5.5d_pa_vs_cnc1_before.pdf",runno);
   c5->Print(fname);
   sprintf(fname,"slew%5.5d_pa_vs_cnc1_after.pdf",runno);
   c6->Print(fname);
   
   if (noofdet == 13)
   {
      //pa cnc2
      TCanvas *c7 = new TCanvas("c7","");
      TCanvas *c8 = new TCanvas("c8","");
      c7->Divide(2,2);
      c8->Divide(2,2);
      char fname[64];
      for (int i = 0; i < 4; i++)
      {
         c7->cd(i+1);
         mean[i+12] = raw_slew_cnc2[i]->GetMean(2);
         //gPad->DrawFrame(0.,mean[i+12]-5000.,5.,mean[i+12]+5000.);
         raw_slew_cnc2[i]->SetXTitle("Charge(E/E_mip)");
         raw_slew_cnc2[i]->SetYTitle("Time(ps)");
         raw_slew_cnc2[i]->GetYaxis()->SetRangeUser(mean[i+12]-10000.,mean[i+12]+10000);
         raw_slew_cnc2[i]->Draw();
         //gPad->RedrawAxis();
         c8->cd(i+1);
         smean[i+12] = check_slew_cnc2[i]->GetMean(2);
         //gPad->DrawFrame(0.,smean[i+12]-5000.,5.,smean[i+12]+5000.);
         check_slew_cnc2[i]->SetXTitle("Charge(E/E_mip)");
         check_slew_cnc2[i]->SetYTitle("Time(ps)");
         check_slew_cnc2[i]->GetYaxis()->SetRangeUser(smean[i+12]-10000.,smean[i+12]+10000);
         check_slew_cnc2[i]->Draw();
         //check_slew_cnc2[i]->GetYaxis()->ZoomOut(0.05,smean[i+12]-5000.);
         //gPad->RedrawAxis();
      }
      sprintf(fname,"slew%5.5d_pa_vs_cnc2_before.pdf",runno);
      c7->Print(fname);
      sprintf(fname,"slew%5.5d_pa_vs_cnc2_after.pdf",runno);
      c8->Print(fname);
   }

   /* //save slew mean
   char smname[64];
   sprintf(smname,"run%4.4d.slew_mean",runno);
   ofstream output_file(smname);
   output_file << fixed << setprecision(2);
   for (int i=0; i<16; i++) 
   {
     output_file << mean[i] << " " << smean[i] << endl;
   }
   output_file.close(); */
} 

void ana_eito::check_ref_slew()
{
   if (noofdet < 15) return;

   //parameter intiarise
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> slew_no[i];
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
         slew_no[i] = 0;
      }
   }
   slew_pass_file.close();

   //make hist
   char hname[64];
   double hbin = 10000.;
   double hmin = -200000.;
   double hmax = 200000.;
   for (int i = 0; i < 2; i++)
   {
      sprintf(hname,"before_pa1,ref_vs_adc[%d]",detno[i+13]);
      raw_slew_ref[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa1,ref_vs_adc[%d]",detno[i+13]);
      check_slew_ref[i] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"before_pa2,ref_vs_adc[%d]",detno[i+13]);
      raw_slew_ref[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
      sprintf(hname,"after_pa2,ref_vs_adc[%d]",detno[i+13]);
      check_slew_ref[i+2] = new TH2F(hname,hname,100,0.,5.,hbin,hmin,hmax);
   }

   //loop
   double rawpa1, rawpa2, rawref;
   double slewpa1, slewpa2, slewref;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit(2))
      {
         rawpa1 = (time[3]+time[4])/2;
         rawpa2 = (time[5]+time[6])/2;
         rawref = (time[13]+time[14])/2;
         slewpa1 = (slew_par0[3] + slew_par1[3]/sqrt(charge[3]) + slew_par2[3]*sqrt(charge[3]) + slew_par0[4] + slew_par1[4]/sqrt(charge[4]) + slew_par2[4]*sqrt(charge[4]))/2;
         slewpa2 = (slew_par0[5] + slew_par1[5]/sqrt(charge[5]) + slew_par2[5]*sqrt(charge[5]) + slew_par0[6] + slew_par1[6]/sqrt(charge[6]) + slew_par2[6]*sqrt(charge[6]))/2;
         slewref = (slew_par0[13] + slew_par1[13]/sqrt(charge[13]) + slew_par2[13]*sqrt(charge[13]) + slew_par0[14] + slew_par1[14]/sqrt(charge[14]) + slew_par2[14]*sqrt(charge[14]))/2;

         //ref pa
         raw_time = rawref - rawpa1;
         slewed_time = raw_time + slewpa1 - slewref;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_ref[i]->Fill(charge[i+13],raw_time);
            check_slew_ref[i]->Fill(charge[i+13],slewed_time);
         }
         raw_time = rawref - rawpa2;
         slewed_time = raw_time + slewpa2 - slewref;
         for (int i = 0; i < 2; i++)
         {
            raw_slew_ref[i+2]->Fill(charge[i+13],raw_time);
            check_slew_ref[i+2]->Fill(charge[i+13],slewed_time);
         }
      }
   }
   
   //pa ref
   TCanvas *c3 = new TCanvas("c3","");
   TCanvas *c4 = new TCanvas("c4","");
   c3->Divide(2,2);
   c4->Divide(2,2);
   double refmean, refsmean;
   for (int i = 0; i < 4; i++)
   {
      c3->cd(i+1);
      refmean = raw_slew_ref[i]->GetMean(2);
      //gPad->DrawFrame(0.,mean[i+4]-5000.,5.,mean[i+4]+5000.);
      raw_slew_ref[i]->GetYaxis()->SetRangeUser(refmean-10000.,refmean+10000.);
      raw_slew_ref[i]->SetXTitle("Charge(E/E_mip)");
      raw_slew_ref[i]->SetYTitle("Time(ps)");
      raw_slew_ref[i]->Draw();
      //gPad->RedrawAxis();
      c4->cd(i+1);
      refsmean = check_slew_ref[i]->GetMean(2);
      //gPad->DrawFrame(0.,smean[i+4]-5000.,5.,smean[i+4]+5000.);
      check_slew_ref[i]->GetYaxis()->SetRangeUser(-10000.,10000.);
      check_slew_ref[i]->SetXTitle("Charge(E/E_mip)");
      check_slew_ref[i]->SetYTitle("Time(ps)");
      check_slew_ref[i]->Draw();
      //gPad->RedrawAxis();
   }
   char fname[64];
   sprintf(fname,"slew%5.5d_pa_vs_ref_before.pdf",runno);
   c3->Print(fname);
   sprintf(fname,"slew%5.5d_pa_vs_ref_after.pdf",runno);
   c4->Print(fname);
}

void ana_eito::Convert() 
{
   for (int i=0; i<noofdet; i++)
  {
   //tdc
   if (leading[detno[i]][0]!=-1){
      rawtdc[i] =  leading[detno[i]][0]-t_offset[i];
      time[i]   = (leading[detno[i]][0]-t_offset[i])*timeConv[i];
   }else{
      rawtdc[i] = -1;
      time[i] = -1;
   }
   
   //adc
   if (7<=i && adc[detno[i]]>3800) 
   {
      charge[i] = -1.;
   } else {
      charge[i] =  (adc[detno[i]]-pedestal[i])/(mip[i]-pedestal[i]);
   }
  }
  return;
}

void ana_eito::dt_fit_ref(int test)
{
   //return;
   /*
      test = 0 cdh
           = 1 cnc1
           = 2 cnc2
           = 3 ref
   */

  //define fit channel name
   int test_r,test_l,ref1_r,ref1_l,ref2_r,ref2_l;
   std::vector<int> vec;
   ref1_r = 3;
   ref1_l = 4;
   ref2_r = 5;
   ref2_l = 6; 
   test_r = 13; 
   test_l = 14;
   vec.push_back(test_r);
   vec.push_back(test_l);
   if (test == 0) {ref1_r = 7; ref1_l = 8; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 1) {ref1_r = 9; ref1_l = 10; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 2) {ref1_r = 11; ref1_l = 12; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 3) {
      ref1_r = 3;
      ref1_l = 4;
      ref2_r = 5;
      ref2_l = 6;
      vec.push_back(ref1_r); vec.push_back(ref1_l);
      vec.push_back(ref2_r); vec.push_back(ref2_l);
   }
   char test_name[8];
   if (test == 0) sprintf(test_name,"cdh");
   if (test == 1) sprintf(test_name,"cnc1");
   if (test == 2) sprintf(test_name,"cnc2"); 
   if (test == 3) sprintf(test_name,"ref"); 

   //check ped_mip
   /* bool x_conv_frag = (pedestal[0]-mip[0])==0.||(pedestal[1]-mip[1])==0.
                        ||(pedestal[2]-mip[2])==0.||(pedestal[3]-mip[3])==0.;
   if(x_conv_frag)
   {
      cout << "do ped_mip() before phtw()" << endl;
      return;
   } */

   //slew parameter intiarise
   int dummy;
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> dummy;
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
      }
   }
   slew_pass_file.close();

   //make hist
   char histname[32];
   /* for (size_t i = 0; i < 3; i++)
   {
      if (i == 0) sprintf(histname,"dt_pa1_vs_pa2");
      if (i == 1) sprintf(histname,"dt_pa1_vs_%s",test_name);
      if (i == 2) sprintf(histname,"dt_pa2_vs_%s",test_name);
      h_dt[i] = new TH1F(histname,histname,20000,-30000,30000);
   } */
   sprintf(histname,"dt_ref_vs_%s",test_name);
   h_dt_ref = new TH1F(histname,histname,20000,-30000,30000);
   if(test == 3){
      for (size_t i = 0; i < 3; i++)
      {
         if (i == 0) sprintf(histname,"dt_pa1_vs_pa2");
         if (i == 1) sprintf(histname,"dt_pa1_vs_%s",test_name);
         if (i == 2) sprintf(histname,"dt_pa2_vs_%s",test_name);
         h_dt[i] = new TH1F(histname,histname,2000,-10000,10000);
      } 
   }

   if (fChain == 0) return;

   //loop fill
   double raw_ref1, raw_ref2, raw_test, slew_pa1, slew_pa2, slew_test;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit2(vec))
      {
         raw_ref1 = (time[ref1_r] + time[ref1_l])/2;
         raw_ref2 = (time[ref2_r] + time[ref2_l])/2;
         raw_test = (time[test_r] + time[test_l])/2;
         slew_pa1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]) 
                        + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2;
         slew_pa2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]) 
                        + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2;
         slew_test = (slew_par0[test_r] + slew_par1[test_r]/sqrt(charge[test_r]) + slew_par2[test_r]*sqrt(charge[test_r]) 
                        + slew_par0[test_l] + slew_par1[test_l]/sqrt(charge[test_l]) + slew_par2[test_l]*sqrt(charge[test_l]))/2;
         if (hit(0))
         {
            if(test == 3) {
               h_dt[0]->Fill(raw_ref1 - raw_ref2 - slew_pa1 + slew_pa2);
               h_dt[1]->Fill(raw_test - raw_ref1 - slew_test + slew_pa1);
               h_dt[2]->Fill(raw_test - raw_ref2 - slew_test + slew_pa2); 
            } 
            h_dt_ref->Fill(raw_test - raw_ref1 - slew_test + slew_pa1);
         }
      }
   }
   double sigma[3], sig, mean;
   double res[3] = {115., 168., 159.,};
   double pa1, pa2, tes;
   TFitResultPtr results1;
   TCanvas *c1 = new TCanvas("c1","");
   if(test == 3){
      c1->Divide(2,2);
      for (int i = 0; i < 3; i++)
      {
         c1->cd(i+1);
         mean = h_dt[i]->GetMean(1);
         h_dt[i]->GetXaxis()->SetRangeUser(mean - 2500.,mean + 2500.);
         h_dt[i]->Draw();
         results1 = h_dt[i]->Fit("gaus","S");
         sigma[i] = results1->Parameter(2);
      }   
   } else {
      mean = h_dt_ref->GetMean(1);
      h_dt_ref->GetXaxis()->SetRangeUser(mean - 2500.,mean + 2500.);
      h_dt_ref->Draw();
      results1 = h_dt_ref->Fit("gaus","S");
      sig = results1->Parameter(2);
      tes = sqrt(sig*sig - res[test]*res[test]);
      cout << "test :    " << tes << endl;
   }

   if(test == 3){
      pa1 = sqrt((sigma[0]*sigma[0] + sigma[1]*sigma[1] - sigma[2]*sigma[2])/2.);
      pa2 = sqrt((sigma[0]*sigma[0] + sigma[2]*sigma[2] - sigma[1]*sigma[1])/2.);
      tes = sqrt((sigma[2]*sigma[2] + sigma[1]*sigma[1] - sigma[0]*sigma[0])/2.);

      std::cout << "pa1 : " << pa1 << std::endl;
      std::cout << "pa2 : " << pa2 << std::endl;
      std::cout << "test : " << tes << std::endl;
   }
}

bool ana_eito::hit2(std::vector<int> vec)
{
   bool frag = false;
   for(int i = 0; i < vec.size(); i++){
      if(i == 0){
         if(time[vec.at(i)]!=-1) frag = true;
         else frag = false;
      } else {
         if(frag && time[vec.at(i)]!=-1) frag = true;
         else frag = false;
      }
   }
   return frag;
}



void ana_eito::dt_fit2(int t0, int t1, int t2)
{
   //return;
   /*
      test = 0 cdh
           = 1 cnc1
           = 2 cnc2
           = 3 ref
   */
  if(t0==t1||t1==t0||t2==t0){
   std::cout << "failed" << endl;
   return;
  }

  //define fit channel name
   int test_r,test_l,ref1_r,ref1_l,ref2_r,ref2_l, test;
   std::vector<int> vec;
   test = t0;
   if (test == 0) {ref1_r = 7; ref1_l = 8; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 1) {ref1_r = 9; ref1_l = 10; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 2) {ref1_r = 11; ref1_l = 12; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   if (test == 3) {ref1_r = 13; ref1_l = 14; vec.push_back(ref1_r); vec.push_back(ref1_l);}
   test = t1;
   if (test == 0) {ref2_r = 7; ref2_l = 8; vec.push_back(ref2_r); vec.push_back(ref2_l);}
   if (test == 1) {ref2_r = 9; ref2_l = 10; vec.push_back(ref2_r); vec.push_back(ref2_l);}
   if (test == 2) {ref2_r = 11; ref2_l = 12; vec.push_back(ref2_r); vec.push_back(ref2_l);}
   if (test == 3) {ref2_r = 13; ref2_l = 14; vec.push_back(ref2_r); vec.push_back(ref2_l);}
   test = t2;
   if (test == 0) {test_r = 7; test_l = 8; vec.push_back(test_r); vec.push_back(test_l);}
   if (test == 1) {test_r = 9; test_l = 10; vec.push_back(test_r); vec.push_back(test_l);}
   if (test == 2) {test_r = 11; test_l = 12; vec.push_back(test_r); vec.push_back(test_l);}
   if (test == 3) {test_r = 13; test_l = 14; vec.push_back(test_r); vec.push_back(test_l);}
   char test_name[8];
   if (test == 0) sprintf(test_name,"cdh");
   if (test == 1) sprintf(test_name,"cnc1");
   if (test == 2) sprintf(test_name,"cnc2"); 
   if (test == 3) sprintf(test_name,"ref"); 

   //check ped_mip
   /* bool x_conv_frag = (pedestal[0]-mip[0])==0.||(pedestal[1]-mip[1])==0.
                        ||(pedestal[2]-mip[2])==0.||(pedestal[3]-mip[3])==0.;
   if(x_conv_frag)
   {
      cout << "do ped_mip() before phtw()" << endl;
      return;
   }  */

   //slew parameter intiarise
   int dummy;
   char slew_pass_file_name[64]; 
   sprintf(slew_pass_file_name,"run%4.4d.slew_pass",runno);
   ifstream slew_pass_file(slew_pass_file_name);
   if (slew_pass_file) 
   {
      cout << "########### " << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> dummy;
      }
   } else {
      cout << "Fail to open " << slew_pass_file_name << endl;
      for (int i=0; i<noofdet; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
      }
   }
   slew_pass_file.close();

   //make hist
   char histname[32];
   for (size_t i = 0; i < 3; i++)
   {
      if (i == 0) sprintf(histname,"dt_1vs2");
      if (i == 1) sprintf(histname,"dt_1vs3");
      if (i == 2) sprintf(histname,"dt_2vs3");
      h_dt[i] = new TH1F(histname,histname,20000,-30000,30000);
   } 

   if (fChain == 0) return;

   //loop fill
   double raw_ref1, raw_ref2, raw_test, slew_pa1, slew_pa2, slew_test;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (hit2(vec))
      {
         raw_ref1 = (time[ref1_r] + time[ref1_l])/2;
         raw_ref2 = (time[ref2_r] + time[ref2_l])/2;
         raw_test = (time[test_r] + time[test_l])/2;
         slew_pa1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]) 
                        + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2;
         slew_pa2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]) 
                        + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2;
         slew_test = (slew_par0[test_r] + slew_par1[test_r]/sqrt(charge[test_r]) + slew_par2[test_r]*sqrt(charge[test_r]) 
                        + slew_par0[test_l] + slew_par1[test_l]/sqrt(charge[test_l]) + slew_par2[test_l]*sqrt(charge[test_l]))/2;
         if (hit(0))
         {
            if(1) {
               h_dt[0]->Fill(raw_ref1 - raw_ref2 - slew_pa1 + slew_pa2);
               h_dt[1]->Fill(raw_test - raw_ref1 - slew_test + slew_pa1);
               h_dt[2]->Fill(raw_test - raw_ref2 - slew_test + slew_pa2); 
            } 
         }
      }
   }
   double sigma[3], sig, mean;
   double res[3] = {115., 168., 159.,};
   double pa1, pa2, tes;
   TFitResultPtr results1;
   TCanvas *c1 = new TCanvas("c1","");
   if(1){
      c1->Divide(2,2);
      for (int i = 0; i < 3; i++)
      {
         c1->cd(i+1);
         mean = h_dt[i]->GetMean(1);
         h_dt[i]->GetXaxis()->SetRangeUser(mean - 2500.,mean + 2500.);
         h_dt[i]->Draw();
         results1 = h_dt[i]->Fit("gaus","S");
         sigma[i] = results1->Parameter(2);
      }   
   } else {
      mean = h_dt_ref->GetMean(1);
      h_dt_ref->GetXaxis()->SetRangeUser(mean - 2500.,mean + 2500.);
      h_dt_ref->Draw();
      results1 = h_dt_ref->Fit("gaus","S");
      sig = results1->Parameter(2);
      tes = sqrt(sig*sig - res[test]*res[test]);
      cout << "test :    " << tes << endl;
   }
   char a[32];
   sprintf(a,"%d_%d_%d_%d.pdf",runno,t0,t1,t2);
   c1->Print(a);

   if(1){
      pa1 = sqrt((sigma[0]*sigma[0] + sigma[1]*sigma[1] - sigma[2]*sigma[2])/2.);
      pa2 = sqrt((sigma[0]*sigma[0] + sigma[2]*sigma[2] - sigma[1]*sigma[1])/2.);
      tes = sqrt((sigma[2]*sigma[2] + sigma[1]*sigma[1] - sigma[0]*sigma[0])/2.);

      std::cout << "1 : " << pa1 << std::endl;
      std::cout << "2 : " << pa2 << std::endl;
      std::cout << "3 : " << tes << std::endl;
   } 
}