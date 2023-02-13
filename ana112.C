#define ana112_cxx
#include "ana112.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ana112::Load_file(int id){
   if(ana_frag==false){return;}

   runno = (double)id + 100000000;
   sprintf(work_dirname,"./ana/%d",id);
   std::filesystem::create_directory(work_dirname);
   sprintf(work_filename,"run%9.0lf",runno);
   sprintf(data_filename,"./data/%s.root",work_filename);
   cout << "load    :" << work_filename << endl;
   cout << "work at : " << work_dirname << endl;
   TTree *tree;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(data_filename);
   if (!f || !f->IsOpen()) 
   { f = new TFile(data_filename); }
   f->GetObject("tree",tree);
   ana_frag = false;
   Init(tree);
}

void ana112::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ana112.C
//      root> ana112 t
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
   if(ana_frag==true){return;}

   TH1D *h0[16], *h1[16];
   int hit[16];
   char histname[32];
   for(int i=0; i<16; i++){
      sprintf(histname,"raw%d",i);
      h0[i] = new TH1D(histname,histname,4095,0,4095);
      sprintf(histname,"hit%d",i);
      h1[i] = new TH1D(histname,histname,4095,0,4095);
      hit[i] = 0;
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
         if(tdc[i]!=0){
            h1[i]->Fill(adc[i]);
            hit[i] += 1;
         }
      }
   }

   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(4,4);
   for (int i = 0; i < 16; i++){
      c1->cd(i+1);
      gPad->SetLogy();
      h0[i]->Draw();
      if(hit[i]!=0){
         h1[i]->SetLineColor(2);
         h1[i]->Draw("same");
      }
   }
   char figname[64];
   sprintf(figname,"%s/%s_Loop.pdf",work_dirname,work_filename);
   c1->Print(figname);
}

void ana112::Phtw2(int ref_r,int ref_l, int cor_r,int cor_l)
{
   if(ana_frag==true){return;}
   if (fChain == 0) return;
   
   //set id
   std::vector<int> vec;
   vec.push_back(ref_r); vec.push_back(ref_l); 
   vec.push_back(cor_r); vec.push_back(cor_l); 


   //swel_no check
   if (slew_no[cor_r]!=slew_no[cor_l]||slew_no[ref_r]!=slew_no[ref_l])
   {
      /* cout << "input id is bat. retry." << endl;
      cout << "cor1 : " << cor1 << " slew_no : " <<slew_no[cor_r] << endl;
      cout << "cor2 : " << cor2 << " slew_no : " <<slew_no[cor_l] << endl;
      cout << "ref1 : " << ref1 << " slew_no : " <<slew_no[ref_r] << endl;
      cout << "ref2 : " << ref2 << " slew_no : " <<slew_no[ref_l] << endl;
      return; */
   }
   

   //hist, fit max min (x axis)
   double hist_max[2], hist_min[2], fit_max[2], fit_min[2];
   for(int i=0; i<2; i++){
      hist_max[i] = 5.0;
      hist_min[i] = 0.;
      fit_max[i] = 4.0;
      fit_min[i] = 0.5;
   }

   //make hist
   char hist_name[64];
   sprintf(hist_name,"ref_vs_adc[%d]",cor_r);
   h_slew[0] = new TH2D(hist_name,hist_name,32,hist_min[0],hist_max[0],1000,-150000,150000);
   sprintf(hist_name,"ref_vs_adc[%d]",cor_l);
   h_slew[1] = new TH2D(hist_name,hist_name,32,hist_min[0],hist_max[0],1000,-150000,150000);

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
      if (Hit(hit_frag))
      {
         t_cor_r = time[cor_r] - (slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r]));
         t_cor_l = time[cor_l] - (slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]));
         t_ref_r = time[ref_r] + (slew_par0[ref_r] + slew_par1[ref_r]/sqrt(charge[ref_r]) + slew_par2[ref_r]*sqrt(charge[ref_r]));
         t_ref_l = time[ref_l] + (slew_par0[ref_l] + slew_par1[ref_l]/sqrt(charge[ref_l]) + slew_par2[ref_l]*sqrt(charge[ref_l]));
         t_cor = (t_cor_r + t_cor_l)/2.0;
         t_ref = (t_ref_r + t_ref_l)/2.0;
         delta_time = t_cor - t_ref;
         if (slew_no[ref_r]==0||slew_no[ref_l]==0){
            if (charge[ref_r]>1.2&&charge[ref_l]>1.2) {
               h_slew[0]->Fill(charge[cor_r],delta_time);
               h_slew[1]->Fill(charge[cor_l],delta_time);
            }
         } else {
            h_slew[0]->Fill(charge[cor_r],delta_time);
            h_slew[1]->Fill(charge[cor_l],delta_time);
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
   h_slew[0]->ProfileX(prof_name);
   prof[0] = (TProfile*)gROOT->FindObject(prof_name);
   result = prof[0]->Fit("func","S","",fit_min[0],fit_max[0]);
   slew_par0[cor_r] += result->Parameter(0);
   slew_par1[cor_r] += result->Parameter(1);
   slew_par2[cor_r] += result->Parameter(2);
   slew_no[cor_r] += 1;
   //left
   c1->cd(2);
   sprintf(prof_name,"slew_left");
   h_slew[1]->ProfileX(prof_name);
   prof[1] = (TProfile*)gROOT->FindObject(prof_name);
   result = prof[1]->Fit("func","S","",fit_min[1],fit_max[1]);
   slew_par0[cor_l] += result->Parameter(0);
   slew_par1[cor_l] += result->Parameter(1);
   slew_par2[cor_l] += result->Parameter(2);
   slew_no[cor_l] += 1;
   
   //save
   //par
   ofstream slew_pass_ofile(slew_pass_filename);
   slew_pass_ofile << fixed << setprecision(2);
   for (int i = 0; i < 16; i++)
   {
      slew_pass_ofile << slew_par0[i] << " " << slew_par1[i] << " " << slew_par2[i] << " " << slew_no[i] << endl;
   }
   slew_pass_ofile.close();
}

void ana112::Phtw(int id)
{
   if(ana_frag==true){return;}
   //set ch number
   int ref1_r, ref1_l, ref2_r, ref2_l;
   int cor_r, cor_l;
   std::string det_name;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;
   if(id==0){ return; cor_r = 6; cor_l = 7; det_name = "ref"; if(detch[cor_r]==-1){return;}}
   else if (id==1){ cor_r = 6; cor_l = 7; det_name = "cdh"; if(detch[cor_r]==-1){return;}}
   else if (id==2){ cor_r = 8; cor_l = 9; det_name = "cnc1"; if(detch[cor_r]==-1){return;}}
   else if (id==3){ cor_r = 10; cor_l = 11; det_name = "cnc2"; if(detch[cor_r]==-1){return;}}
   //else if (id==4){ cor_r = 8; cor_l = 9; det_name = "test"; if(detch[cor_r]==-1){return;}}
   else if (id==5){ cor_r = 13; cor_l = 14; det_name = "mppc"; if(detch[cor_r]==-1){return;}}
   else {return;}

   //max min define
   double hist_min[4];
   double hist_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      hist_min[i] = 0.;
      //hist_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      hist_max[i] = 5.;
   }
   double fit_min[4];
   double fit_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.4;
      //fit_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      fit_max[i] = 4.0;
   }

   string hname;
   for (int i = 0; i < 4; i++)
   {
      hname = "ref1_vs_" + det_name + std::to_string(i);
      h_slew[i] = new TH2D(hname.c_str(),hname.c_str(),16,0.,5.,10000,-300000,300000);
   }
   
   func = new TF1("func","[0] + [1]/sqrt(x) + [2]*sqrt(x)",0.,7.);

   int slew_det[4];
   slew_det[0] = cor_l;
   slew_det[1] = cor_r;
   slew_det[2] = cor_l;
   slew_det[3] = cor_r;
   
   TProfile *prof[4];
   double t_ref1_r, t_ref1_l, t_ref2_r, t_ref2_l, t_cor_r, t_cor_l, t_cor, t_ref1, t_ref2, delta_time1, delta_time2;
   double mean[4];
   if (slew_no[cor_r]==0||slew_no[cor_l]==0)
   {
      Long64_t nentries = fChain->GetEntriesFast();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
      {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         Convert();
         if (Hit(hit_frag))
         {
            t_cor_r = time[cor_r] - (slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r]));
            t_cor_l = time[cor_l] - (slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]));
            t_ref1_r = time[ref1_r] + (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]));
            t_ref1_l = time[ref1_l] + (slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]));
            t_ref2_r = time[ref2_r] + (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]));
            t_ref2_l = time[ref2_l] + (slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]));
            t_cor = (t_cor_r + t_cor_l)/2.0;
            t_ref1 = (t_ref1_r + t_ref1_l)/2.0;
            t_ref2 = (t_ref2_r + t_ref2_l)/2.0;
            delta_time1 = t_cor - t_ref1;
            delta_time2 = t_cor - t_ref2;
            /* double t_ref1 = (time[ref1_r]+time[ref1_l])/2;
            double t_ref2 = (time[ref2_r]+time[ref2_l])/2;
            double t_test = (time[cor_r]+time[cor_l])/2;
            double delta_ref1 = slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                                 + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]);
            double delta_ref2 = slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                                 + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]);
            double delta_test = slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r])
                                 + slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]);                     
            double delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
            double delta_time1 = (t_test - delta_test) - (t_ref1 - delta_ref1);
            double delta_time2 = (t_test - delta_test) - (t_ref2 - delta_ref2); */
            if (charge[cor_l]>1.1&&charge[cor_r]>1.1) {
               h_slew[0]->Fill(charge[cor_r],delta_time1);
               h_slew[1]->Fill(charge[cor_l],delta_time1);
               h_slew[2]->Fill(charge[cor_r],delta_time2);
               h_slew[3]->Fill(charge[cor_l],delta_time2);
            }
         }
      }

      TCanvas *c1 = new TCanvas("c1","");
      c1->Divide(2,2);
      char prof_name[32];
      for (int i = 0; i < 4; i++)
      {
         c1->cd(i+1);
         sprintf(prof_name,"prof_hist%d",i);
         h_slew[i]->ProfileX(prof_name);
         prof[i] = (TProfile*)gROOT->FindObject(prof_name);
         //dt[i]->Draw();
         //prof[i]->Draw("same");
         TFitResultPtr result1 = prof[i]->Fit("func","S","",fit_min[i],fit_max[i]);
         prof[i]->GetYaxis()->SetRangeUser(mean[i]-200.,mean[i]+200.);
         //h_slew[i]->Draw("same");
         slew_par0[slew_det[i]] += result1->Parameter(0);
         slew_par1[slew_det[i]] += result1->Parameter(1);
         slew_par2[slew_det[i]] += result1->Parameter(2);
         slew_no[slew_det[i]] += 1;
      }
      char phtw_1st[64];
      sprintf(phtw_1st,"%s/%s_1stphtw_%s.pdf",work_dirname,work_filename,det_name.c_str());
      c1->Print(phtw_1st);

   //2nd or more slew
   } else {
      Long64_t nentries = fChain->GetEntriesFast();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
      {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         Convert();
         if (Hit(hit_frag))
         {
            t_cor_r = time[cor_r] - (slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r]));
            t_cor_l = time[cor_l] - (slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]));
            t_ref1_r = time[ref1_r] + (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]));
            t_ref1_l = time[ref1_l] + (slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]));
            t_ref2_r = time[ref2_r] + (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]));
            t_ref2_l = time[ref2_l] + (slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]));
            t_cor = (t_cor_r + t_cor_l)/2.0;
            t_ref1 = (t_ref1_r + t_ref1_l)/2.0;
            t_ref2 = (t_ref2_r + t_ref2_l)/2.0;
            delta_time1 = t_cor - t_ref1;
            delta_time2 = t_cor - t_ref2;
            /* double t_ref1 = (time[ref1_r]+time[ref1_l])/2;
            double t_ref2 = (time[ref2_r]+time[ref2_l])/2;
            double t_test = (time[cor_r]+time[cor_l])/2;
            double delta_ref1 = slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                                 + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]);
            double delta_ref2 = slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                                 + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]);
            double delta_test = slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r])
                                 + slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]);                     
            double delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
            double delta_time1 = (t_test - delta_test) - (t_ref1 - delta_ref1);
            double delta_time2 = (t_test - delta_test) - (t_ref2 - delta_ref2); */
            if (1) {
               h_slew[0]->Fill(charge[cor_l],delta_time1);
               h_slew[1]->Fill(charge[cor_r],delta_time1);
               h_slew[2]->Fill(charge[cor_l],delta_time2);
               h_slew[3]->Fill(charge[cor_r],delta_time2);
            }
         }
      }
      TCanvas *c1 = new TCanvas("c1","");
      c1->Divide(2,2);
      char prof_name[32];
      
      for (int i = 0; i < 4; i++)
      {
         c1->cd(i+1);
         sprintf(prof_name,"prof_hist%d",i+6);
         h_slew[i]->ProfileX(prof_name);
         prof[i] = (TProfile*)gROOT->FindObject(prof_name);
         //dt[i]->Draw();
         //prof[i]->Draw("same");
         TFitResultPtr result1 = prof[i]->Fit("func","S","",fit_min[i],fit_max[i]);
         prof[i]->GetYaxis()->SetRangeUser(mean[i]-200.,mean[i]+200.);
         //h_slew[i]->Draw("same");
         slew_par0[slew_det[i]] += result1->Parameter(0);
         slew_par1[slew_det[i]] += result1->Parameter(1);
         slew_par2[slew_det[i]] += result1->Parameter(2);
         slew_no[slew_det[i]] += 1;
      }
      char phtw[64];
      sprintf(phtw,"%s/%s_phtw_%s.pdf",work_dirname,work_filename,det_name.c_str());
      c1->Print(phtw);
   }

   //save parameter
   ofstream slew_pass_file(slew_pass_filename);
   slew_pass_file << fixed << setprecision(2);
   for (int i = 0; i < 16; i++)
   {
      slew_pass_file << slew_par0[i] << " " << slew_par1[i] << " " << slew_par2[i] << " " << slew_no[i] << endl;
   }
   slew_pass_file.close();
}

void ana112::Phtw_ref()
{
   if(ana_frag==true){return;}

   int ref1_r, ref1_l, ref2_r, ref2_l;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;

   //max min define
   double hist_min[4];
   double hist_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      hist_min[i] = 0.;
      //hist_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      hist_max[i] = 5.;
   }
   double fit_min[4];
   double fit_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.3;
      //fit_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      fit_max[i] = 4.0;
   }

   string hname;
   for (int i = 0; i < 4; i++)
   {
      hname = "ref" + std::to_string(i+2);
      h_ref_swel[i] = new TH2D(hname.c_str(),hname.c_str(),64,0.,5.,100000,-500000,500000);
   }

   func = new TF1("func","[0] + [1]/sqrt(x) + [2]*sqrt(x)",0.,7.);

   int slew_det[4];
   slew_det[0] = ref1_r;
   slew_det[1] = ref1_l;
   slew_det[2] = ref2_r;
   slew_det[3] = ref2_l;
   
   TProfile *prof[4];
   if (slew_no[ref1_r]==0||slew_no[ref2_l]==0)
   {
      Long64_t nentries = fChain->GetEntriesFast();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
      {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         Convert();
         if (Hit(2))
         {
            double t_ref1 = (time[ref1_r]+time[ref1_l])/2;
            double t_ref2 = (time[ref2_r]+time[ref2_l])/2;
            double delta_ref1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                                 + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2.;
            double delta_ref2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                                 + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2.;
            double delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
            if (1) {
               h_ref_swel[0]->Fill(charge[ref1_r],delta_ref);
               h_ref_swel[1]->Fill(charge[ref1_l],delta_ref);
               h_ref_swel[2]->Fill(charge[ref2_r],-1*delta_ref);
               h_ref_swel[3]->Fill(charge[ref2_l],-1*delta_ref);
            }
         }
      }

      TCanvas *c1 = new TCanvas("c1","");
      c1->Divide(2,2);
      char prof_name[32];
      double mean[4];
      for (int i = 0; i < 4; i++)
      {
         c1->cd(i+1);
         sprintf(prof_name,"prof_hist%d",i);
         h_ref_swel[i]->ProfileX(prof_name);
         prof[i] = (TProfile*)gROOT->FindObject(prof_name);
         mean[i] = prof[i]->GetMean(2);
         //prof[i]->GetYaxis()->SetRangeUser(mean[i]-5000.,mean[i]+5000.);
         //[i]->Draw();
         //prof[i]->Draw("same");
         TFitResultPtr result1 = prof[i]->Fit("func","S","",fit_min[i],fit_max[i]);
         prof[i]->GetYaxis()->SetRangeUser(mean[i]-500.,mean[i]+500.);
         //h_slew[i]->Draw("same");
         slew_par0[slew_det[i]] += result1->Parameter(0);
         slew_par1[slew_det[i]] += result1->Parameter(1);
         slew_par2[slew_det[i]] += result1->Parameter(2);
         slew_no[slew_det[i]] += 1;
      }
      char phtw_1st[64];
      sprintf(phtw_1st,"%s/%s_1stphtw_ref.pdf",work_dirname,work_filename);
      c1->Print(phtw_1st);

   //2nd or more slew
   } else {
      Long64_t nentries = fChain->GetEntriesFast();

      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
      {
         Long64_t ientry = LoadTree(jentry);
         if (ientry < 0) break;
         nb = fChain->GetEntry(jentry);   nbytes += nb;
         Convert();
         if (Hit(hit_frag))
         {
            double t_ref1 = (time[ref1_r]+time[ref1_l])/2;
            double t_ref2 = (time[ref2_r]+time[ref2_l])/2;
            double delta_ref1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                                 + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2.;
            double delta_ref2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                                 + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2.;
            double delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
            if (1) {
               h_ref_swel[0]->Fill(charge[ref1_r],delta_ref);
               h_ref_swel[1]->Fill(charge[ref1_l],delta_ref);
               h_ref_swel[2]->Fill(charge[ref2_r],-1*delta_ref);
               h_ref_swel[3]->Fill(charge[ref2_l],-1*delta_ref);
            }
         }
      }
      TCanvas *c1 = new TCanvas("c1","");
      c1->Divide(2,2);
      char prof_name[32];
      double mean[4];
      for (int i = 0; i < 4; i++)
      {
         c1->cd(i+1);
         sprintf(prof_name,"prof_hist%d",i+6);
         h_ref_swel[i]->ProfileX(prof_name);
         prof[i] = (TProfile*)gROOT->FindObject(prof_name);
         mean[i] = prof[i]->GetMean(2);
         //prof[i]->GetYaxis()->SetRangeUser(mean[i]-5000.,mean[i]+5000.);
         //dt[i]->Draw();
         //prof[i]->Draw("same");
         TFitResultPtr result1 = prof[i]->Fit("func","S","",fit_min[i],fit_max[i]);
         prof[i]->GetYaxis()->SetRangeUser(mean[i]-500.,mean[i]+500.);
         //h_slew[i]->Draw("same");
         slew_par0[slew_det[i]] += result1->Parameter(0);
         slew_par1[slew_det[i]] += result1->Parameter(1);
         slew_par2[slew_det[i]] += result1->Parameter(2);
         slew_no[slew_det[i]] += 1;
      }
      char phtw[64];
      sprintf(phtw,"%s/%s_phtw_ref.pdf",work_dirname,work_filename);
      c1->Print(phtw);
   }

  //save parameter
   ofstream slew_pass_file(slew_pass_filename);
   slew_pass_file << fixed << setprecision(2);
   for (int i = 0; i < 16; i++)
   {
      slew_pass_file << slew_par0[i] << " " << slew_par1[i] << " " << slew_par2[i] << " " << slew_no[i] << endl;
   }
   slew_pass_file.close();

}

void ana112::Ped_mip()
{
   if(ana_frag==true){return;}

   //return;
   //initiarise
   for (int i=0; i<16; i++){
      pedestal[i]=0.;
      ped_sigma[i]=0.;
      mip[i]=0.;
      t_offset[i]=0.;
      t_width[i]=0.;
   }

   //make hist
   char chnum[32];
   for (int i = 0; i < 16; i++){
      sprintf(chnum,"h_adc_%d",i);
      h_adc[i] = new TH1D(chnum,chnum,4095,0.,4095.);
      sprintf(chnum,"h_adcmip_%d",i);
      h_adcmip[i] = new TH1D(chnum,chnum,4095,0.,4095.);
      sprintf(chnum,"h_adcped_%d",i);
      h_adcped[i] = new TH1D(chnum,chnum,500,0.,500.);
      sprintf(chnum,"h_tdc_%d",i);
      h_tdc[i] = new TH1D(chnum,chnum,4095,0.,4095.);
   }

   if (fChain == 0) return;

   //loop
   int nhit[16];
   for(int i=0; i<16; i++){nhit[i]=0;}
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //Convert();
      for (int i=0; i<16; i++) {
        h_adc[i]->Fill(adc[i]);
        if(tdc[i]==0){ h_adcped[i]->Fill(adc[i]); }
        if(tdc[i]!=0){ h_adcmip[i]->Fill(adc[i]);
                       h_tdc[i]->Fill(tdc[i]); 
                       nhit[i]++;}
      }
   }
   
   //ped
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(4,4);
   TFitResultPtr results1;
   double maxbin_ped[16];
   for (int i = 0; i < 16; i++){
      if(1){
         c1->cd(i+1);
         if (nhit[i]>0)
         {
         if(i>1 && i!=3 && i!=5){
            maxbin_ped[i] = h_adcped[i]->GetMaximumBin();
            results1 = h_adcped[i]->Fit("gaus","S","",maxbin_ped[i]-50,maxbin_ped[i]+50);
            pedestal[i] = results1->Parameter(1);
            ped_sigma[i] = results1->Parameter(2);
         } else if(i==3 || i==5){
            maxbin_ped[i] = h_adcped[i]->GetMaximumBin();
            results1 = h_adcped[i]->Fit("gaus","S","",maxbin_ped[i]-10,maxbin_ped[i]+10);
            pedestal[i] = results1->Parameter(1);
            ped_sigma[i] = results1->Parameter(2);
         } else {
            pedestal[i] = 0.;
            ped_sigma[i] = 0.;
         }
         }
       else {
         pedestal[i] = 0.;
         ped_sigma[i] = 0.;
      }
      }
   }
   char ped_figname[64];
   sprintf(ped_figname,"%s/%s_ped.pdf",work_dirname,work_filename);
   c1->Print(ped_figname);

   //mip
   TCanvas *c2 = new TCanvas("c2","c2");
   c2->Divide(4,4);
   TFitResultPtr results2;
   for (int i = 0; i < 16; i++){

      if (nhit[i]>0)
      {
         c2->cd(i+1);
         results2 = h_adcmip[i]->Fit("landau","S","",400,3500);
         mip[i] = results2->Parameter(1);
      }
      else {
      mip[i] = 0.;
   }
   }
   char mip_figname[64];
   sprintf(mip_figname,"%s/%s_mip.pdf",work_dirname,work_filename);
   c2->Print(mip_figname);

   //tdc offset
   TCanvas *c3 = new TCanvas("c3","c3");
   c3->Divide(4,4);   
   TFitResultPtr results3;
   for (int i=0; i<16; i++){
      if(1){
         c3->cd(i+1);
         if (nhit[i]>0)
         {
         results3  = h_tdc[i]  ->Fit("gaus","S");
         t_offset[i]   = results3->Parameter(1);
         t_width[i] = results3->Parameter(2);
         }
       else {
         t_offset[i] = 0.;
         t_width[i] = 0.;
      }
      }
   }
   char tdc_figname[64]; 
   sprintf(tdc_figname,"%s/%s_tdc.pdf",work_dirname,work_filename);
   c3->Print(tdc_figname);

   //write file
   ofstream pedmip_file(pede_filename);
   pedmip_file << fixed << setprecision(2);
   for (int i=0; i<16; i++) 
   {
     pedmip_file << pedestal[i]  << " " << ped_sigma[i] << " "<< mip[i] << endl;
   }
   pedmip_file.close();

   ofstream offset_file(offset_filename);
   offset_file << fixed << setprecision(2);
   for (int i=0; i<16; i++) 
   {
     offset_file << t_offset[i] << " " << t_width[i] << endl;
   }
   offset_file.close();
}

void ana112::Slewing(int id,int num){
   if(ana_frag==true){return;}
   //return
   //slewing
   for(int i=0; i<num; i++){
      if (id==0) {
         Phtw_ref();
      } else {
         Phtw(id);
      }
   }
   if (id==0) {
      Check_ref_slew();
   } else {
      Check_slew(id);
   }
}

void ana112::Resolution(int id){
   if(ana_frag==true){return;}
   //return;
   // i = 0 ref vs ref
   // i = 1 ref vs cdh
   // i = 2 ref vs cnc1
   // i = 3 ref vs cnc2
   // i = 4 ref vs test

   //set ch number
   int ref1_r, ref1_l, ref2_r, ref2_l;
   int test_r, test_l;
   std::string det_name;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;
   if(id==0){ Reso_ref(); return; det_name = "ref"; }
   else if (id==1){ test_r = 6;  test_l = 7;  det_name = "cdh"; if(detch[test_r]==-1){return;}}
   else if (id==2){ test_r = 8;  test_l = 9;  det_name = "cnc1"; if(detch[test_r]==-1){return;}}
   else if (id==3){ test_r = 10; test_l = 11; det_name = "cnc2"; if(detch[test_r]==-1){return;}}
   else if (id==4){ test_r = 8;  test_l = 9;  det_name = "test"; if(detch[test_r]==-1){return;}}
   else if (id==5){ test_r = 13;  test_l = 14;  det_name = "mppc"; if(detch[test_r]==-1){return;}}
   else { cout << "id=0~5" << endl; return; }

   //hist, fit max min (x axis)
   double hist_max[3], hist_min[3], fit_max[3], fit_min[3];
   for(int i=0; i<3; i++){
      hist_max[i] = 100000.;
      hist_min[i] = -100000.;
      fit_max[i] = 5000.;
      fit_min[i] = -5000;
   }

   //make hist
   char hist_name[64];
   sprintf(hist_name,"ref1_vs_ref2");
   h_reso[0] = new TH1D(hist_name, hist_name, 100000, hist_min[0], hist_max[0]);
   sprintf(hist_name,"ref1_vs_%s",det_name.c_str());
   h_reso[1] = new TH1D(hist_name, hist_name, 100000, hist_min[1], hist_max[1]);
   sprintf(hist_name,"ref2_vs_%s",det_name.c_str());
   h_reso[2] = new TH1D(hist_name, hist_name, 100000, hist_min[2], hist_max[2]);

   //Loop
   if (fChain == 0) return;
   double raw_ref1, raw_ref2, raw_test, slew_pa1, slew_pa2, slew_test;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (1){
         raw_ref1 = (time[ref1_r] + time[ref1_l])/2;
         raw_ref2 = (time[ref2_r] + time[ref2_l])/2;
         raw_test = (time[test_r] + time[test_l])/2;
         slew_pa1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]) 
                        + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2;
         slew_pa2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]) 
                        + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2;
         slew_test = (slew_par0[test_r] + slew_par1[test_r]/sqrt(charge[test_r]) + slew_par2[test_r]*sqrt(charge[test_r]) 
                        + slew_par0[test_l] + slew_par1[test_l]/sqrt(charge[test_l]) + slew_par2[test_l]*sqrt(charge[test_l]))/2;
         
         if (Hit(hit_frag))
         {
            h_reso[0]->Fill(raw_ref1 - raw_ref2 - slew_pa1 + slew_pa2);
            h_reso[1]->Fill(raw_test - raw_ref1 - slew_test + slew_pa1);
            h_reso[2]->Fill(raw_test - raw_ref2 - slew_test + slew_pa2);
         }
      }
   }

   //fit
   //draw fit
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->Divide(2,2);
   TFitResultPtr result;
   double mean[3];
   double sigma[3];
   for (int i = 0; i < 3; i++){
      if(id==0&&i>0)continue;
      c1->cd(i+1);
      mean[i] = h_reso[i]->GetMean(1);
      h_reso[i]->GetXaxis()->SetRangeUser(mean[i] - 1000.,mean[i] + 1000.);
      h_reso[i]->Draw();
      result = h_reso[i]->Fit("gaus","S","",mean[i] - 1000.,mean[i] + 1000.);
      sigma[i] = result->Parameter(2);
   }
   char figname[64];
   sprintf(figname,"%s/%s_reso_%s.pdf",work_dirname,work_filename,det_name.c_str());
   c1->Print(figname);

   double pa1_res, pa2_res, test_res, test_res1, test_res2;
      pa1_res = sigma[0]/sqrt(2);
      pa2_res = sigma[0]/sqrt(2);
      //pa1_res = sqrt((sigma[0]*sigma[0] + sigma[1]*sigma[1] - sigma[2]*sigma[2])/2);
      //pa2_res = sqrt((sigma[0]*sigma[0] + sigma[2]*sigma[2] - sigma[1]*sigma[1])/2); 
      //test_res = sqrt((sigma[1]*sigma[1] + sigma[2]*sigma[2] - sigma[0]*sigma[0])/2);
      test_res1 = sqrt(sigma[1]*sigma[1] - pa1_res*pa1_res);
      test_res2 = sqrt(sigma[2]*sigma[2] - pa2_res*pa2_res);
      if(test_res1 > test_res2){
         resolution[test_l] = test_res2;
         resolution[test_r] = test_res2;
      } else {
         resolution[test_l] = test_res1;
         resolution[test_r] = test_res1;
      }
   
   

   //save
   ofstream reso_file(reso_filename);
   reso_file << fixed << setprecision(2);
   for (int i=0; i<16; i++) {
     reso_file << i << " " << resolution[i] << endl;
   }
   reso_file.close();

   //print result
   if(id==0){
      cout << endl;
      cout << "####################################" << endl;
      cout << "          time resolution           " << endl;
      cout << "       " << det_name << "   = " << test_res << " ps" << endl;
      //cout << "       test   = " << test_res << " ps" << endl;
      cout << "####################################" << endl;
   } else {
      cout << endl;
      cout << "####################################" << endl;
      cout << "          time resolution           " << endl;
      cout << "       pa1    = " << pa1_res << " ps" << endl;
      cout << "       pa2    = " << pa2_res << " ps" << endl;
      //cout << "       " << det_name << "   = " << test_res << " ps" << endl;
      cout << "       " << det_name << "   = " << test_res1 << " ps" << endl;
      cout << "       " << det_name << "   = " << test_res2 << " ps" << endl;
      //cout << "       test   = " << test_res << " ps" << endl;
      cout << "####################################" << endl;
   }
}

void ana112::Reso_ref(){
   if(ana_frag==true){return;}

   int ref1_r, ref1_l, ref2_r, ref2_l;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;
   double hmax, hmin;
   hmax = 150000.;
   hmin = -150000.;

   h_ref_reso = new TH1D("ref_reso","ref_reso",(int)((hmax-hmin)/16),hmin,hmax);

   double t_ref1, t_ref2, delta_ref1, delta_ref2, delta_ref;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (Hit(hit_frag))
      {
         t_ref1 = (time[ref1_r]+time[ref1_l])/2;
         t_ref2 = (time[ref2_r]+time[ref2_l])/2;
         delta_ref1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                              + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2.;
         delta_ref2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                              + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2.;
         delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
         if (1) {
            h_ref_reso->Fill(delta_ref);
         }
      }
   }
   TCanvas *c1 = new TCanvas("c1", "c1");
   double mean,reso;
   TFitResultPtr result;
   mean = h_ref_reso->GetMean();
   h_ref_reso->GetXaxis()->SetRangeUser(mean-1000.,mean+1000.);
   h_ref_reso->Draw();
   result = h_ref_reso->Fit("gaus","S","",mean-1000.,mean+1000.);
   reso = (result->Parameter(2))/sqrt(2);
   char figname[64];
   sprintf(figname,"%s/%s_reso_ref.pdf",work_dirname,work_filename);
   c1->Print(figname);

   resolution[ref1_l] = reso;
   resolution[ref1_r] = reso;
   resolution[ref2_l] = reso;
   resolution[ref2_r] = reso;

   ofstream reso_file(reso_filename);
   reso_file << fixed << setprecision(2);
   for (int i=0; i<16; i++) {
     reso_file << i << " " << resolution[i] << endl;
   }
   reso_file.close();

   cout << endl;
   cout << "####################################" << endl;
   cout << "          time resolution           " << endl;
   cout << "       ref     = " << reso << " ps" << endl;
   //cout << "       test   = " << test_res << " ps" << endl;
   cout << "####################################" << endl;


}

void ana112::Check_slew(int id){
   if(ana_frag==true){return;}
   //set ch number
   int ref1_r, ref1_l, ref2_r, ref2_l;
   int cor_r, cor_l;
   std::string det_name;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;
   if(id==0){ return; cor_r = 6; cor_l = 7; det_name = "ref"; }
   else if (id==1){ cor_r = 6; cor_l = 7; det_name = "cdh"; }
   else if (id==2){ cor_r = 8; cor_l = 9; det_name = "cnc1"; }
   else if (id==3){ cor_r = 10; cor_l = 11; det_name = "cnc2"; }
   else if (id==4){ cor_r = 8; cor_l = 9; det_name = "test"; }
   else if (id==5){ cor_r = 13; cor_l = 14; det_name = "test"; }
   else {return;}

   //max min define
   double hist_min[4];
   double hist_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      hist_min[i] = 0.;
      //hist_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      hist_max[i] = 5.;
   }
   double fit_min[4];
   double fit_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.1;
      //fit_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      fit_max[i] = 4.0;
   }

   string hname;
   for (int i = 0; i < 4; i++)
   {
      hname = "before_bref1_vs_" + det_name + std::to_string(i);
      h_slew_before[i] = new TH2D(hname.c_str(),hname.c_str(),2000,0.,5.,100000,-300000,300000);
      hname = "after_bref1_vs_" + det_name + std::to_string(i);
      h_slew_after[i] = new TH2D(hname.c_str(),hname.c_str(),2000,0.,5.,100000,-300000,300000);
   }

   int slew_det[4];
   slew_det[0] = cor_l;
   slew_det[1] = cor_r;
   slew_det[2] = cor_l;
   slew_det[3] = cor_r;
   
   func = new TF1("func","[0] + [1]/sqrt(x) + [2]*sqrt(x)",0.,7.);
   
   double t_ref1_r, t_ref1_l, t_ref2_r, t_ref2_l, t_cor_r, t_cor_l, t_cor, t_ref1, t_ref2, delta_time1, delta_time2;
   double raw_ref1, raw_ref2, raw_cor, time1, time2;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (Hit(hit_frag))
      {
         t_cor_r = time[cor_r] - (slew_par0[cor_r] + slew_par1[cor_r]/sqrt(charge[cor_r]) + slew_par2[cor_r]*sqrt(charge[cor_r]));
         t_cor_l = time[cor_l] - (slew_par0[cor_l] + slew_par1[cor_l]/sqrt(charge[cor_l]) + slew_par2[cor_l]*sqrt(charge[cor_l]));
         t_ref1_r = time[ref1_r] + (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r]));
         t_ref1_l = time[ref1_l] + (slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]));
         t_ref2_r = time[ref2_r] + (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r]));
         t_ref2_l = time[ref2_l] + (slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]));
         t_cor = (t_cor_r + t_cor_l)/2.0;
         t_ref1 = (t_ref1_r + t_ref1_l)/2.0;
         t_ref2 = (t_ref2_r + t_ref2_l)/2.0;
         delta_time1 = t_cor - t_ref1;
         delta_time2 = t_cor - t_ref2;
         raw_ref1 = (time[ref1_r] + time[ref1_l])/2.;
         raw_ref2 = (time[ref2_r] + time[ref2_l])/2.;
         raw_cor = (time[cor_r] + time[cor_l])/2.;
         time1 = raw_cor - t_ref1;
         time2 = raw_cor - t_ref2;
         if (1) {
            h_slew_before[0]->Fill(charge[cor_l],time1);
            h_slew_before[1]->Fill(charge[cor_r],time1);
            h_slew_before[2]->Fill(charge[cor_l],time2);
            h_slew_before[3]->Fill(charge[cor_r],time2);
            h_slew_after[0]->Fill(charge[cor_l],delta_time1);
            h_slew_after[1]->Fill(charge[cor_r],delta_time1);
            h_slew_after[2]->Fill(charge[cor_l],delta_time2);
            h_slew_after[3]->Fill(charge[cor_r],delta_time2);
         }
      }
   }
   TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,2);
   double mean[4];
   for (int i = 0; i < 4; i++)
   {
      c1->cd(i+1);
      mean[i] = h_slew_before[i]->GetMean(2);
      h_slew_before[i]->GetYaxis()->SetRangeUser(mean[i]-7000.,mean[i]+7000.);
      h_slew_before[i]->Draw("colz");
   }
   char phtw[64];
   sprintf(phtw,"%s/%s_before_%s.pdf",work_dirname,work_filename,det_name.c_str());
   c1->Print(phtw);

   TCanvas *c2 = new TCanvas("c2","");
   c2->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c2->cd(i+1);
      mean[i] = h_slew_after[i]->GetMean(2);
      h_slew_after[i]->GetYaxis()->SetRangeUser(mean[i]-7000.,mean[i]+7000.);
      h_slew_after[i]->Draw("colz");
   }
   sprintf(phtw,"%s/%s_after_%s.pdf",work_dirname,work_filename,det_name.c_str());
   c2->Print(phtw);
}

void ana112::Check_ref_slew(){
   if(ana_frag==true){return;}
   //return;
   int ref1_r, ref1_l, ref2_r, ref2_l;
   ref1_r = 2; ref1_l = 3;
   ref2_r = 4; ref2_l = 5;

   //max min define
   double hist_min[4];
   double hist_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      hist_min[i] = 0.;
      //hist_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      hist_max[i] = 5.;
   }
   double fit_min[4];
   double fit_max[4];
   for (size_t i = 0; i < 4; i++)
   {
      //fit_min[i] = 10.0*ped_sigma[i]/(mip[i] - pedestal[i]);
      fit_min[i] = 0.4;
      //fit_max[i] = (4095. - pedestal[i])/(mip[i] - pedestal[i]);
      fit_max[i] = 4.0;
   }

   string hname;
   for (int i = 0; i < 4; i++)
   {
      hname = "ref_before" + std::to_string(i+2);
      h_ref_before[i] = new TH2D(hname.c_str(),hname.c_str(),2000,0.,5.,100000,-300000,300000);
      hname = "ref_after" + std::to_string(i+2);
      h_ref_after[i] = new TH2D(hname.c_str(),hname.c_str(),2000,0.,5.,100000,-300000,300000);
   }

   func = new TF1("func","[0] + [1]/sqrt(x) + [2]*sqrt(x)",0.,7.);

   int slew_det[4];
   slew_det[0] = ref1_r;
   slew_det[1] = ref1_l;
   slew_det[2] = ref2_r;
   slew_det[3] = ref2_l;
   
   TProfile *prof[4];
   double t_ref1, t_ref2, delta_ref1, delta_ref2, delta_ref, delta_ref0;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      Convert();
      if (Hit(hit_frag))
      {
         t_ref1 = (time[ref1_r]+time[ref1_l])/2;
         t_ref2 = (time[ref2_r]+time[ref2_l])/2;
         delta_ref1 = (slew_par0[ref1_r] + slew_par1[ref1_r]/sqrt(charge[ref1_r]) + slew_par2[ref1_r]*sqrt(charge[ref1_r])
                              + slew_par0[ref1_l] + slew_par1[ref1_l]/sqrt(charge[ref1_l]) + slew_par2[ref1_l]*sqrt(charge[ref1_l]))/2.;
         delta_ref2 = (slew_par0[ref2_r] + slew_par1[ref2_r]/sqrt(charge[ref2_r]) + slew_par2[ref2_r]*sqrt(charge[ref2_r])
                              + slew_par0[ref2_l] + slew_par1[ref2_l]/sqrt(charge[ref2_l]) + slew_par2[ref2_l]*sqrt(charge[ref2_l]))/2.;
         delta_ref =  (t_ref1 - delta_ref1) - (t_ref2 - delta_ref2);
         delta_ref0 = t_ref1 - t_ref2;
         if (1) {
            h_ref_before[0]->Fill(charge[ref1_r],delta_ref0);
            h_ref_before[1]->Fill(charge[ref1_l],delta_ref0);
            h_ref_before[2]->Fill(charge[ref2_r],delta_ref0);
            h_ref_before[3]->Fill(charge[ref2_l],delta_ref0);
            h_ref_after[0]->Fill(charge[ref1_r],delta_ref);
            h_ref_after[1]->Fill(charge[ref1_l],delta_ref);
            h_ref_after[2]->Fill(charge[ref2_r],delta_ref);
            h_ref_after[3]->Fill(charge[ref2_l],delta_ref);
         }
      }
   }
   TCanvas *c1 = new TCanvas("c1","");
   c1->Divide(2,2);
   char prof_name[32];
   double mean[4];
   for (int i = 0; i < 4; i++)
   {
      c1->cd(i+1);
      mean[i] = h_ref_before[i]->GetMean(2);
      h_ref_before[i]->GetYaxis()->SetRangeUser(mean[i]-2000.,mean[i]+2000.);
      h_ref_before[i]->Draw("colz");
   }
   char phtw[64];
   sprintf(phtw,"%s/%s_before_ref.pdf",work_dirname,work_filename);
   c1->Print(phtw);
   TCanvas *c2 = new TCanvas("c2","");
   c2->Divide(2,2);
   for (int i = 0; i < 4; i++)
   {
      c2->cd(i+1);
      mean[i] = h_ref_after[i]->GetMean(2);
      h_ref_after[i]->GetYaxis()->SetRangeUser(mean[i]-2000.,mean[i]+2000.);
      h_ref_after[i]->Draw("colz");
   }
   sprintf(phtw,"%s/%s_after_ref.pdf",work_dirname,work_filename);
   c1->Print(phtw);
}

void ana112::Convert(){
   //return;
   //adc
   for(int i=0; i<16; i++){
      if(mip[i]!=0. && adc[i]<3500){
         charge[i] = (adc[i]-pedestal[i])/(mip[i]-pedestal[i]);
      } else {
         charge[i] = 0.;
      }
   }

   //tdc
   for(int i=0; i<16; i++){
      if(t_offset[i]!=0.){
         time[i] = timeConv*(tdc[i] - t_offset[i]);
      } else {
         time[i] = 0.;
      }
   }
}

bool ana112::Hit(int id){
   //return;
   bool frag = true;
   for(int i=0; i<16; i++){
      if(detch[i]==-1) continue;
      if(tdc[i]==0){frag = false; break;}
   }
   if(id==6||id==7){
      if(frag && adc[1]<1200) frag = true;
      else frag = false;
      if(frag && adc[8]<1200 && adc[9]<1200) frag = true;
      else frag = false;
      if(frag && adc[10]<1350 && adc[11]<1350) frag = true;
      else frag = false;
      if(frag && adc[13]<1500 && adc[14]<2000) frag = true;
      else frag = false;
   } else if (id==8||id==9)
   {
      if(frag && adc[1]<1200) frag = true;
      else frag = false;
      if(frag && adc[6]<1500 && adc[7]<1400) frag = true;
      else frag = false;
      if(frag && adc[10]<1350 && adc[11]<1350) frag = true;
      else frag = false;
      if(frag && adc[13]<1500 && adc[14]<2000) frag = true;
      else frag = false;
   }
   
}

void ana112::Hit(std::vector<int> vec){
   return;
}

void ana112::Set_exist_det(){
   //return;
   string num;
   cout << endl;
   cout << "0:  cdh:yes, cnc1:yes, cnc2:yes, trig3:no,  mppc:no " << endl;
   cout << "1:  cdh:no,  cnc1:yes, cnc2:yes, trig3:no,  mppc:no " << endl;
   cout << "2:  cdh:no,  cnc1:yes, cnc2:yes, trig3:yes, mppc:no " << endl;
   cout << "3:  cdh:yes, cnc1:yes, cnc2:yes, trig3:yes, mppc:no " << endl;
   cout << "4:  cdh:yes, cnc1:yes, cnc2:no,  trig3:yes, mppc:no " << endl;
   cout << "5:  cdh:yes, cnc1:yes, cnc2:no,  trig3:yes, mppc:yes" << endl;
   cout << "6:  cdh:yes, cnc1:yes, cnc2:no,  trig3:no,  mppc:yes" << endl;
   cout << "input >> ";
   cin >> num;

   int id;
   id = atoi(num.c_str());
   char filename[64], dummy[32];
   sprintf(filename,"./ana/det%d.txt",id);
   ifstream file(filename);
   for (int i = 0; i < 16; i++)
   {
      file >> dummy >>detch[i];
      cout << dummy << " " << detch[i] << endl;
   } 
}

void ana112::Set_exist_det(int id){
   //return;
   char filename[64], dummy[32];
   sprintf(filename,"./ana/det%d.txt",id);
   ifstream file(filename);
   for (int i = 0; i < 16; i++)
   {
      file >> dummy >>detch[i];
      cout << dummy << " " << detch[i] << endl;
   } 
}

void ana112::Set_cut_type(){
   //return;
   string num;
   cout << endl;
   cout << "0:  no cut" << endl;
   cout << "1:  all tdc != 0" << endl;
   cout << "2:  all adc has good value" << endl;
   cout << "input >> ";
   cin >> num;
   
   hit_frag = atoi(num.c_str());
}

void ana112::Save_resolution(){
   if(ana_frag==true){return;}
   //return;
   bool frag = true;
   for(int i=-1; i<256; i++){
      if((reso_value[i][0]==-1 || reso_value[i][0]==(double)runno) && frag){
         reso_value[i][0] = runno - 100000000;
         reso_value[i][1] = resolution[2]; //ref1
         reso_value[i][2] = resolution[4]; //ref2
         reso_value[i][3] = resolution[6]; //cdh
         reso_value[i][4] = resolution[8]; //cnc1
         reso_value[i][5] = resolution[10]; //cnc2
         reso_value[i][6] = resolution[8]; //test
         reso_value[i][7] = resolution[13]; //mppc
         frag = false;
      } else {
         for(int j=0; j<8; j++){
            reso_value[i][j] = -1;
         }
      }
   }
   
   ofstream reso_file2(reso_filename2);
   reso_file2 << fixed << setprecision(2);
   for (int i=-1; i<256; i++) {
      if(i==-1){
         reso_file2 << "runno  ref1  ref2  cdh   cnc1  cnc2  test  mppc" << endl;
      } else {
         for(int j=0; j<8; j++){
            reso_file2 << reso_value[i][j] << " ";
         }
         reso_file2 << endl;
      }
   }
   reso_file2.close();
   cout << "resolution save complete !!!" << endl;
}

void ana112::Ana_gain_cnc(int id){
   return;
}