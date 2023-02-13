//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 14 10:16:49 2022 by ROOT version 6.24/06
// from TTree tree/tree
// found on file: ana_0229.root
//////////////////////////////////////////////////////////

#ifndef ana_eito_h
#define ana_eito_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define noofdet 15

// Header file for the classes stored in the TTree if any.

class ana_eito {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           adc[16];
   Int_t           leading[16][4];
   Int_t           trailing[16][4];
   Int_t           leading2[16][4];
   Int_t           trailing2[16][4];

   // List of branches
   TBranch        *b_adc;   //!
   TBranch        *b_leading;   //!
   TBranch        *b_trailing;   //!
   TBranch        *b_leading2;   //!
   TBranch        *b_trailing2;   //!


   //grobal
   Int_t runno;
   double timeConv[16];
   //Int_t noofdet;
   Int_t detno[noofdet];
   

   //ped_mip
   double pedestal[noofdet];
   double ped_sigma[noofdet];
   double mip[noofdet];
   double t_offset[noofdet];
   double t_width[noofdet];
   TH1F *h_adcmip[noofdet];
   TH1F *h_adcped[noofdet];
   TH1F *h_rawtdc[noofdet];

   //phtw
   double slew_par0[noofdet];
   double slew_par1[noofdet];
   double slew_par2[noofdet];
   int slew_no[noofdet];
   TH2F *dt[2];
   TF1 *func;
   double rawtdc[noofdet];
   double charge[noofdet];
   double time[noofdet];
   int cor_r;
   int cor_l;
   int ref_r;
   int ref_l;

   //dt_fit
   TH1F *h_dt[3];
   double fit_par[3];
   double fit_error[3];

   //changed_data
   TH1F *changed_adc[noofdet];
   TH1F *changed_tdc[noofdet];

   //check
   TH1F *raw_adc[noofdet];
   TH1F *raw_tdc[noofdet];
   TH1F *hit_adc[noofdet];

   //check_swel
   TH2F *raw_slew_pa[4];
   TH2F *raw_slew_cdh[4];
   TH2F *raw_slew_cnc1[4];
   TH2F *raw_slew_cnc2[4];
   TH2F *raw_slew_ref[4];
   TH2F *check_slew_pa[4];
   TH2F *check_slew_cdh[4];
   TH2F *check_slew_cnc1[4];
   TH2F *check_slew_cnc2[4];
   TH2F *check_slew_ref[4];
   double mean[16];
   double smean[16];
   double raw_time;
   double slewed_time;

   //ref
   TH1F *h_dt_ref;

   ana_eito(Int_t id = 96);
   virtual ~ana_eito();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     ped_mip();
   virtual void     phtw(int cor1,int cor2,int ref1,int ref2);
   virtual void     dt_fit(int test);
   virtual bool     hit(int t);
   virtual void     changed_data();
   virtual void     check();
   //virtual void     clock_test();
   virtual void     check_slew();
   virtual void     check_ref_slew();
   virtual void     Convert();
   virtual void     dt_fit_ref(int test);
   virtual bool     hit2(std::vector<int> vec);
   virtual void     dt_fit2(int t0, int t1, int t2);
};

#endif

#ifdef ana_eito_cxx
ana_eito::ana_eito(Int_t id) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   runno = id;
   TTree* tree;
   char fname[64];
   sprintf(fname,"test_%5.5d.root",runno);
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
   if (!f || !f->IsOpen()) 
   {
      f = new TFile(fname);
   }
   f->GetObject("raw",tree);
   Init(tree);
}

ana_eito::~ana_eito()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ana_eito::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ana_eito::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ana_eito::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("leading", leading, &b_leading);
   fChain->SetBranchAddress("trailing", trailing, &b_trailing);
   fChain->SetBranchAddress("leading2", leading2, &b_leading2);
   fChain->SetBranchAddress("trailing2", trailing2, &b_trailing2);
   Notify();

   //define number of det
   //noofdet = 11;
   //intiarise  
   char pede_file_name[64];   
   sprintf(pede_file_name,  "run%4.4d.ped_mip", runno);
   ifstream pede_file(pede_file_name);
   if (pede_file) {
     for (int i=0; i<noofdet; i++) {
       pede_file >> pedestal[i] >> ped_sigma[i] >> mip[i];
     }
   } else {
     for (int i=0; i<noofdet; i++) {
       pedestal[i]  = 0.0;
       ped_sigma[i] = 0.0;
       mip[i] = 0.0;
     }
   }
   pede_file.close();

   char offset_file_name[64]; 
   sprintf(offset_file_name,"run%4.4d.t_offset",runno);
   ifstream offset_file(offset_file_name);
   if (offset_file) {
     for (int i=0; i<noofdet; i++) {
       offset_file >> t_offset[i] >> t_width[i];
     }
   } else {
     for (int i=0; i<noofdet; i++) {
       t_offset[i] = 0.0;
       t_width[i] = 0.0;
     }
   }
   offset_file.close();

   ifstream time_file("timeConv.dat");
   for (int i=0; i<16; i++) {
     time_file >> timeConv[i];
   }
   time_file.close();

   ifstream config_file("Run.conf");
   char dummy[32];
   int ith = 0;
   for (int i=0; i<noofdet+1; i++) {
     if (i == 0)
     {
      config_file >> dummy >>dummy;
      continue;
     }
     config_file >> dummy >> detno[ith];
     std::cout << detno[ith] << std::endl;
     ith++;
   }
   config_file.close();

   char smname[64];
   sprintf(smname,"run%4.4d.slew_mean",runno);
   ifstream mean_file(smname);
   if (mean_file){
   for (int i=0; i<16; i++) {
     mean_file >> mean[i] >> smean[i];
   }
   }else{
      for (int i = 0; i < 16; i++)
      {
         mean[i] = 0.;
         smean[i] = 0.;
      }
   }
   mean_file.close();
}

Bool_t ana_eito::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ana_eito::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ana_eito::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ana_eito_cxx
