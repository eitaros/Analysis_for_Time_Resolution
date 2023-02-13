//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 19 00:07:44 2023 by ROOT version 6.26/06
// from TTree tree/tree
// found on file: ./data/run100000097.root
//////////////////////////////////////////////////////////

#ifndef aiueo_h
#define aiueo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <ctime>

// Header file for the classes stored in the TTree if any.

class aiueo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           adc[16];
   Int_t           tdc[16];

   // List of branches
   TBranch        *b_adc;   //!
   TBranch        *b_tdc;   //!


   double runno;
   int detch[16];
   string anadir_name;
   int ananum;
   double mean_tdc[16], sigma_tdc[16];
   double tmax[16], tmin[16];
   double mip_mean[16], mip_sigma[16];

   aiueo(int id, int id2);
   aiueo(int id);
   virtual ~aiueo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int id);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual bool     Hit(int id);
   virtual void     Set_exist_det(int id);
   virtual void     Search_tdc_garbage();
   virtual bool     Hit2(int id);
   virtual void     Check();
   virtual void     MipTest();
   virtual void     MipCheck();
};

#endif

#ifdef aiueo_cxx
aiueo::aiueo(int id) : fChain(0)
{
   runno = (double)id + 100000000;
   char fname[128];
   snprintf(fname,sizeof(fname),"./data/run%9.0lf.root",runno);
   TTree *tree = 0;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
      if (!f || !f->IsOpen()) {
         f = new TFile(fname);
      }
      f->GetObject("tree",tree);

   }
   anadir_name = "./anachan/";
   Init(tree);
}

aiueo::aiueo(int id, int id2) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   runno = (double)id + 100000000;
   char fname[128];
   snprintf(fname,sizeof(fname),"./data/run%9.0lf.root",runno);
   TTree *tree = 0;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
      if (!f || !f->IsOpen()) {
         f = new TFile(fname);
      }
      f->GetObject("tree",tree);

   }
   anadir_name = "./aiueo/";
   Init(tree);

   //Loop(id2);
}


aiueo::~aiueo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t aiueo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t aiueo::LoadTree(Long64_t entry)
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

void aiueo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   string ananumf_name;
   ananumf_name = "./ana/cut/ananum.txt";
   ifstream ananumf(ananumf_name.c_str());
   if(ananumf){
      ananumf >> ananum;
   } else {
      ananum = 999;
   }
   ananumf.close();
   int ananum_save;
   ananum_save = ananum + 1;
   ofstream ananumfs(ananumf_name.c_str(), std::ios::trunc);
   if(ananumfs){
      ananumfs << ananum_save;
   }
   ananumfs.close();

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("tdc", tdc, &b_tdc);
   Notify();

   int num = (int)(runno - 100000000);
   if(num==26 || num ==27){Set_exist_det(0);}
   else if(num==30 || num==32){Set_exist_det(1);}
   else if(num==35 || num==37){Set_exist_det(3);}
   else if(num>50 && num<94){Set_exist_det(5);}
   else if(num>94){Set_exist_det(7);}


   // get tdc mean, sigma
   string fname;
   fname = "./ana/cut/tdc_mean.txt";
   int buf;
   ifstream ifs1(fname.c_str());
   if(ifs1){
      ifs1 >> buf;
      for(int i=0; i<16; i++){
         ifs1 >> mean_tdc[i];
      }
      cout << "use " << buf << "ana tdc mean" << endl;
   } else {
      for(int i=0; i<16; i++){
         mean_tdc[i] = -1;
      }
   }
   ifs1.close();
   fname = "./ana/cut/tdc_sigma.txt";
   ifstream ifs2(fname.c_str());
   if(ifs2){
      ifs2 >> buf;
      for(int i=0; i<16; i++){
         ifs2 >> sigma_tdc[i];
      }
      cout << "use " << buf << "ana tdc sigma" << endl;
   } else {
      for(int i=0; i<16; i++){
         sigma_tdc[i] = -1;
      }
   }
   ifs2.close();
}

Bool_t aiueo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void aiueo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t aiueo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef aiueo_cxx
