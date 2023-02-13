//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan  6 16:22:50 2023 by ROOT version 6.26/06
// from TTree tree/tree
// found on file: ./data/run132243516.root
//////////////////////////////////////////////////////////

#ifndef ana112_h
#define ana112_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

// Header file for the classes stored in the TTree if any.

class ana112 {
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

   //grobal
   double runno;
   double timeConv;
   char   work_filename[64];
   char   work_dirname[64];
   char   runnum_filename[64];
   char   offset_filename[64];
   char   data_filename[64];
   char   pede_filename[64];
   char   slew_pass_filename[64];
   char   reso_filename[64];
   char   reso_filename2[64];
   double pedestal[16];
   double ped_sigma[16];
   double mip[16];
   double t_offset[16];
   double t_width[16];
   double charge[16];
   double time[16];
   double slew_par0[16];
   double slew_par1[16];
   double slew_par2[16];
   int    slew_no[16];
   double resolution[16];
   int    detch[16];
   int    hit_frag;
   double reso_value[256][8];
   bool   ana_frag;

   //hist
   TH1D *h_adc[16], *h_adcmip[16], *h_adcped[16], *h_tdc[16];
   TH2D *h_slew[4];
   TH1D *h_reso[3];
   TH2D *h_slew_before[4], *h_slew_after[4];
   TH2D *h_ref_swel[4];
   TH2D *h_ref_before[4], *h_ref_after[4];
   TH1D *h_ref_reso;

   //function
   TF1 *func;



   ana112(double id = -1);
   virtual ~ana112();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   virtual void     Ped_mip();  //get pedestal, mip, tdc mean
   virtual void     Slewing(int id, int num);  //slewing
   virtual void     Phtw2(int ref_r,int ref_l, int cor_r,int cor_l); 
   virtual void     Phtw(int id);
   virtual void     Phtw_ref();
   virtual void     Resolution(int id);  //check resolution
   virtual void     Reso_ref();
   virtual void     Check_slew(int id);  //check slewing
   virtual void     Check_ref_slew();
   virtual void     Convert();     //convert adc channel
   virtual bool     Hit(int id);         //define hit
   virtual void     Hit(std::vector<int> vec);         //define hit
   virtual void     Set_exist_det();  //define detector
   virtual void     Set_exist_det(int id);  //define detector
   virtual void     Set_cut_type();        //set cut type
   virtual void     Save_resolution();
   virtual void     Load_file(int id);
   virtual void     Ana_gain_cnc(int id);
};

#endif

#ifdef ana112_cxx
ana112::ana112(double id = 26) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   int j = 0;
   j = (int)id;
   runno = id + 100000000;
   if(id != -1){
      sprintf(work_dirname,"./ana/%d",j);
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
   } else {
      cout << "Analysis mode" << endl;
      sprintf(work_dirname,"./ana/ana");
      cout << "work at : " << work_dirname << endl;
      ana_frag = true;
   }
}

ana112::~ana112()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ana112::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ana112::LoadTree(Long64_t entry)
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

void ana112::Init(TTree *tree)
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
   fChain->SetBranchAddress("tdc", tdc, &b_tdc);
   Notify();

   //load ped_mip
   sprintf(pede_filename, "%s/%s.ped_mip", work_dirname, work_filename);
   ifstream pede_file(pede_filename);
   if (pede_file) {
     for (int i=0; i<16; i++) {
       pede_file >> pedestal[i] >> ped_sigma[i] >> mip[i];
     }
     cout << "load pedestal,mip file" << endl;
   } else {
     for (int i=0; i<16; i++) {
       pedestal[i]  = 0.0;
       ped_sigma[i] = 0.0;
       mip[i] = 0.0;
     }
     cout << "cant load pedestal,mip file" << endl;
   }
   pede_file.close();

   //load tdc offset 
   sprintf(offset_filename, "%s/%s.t_offset", work_dirname, work_filename);
   ifstream offset_file(offset_filename);
   if (offset_file) {
     for (int i=0; i<16; i++) {
       offset_file >> t_offset[i] >> t_width[i];
     }
     cout << "load tdc offdet file" << endl;
   } else {
     for (int i=0; i<16; i++) {
       t_offset[i] = 0.0;
       t_width[i] = 0.0;
     }
     cout << "cant load tdc offdet file" << endl;
   }
   offset_file.close();

   //tdc time(ps)/ch
   ifstream time_file("./ana/timeConv.dat");
   if(time_file){
      time_file >> timeConv;
      cout << "load timeConv file" << endl;
   } else {
      timeConv = 35.;
      cout << "cant load timeConv file" << endl;
   }
   time_file.close();

   //slew pass
   sprintf(slew_pass_filename, "%s/%s.slew_pass", work_dirname, work_filename);
   ifstream slew_pass_file(slew_pass_filename);
   if(slew_pass_file){
      for(int i=0; i<16; i++){
         slew_pass_file >> slew_par0[i] >> slew_par1[i] >> slew_par2[i] >> slew_no[i];
      }
      cout << "load slew pass file" << endl;
   } else {
      for (int i=0; i<16; i++) 
      {
         slew_par0[i] = 0.0; 
         slew_par1[i] = 0.0;
         slew_par2[i] = 0.0; 
         slew_no[i] = 0;
      }
      cout << "cant load slew pass file" << endl;
   }
   slew_pass_file.close();

   //resolution
   sprintf(reso_filename, "%s/%s.reso", work_dirname, work_filename);
   ifstream reso_file(reso_filename);
   int buf;
   if(reso_file){
      for(int i=0; i<16; i++){
         reso_file >> buf >>resolution[i];
      }
      cout << "load reso file" << endl;
   } else {
      for (int i=0; i<16; i++) {
         resolution[i] = 0.;
      }
      cout << "cant load reso file" << endl;
   }
   reso_file.close();
   
   //detch
   //Set_exist_det();
   int num = (int)(runno - 100000000);
   if(num==26 || num ==27){Set_exist_det(0);}
   else if(num==30 || num==32){Set_exist_det(1);}
   else if(num==35 || num==37){Set_exist_det(3);}
   else if(num>50 && num<94){Set_exist_det(5);}
   else if(num>94){Set_exist_det(7);}

   //hit_frag
   hit_frag = 1;
   //Set_cut_type();

   //resolution save
   sprintf(reso_filename2,"./resolution.txt");
   ifstream reso_file2(reso_filename2);
   char buf2[256];
   if(reso_file2){
      for(int i=-1; i<256; i++){
         for(int j=0; j<8; j++){
            if(i==-1){
               reso_file2 >> buf2;
            } else {
               reso_file2 >> reso_value[i][j];
            }
         }
      }
      cout << "load reso2 file" << endl;
   } else {
      for(int i=0; i<256; i++){
         for(int j=0; j<8; j++){
            reso_value[i][j] = -1;
         }
      }
      cout << "cant load reso2 file" << endl;
   }
}

Bool_t ana112::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ana112::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ana112::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ana112_cxx
