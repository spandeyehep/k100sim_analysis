//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 16 18:22:34 2018 by ROOT version 6.06/01
// from TTree hits/HGC rechits
// found on file: muon_v10.root
//////////////////////////////////////////////////////////

#ifndef NtupleVariables_h
#define NtupleVariables_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class NtupleVariables {
public :

   NtupleVariables(TTree * /*tree*/ =0) : fChain(0) { }
   ~NtupleVariables() { }
   //void    Init(TTree *tree);
   void    Init(TTree *tree);
   Bool_t  Notify();
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
  
   bool IsUnique(vector<long> tracks, long thisTrack);  
   vector<double> sort(vector<double> vec);
   // Fixed size dimensions of array or collections stored in the TTree if any.

   Long64_t           EV;
   vector<int>     *DT;
   vector<long>    *TS;
   vector<long>    *P;
   vector<int>     *Type;
   vector<float>   *E1;
   vector<float>   *D3;
   vector<float>   *X1;
   vector<float>   *Y1;
   vector<float>   *Z1;
   vector<float>   *X3;
   vector<float>   *Y3;
   vector<float>   *Z3;
   vector<float>   *PX1;
   vector<float>   *PY1;
   vector<float>   *PZ1;
   vector<float>   *PX3;
   vector<float>   *PY3;
   vector<float>   *PZ3;
   vector<float>   *time1;
   vector<float>   *time3;
   vector<int>     *nCap;

   // List of branches
   TBranch        *b_EV;   //!
   TBranch        *b_DT;   //!
   TBranch        *b_TS;   //!
   TBranch        *b_P;   //!
   TBranch        *b_Type;   //!
   TBranch        *b_E1;   //!
   TBranch        *b_D3;   //!
   TBranch        *b_X1;   //!
   TBranch        *b_Y1;   //!
   TBranch        *b_Z1;   //!
   TBranch        *b_X3;   //!
   TBranch        *b_Y3;   //!
   TBranch        *b_Z3;   //!
   TBranch        *b_PX1;   //!
   TBranch        *b_PY1;   //!
   TBranch        *b_PZ1;   //!
   TBranch        *b_PX3;   //!
   TBranch        *b_PY3;   //!
   TBranch        *b_PZ3;   //!
   TBranch        *b_time1;   //!
   TBranch        *b_time3;   //!
   TBranch        *b_nCap;   //!



};

#endif

#ifdef NtupleVariables_cxx

void NtupleVariables::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //DT->clear();
   
   TS = 0;
   DT = 0;
   P = 0;
   Type = 0;
   E1 = 0;
   D3 = 0;
   X1 = 0;
   Y1 = 0;
   Z1 = 0;
   X3 = 0;
   Y3 = 0;
   Z3 = 0;
   PX1 = 0;
   PY1 = 0;
   PZ1 = 0;
   PX3 = 0;
   PY3 = 0;
   PZ3 = 0;
   time1 = 0;
   time3 = 0;
   nCap = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EV", &EV, &b_EV);
   fChain->SetBranchAddress("DT", &DT, &b_DT);
   fChain->SetBranchAddress("TS", &TS, &b_TS);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("Type", &Type, &b_Type);
   fChain->SetBranchAddress("E1", &E1, &b_E1);
   fChain->SetBranchAddress("D3", &D3, &b_D3);
   fChain->SetBranchAddress("X1", &X1, &b_X1);
   fChain->SetBranchAddress("Y1", &Y1, &b_Y1);
   fChain->SetBranchAddress("Z1", &Z1, &b_Z1);
   fChain->SetBranchAddress("X3", &X3, &b_X3);
   fChain->SetBranchAddress("Y3", &Y3, &b_Y3);
   fChain->SetBranchAddress("Z3", &Z3, &b_Z3);
   fChain->SetBranchAddress("PX1", &PX1, &b_PX1);
   fChain->SetBranchAddress("PY1", &PY1, &b_PY1);
   fChain->SetBranchAddress("PZ1", &PZ1, &b_PZ1);
   fChain->SetBranchAddress("PX3", &PX3, &b_PX3);
   fChain->SetBranchAddress("PY3", &PY3, &b_PY3);
   fChain->SetBranchAddress("PZ3", &PZ3, &b_PZ3);
   fChain->SetBranchAddress("time1", &time1, &b_time1);
   fChain->SetBranchAddress("time3", &time3, &b_time3);
   fChain->SetBranchAddress("nCap", &nCap, &b_nCap);
   Notify();

   
}

Bool_t NtupleVariables::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef NtupleVariables_cxx
