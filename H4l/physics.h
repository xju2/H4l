//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  4 21:03:44 2015 by ROOT version 5.34/32
// from TTree tree_incl_all/tree_incl_all
// found on file: root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/monoH4l/minitrees/mc15_13TeV_v1/mc_341748_zphxx_mzp200_mx1_fullsim_mono2.root
//////////////////////////////////////////////////////////

#ifndef __H4L_PHYSICS_H__
#define __H4L_PHYSICS_H__ 

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class physics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           event_type;
   Float_t         m4l_unconstrained;
   Float_t         mZ1_unconstrained;
   Float_t         mZ2_unconstrained;
   Float_t         m4l_fsr;
   Float_t         mZ1_fsr;
   Float_t         mZ2_fsr;
   Float_t         m4l_constrained;
   Float_t         mZ1_constrained;
   Float_t         mZ2_constrained;
   Float_t         weight_corr;
   Float_t         weight_lumi;
   Float_t         weight_sampleoverlap;
   Float_t         weight;
   Float_t         Z1_lepplus_pt;
   Float_t         Z1_lepminus_pt;
   Float_t         Z2_lepplus_pt;
   Float_t         Z2_lepminus_pt;
   Float_t         Z1_lepplus_eta;
   Float_t         Z1_lepminus_eta;
   Float_t         Z2_lepplus_eta;
   Float_t         Z2_lepminus_eta;
   Float_t         Z1_lepplus_phi;
   Float_t         Z1_lepminus_phi;
   Float_t         Z2_lepplus_phi;
   Float_t         Z2_lepminus_phi;
   Float_t         Z1_lepplus_m;
   Float_t         Z1_lepminus_m;
   Float_t         Z2_lepplus_m;
   Float_t         Z2_lepminus_m;
   Int_t           n_jets;
   Int_t           n_good_jets;
   Float_t         dijet_invmass;
   Float_t         leading_jet_pt;
   Float_t         subleading_jet_pt;
   Float_t         BDT_discriminant;
   Float_t         BDT_discriminant_VBF;
   Float_t         BDT_discriminant_HadVH;
   Float_t         MET;
   Float_t         MET_Final;
   Float_t         MET_phi;
   Float_t         Hpt;
   Float_t         dphi_met_higgs;
   Float_t         MET_noHiggs;
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           mc_channel_number;
   Float_t         MCWeight;
   Float_t         actualIPC;
   Float_t         averageIPC;
   UInt_t          bcid;
   UInt_t          lumiblock;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_event_type;   //!
   TBranch        *b_m4l_unconstrained;   //!
   TBranch        *b_mZ1_unconstrained;   //!
   TBranch        *b_mZ2_unconstrained;   //!
   TBranch        *b_m4l_fsr;   //!
   TBranch        *b_mZ1_fsr;   //!
   TBranch        *b_mZ2_fsr;   //!
   TBranch        *b_m4l_constrained;   //!
   TBranch        *b_mZ1_constrained;   //!
   TBranch        *b_mZ2_constrained;   //!
   TBranch        *b_weight_corr;   //!
   TBranch        *b_weight_lumi;   //!
   TBranch        *b_weight_sampleoverlap;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_Z1_lepplus_pt;   //!
   TBranch        *b_Z1_lepminus_pt;   //!
   TBranch        *b_Z2_lepplus_pt;   //!
   TBranch        *b_Z2_lepminus_pt;   //!
   TBranch        *b_Z1_lepplus_eta;   //!
   TBranch        *b_Z1_lepminus_eta;   //!
   TBranch        *b_Z2_lepplus_eta;   //!
   TBranch        *b_Z2_lepminus_eta;   //!
   TBranch        *b_Z1_lepplus_phi;   //!
   TBranch        *b_Z1_lepminus_phi;   //!
   TBranch        *b_Z2_lepplus_phi;   //!
   TBranch        *b_Z2_lepminus_phi;   //!
   TBranch        *b_Z1_lepplus_m;   //!
   TBranch        *b_Z1_lepminus_m;   //!
   TBranch        *b_Z2_lepplus_m;   //!
   TBranch        *b_Z2_lepminus_m;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_n_good_jets;   //!
   TBranch        *b_dijet_invmass;   //!
   TBranch        *b_leading_jet_pt;   //!
   TBranch        *b_subleading_jet_pt;   //!
   TBranch        *b_BDT_discriminant;   //!
   TBranch        *b_BDT_discriminant_VBF;   //!
   TBranch        *b_BDT_discriminant_HadVH;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_Final;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_Hpt;   //!
   TBranch        *b_dphi_met_higgs;   //!
   TBranch        *b_MET_noHiggs;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_mc_channel_number;   //!
   TBranch        *b_MCWeight;   //!
   TBranch        *b_actualIPC;   //!
   TBranch        *b_averageIPC;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_lumiblock;   //!

   physics(TTree *tree=0);
   virtual ~physics();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef h4l_physics_cxx_
physics::physics(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/monoH4l/minitrees/mc15_13TeV_v1/mc_341748_zphxx_mzp200_mx1_fullsim_mono2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/monoH4l/minitrees/mc15_13TeV_v1/mc_341748_zphxx_mzp200_mx1_fullsim_mono2.root");
      }
      f->GetObject("tree_incl_all",tree);

   }
   Init(tree);
}

physics::~physics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t physics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t physics::LoadTree(Long64_t entry)
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

void physics::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("event_type", &event_type, &b_event_type);
   fChain->SetBranchAddress("m4l_unconstrained", &m4l_unconstrained, &b_m4l_unconstrained);
   fChain->SetBranchAddress("mZ1_unconstrained", &mZ1_unconstrained, &b_mZ1_unconstrained);
   fChain->SetBranchAddress("mZ2_unconstrained", &mZ2_unconstrained, &b_mZ2_unconstrained);
   fChain->SetBranchAddress("m4l_fsr", &m4l_fsr, &b_m4l_fsr);
   fChain->SetBranchAddress("mZ1_fsr", &mZ1_fsr, &b_mZ1_fsr);
   fChain->SetBranchAddress("mZ2_fsr", &mZ2_fsr, &b_mZ2_fsr);
   fChain->SetBranchAddress("m4l_constrained", &m4l_constrained, &b_m4l_constrained);
   fChain->SetBranchAddress("mZ1_constrained", &mZ1_constrained, &b_mZ1_constrained);
   fChain->SetBranchAddress("mZ2_constrained", &mZ2_constrained, &b_mZ2_constrained);
   fChain->SetBranchAddress("weight_corr", &weight_corr, &b_weight_corr);
   fChain->SetBranchAddress("weight_lumi", &weight_lumi, &b_weight_lumi);
   fChain->SetBranchAddress("weight_sampleoverlap", &weight_sampleoverlap, &b_weight_sampleoverlap);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("Z1_lepplus_pt", &Z1_lepplus_pt, &b_Z1_lepplus_pt);
   fChain->SetBranchAddress("Z1_lepminus_pt", &Z1_lepminus_pt, &b_Z1_lepminus_pt);
   fChain->SetBranchAddress("Z2_lepplus_pt", &Z2_lepplus_pt, &b_Z2_lepplus_pt);
   fChain->SetBranchAddress("Z2_lepminus_pt", &Z2_lepminus_pt, &b_Z2_lepminus_pt);
   fChain->SetBranchAddress("Z1_lepplus_eta", &Z1_lepplus_eta, &b_Z1_lepplus_eta);
   fChain->SetBranchAddress("Z1_lepminus_eta", &Z1_lepminus_eta, &b_Z1_lepminus_eta);
   fChain->SetBranchAddress("Z2_lepplus_eta", &Z2_lepplus_eta, &b_Z2_lepplus_eta);
   fChain->SetBranchAddress("Z2_lepminus_eta", &Z2_lepminus_eta, &b_Z2_lepminus_eta);
   fChain->SetBranchAddress("Z1_lepplus_phi", &Z1_lepplus_phi, &b_Z1_lepplus_phi);
   fChain->SetBranchAddress("Z1_lepminus_phi", &Z1_lepminus_phi, &b_Z1_lepminus_phi);
   fChain->SetBranchAddress("Z2_lepplus_phi", &Z2_lepplus_phi, &b_Z2_lepplus_phi);
   fChain->SetBranchAddress("Z2_lepminus_phi", &Z2_lepminus_phi, &b_Z2_lepminus_phi);
   fChain->SetBranchAddress("Z1_lepplus_m", &Z1_lepplus_m, &b_Z1_lepplus_m);
   fChain->SetBranchAddress("Z1_lepminus_m", &Z1_lepminus_m, &b_Z1_lepminus_m);
   fChain->SetBranchAddress("Z2_lepplus_m", &Z2_lepplus_m, &b_Z2_lepplus_m);
   fChain->SetBranchAddress("Z2_lepminus_m", &Z2_lepminus_m, &b_Z2_lepminus_m);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("n_good_jets", &n_good_jets, &b_n_good_jets);
   fChain->SetBranchAddress("dijet_invmass", &dijet_invmass, &b_dijet_invmass);
   fChain->SetBranchAddress("leading_jet_pt", &leading_jet_pt, &b_leading_jet_pt);
   fChain->SetBranchAddress("subleading_jet_pt", &subleading_jet_pt, &b_subleading_jet_pt);
   fChain->SetBranchAddress("BDT_discriminant", &BDT_discriminant, &b_BDT_discriminant);
   fChain->SetBranchAddress("BDT_discriminant_VBF", &BDT_discriminant_VBF, &b_BDT_discriminant_VBF);
   fChain->SetBranchAddress("BDT_discriminant_HadVH", &BDT_discriminant_HadVH, &b_BDT_discriminant_HadVH);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_Final", &MET_Final, &b_MET_Final);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("Hpt", &Hpt, &b_Hpt);
   fChain->SetBranchAddress("dphi_met_higgs", &dphi_met_higgs, &b_dphi_met_higgs);
   fChain->SetBranchAddress("MET_noHiggs", &MET_noHiggs, &b_MET_noHiggs);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
   fChain->SetBranchAddress("MCWeight", &MCWeight, &b_MCWeight);
   fChain->SetBranchAddress("actualIPC", &actualIPC, &b_actualIPC);
   fChain->SetBranchAddress("averageIPC", &averageIPC, &b_averageIPC);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   Notify();
}

Bool_t physics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void physics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t physics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef physics_cxx
