#ifndef H4l_Minitree_H
#define H4l_Minitree_H

#ifndef H4l_QuadrupletSelection_H
#include <H4l/QuadrupletSelection.h>
#else
struct quadruplet;
#endif

#include "TFile.h"
#include "TTree.h"

class Minitree {
    public:
        Minitree(bool debug = false);
        ~Minitree();
        
        void fillTree(int runNumber, int eventNumber, quadruplet *quad, quadruplet *quad_fsr);
        void setDirectory(TFile *outFile);
        void cleanTree();

    public:
        void makeBranches(TTree *tree);

        bool m_debug; //!

        // trees
        TTree *m_tree_incl_all; //!

        // branches
        int run; //!
        int event; //!
        int mc_channel_number_; //!
        int event_type; //!
        float m4l_unconstrained; //!
        float mZ1_unconstrained; //!
        float mZ2_unconstrained; //!
        float m4l_fsr; //!
        float mZ1_fsr; //!
        float mZ2_fsr; //!
        float m4l_constrained; //!
        float mZ1_constrained; //!
        float mZ2_constrained; //!
        float weight_corr; //!
        float weight_lumi; //!
        float weight_sampleoverlap; //!
        float weight; //!
        float Z1_lepplus_pt; //!
        float Z1_lepminus_pt; //!
        float Z2_lepplus_pt; //!
        float Z2_lepminus_pt; //!
        float Z1_lepplus_eta; //!
        float Z1_lepminus_eta; //!
        float Z2_lepplus_eta; //!
        float Z2_lepminus_eta; //!
        float Z1_lepplus_phi; //!
        float Z1_lepminus_phi; //!
        float Z2_lepplus_phi; //!
        float Z2_lepminus_phi; //!
        float Z1_lepplus_m; //!
        float Z1_lepminus_m; //!
        float Z2_lepplus_m; //!
        float Z2_lepminus_m; //!
        int n_jets; //!
        int n_good_jets; //!
        float dijet_invmass; //!
        float leading_jet_pt; //!
        float subleading_jet_pt; //!
        float BDT_discriminant; //!
        float BDT_discriminant_VBF; //!
        float BDT_discriminant_HadVH; //!
        float MET; //!
        float MET_final_; //!
        float Hpt_; //!
        float MET_phi_; //!
        float dphi_met_higgs_; //!
        float MET_noHiggs_;

};

#endif
