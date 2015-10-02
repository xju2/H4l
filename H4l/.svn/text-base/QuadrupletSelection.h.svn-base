#ifndef H4l_QuadrupletSelection_H
#define H4l_QuadrupletSelection_H

#include <math.h>
#include <algorithm>

#include "TLorentzVector.h"

#include "xAODRootAccess/Init.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "MyXAODTools/CPToolsHelper.h"

const float mZ = 91187.6;

enum quadType {
    _4mu,
    _4e,
    _2mu2e,
    _2e2mu
};

typedef enum {
    Electron = 0,
    Muon_combined,
    Muon_CALO,
    Muon_SA,
    Muon_ST,
    NONOTYPE
} ParticleType;


typedef struct {
    TLorentzVector l1, l2, l3, l4, z1, z2, q;
    int type;
    float met;
    float met_final;
    float met_phi;
} Quadruplet ;

struct quadruplet {
    TLorentzVector l1, l2, l3, l4, z1, z2, q;
    int type;
    float met;
    float met_final;
    float met_phi;
    float mpx;
    float mpy;
};

struct quadruplet_4mu : quadruplet {
    const xAOD::Muon *m1, *m2, *m3, *m4;
};

struct quadruplet_2e2mu : quadruplet {
    const xAOD::Electron *m1, *m2;
    const xAOD::Muon *m3, *m4;
};

struct quadruplet_2mu2e : quadruplet {
    const xAOD::Muon *m1, *m2;
    const xAOD::Electron *m3, *m4;
};

struct quadruplet_4e : quadruplet {
    const xAOD::Electron *m1, *m2, *m3, *m4;
};

class QuadrupletSelection {
    public:
        QuadrupletSelection(xAOD::ElectronContainer *in_el, xAOD::MuonContainer *in_mu, bool debug = false);
        ~QuadrupletSelection();

        std::vector<int> applyCuts_Inclusive();
        std::vector<int> applyCuts_Type();
        void fillVecs(quadruplet *quad, std::vector<const xAOD::Electron*> &els, std::vector<const xAOD::Muon*> &mus);
        bool pass() {return m_pass;};
        quadruplet* getQuadruplet() {return m_quad;};
        void doCR(); 
        void SetCPTool(CPToolsHelper* cp) { cp_tools_ = cp;}
        void SetEventInfo(const xAOD::EventInfo* _ei){ ei = _ei; }

    private:
        int applyCuts(quadruplet *quad);
        
        bool doCR_;
        bool m_pass; //!
        bool m_debug; //!
        xAOD::ElectronContainer *m_el; //!
        xAOD::MuonContainer *m_mu; //!
        quadruplet *m_quad; //!
        quadruplet_4mu *m_quad_4mu; //!
        quadruplet_2e2mu *m_quad_2e2mu; //!
        quadruplet_2mu2e *m_quad_2mu2e; //!
        quadruplet_4e *m_quad_4e; //!
        std::vector<quadruplet*> m_quads, m_quads_kine, m_quads_trig, m_quads_z1; //!
        std::vector<quadruplet_4mu*> m_quads_4mu, m_quads_4mu_kine, m_quads_4mu_trig, m_quads_4mu_z1; //!
        std::vector<quadruplet_2e2mu*> m_quads_2e2mu, m_quads_2e2mu_kine, m_quads_2e2mu_trig, m_quads_2e2mu_z1; //!
        std::vector<quadruplet_2mu2e*> m_quads_2mu2e, m_quads_2mu2e_kine, m_quads_2mu2e_trig, m_quads_2mu2e_z1; //!
        std::vector<quadruplet_4e*> m_quads_4e, m_quads_4e_kine, m_quads_4e_trig, m_quads_4e_z1; //!

        CPToolsHelper* cp_tools_;
        const xAOD::EventInfo*  ei;
};

#endif
