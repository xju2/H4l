#ifndef H4l_OverlapRemoval_H
#define H4l_OverlapRemoval_H

#include <math.h>

#include "TMath.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TStore.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"

static SG::AuxElement::Decorator<char> dec_passOR("passOR");
static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");

namespace OverlapRemoval {
    void applyOverlapRemoval(xAOD::ElectronContainer *m_el_pass, 
            xAOD::MuonContainer *m_mu_pass, 
            xAOD::JetContainer *m_j_pass, bool debug = false);
}

#endif
