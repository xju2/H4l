#ifndef H4l_FSRCORRECTION_H
#define H4l_FSRCORRECTION_H

#ifndef H4l_QuadrupletSelection_H
#include <H4l/QuadrupletSelection.h>
#else
struct quadruplet;
#endif

#include <math.h>

#include "TLorentzVector.h"

#include "xAODRootAccess/Init.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/PhotonContainer.h"

#include "FsrUtils/FsrPhotonTool.h"

class FSRCorrection {
    public:
        FSRCorrection(quadruplet *in_quad, bool debug = false);
        ~FSRCorrection();

        void applyFSR(FSR::FsrPhotonTool *fsrTool, const xAOD::ElectronContainer *in_el, const xAOD::PhotonContainer *in_ph);
        quadruplet* getQuadruplet() {return m_quad_fsr;};

    private:
        bool m_debug; //!
        quadruplet *m_quad; //!
        quadruplet *m_quad_fsr; //!
};

#endif
