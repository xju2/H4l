#ifndef H4l_ElectronSelection_H
#define H4l_ElectronSelection_H

#include <math.h>

#include "xAODRootAccess/Init.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODTracking/VertexContainer.h"

class AsgElectronLikelihoodTool;

class ElectronSelection {
    public:
        ElectronSelection(bool debug = false);
        ~ElectronSelection();

        int applyCuts(const xAOD::Electron* m_el, const xAOD::Vertex* m_pvx);
        bool pass() {return m_pass;};

    private:
        xAOD::Electron *m_el; //!
        const xAOD::Vertex *m_pvx; //!
        bool m_pass; //!
        bool m_debug; //!
        AsgElectronLikelihoodTool *m_likelihoodCut; //!
};

#endif
