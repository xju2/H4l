#ifndef H4l_JetSelection_H
#define H4l_JetSelection_H

#include <math.h>

#include "xAODRootAccess/Init.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/Jet.h"

class JetSelection {
    public:
        JetSelection(bool debug = false);
        ~JetSelection();

        int applyCuts(const xAOD::Jet* in_j, uint32_t runNumber = 0);
        bool pass() {return m_pass;};

    private:
        bool m_pass; //!
        bool m_debug; //!
};

#endif
