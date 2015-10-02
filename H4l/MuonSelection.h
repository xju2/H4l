#ifndef H4l_MuonSelection_H
#define H4l_MuonSelection_H

#include <math.h>

#include "xAODRootAccess/Init.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTracking/VertexContainer.h"

namespace CP {
    class MuonSelectionTool;
}

class MuonSelection {
    public:
        MuonSelection(bool debug = false);
        ~MuonSelection();

        int applyCuts(const xAOD::Muon* m_mu, const xAOD::Vertex* m_pvx);
        bool pass() {return m_pass;};

    private:
        bool m_pass; //!
        bool m_debug; //!
        CP::MuonSelectionTool* m_selectionTool;
};
#endif
