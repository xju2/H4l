#include <H4l/ElectronSelection.h>
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

using namespace std;

ElectronSelection::ElectronSelection(bool debug)
{
    m_debug = debug;
    m_pass = false;

    // Electron likelihood selector
    m_likelihoodCut = new AsgElectronLikelihoodTool("likelihoodCut");
    // Loose config
    m_likelihoodCut->setProperty("primaryVertexContainer", "PrimaryVertices");
    m_likelihoodCut->setProperty("ConfigFile", "ElectronPhotonSelectorTools/offline/mc15_20150429/ElectronLikelihoodLooseOfflineConfig2015.conf");
    if(!m_likelihoodCut->initialize().isSuccess()){
        cerr << "[Error] Cannnot initialize Electron Likelihood Tool" << endl;
    } else {
        cout << "Electron Likelihood tool is initialized" << endl;
    }
}

ElectronSelection::~ElectronSelection()
{
    delete m_likelihoodCut;
}

int ElectronSelection::applyCuts(const xAOD::Electron* m_el, const xAOD::Vertex* m_pvx)
{
    m_pass = false;
    float et = (float) m_el->caloCluster()->e() / (float) cosh(m_el->trackParticle()->eta());
    float eta = (float) m_el->caloCluster()->eta();
    int oq = (int) m_el->auxdata<uint32_t>("OQ") & 1446;
    float z0 = (float) m_el->trackParticle()->z0() + (float) m_el->trackParticle()->vz() - (float) m_pvx->z();
    float z0_sintheta = fabs(z0*sin(m_el->trackParticle()->theta()));
    
    if (m_debug) {
        cout << "Applying cuts to electron." << endl;
        cout << "Et:\t\t\t" << et << endl;
        cout << "eta_cl:\t\t\t" << eta << endl;
        cout << "OQ:\t\t\t" << oq << endl;
        cout << "z0:\t\t\t" << z0 << endl;
    }

    if(m_el->author() != 1 && m_el->author() != 16) {
        return 3;
    }

    // Add nBL Hits cut
    if(el->trackParticleSummaryIntValue(xAOD::numberOfInnermostPixelLayerHits) < 1){
        return 9;
    }

    // Loose Likelihood cut
    if (!(m_likelihoodCut->accept(m_el))) {
        if (m_debug) cout << "Electron failed Loose Likelihood cut." << endl;
        return 4;
    }

    // eta cut
    if (!(fabs(eta) < 2.47)) {
        if (m_debug) cout << "Electron failed eta cut." << endl;
        return 5;
    }
    
    // Et cut
    if (!(et > 7000.0)) {
        if (m_debug) cout << "Electron failed Et cut." << endl;
        return 6;
    }

    // OQ cut
    if (!(oq == 0)) {
        if (m_debug) cout << "Electron failed OQ cut." << endl;
        return 7;
    }

    // z0 cut
    if (!(fabs(z0_sintheta) < .5)) {
        if (m_debug) cout << "Electron failed z0-sintheta cut." << endl;
        return 8;
    }

    // electron passed all cuts
    if (m_debug) cout << "Electron passed all cuts." << endl;
    m_pass = true;

    return 0;
}
