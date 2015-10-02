#include <H4l/JetSelection.h>

using namespace std;

JetSelection::JetSelection(bool debug) {
    m_debug = debug;
    m_pass = false;
}

JetSelection::~JetSelection(){

}

int JetSelection::applyCuts(const xAOD::Jet* m_j, uint32_t runNumber) 
{
    m_pass = false;
    float pt = (float) m_j->pt();
    float eta = (float) m_j->jetP4(xAOD::JetEMScaleMomentum).eta();

    if (m_debug) {
        cout << "Applying cuts to jet." << endl;
        cout << "pt:\t\t\t" << pt << endl;
        cout << "eta:\t\t\t" << eta << endl;
        cout << "Run Number:\t\t\t" << runNumber << endl;
    }

    // pt cut
    if (!(pt > 25000.0)) {
        if (m_debug) cout << "Jet failed pt cut." << endl;
        return 3;
    }

    // eta cut
    if (!((pt <= 30000.0 && fabs(eta) < 2.4) || (pt > 30000.0 && fabs(eta) < 4.5))) {
        if (m_debug) cout << "Jet failed eta cut." << endl;
        return 4;
    }

    // Pile Up Removal cut
    if (pt < 50000.0 && fabs(eta) < 2.4) {
        vector<float> jvfVal;
        m_j->getAttribute(xAOD::JetAttribute::JVF, jvfVal);
        if (!(fabs(jvfVal[0]) > 0.50)) {
            if (m_debug) cout << "Jet failed Pile Up Removal cut." << endl;
            return 5;
        }
    }

    // Jet Cleaning cut
    bool b_BadLooseMinus = false, b_HotTileCell = false;
    Bool_t isBadLoose = 0;
    m_j->getAttribute(xAOD::JetAttribute::isBadLoose, isBadLoose);
    if (isBadLoose != 0) b_BadLooseMinus = true;
    float fmax = 0, smax = 0;
    m_j->getAttribute(xAOD::JetAttribute::FracSamplingMax, fmax);
    m_j->getAttribute(xAOD::JetAttribute::SamplingMax, smax);
    float jc_eta = (float) m_j->eta();
    float jc_phi = (float) m_j->phi();
    bool etaphi28 = (jc_eta > -0.2 && jc_eta < -0.1 && jc_phi > 2.65 && jc_phi < 2.75);
    bool badrun = (runNumber == 202660 || runNumber == 202668 || runNumber == 202712 || runNumber == 202740 || runNumber == 202965 || runNumber == 202987 || runNumber == 202991 || runNumber == 203027);
    if (fmax > 0.6 && smax == 13 && etaphi28 && badrun) b_HotTileCell = true;
    if (b_BadLooseMinus || b_HotTileCell) {
        if (m_debug) cout << "Jet failed Jet Cleaning cut." << endl;
        return 6;
    }

    // jet passed all cuts
    if (m_debug) cout << "Jet passed all cuts." << endl;
    m_pass = true;

    return 0;
}
