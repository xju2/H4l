#include <H4l/MuonSelection.h>
#include "MuonSelectorTools/MuonSelectionTool.h"

using namespace std;

MuonSelection::MuonSelection(bool debug) {
    m_debug = debug;
    m_pass = false;

    // Muon selection tool
    m_selectionTool = new CP::MuonSelectionTool("MuonSelection");
    m_selectionTool->setProperty("MaxEta", 2.7);
    m_selectionTool->setProperty("MuQuality",(int) xAOD::Muon::Loose).ignore();
    if(! m_selectionTool->initialize().isSuccess()){
        cerr << "Muon Selector Tool failed initialization" << endl;         
    } else {
        cout << "Muon Selector Tools initialized" << endl;
    }
}

MuonSelection::~MuonSelection(){
    delete m_selectionTool;
}

int MuonSelection::applyCuts(const xAOD::Muon* m_mu, const xAOD::Vertex* m_pvx)
{
    m_pass = false;
    int author = (int) m_mu->author();
    int type = (int) m_mu->muonType();
    float pt = (float) m_mu->pt();
    float eta = (float) m_mu->eta();
    int nPixHits = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfPixelHits");
    int nPixelDeadSensors = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfPixelDeadSensors");
    int nSCTHits = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfSCTHits");
    int nSCTDeadSensors = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfSCTDeadSensors");
    int nPixHoles = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfPixelHoles");
    int nSCTHoles = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfSCTHoles");
    int nTRTHits = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfTRTHits");
    int nTRTOutliers = (type == xAOD::Muon::MuonType::MuonStandAlone) ? 0 : (int) m_mu->primaryTrackParticle()->auxdata<uint8_t>("numberOfTRTOutliers");
    float d0 = (float) m_mu->primaryTrackParticle()->d0();
    float z0 = (float) m_mu->primaryTrackParticle()->z0() + (float) m_mu->primaryTrackParticle()->vz() - (float) m_pvx->z();
    float z0_sintheta = z0* sin(m_mu->primaryTrackParticle()->theta());

    if (m_debug) {
        cout << "Applying cuts to muon." << endl;
        // cout << "author:\t\t\t" << author << endl;
        cout << "type:\t\t\t" << type << endl;
        cout << "pt:\t\t\t" << pt << endl;
        cout << "eta:\t\t\t" << eta << endl;
        cout << "nPixHits:\t\t\t" << nPixHits << endl;
        cout << "nPixelDeadSensors:\t\t\t" << nPixelDeadSensors << endl;
        cout << "nSCTHits:\t\t\t" << nSCTHits << endl;
        cout << "nSCTDeadSensors:\t\t\t" << nSCTDeadSensors << endl;
        cout << "nPixHoles:\t\t\t" << nPixHoles << endl;
        cout << "nSCTHoles:\t\t\t" << nSCTHoles << endl;
        cout << "nTRTHits:\t\t\t" << nTRTHits << endl;
        cout << "nTRTOutliers:\t\t\t" << nTRTOutliers << endl;
        cout << "d0:\t\t\t" << d0 << endl;
        cout << "z0:\t\t\t" << z0 << endl;
    }

    // split by type
    if (type == xAOD::Muon::MuonType::Combined || type == xAOD::Muon::MuonType::SegmentTagged) {
        
        // pt cut
        if (!(pt > 6000)) {
            if (m_debug) cout << "Muon failed pt cut." << endl;
            return 4;
        }

        // eta cut
        if (!(fabs(eta) < 2.7)) {
            if (m_debug) cout << "Muon failed eta cut." << endl;
            return 5;
        }

        // ID cuts w/ MuonSelectionTool
        if (!(m_selectionTool->accept(*m_mu))) {
            if (m_debug) cout << "Muon failed ID cut." << endl;
            return 10;
        }

        // d0 cut
        if (!(fabs(d0) < 1.0 && fabs(z0_sintheta) < 0.5)) {
            if (m_debug) cout << "Muon failed d0 or z0*sinthea cut." << endl;
            return 11;
        }
    }
    else if (type == xAOD::Muon::MuonType::CaloTagged) {
        
        // pt cut
        if (!(pt > 15000)) {
            if (m_debug) cout << "Muon failed pt cut." << endl;
            return 4;
        }

        // eta cut
        //if (!(fabs(eta) < 0.1)) 
        if (!(fabs(eta) < 2.7)) 
        {
            if (m_debug) cout << "Muon failed eta cut." << endl;
            return 5;
        }

        // ID cuts w/ MuonSelectionTool
        if (!(m_selectionTool->accept(*m_mu))) {
            if (m_debug) cout << "Muon failed ID cut." << endl;
            return 10;
        }

        // d0 cut
        if (!(fabs(d0) < 1.0 && fabs(z0_sintheta) < 0.5)) {
            if (m_debug) cout << "Muon failed d0/z0 cut." << endl;
            return 11;
        }
    }
    else if (type == xAOD::Muon::MuonType::MuonStandAlone) {
        // author cut
//        if (!(author == 6)) {
//            if (m_debug) cout << "Muon failed author cut." << endl;
//            return 3;
//        }
        
        // pt cut
        if (!(pt > 6000)) {
            if (m_debug) cout << "Muon failed pt cut." << endl;
            return 4;
        }

        // eta cut
        //if (!((fabs(eta) > 2.5) && (fabs(eta) < 2.7))) 
        if (!(fabs(eta) < 2.7)) 
        {
            if (m_debug) cout << "Muon failed eta cut." << endl;
            return 5;
        }

        // ID cuts w/ MuonSelectionTool
        if (!(m_selectionTool->accept(*m_mu))) {
            if (m_debug) cout << "Muon failed ID cut." << endl;
            return 10;
        }
    }
    else {
        cout << "Muon failed type cut." << endl;
        return 3;
    }

    // muon passed all cuts
    if (m_debug) cout << "Muon passed all cuts." << endl;
    m_pass = true;

    return 0;
}
