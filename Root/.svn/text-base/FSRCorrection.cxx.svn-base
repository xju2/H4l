#include <H4l/FSRCorrection.h>

#include "FsrUtils/FsrPhotonTool.h"

using namespace std;

FSRCorrection::FSRCorrection(quadruplet *in_quad, bool debug) {
    m_quad = in_quad;
    m_debug = debug;

    m_quad_fsr = new quadruplet;
}

FSRCorrection::~FSRCorrection() {
    delete m_quad_fsr;
}

void FSRCorrection::applyFSR(FSR::FsrPhotonTool *fsrTool, const xAOD::ElectronContainer *in_el, const xAOD::PhotonContainer *in_ph) {
    FSR::FsrCandidate cand1, cand2, cand3, cand4;
    vector<FSR::FsrCandidate> cands;
    if (m_quad->type == _4mu) {
        if ((static_cast<quadruplet_4mu*>(m_quad))->m1->muonType() != xAOD::Muon::MuonType::CaloTagged || 
                (static_cast<quadruplet_4mu*>(m_quad))->m1->muonType() != xAOD::Muon::MuonType::MuonStandAlone) 
            fsrTool->getFsrPhoton((static_cast<quadruplet_4mu*>(m_quad))->m1, cand1, in_ph, in_el);

        if ((static_cast<quadruplet_4mu*>(m_quad))->m2->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_4mu*>(m_quad))->m2->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_4mu*>(m_quad))->m2, cand2, in_ph, in_el);
        if ((static_cast<quadruplet_4mu*>(m_quad))->m3->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_4mu*>(m_quad))->m3->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_4mu*>(m_quad))->m3, cand3, in_ph, in_el);
        if ((static_cast<quadruplet_4mu*>(m_quad))->m4->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_4mu*>(m_quad))->m4->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_4mu*>(m_quad))->m4, cand4, in_ph, in_el);
    }
    else if (m_quad->type == _2mu2e) {
        if ((static_cast<quadruplet_2mu2e*>(m_quad))->m1->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_2mu2e*>(m_quad))->m1->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_2mu2e*>(m_quad))->m1, cand1, in_ph, in_el);
        if ((static_cast<quadruplet_2mu2e*>(m_quad))->m2->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_2mu2e*>(m_quad))->m2->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_2mu2e*>(m_quad))->m2, cand2, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_2mu2e*>(m_quad))->m3, cand3, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_2mu2e*>(m_quad))->m4, cand4, in_ph, in_el);
    }
    else if (m_quad->type == _2e2mu) {
        fsrTool->getFsrPhoton((static_cast<quadruplet_2e2mu*>(m_quad))->m1, cand1, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_2e2mu*>(m_quad))->m2, cand2, in_ph, in_el);
        if ((static_cast<quadruplet_2e2mu*>(m_quad))->m3->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_2e2mu*>(m_quad))->m3->muonType() != xAOD::Muon::MuonType::MuonStandAlone) 
            fsrTool->getFsrPhoton((static_cast<quadruplet_2e2mu*>(m_quad))->m3, cand3, in_ph, in_el);
        if ((static_cast<quadruplet_2e2mu*>(m_quad))->m4->muonType() != xAOD::Muon::MuonType::CaloTagged || (static_cast<quadruplet_2e2mu*>(m_quad))->m4->muonType() != xAOD::Muon::MuonType::MuonStandAlone) fsrTool->getFsrPhoton((static_cast<quadruplet_2e2mu*>(m_quad))->m4, cand4, in_ph, in_el);
    }
    else {
        fsrTool->getFsrPhoton((static_cast<quadruplet_4e*>(m_quad))->m1, cand1, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_4e*>(m_quad))->m2, cand2, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_4e*>(m_quad))->m3, cand3, in_ph, in_el);
        fsrTool->getFsrPhoton((static_cast<quadruplet_4e*>(m_quad))->m4, cand4, in_ph, in_el);
    }
    cands.push_back(cand1);
    cands.push_back(cand2);
    cands.push_back(cand3);
    cands.push_back(cand4);

    // Collinear FSR
    if (m_quad->type == _4mu || m_quad->type == _2mu2e) {
        FSR::FsrCandidate cand;
        int icand = -1;
        float high_et = 0;
        for (int i = 0; i < 2; i++) {
            if (cands[i].Et != -1 && cands[i].Et > high_et) {
                icand = i;
                high_et = cands[i].Et;
            }
        }
        if (icand != -1) {
            cand = cands[icand];
            TLorentzVector z = m_quad->z1;
            TLorentzVector c = cand.particle->p4();
            TLorentzVector zc = z + c;
            if ((z.M() > 66000. && z.M() < 89000.) && zc.M() < 100000.) {
                m_quad_fsr->type = m_quad->type;
                m_quad_fsr->l1 = m_quad->l1;
                m_quad_fsr->l2 = m_quad->l2;
                m_quad_fsr->l3 = m_quad->l3;
                m_quad_fsr->l4 = m_quad->l4;
                m_quad_fsr->z1 = zc;
                m_quad_fsr->z2 = m_quad->z2;
                m_quad_fsr->q = m_quad_fsr->z1 + m_quad_fsr->z2;
                return;
            }
        }
    }

    // Far FSR
    FSR::FsrCandidate cand;
    int icand = -1;
    float high_et = 0;
    for (int i = 0; i < 2; i++) {
        if (cands[i].Et != -1 && cands[i].Et > high_et) {
            icand = i;
            high_et = cands[i].Et;
        }
    }
    if (icand == -1) {m_quad_fsr = m_quad; return;}
    cand = cands[icand];
    TLorentzVector z = (icand < 2) ? m_quad->z1 : m_quad->z2;
    TLorentzVector c = cand.particle->p4();
    TLorentzVector zc = z + c;
    if (!(c.DeltaR(m_quad->l1) > 0.15 && c.DeltaR(m_quad->l2) > 0.15 && c.DeltaR(m_quad->l3) > 0.15 && c.DeltaR(m_quad->l4) > 0.15)) {m_quad_fsr = m_quad; return;}
    if (z.M() < 81000. && zc.M() < 100000.) {
        m_quad_fsr->type = m_quad->type;
        m_quad_fsr->l1 = m_quad->l1;
        m_quad_fsr->l2 = m_quad->l2;
        m_quad_fsr->l3 = m_quad->l3;
        m_quad_fsr->l4 = m_quad->l4;
        if (icand < 2) {
            m_quad_fsr->z1 = zc;
            m_quad_fsr->z2 = m_quad->z2;
        }
        else {
            m_quad_fsr->z1 = m_quad->z1;
            m_quad_fsr->z2 = zc;
        }
        m_quad_fsr->q = m_quad_fsr->z1 + m_quad_fsr->z2;
    }

    m_quad_fsr = m_quad;
    return;
}
