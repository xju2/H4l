#include <H4l/QuadrupletSelection.h>
#include "xAODTracking/TrackParticlexAODHelpers.h"

using namespace std;

QuadrupletSelection::QuadrupletSelection(xAOD::ElectronContainer *in_el, xAOD::MuonContainer *in_mu, bool debug) {
    m_el = in_el;
    m_mu = in_mu;
    m_debug = debug;
    m_pass = false;
    doCR_ = false;

    m_quad = new quadruplet;
    m_quad_4mu = new quadruplet_4mu;
    m_quad_2e2mu = new quadruplet_2e2mu;
    m_quad_2mu2e = new quadruplet_2mu2e;
    m_quad_4e = new quadruplet_4e;
}

QuadrupletSelection::~QuadrupletSelection(){
    delete m_quad;
    delete m_quad_4mu;
    delete m_quad_2e2mu;
    delete m_quad_2mu2e;
    delete m_quad_4e;
}

vector<int> QuadrupletSelection::applyCuts_Inclusive() {
    vector<int> result;
    for (int i = 0; i < 4; i++) result.push_back(0);

    return result;
}

vector<int> QuadrupletSelection::applyCuts_Type() {
    vector<int> result;
    for (int i = 0; i < 4; i++) result.push_back(0);

    // Selected Leptons cut
    if (!((int) m_mu->size() >= 4)) {
        if (m_debug) cout << "Failed Selected Leptons (4mu) cut." << endl;
        result[_4mu] = 5;
    }
    if (!((int) m_el->size() >= 4)) {
        if (m_debug) cout << "Failed Selected Leptons (4e) cut." << endl;
        result[_4e] = 5;
    }
    if (!((int) m_el->size() >= 2 && (int) m_mu->size() >= 2)) {
        if (m_debug) cout << "Failed Selected Leptons (2e2mu/2mu2e) cut." << endl;
        result[_2e2mu] = 5;
        result[_2mu2e] = 5;
    }

    // SFOS cut (form quadruplets)
    // 4mu
    if (result[_4mu] == 0) {
        for (int imu1 = 0; imu1 < (int) m_mu->size() - 3; imu1++) {
            for (int imu2 = imu1 + 1; imu2 < (int) m_mu->size() - 2; imu2++) {
                for (int imu3 = imu2 + 1; imu3 < (int) m_mu->size() - 1; imu3++) {
                    for (int imu4 = imu3 + 1; imu4 < (int) m_mu->size(); imu4++) {
                        vector<const xAOD::Muon*> muons;
                        muons.push_back((*m_mu)[imu1]);
                        muons.push_back((*m_mu)[imu2]);
                        muons.push_back((*m_mu)[imu3]);
                        muons.push_back((*m_mu)[imu4]);
                        if (!((muons[0]->charge() == -1 * muons[1]->charge() && muons[2]->charge() == -1 * muons[3]->charge()) ||
                              (muons[0]->charge() == -1 * muons[2]->charge() && muons[1]->charge() == -1 * muons[3]->charge()) ||
                              (muons[0]->charge() == -1 * muons[3]->charge() && muons[1]->charge() == -1 * muons[2]->charge()))) continue;
                        int m1_sacalo = (muons[0]->muonType() == xAOD::Muon::MuonStandAlone || muons[0]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                        int m2_sacalo = (muons[1]->muonType() == xAOD::Muon::MuonStandAlone || muons[1]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                        int m3_sacalo = (muons[2]->muonType() == xAOD::Muon::MuonStandAlone || muons[2]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                        int m4_sacalo = (muons[3]->muonType() == xAOD::Muon::MuonStandAlone || muons[3]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                        if (m1_sacalo + m2_sacalo + m3_sacalo + m4_sacalo > 1) continue;
                        float mZ_diff = 100000.;
                        int il1 = 0, il2 = 1;
                        for (int i1 = 0; i1 < 3; i1++) {
                            for (int i2 = i1 + 1; i2 < 4; i2++) {
                                if (!(muons[i1]->charge() == -1 * muons[i2]->charge())) continue;
                                TLorentzVector l1 = muons[i1]->p4();
                                TLorentzVector l2 = muons[i2]->p4();
                                TLorentzVector z1 = l1 + l2;
                                float dM = fabs(z1.M() - mZ);
                                if (dM < mZ_diff) {
                                    il1 = i1;
                                    il2 = i2;
                                    mZ_diff = dM;
                                }
                            }
                        }
                        vector<const xAOD::Muon*> muons_z1, muons_z2;
                        for (int i = 0; i < 4; i++) {
                            if (i == il1 || i == il2) muons_z1.push_back(muons[i]);
                            else muons_z2.push_back(muons[i]);
                        }
                        quadruplet_4mu *quad_4mu = new quadruplet_4mu;
                        quad_4mu->type = _4mu;
                        if (muons_z1[0]->charge() > 0) {
                            quad_4mu->m1 = muons_z1[0];
                            quad_4mu->m2 = muons_z1[1];
                        }
                        else {
                            quad_4mu->m1 = muons_z1[1];
                            quad_4mu->m2 = muons_z1[0];
                        }
                        if (muons_z2[0]->charge() > 0) {
                            quad_4mu->m3 = muons_z2[0];
                            quad_4mu->m4 = muons_z2[1];
                        }
                        else {
                            quad_4mu->m3 = muons_z2[1];
                            quad_4mu->m4 = muons_z2[0];
                        }
                        quad_4mu->l1 = quad_4mu->m1->p4();
                        quad_4mu->l2 = quad_4mu->m2->p4();
                        quad_4mu->l3 = quad_4mu->m3->p4();
                        quad_4mu->l4 = quad_4mu->m4->p4();
                        quad_4mu->z1 = quad_4mu->l1 + quad_4mu->l2;
                        quad_4mu->z2 = quad_4mu->l3 + quad_4mu->l4;
                        quad_4mu->q = quad_4mu->z1 + quad_4mu->z2;
                        m_quads_4mu.push_back(quad_4mu);
                        m_quads.push_back(quad_4mu);
                    }
                }
            }
        }
        if (m_quads_4mu.size() == 0) {
            if (m_debug) cout << "Failed SFOS (4mu) cut." << endl;
            result[_4mu] = 6;
        }
    }
    // 4e
    if (result[_4e] == 0) {
        for (int ie1 = 0; ie1 < (int) m_el->size() - 3; ie1++) {
            for (int ie2 = ie1 + 1; ie2 < (int) m_el->size() - 2; ie2++) {
                for (int ie3 = ie2 + 1; ie3 < (int) m_el->size() - 1; ie3++) {
                    for (int ie4 = ie3 + 1; ie4 < (int) m_el->size(); ie4++) {
                        vector<const xAOD::Electron*> electrons;
                        electrons.push_back((*m_el)[ie1]);
                        electrons.push_back((*m_el)[ie2]);
                        electrons.push_back((*m_el)[ie3]);
                        electrons.push_back((*m_el)[ie4]);
                        if (!((electrons[0]->charge() == -1 * electrons[1]->charge() && electrons[2]->charge() == -1 * electrons[3]->charge()) ||
                              (electrons[0]->charge() == -1 * electrons[2]->charge() && electrons[1]->charge() == -1 * electrons[3]->charge()) ||
                              (electrons[0]->charge() == -1 * electrons[3]->charge() && electrons[1]->charge() == -1 * electrons[2]->charge()))) continue;
                        float mZ_diff = 100000.;
                        int il1 = 0, il2 = 1;
                        for (int i1 = 0; i1 < 3; i1++) {
                            for (int i2 = i1 + 1; i2 < 4; i2++) {
                                if (!(electrons[i1]->charge() == -1 * electrons[i2]->charge())) continue;
                                TLorentzVector l1 = electrons[i1]->p4();
                                TLorentzVector l2 = electrons[i2]->p4();
                                TLorentzVector z1 = l1 + l2;
                                float dM = fabs(z1.M() - mZ);
                                if (dM < mZ_diff) {
                                    il1 = i1;
                                    il2 = i2;
                                    mZ_diff = dM;
                                }
                            }
                        }
                        vector<const xAOD::Electron*> electrons_z1, electrons_z2;
                        for (int i = 0; i < 4; i++) {
                            if (i == il1 || i == il2) electrons_z1.push_back(electrons[i]);
                            else electrons_z2.push_back(electrons[i]);
                        }
                        quadruplet_4e *quad_4e = new quadruplet_4e;
                        quad_4e->type = _4e;
                        if (electrons_z1[0]->charge() > 0) {
                            quad_4e->m1 = electrons_z1[0];
                            quad_4e->m2 = electrons_z1[1];
                        }
                        else {
                            quad_4e->m1 = electrons_z1[1];
                            quad_4e->m2 = electrons_z1[0];
                        }
                        if (electrons_z2[0]->charge() > 0) {
                            quad_4e->m3 = electrons_z2[0];
                            quad_4e->m4 = electrons_z2[1];
                        }
                        else {
                            quad_4e->m3 = electrons_z2[1];
                            quad_4e->m4 = electrons_z2[0];
                        }
                        quad_4e->l1 = quad_4e->m1->p4();
                        quad_4e->l2 = quad_4e->m2->p4();
                        quad_4e->l3 = quad_4e->m3->p4();
                        quad_4e->l4 = quad_4e->m4->p4();
                        quad_4e->z1 = quad_4e->l1 + quad_4e->l2;
                        quad_4e->z2 = quad_4e->l3 + quad_4e->l4;
                        quad_4e->q = quad_4e->z1 + quad_4e->z2;
                        m_quads_4e.push_back(quad_4e);
                        m_quads.push_back(quad_4e);
                    }
                }
            }
        }
        if (m_quads_4e.size() == 0) {
            if (m_debug)  cout << "Failed SFOS (4e) cut." << endl;
            result[_4e] = 6;
        }
    }
    // 2e2mu/2mu2e
    if (result[_2e2mu] == 0) {
        for (int imu1 = 0; imu1 < (int) m_mu->size() - 1; imu1++) {
            for (int imu2 = imu1 + 1; imu2 < (int) m_mu->size(); imu2++) {
                vector<const xAOD::Muon*> muons;
                muons.push_back((*m_mu)[imu1]);
                muons.push_back((*m_mu)[imu2]);
                if (!(muons[0]->charge() == -1 * muons[1]->charge())) continue;
                int m1_sacalo = (muons[0]->muonType() == xAOD::Muon::MuonStandAlone || muons[0]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                int m2_sacalo = (muons[1]->muonType() == xAOD::Muon::MuonStandAlone || muons[1]->muonType() == xAOD::Muon::CaloTagged) ? 1 : 0;
                if (m1_sacalo + m2_sacalo > 1) continue;
                for (int ie1 = 0; ie1 < (int) m_el->size() - 1; ie1++) {
                    for (int ie2 = ie1 + 1; ie2 < (int) m_el->size(); ie2++) {
                        vector<const xAOD::Electron*> electrons;
                        electrons.push_back((*m_el)[ie1]);
                        electrons.push_back((*m_el)[ie2]);
                        if (!(electrons[0]->charge() == -1 * electrons[1]->charge())) continue;
                        TLorentzVector m1 = muons[0]->p4();
                        TLorentzVector m2 = muons[1]->p4();
                        TLorentzVector m = m1 + m2;
                        TLorentzVector e1 = electrons[0]->p4();
                        TLorentzVector e2 = electrons[1]->p4();
                        TLorentzVector e = e1 + e2;
                        float dM_m = fabs(m.M() - mZ);
                        float dM_e = fabs(e.M() - mZ);
                        if (dM_m <= dM_e) {
                            quadruplet_2mu2e *quad_2mu2e = new quadruplet_2mu2e;
                            quad_2mu2e->type = _2mu2e;
                            if (muons[0]->charge() > 0) {
                                quad_2mu2e->m1 = muons[0];
                                quad_2mu2e->m2 = muons[1];
                            }
                            else {
                                quad_2mu2e->m1 = muons[1];
                                quad_2mu2e->m2 = muons[0];
                            }
                            if (electrons[0]->charge() > 0) {
                                quad_2mu2e->m3 = electrons[0];
                                quad_2mu2e->m4 = electrons[1];
                            }
                            else {
                                quad_2mu2e->m3 = electrons[1];
                                quad_2mu2e->m4 = electrons[0];
                            }
                            quad_2mu2e->l1 = quad_2mu2e->m1->p4();
                            quad_2mu2e->l2 = quad_2mu2e->m2->p4();
                            quad_2mu2e->l3 = quad_2mu2e->m3->p4();
                            quad_2mu2e->l4 = quad_2mu2e->m4->p4();
                            quad_2mu2e->z1 = quad_2mu2e->l1 + quad_2mu2e->l2;
                            quad_2mu2e->z2 = quad_2mu2e->l3 + quad_2mu2e->l4;
                            quad_2mu2e->q = quad_2mu2e->z1 + quad_2mu2e->z2;
                            m_quads_2mu2e.push_back(quad_2mu2e);
                            m_quads.push_back(quad_2mu2e);
                        }
                        else {
                            quadruplet_2e2mu *quad_2e2mu = new quadruplet_2e2mu;
                            quad_2e2mu->type = _2e2mu;
                            if (muons[0]->charge() > 0) {
                                quad_2e2mu->m3 = muons[0];
                                quad_2e2mu->m4 = muons[1];
                            }
                            else {
                                quad_2e2mu->m3 = muons[1];
                                quad_2e2mu->m4 = muons[0];
                            }
                            if (electrons[0]->charge() > 0) {
                                quad_2e2mu->m1 = electrons[0];
                                quad_2e2mu->m2 = electrons[1];
                            }
                            else {
                                quad_2e2mu->m1 = electrons[1];
                                quad_2e2mu->m2 = electrons[0];
                            }
                            quad_2e2mu->l1 = quad_2e2mu->m1->p4();
                            quad_2e2mu->l2 = quad_2e2mu->m2->p4();
                            quad_2e2mu->l3 = quad_2e2mu->m3->p4();
                            quad_2e2mu->l4 = quad_2e2mu->m4->p4();
                            quad_2e2mu->z1 = quad_2e2mu->l1 + quad_2e2mu->l2;
                            quad_2e2mu->z2 = quad_2e2mu->l3 + quad_2e2mu->l4;
                            quad_2e2mu->q = quad_2e2mu->z1 + quad_2e2mu->z2;
                            m_quads_2e2mu.push_back(quad_2e2mu);
                            m_quads.push_back(quad_2e2mu);
                        }
                    }
                }
            }
        }
        if (m_quads_2mu2e.size() == 0) {
            if (m_debug) cout << "Failed SFOS (2mu2e) cut." << endl;
            result[_2mu2e] = 6;
        }
        if (m_quads_2e2mu.size() == 0) {
            if (m_debug) cout << "Failed SFOS (2e2mu) cut." << endl;
            result[_2e2mu] = 6;
        }
    }

    // Kinematics cut
    for (quadruplet *quad : m_quads) {
        vector<float> pts;
        pts.push_back(quad->l1.Pt());
        pts.push_back(quad->l2.Pt());
        pts.push_back(quad->l3.Pt());
        pts.push_back(quad->l4.Pt());

        sort(pts.begin(), pts.end());
        if (!(pts[3] > 20000.0 && pts[2] > 15000.0 && pts[1] > 10000.0)) continue;
        m_quads_kine.push_back(quad);
        if (quad->type == _4mu) m_quads_4mu_kine.push_back(static_cast<quadruplet_4mu*>(quad));
        if (quad->type == _2e2mu) m_quads_2e2mu_kine.push_back(static_cast<quadruplet_2e2mu*>(quad));
        if (quad->type == _2mu2e) m_quads_2mu2e_kine.push_back(static_cast<quadruplet_2mu2e*>(quad));
        if (quad->type == _4e) m_quads_4e_kine.push_back(static_cast<quadruplet_4e*>(quad));
    }
    if (m_quads_4mu_kine.size() == 0 && result[_4mu] == 0) {
        if (m_debug) cout << "Failed Kinematics (4mu) cut." << endl;
        result[_4mu] = 7;
    }
    if (m_quads_2e2mu_kine.size() == 0 && result[_2e2mu] == 0) {
        if (m_debug) cout << "Failed Kinematics (2e2mu) cut." << endl;
        result[_2e2mu] = 7;
    }
    if (m_quads_2mu2e_kine.size() == 0 && result[_2mu2e] == 0) {
        if (m_debug) cout << "Failed Kinematics (2mu2e) cut." << endl;
        result[_2mu2e] = 7;
    }
    if (m_quads_4e_kine.size() == 0 && result[_4e] == 0) {
        if (m_debug) cout << "Failed Kinematics (4e) cut." << endl;
        result[_4e] = 7;
    }

    // Trigger Match cut -- not implemented yet
    for (quadruplet *quad : m_quads_kine) {
        // cut to go here
        m_quads_trig.push_back(quad);
        if (quad->type == _4mu) m_quads_4mu_trig.push_back(static_cast<quadruplet_4mu*>(quad));
        if (quad->type == _2e2mu) m_quads_2e2mu_trig.push_back(static_cast<quadruplet_2e2mu*>(quad));
        if (quad->type == _2mu2e) m_quads_2mu2e_trig.push_back(static_cast<quadruplet_2mu2e*>(quad));
        if (quad->type == _4e) m_quads_4e_trig.push_back(static_cast<quadruplet_4e*>(quad));
    }
    if (m_quads_4mu_trig.size() == 0 && result[_4mu] == 0) {
        if (m_debug) cout << "Failed Trigger Match (4mu) cut." << endl;
        result[_4mu] = 8;
    }
    if (m_quads_2e2mu_trig.size() == 0 && result[_2e2mu] == 0) {
        if (m_debug) cout << "Failed Trigger Match (2e2mu) cut." << endl;
        result[_2e2mu] = 8;
    }
    if (m_quads_2mu2e_trig.size() == 0 && result[_2mu2e] == 0) {
        if (m_debug) cout << "Failed Trigger Match (2mu2e) cut." << endl;
        result[_2mu2e] = 8;
    }
    if (m_quads_4e_trig.size() == 0 && result[_4e] == 0) {
        if (m_debug) cout << "Failed Trigger Match (4e) cut." << endl;
        result[_4e] = 8;
    }

    // single quadruplet selection per type
    if (result[_4mu] == 0) {
        if (m_quads_4mu_trig.size() == 1) {
            m_quad_4mu = m_quads_4mu_trig[0];
        }
        else {
            float mZ1_diff = 100000.;
            for (int i = 0; i < (int) m_quads_4mu_trig.size(); i++) {
                float dM = fabs(m_quads_4mu_trig[i]->z1.M() - mZ);
                if (dM < mZ1_diff) {
                    mZ1_diff = dM;
                }
            }
            for (int i = 0; i < (int) m_quads_4mu_trig.size(); i++) {
                float dM = fabs(m_quads_4mu_trig[i]->z1.M() - mZ);
                if (dM == mZ1_diff) {
                    m_quads_4mu_z1.push_back(m_quads_4mu_trig[i]);
                }
            }
            float mZ2_diff = 100000.;
            int index = 0;
            for (int i = 0; i < (int) m_quads_4mu_z1.size(); i++) {
                float dM = fabs(m_quads_4mu_z1[i]->z2.M() - mZ);
                if (dM < mZ2_diff) {
                    mZ2_diff = dM;
                    index = i;
                }
            }
            m_quad_4mu = m_quads_4mu_z1[index];
        }
    }
    if (result[_2e2mu] == 0) {
        if (m_quads_2e2mu_trig.size() == 1) {
            m_quad_2e2mu = m_quads_2e2mu_trig[0];
        }
        else {
            float mZ1_diff = 100000.;
            for (int i = 0; i < (int) m_quads_2e2mu_trig.size(); i++) {
                float dM = fabs(m_quads_2e2mu_trig[i]->z1.M() - mZ);
                if (dM < mZ1_diff) {
                    mZ1_diff = dM;
                }
            }
            for (int i = 0; i < (int) m_quads_2e2mu_trig.size(); i++) {
                float dM = fabs(m_quads_2e2mu_trig[i]->z1.M() - mZ);
                if (dM == mZ1_diff) {
                    m_quads_2e2mu_z1.push_back(m_quads_2e2mu_trig[i]);
                }
            }
            float mZ2_diff = 100000.;
            int index = 0;
            for (int i = 0; i < (int) m_quads_2e2mu_z1.size(); i++) {
                float dM = fabs(m_quads_2e2mu_z1[i]->z2.M() - mZ);
                if (dM < mZ2_diff) {
                    mZ2_diff = dM;
                    index = i;
                }
            }
            m_quad_2e2mu = m_quads_2e2mu_z1[index];
        }
    }
    if (result[_2mu2e] == 0) {
        if (m_quads_2mu2e_trig.size() == 1) {
            m_quad_2mu2e = m_quads_2mu2e_trig[0];
        }
        else {
            float mZ1_diff = 100000.;
            for (int i = 0; i < (int) m_quads_2mu2e_trig.size(); i++) {
                float dM = fabs(m_quads_2mu2e_trig[i]->z1.M() - mZ);
                if (dM < mZ1_diff) {
                    mZ1_diff = dM;
                }
            }
            for (int i = 0; i < (int) m_quads_2mu2e_trig.size(); i++) {
                float dM = fabs(m_quads_2mu2e_trig[i]->z1.M() - mZ);
                if (dM == mZ1_diff) {
                    m_quads_2mu2e_z1.push_back(m_quads_2mu2e_trig[i]);
                }
            }
            float mZ2_diff = 100000.;
            int index = 0;
            for (int i = 0; i < (int) m_quads_2mu2e_z1.size(); i++) {
                float dM = fabs(m_quads_2mu2e_z1[i]->z2.M() - mZ);
                if (dM < mZ2_diff) {
                    mZ2_diff = dM;
                    index = i;
                }
            }
            m_quad_2mu2e = m_quads_2mu2e_z1[index];
        }
    }
    if (result[_4e] == 0) {
        if (m_quads_4e_trig.size() == 1) {
            m_quad_4e = m_quads_4e_trig[0];
        }
        else {
            float mZ1_diff = 100000.;
            for (int i = 0; i < (int) m_quads_4e_trig.size(); i++) {
                float dM = fabs(m_quads_4e_trig[i]->z1.M() - mZ);
                if (dM < mZ1_diff) {
                    mZ1_diff = dM;
                }
            }
            for (int i = 0; i < (int) m_quads_4e_trig.size(); i++) {
                float dM = fabs(m_quads_4e_trig[i]->z1.M() - mZ);
                if (dM == mZ1_diff) {
                    m_quads_4e_z1.push_back(m_quads_4e_trig[i]);
                }
            }
            float mZ2_diff = 100000.;
            int index = 0;
            for (int i = 0; i < (int) m_quads_4e_z1.size(); i++) {
                float dM = fabs(m_quads_4e_z1[i]->z2.M() - mZ);
                if (dM < mZ2_diff) {
                    mZ2_diff = dM;
                    index = i;
                }
            }
            m_quad_4e = m_quads_4e_z1[index];
        }
    }

    // apply additional cuts per quadruplet in order: 4mu, 2e2mu, 2mu2e, 4e
    if (result[_4mu] == 0) {
        int code = applyCuts(m_quad_4mu);
        result[_4mu] = code;
        if (code == 0) {
            if (result[_2e2mu] == 0) result[_2e2mu] = 8;
            if (result[_2mu2e] == 0) result[_2mu2e] = 8;
            if (result[_4e] == 0) result[_4e] = 8;
            m_quad = m_quad_4mu;
            m_pass = true;
        }
    }
    if (result[_2e2mu] == 0 && result[_4mu] != 0) {
        int code = applyCuts(m_quad_2e2mu);
        result[_2e2mu] = code;
        if (code == 0) {
            if (result[_2mu2e] == 0) result[_2mu2e] = 8;
            if (result[_4e] == 0) result[_4e] = 8;
            m_quad = m_quad_2e2mu;
            m_pass = true;
        }
    }
    if (result[_2mu2e] == 0 && result[_4mu] != 0 && result[_2e2mu] != 0) {
        int code = applyCuts(m_quad_2mu2e);
        result[_2mu2e] = code;
        if (code == 0) {
            if (result[_4e] == 0) result[_4e] = 8;
            m_quad = m_quad_2mu2e;
            m_pass = true;
        }
    }
    if (result[_4e] == 0 && result[_4mu] != 0 && result[_2e2mu] != 0 && result[_2mu2e] != 0) {
        int code = applyCuts(m_quad_4e);
        result[_4e] = code;
        if (code == 0) {
            m_quad = m_quad_4e;
            m_pass = true;
        }
    }

    return result;
}

int QuadrupletSelection::applyCuts(quadruplet *quad) 
{
    // Z1 Mass cut
    if (!(quad->z1.M() > 50000.0 && quad->z1.M() < 106000.0)) {
        if (m_debug) cout << "Quadruplet failed Z1 Mass cut." << endl;
        return 9;
    }

    // Z2 Mass cut
    float z2_low;
    if (quad->q.M() < 140000.0) z2_low = 12000.0;
    else if (quad->q.M() > 190000.0) z2_low = 50000.0;
    else z2_low = 12000.0 + 0.76 * (quad->q.M() - 140000.0);
    if (!(quad->z2.M() > z2_low && quad->z2.M() < 115000.0)) {
        if (m_debug) cout << "Quadruplet failed Z2 Mass cut." << endl;
        return 10;
    }

    vector<const xAOD::Electron*> els;
    vector<const xAOD::Muon*> mus;
    fillVecs(quad, els, mus);

    // DeltaR / J/Psi cut
    if (quad->type == _4mu) {
        bool b_dR = true, b_jpsi = true;
        for (int imu1 = 0; imu1 < 3; imu1++) {
            for (int imu2 = imu1 + 1; imu2 < 4; imu2++) {
                if (b_dR == false) break;
                float phi1 = (float) mus[imu1]->phi();
                float phi2 = (float) mus[imu2]->phi();
                float eta1 = (float) mus[imu1]->eta();
                float eta2 = (float) mus[imu2]->eta();
                float dPhi = (fabs(phi1 - phi2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi1 - phi2) : fabs(phi1 - phi2);
                float dEta = fabs(eta1 - eta2);
                float dR = sqrt(dPhi*dPhi + dEta*dEta);
                if (!(dR > 0.10)) b_dR = false;
            }
        }
        TLorentzVector cross1 = quad->l1 + quad->l4;
        TLorentzVector cross2 = quad->l2 + quad->l3;
        if (cross1.M() < 5000.0 || cross2.M() < 5000.0) b_jpsi = false;
        if (!(b_dR && b_jpsi)) {
            if (m_debug) cout << "Quadruplet failed DeltaR / J/Psi cut." << endl;
            return 11;
        }
    }
    else if (quad->type == _4e) {
        bool b_dR = true, b_jpsi = true;
        for (int iel1 = 0; iel1 < 3; iel1++) {
            for (int iel2 = iel1 + 1; iel2 < 4; iel2++) {
                if (b_dR == false) break;
                float phi1 = (float) els[iel1]->trackParticle()->phi();
                float phi2 = (float) els[iel2]->trackParticle()->phi();
                float eta1 = (float) els[iel1]->trackParticle()->eta();
                float eta2 = (float) els[iel2]->trackParticle()->eta();
                float dPhi = (fabs(phi1 - phi2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi1 - phi2) : fabs(phi1 - phi2);
                float dEta = fabs(eta1 - eta2);
                float dR = sqrt(dPhi*dPhi + dEta*dEta);
                if (!(dR > 0.10)) b_dR = false;
            }
        }
        TLorentzVector cross1 = quad->l1 + quad->l4;
        TLorentzVector cross2 = quad->l2 + quad->l3;
        if (cross1.M() < 5000.0 || cross2.M() < 5000.0) b_jpsi = false;
        if (!(b_dR && b_jpsi)) {
            if (m_debug) cout << "Quadruplet failed DeltaR / J/Psi cut." << endl;
            return 11;
        }
    }
    else {
        float phi_e1 = (float) els[0]->trackParticle()->phi();
        float phi_e2 = (float) els[1]->trackParticle()->phi();
        float phi_m1 = (float) mus[0]->phi();
        float phi_m2 = (float) mus[1]->phi();
        float eta_e1 = (float) els[0]->trackParticle()->eta();
        float eta_e2 = (float) els[1]->trackParticle()->eta();
        float eta_m1 = (float) mus[0]->eta();
        float eta_m2 = (float) mus[1]->eta();
        // two electrons
        float dPhi_e = (fabs(phi_e1 - phi_e2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_e1 - phi_e2) : fabs(phi_e1 - phi_e2);
        float dEta_e = fabs(eta_e1 - eta_e2);
        float dR_e = sqrt(dPhi_e*dPhi_e + dEta_e*dEta_e);
        // two muons
        float dPhi_m = (fabs(phi_m1 - phi_m2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_m1 - phi_m2) : fabs(phi_m1 - phi_m2);
        float dEta_m = fabs(eta_m1 - eta_m2);
        float dR_m = sqrt(dPhi_m*dPhi_m + dEta_m*dEta_m);
        // e1 m1
        float dPhi_e1m1 = (fabs(phi_e1 - phi_m1) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_e1 - phi_m1) : fabs(phi_e1 - phi_m1);
        float dEta_e1m1 = fabs(eta_e1 - eta_m1);
        float dR_e1m1 = sqrt(dPhi_e1m1*dPhi_e1m1 + dEta_e1m1*dEta_e1m1);
        // e1 m2
        float dPhi_e1m2 = (fabs(phi_e1 - phi_m2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_e1 - phi_m2) : fabs(phi_e1 - phi_m2);
        float dEta_e1m2 = fabs(eta_e1 - eta_m2);
        float dR_e1m2 = sqrt(dPhi_e1m2*dPhi_e1m2 + dEta_e1m2*dEta_e1m2);
        // e2 m1
        float dPhi_e2m1 = (fabs(phi_e2 - phi_m1) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_e2 - phi_m1) : fabs(phi_e2 - phi_m1);
        float dEta_e2m1 = fabs(eta_e2 - eta_m1);
        float dR_e2m1 = sqrt(dPhi_e2m1*dPhi_e2m1 + dEta_e2m1*dEta_e2m1);
        // e2 m2
        float dPhi_e2m2 = (fabs(phi_e2 - phi_m2) > TMath::Pi()) ? 2 * TMath::Pi() - fabs(phi_e2 - phi_m2) : fabs(phi_e2 - phi_m2);
        float dEta_e2m2 = fabs(eta_e2 - eta_m2);
        float dR_e2m2 = sqrt(dPhi_e2m2*dPhi_e2m2 + dEta_e2m2*dEta_e2m2);
        if (!(dR_e > 0.10 && dR_m > 0.10 && dR_e1m1 > 0.20 && dR_e1m2 > 0.20 && dR_e2m1 > 0.20 && dR_e2m2 > 0.20)) {
            if (m_debug) cout << "Quadruplet failed DeltaR / J/Psi cut." << endl;
            return 11;
        }
    }

    int npart = 4;

    bool b_trackiso = true;
    bool b_d0 = true;
    if (quad->type == _4mu) {
        for (int imu = 0; imu < npart; imu++) {
            float d0sig = xAOD::TrackingHelpers::d0significance(mus.at(imu)->primaryTrackParticle(),
                    ei->beamPosSigmaX(), ei->beamPosSigmaY(), ei->beamPosSigmaXY());
            if(!cp_tools_->PassIsolation(*(mus.at(imu)))) b_trackiso = false;
            if (!(d0sig < 3.5)) b_d0 = false;
        }
    } else if (quad->type == _4e) {
        for (int iel = 0; iel < npart; iel++) {
            float d0sig = xAOD::TrackingHelpers::d0significance(els.at(iel)->trackParticle(),
                    ei->beamPosSigmaX(), ei->beamPosSigmaY(), ei->beamPosSigmaXY());
            if(!cp_tools_->PassIsolation(*(els.at(iel)))) b_trackiso = false;
            if (!(d0sig < 6.5)) b_d0 = false;
        }
    } else {
        for (int imu = 0; imu < 2; imu++) {
            float d0sig = xAOD::TrackingHelpers::d0significance(mus.at(imu)->primaryTrackParticle(),
                    ei->beamPosSigmaX(), ei->beamPosSigmaY(), ei->beamPosSigmaXY());
            if(!cp_tools_->PassIsolation(*(mus.at(imu)))) b_trackiso = false;
            if (!(d0sig < 3.5)) b_d0 = false;
        }
        for (int iel = 0; iel < 2; iel++) {
            float d0sig = xAOD::TrackingHelpers::d0significance(els.at(iel)->trackParticle(),
                    ei->beamPosSigmaX(), ei->beamPosSigmaY(), ei->beamPosSigmaXY());
            if(!cp_tools_->PassIsolation(*(els.at(iel)))) b_trackiso = false;
            if (!(d0sig < 6.5)) b_d0 = false;
        }
    }

    // Isolation cut
    // d0 Significance cut
    bool b_caloiso = b_trackiso;
    if(!doCR_){
        if (!(b_trackiso)) {
            if (m_debug) cout << "Quadruplet failed Track Isolation cut." << endl;
            return 12;
        }
        if (!(b_caloiso)) {
            if (m_debug) cout << "Quadruplet failed Calorimeter Isolation cut." << endl;
            return 13;
        }

        if (!(b_d0)) {
            if (m_debug) cout << "Quadruplet failed d0 Significance cut." << endl;
            return 14;
        }
    } else {
        /* CR defined as if passed isolation and ip, reject the event.
         * which means survival events fail at least one cut.
         * */
        bool all_passed = b_trackiso && b_caloiso && b_d0;
        if(all_passed) return 12;
    }

    if (m_debug) cout << "Quadruplet passed all cuts." << endl;
    return 0;
}

void QuadrupletSelection::fillVecs(quadruplet *quad, vector<const xAOD::Electron*> &els, vector<const xAOD::Muon*> &mus) {
    if (quad->type == _4mu) {
        quadruplet_4mu *quad_4mu = static_cast<quadruplet_4mu*>(quad);
        mus.push_back(quad_4mu->m1);
        mus.push_back(quad_4mu->m2);
        mus.push_back(quad_4mu->m3);
        mus.push_back(quad_4mu->m4);
    }
    else if (quad->type == _4e) {
        quadruplet_4e *quad_4e = static_cast<quadruplet_4e*>(quad);
        els.push_back(quad_4e->m1);
        els.push_back(quad_4e->m2);
        els.push_back(quad_4e->m3);
        els.push_back(quad_4e->m4);
    }
    else if (quad->type == _2e2mu) {
        quadruplet_2e2mu *quad_2e2mu = static_cast<quadruplet_2e2mu*>(quad);
        els.push_back(quad_2e2mu->m1);
        els.push_back(quad_2e2mu->m2);
        mus.push_back(quad_2e2mu->m3);
        mus.push_back(quad_2e2mu->m4);
    }
    else if (quad->type == _2mu2e) {
        quadruplet_2mu2e *quad_2mu2e = static_cast<quadruplet_2mu2e*>(quad);
        mus.push_back(quad_2mu2e->m1);
        mus.push_back(quad_2mu2e->m2);
        els.push_back(quad_2mu2e->m3);
        els.push_back(quad_2mu2e->m4);
    }
}

void QuadrupletSelection::doCR()
{ 
    // cout <<" In CR mode" << endl;
    doCR_ = true; 
}
