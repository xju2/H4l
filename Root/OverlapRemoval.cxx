#include <H4l/OverlapRemoval.h>
#include "xAODEgamma/ElectronxAODHelpers.h"
#include <iostream>
#include <TVector2.h>

using namespace std;

void OverlapRemoval::applyOverlapRemoval(
        xAOD::ElectronContainer *m_el_pass, 
        xAOD::MuonContainer *m_mu_pass, 
        xAOD::JetContainer *m_j_pass, bool debug) 
{
    // set true for all
    for (const auto& electron: *m_el_pass) {
        dec_passOR(*electron) = dec_signal(*electron);
    }
    for (const auto& muon: *m_mu_pass) {
        dec_passOR(*muon) = dec_signal(*muon);
    }
    for (const auto& jet: *m_j_pass) {
        dec_passOR(*jet) = dec_signal(*jet);
    }
    // e-e overlap
    for (int ie1 = 0; ie1 < (int) m_el_pass->size() - 1; ie1++) {
        if (! dec_passOR(*((*m_el_pass)[ie1]))) continue;
        const xAOD::TrackParticle* track_1 = 
            xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(m_el_pass->at(ie1)->trackParticle());
        float e1_d0 = (float) track_1->d0();
        float e1_z0 = (float) track_1->z0();
        float e1_theta = (float) track_1->theta();
        float e1_phi = (float)   track_1->phi();
        float e1_qOverP = (float) track_1->qOverP();
        float e1_et = (float) (*m_el_pass)[ie1]->caloCluster()->et();
        float e1_eta = (float) (*m_el_pass)[ie1]->caloCluster()->eta();
        for (int ie2 = ie1 + 1;  ie2 < (int) m_el_pass->size(); ie2++) {
            if (! dec_passOR(*((*m_el_pass)[ie2]))) continue;
            const xAOD::TrackParticle* track_2 = 
                xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(m_el_pass->at(ie2)->trackParticle());
            float e2_d0 = (float) track_2->d0();
            float e2_z0 = (float) track_2->z0();
            float e2_theta = (float) track_2->theta();
            float e2_phi = (float) track_2->phi();
            float e2_qOverP = (float) track_2->qOverP();
            float e2_et = (float) (*m_el_pass)[ie2]->caloCluster()->et();
            float e2_eta = (float) (*m_el_pass)[ie2]->caloCluster()->eta();

            float d_eta = fabs(e1_eta - e2_eta);
            float d_phi = TVector2::Phi_mpi_pi(e1_phi - e2_phi);
            if ((e1_d0 == e2_d0 && e1_z0 == e2_z0 && e1_theta == e2_theta && e1_phi == e2_phi && e1_qOverP == e2_qOverP) || (d_eta < 3 * 0.025 && d_phi < 5 * 0.025)) {
                if (e1_et < e2_et){ 
                    dec_passOR(*((*m_el_pass)[ie1])) = false;
                } else {
                    dec_passOR(*((*m_el_pass)[ie2])) = false;
                }
            }
        }
    }

    // e-mu overlap
    for (int ie = 0; ie < (int) m_el_pass->size(); ie++) {
        if (! dec_passOR(*((*m_el_pass)[ie]))) continue;
        const xAOD::TrackParticle* track_ele = 
            xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(m_el_pass->at(ie)->trackParticle());
        float e_phi = (float) track_ele->phi();
        float e_eta = (float) track_ele->eta();
        for (auto mu_iter = m_mu_pass->begin(); 
                mu_iter != m_mu_pass->end(); ++mu_iter)
        {
            if(! dec_passOR(**mu_iter) ) continue;
            if((*mu_iter)->muonType() == xAOD::Muon::MuonType::MuonStandAlone) continue;
            float cb_phi = (float)(*mu_iter)->trackParticle(xAOD::Muon::InnerDetectorTrackParticle)->phi();
            float cb_theta = (float)(*mu_iter)->trackParticle(xAOD::Muon::InnerDetectorTrackParticle)->theta();
            float cb_eta = (float) -1.0 * log(tan(cb_theta/2.0));
            float d_eta = fabs(e_eta - cb_eta);
            float d_phi =TVector2::Phi_mpi_pi(e_phi - cb_phi); 
            float d_R = sqrt(d_eta*d_eta + d_phi*d_phi);
            if (d_R < 0.02){ 
                if((*mu_iter)->muonType() != xAOD::Muon::MuonType::Combined && 
                   (*mu_iter)->muonType() != xAOD::Muon::MuonType::SegmentTagged){
                    dec_passOR(*((*m_el_pass)[ie])) = false;
                }
                if((*mu_iter)->muonType() == xAOD::Muon::MuonType::CaloTagged){
                    dec_passOR(**mu_iter) = false;
                }
            }
        }
    }

    // j-e overlap
    for (int ij = 0; ij < (int) m_j_pass->size(); ij++) {
        if (! dec_passOR(*((*m_j_pass)[ij])) ) continue;
        float j_eta = (*m_j_pass)[ij]->jetP4(xAOD::JetEMScaleMomentum).eta();
        float j_phi = (*m_j_pass)[ij]->jetP4(xAOD::JetEMScaleMomentum).phi();
        for (int ie = 0; ie < (int) m_el_pass->size(); ie++) {
            if (! dec_passOR(*((*m_el_pass)[ie])) ||
                ! dec_signal(*((*m_el_pass)[ie]))) continue;
            float e_phi = (float) (*m_el_pass)[ie]->trackParticle()->phi();
            float e_eta = (float) (*m_el_pass)[ie]->trackParticle()->eta();
            float d_eta = fabs(e_eta - j_eta);
            float d_phi = TVector2::Phi_mpi_pi(e_phi - j_phi);
            float d_R = sqrt(d_eta*d_eta + d_phi*d_phi);
            if (d_R < 0.2){ 
                dec_passOR(*((*m_j_pass)[ij])) = false;
            }
        }
    }
}
