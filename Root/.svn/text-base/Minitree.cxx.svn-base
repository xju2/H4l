#include <H4l/Minitree.h>
#include <TVector2.h>
using namespace std;

Minitree::Minitree(bool debug) {
    m_debug = debug;

    m_tree_incl_all = new TTree("tree_incl_all", "tree_incl_all");

    makeBranches(m_tree_incl_all);
}

Minitree::~Minitree(){
    // delete m_tree_incl_all;
}

void Minitree::cleanTree(){
    m4l_constrained = -999;
    mZ1_constrained = -999;
    mZ2_constrained = -999;
    weight_corr = -999;
    weight_lumi = -999;
    weight_sampleoverlap = -999;
    weight = -999;
    n_jets = -999;
    n_good_jets = -999;
    dijet_invmass = -999;
    leading_jet_pt = -999;
    subleading_jet_pt = -999;
    BDT_discriminant = -999;
    BDT_discriminant_VBF = -999;
    BDT_discriminant_HadVH = -999;
}

void Minitree::fillTree(int runNumber, int eventNumber,
        quadruplet *quad, quadruplet *quad_fsr) 
{
    run = runNumber;
    event = eventNumber;
    event_type = quad->type;
    m4l_unconstrained = quad->q.M() / 1000.0;
    mZ1_unconstrained = quad->z1.M() / 1000.0;
    mZ2_unconstrained = quad->z2.M() / 1000.0;
    m4l_fsr = quad_fsr->q.M() / 1000.0;
    mZ1_fsr = quad_fsr->z1.M() / 1000.0;
    mZ2_fsr = quad_fsr->z2.M() / 1000.0;
    Z1_lepplus_pt = quad->l1.Pt() / 1000.0;
    Z1_lepminus_pt = quad->l2.Pt() / 1000.0;
    Z2_lepplus_pt = quad->l3.Pt() / 1000.0;
    Z2_lepminus_pt = quad->l4.Pt() / 1000.0;
    Z1_lepplus_eta = quad->l1.Eta();
    Z1_lepminus_eta = quad->l2.Eta();
    Z2_lepplus_eta = quad->l3.Eta();
    Z2_lepminus_eta = quad->l4.Eta();
    Z1_lepplus_phi = quad->l1.Phi();
    Z1_lepminus_phi = quad->l2.Phi();
    Z2_lepplus_phi = quad->l3.Phi();
    Z2_lepminus_phi = quad->l4.Phi();
    Z1_lepplus_m = quad->l1.M() / 1000.0;
    Z1_lepminus_m = quad->l2.M() / 1000.0;
    Z2_lepplus_m = quad->l3.M() / 1000.0;
    Z2_lepminus_m = quad->l4.M() / 1000.0;
    MET = quad->met / 1000.0;
    MET_final_ = quad->met_final /1e3;
    MET_phi_ = quad->met_phi;
    TLorentzVector higgs_tlv(quad->l1+quad->l2+quad->l3+quad->l4);
    Hpt_ = higgs_tlv.Pt()/1e3;
    dphi_met_higgs_ = fabs(TMath::ACos(TMath::Cos(MET_phi_ - higgs_tlv.Phi())));
    TVector2 final_met(quad->mpx, quad->mpy);
    TVector2 higgs_xy(higgs_tlv.Px(), higgs_tlv.Py());
    MET_noHiggs_ = (final_met - higgs_xy).Mod();

    m_tree_incl_all->Fill();

    return;
}

void Minitree::makeBranches(TTree *tree) {
    tree->Branch("run", &run, "run/I");
    tree->Branch("event", &event, "event/I");
    tree->Branch("event_type", &event_type, "event_type/I");
    tree->Branch("m4l_unconstrained", &m4l_unconstrained, "m4l_unconstrained/F");
    tree->Branch("mZ1_unconstrained", &mZ1_unconstrained, "mZ1_unconstrained/F");
    tree->Branch("mZ2_unconstrained", &mZ2_unconstrained, "mZ2_unconstrained/F");
    tree->Branch("m4l_fsr", &m4l_fsr, "m4l_fsr/F");
    tree->Branch("mZ1_fsr", &mZ1_fsr, "mZ1_fsr/F");
    tree->Branch("mZ2_fsr", &mZ2_fsr, "mZ2_fsr/F");
    tree->Branch("m4l_constrained", &m4l_constrained, "m4l_constrained/F");
    tree->Branch("mZ1_constrained", &mZ1_constrained, "mZ1_constrained/F");
    tree->Branch("mZ2_constrained", &mZ2_constrained, "mZ2_constrained/F");
    tree->Branch("weight_corr", &weight_corr, "weight_corr/F");
    tree->Branch("weight_lumi", &weight_lumi, "weight_lumi/F");
    tree->Branch("weight_sampleoverlap", &weight_sampleoverlap, "weight_sampleoverlap/F");
    tree->Branch("weight", &weight, "weight/F");
    tree->Branch("Z1_lepplus_pt", &Z1_lepplus_pt, "Z1_lepplus_pt/F");
    tree->Branch("Z1_lepminus_pt", &Z1_lepminus_pt, "Z1_lepminus_pt/F");
    tree->Branch("Z2_lepplus_pt", &Z2_lepplus_pt, "Z2_lepplus_pt/F");
    tree->Branch("Z2_lepminus_pt", &Z2_lepminus_pt, "Z2_lepminus_pt/F");
    tree->Branch("Z1_lepplus_eta", &Z1_lepplus_eta, "Z1_lepplus_eta/F");
    tree->Branch("Z1_lepminus_eta", &Z1_lepminus_eta, "Z1_lepminus_eta/F");
    tree->Branch("Z2_lepplus_eta", &Z2_lepplus_eta, "Z2_lepplus_eta/F");
    tree->Branch("Z2_lepminus_eta", &Z2_lepminus_eta, "Z2_lepminus_eta/F");
    tree->Branch("Z1_lepplus_phi", &Z1_lepplus_phi, "Z1_lepplus_phi/F");
    tree->Branch("Z1_lepminus_phi", &Z1_lepminus_phi, "Z1_lepminus_phi/F");
    tree->Branch("Z2_lepplus_phi", &Z2_lepplus_phi, "Z2_lepplus_phi/F");
    tree->Branch("Z2_lepminus_phi", &Z2_lepminus_phi, "Z2_lepminus_phi/F");
    tree->Branch("Z1_lepplus_m", &Z1_lepplus_m, "Z1_lepplus_m/F");
    tree->Branch("Z1_lepminus_m", &Z1_lepminus_m, "Z1_lepminus_m/F");
    tree->Branch("Z2_lepplus_m", &Z2_lepplus_m, "Z2_lepplus_m/F");
    tree->Branch("Z2_lepminus_m", &Z2_lepminus_m, "Z2_lepminus_m/F");
    tree->Branch("n_jets", &n_jets, "n_jets/I");
    tree->Branch("n_good_jets", &n_good_jets, "n_good_jets/I");
    tree->Branch("dijet_invmass", &dijet_invmass, "dijet_invmass/F");
    tree->Branch("leading_jet_pt", &leading_jet_pt, "leading_jet_pt/F");
    tree->Branch("subleading_jet_pt", &subleading_jet_pt, "subleading_jet_pt/F");
    tree->Branch("BDT_discriminant", &BDT_discriminant, "BDT_discriminant/F");
    tree->Branch("BDT_discriminant_VBF", &BDT_discriminant_VBF, "BDT_discriminant_VBF/F");
    tree->Branch("BDT_discriminant_HadVH", &BDT_discriminant_HadVH, "BDT_discriminant_HadVH/F");
    tree->Branch("MET", &MET, "MET/F");
    tree->Branch("MET_Final", &MET_final_, "MET_Final/F");
    tree->Branch("MET_phi", &MET_phi_, "MET_phi/F");
    tree->Branch("Hpt", &Hpt_, "Hpt/F");
    tree->Branch("dphi_met_higgs", &dphi_met_higgs_, "dphi_met_higgs/F");
    tree->Branch("MET_noHiggs", &MET_noHiggs_, "MET_noHiggs/F");

    return;
}

void Minitree::setDirectory(TFile *outFile) {
    m_tree_incl_all->SetDirectory(outFile);

    return;
}
