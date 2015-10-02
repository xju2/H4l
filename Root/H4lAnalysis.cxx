#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include <H4l/H4lAnalysis.h>
#include <H4l/ElectronSelection.h>
#include <H4l/MuonSelection.h>
#include <H4l/JetSelection.h>
#include <H4l/OverlapRemoval.h>
#include <H4l/QuadrupletSelection.h>
#include <H4l/FSRCorrection.h>
#include <H4l/Minitree.h>

#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODCore/ShallowCopy.h"
#include "CPAnalysisExamples/errorcheck.h"


#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "FsrUtils/FsrPhotonTool.h"

// this is needed to distribute the algorithm to the workers
ClassImp(H4lAnalysis)
using namespace std;

const char* H4lAnalysis::APP_NAME = "H4lAnalysis";

H4lAnalysis :: H4lAnalysis (bool inclusive, bool smearing, bool weights, 
        bool doCR, bool is_data, bool atlfast, bool debug)
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

    // Job options
    m_inclusive = inclusive;
    m_smearing = smearing;
    m_weights = weights;
    m_doCR = doCR;
    m_is_data = is_data;
    m_debug = debug;
    m_atlfast = atlfast;

    if(m_doCR)  cout <<"doing Control Region" << endl;
    if(m_debug) cout << "In debug mode!!" << endl;

    m_jetInputType = xAOD::JetInput::EMTopo;

    jet_base_name_ = "AntiKt4"+xAOD::JetInput::typeName(m_jetInputType);
    jet_con_name_ = jet_base_name_+"Jets";
    met_core_name_ = "MET_Core_"+jet_base_name_;
    met_map_name_ = "METAssoc_"+jet_base_name_;

}

H4lAnalysis::~H4lAnalysis(){
}

EL::StatusCode H4lAnalysis :: setupJob (EL::Job& job)
{
    // Here you put code that sets up the job on the submission object
    // so that it is ready to work with your algorithm, e.g. you can
    // request the D3PDReader service or add output files.  Any code you
    // put here could instead also go into the submission script.  The
    // sole advantage of putting it here is that it gets automatically
    // activated/deactivated when you add/remove the algorithm from your
    // job, which may or may not be of value to you.

    job.useXAOD();

    xAOD::Init("H4LAnalysis").ignore();

    EL::OutputStream out("outputLabel");
    job.outputAdd(out);

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.
    
    // Histograms for cutflows
    TFile *outFile = wk()->getOutputFile("outputLabel");
    outFile->cd();
    h_es_fail = new TH1I("h_es_fail", "h_es_fail", 10, -0.5, 9.5);
    h_ms_fail = new TH1I("h_ms_fail", "h_ms_fail", 13, -0.5, 12.5);
    h_js_fail = new TH1I("h_js_fail", "h_js_fail", 8, -0.5, 7.5);
    h_qs_fail_4mu = new TH1I("h_qs_fail_4mu", "h_qs_fail_4mu", 15, -0.5, 14.5);
    h_qs_fail_4e = new TH1I("h_qs_fail_4e", "h_qs_fail_4e", 15, -0.5, 14.5);
    h_qs_fail_2mu2e = new TH1I("h_qs_fail_2mu2e", "h_qs_fail_2mu2e", 15, -0.5, 14.5);
    h_qs_fail_2e2mu = new TH1I("h_qs_fail_2e2mu", "h_qs_fail_2e2mu", 15, -0.5, 14.5);
    // wk()->addOutput(h_es_fail);
    // wk()->addOutput(h_ms_fail);
    // wk()->addOutput(h_js_fail);
    // wk()->addOutput(h_qs_fail_4mu);
    // wk()->addOutput(h_qs_fail_4e);
    // wk()->addOutput(h_qs_fail_2mu2e);
    // wk()->addOutput(h_qs_fail_2e2mu);

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    xAOD::TEvent* event = wk()->xaodEvent();
    uint64_t n_events_process = 0;
    double sum_of_evt_weights = 0;
    double sum_of_evt_weight_sqd = 0;
    CHECK(CPToolsHelper::GetProcessEventsInfo(*event,
                n_events_process,
                sum_of_evt_weights,
                sum_of_evt_weight_sqd));
    // cout << "event process: " << n_events_process << " " << total_evts_pro_ << endl;
    total_evts_pro_ += n_events_process;
    sum_of_evt_w_ += sum_of_evt_weights;
    sum_of_evt_w_sq_ += sum_of_evt_weight_sqd;
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: changeInput (bool firstFile)
{
    // Here you do everything you need to do when we change input files,
    // e.g. resetting branch addresses on trees.  If you are using
    // D3PDReader or a similar service this method is not needed.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.

    m_event = wk()->xaodEvent();
    m_store = wk()->xaodStore();

    TFile *outFile = wk()->getOutputFile("outputLabel");
    mt = new Minitree(m_debug);
    mt->setDirectory(outFile);


    Info("initialize()", "Number of events = %lli", m_event->getEntries());

    m_eventCounter = 0;


    // FSR correction tool
    m_fsrTool = new FSR::FsrPhotonTool("FSRTool");
    m_fsrTool->initialize();
    CHECK(m_fsrTool->initialize());

    // CP tools
    m_muonCorrectionTool = new CP::MuonCalibrationAndSmearingTool("MuonCorrectionTool");
    if (!m_muonCorrectionTool->initialize().isSuccess()) 
    {
        Error("initialize()","Failed to properly initialize the MuonCorrectionTool.");
        return EL::StatusCode::FAILURE;
    }

    m_egammaCorrectionTool = new CP::EgammaCalibrationAndSmearingTool("EgammaCalibrationAndSmearingTool");
    m_egammaCorrectionTool->setProperty("ESModel", "es2015PRE");
    m_egammaCorrectionTool->setProperty("decorrelationModel", "1NP_v1");
    if(m_atlfast) m_egammaCorrectionTool->setProperty("useAFII", true);
    if (!m_egammaCorrectionTool->initialize().isSuccess()) {
        Error("initialize()","Failed to properly initialize the EgammaCalibrationAndSmearingTool.");
        return EL::StatusCode::FAILURE;
    }
    /* config jet calibartion tools
     * add AtlFast II configuration
     * */ 
    string jes_config_file("JES_MC15Prerecommendation_April2015.config");
    if(m_atlfast) jes_config_file = "JES_MC15Prerecommendation_AFII_June2015.config";
    string calibseq = "JetArea_Origin_Residual_EtaJES_GSC";
    if(m_is_data) calibseq += "_Insitu";
    m_jetCorrectionTool = new JetCalibrationTool("JetCalibTool", 
            jet_base_name_.c_str(), 
            jes_config_file,
            calibseq,
            m_is_data);
    if (!m_jetCorrectionTool->initializeTool("JetCalibTool").isSuccess()) {
        Error("initialize()","Failed to properly initialize jet calibration tool.");
        return EL::StatusCode::FAILURE;
    }
    /* MET builder */
    met_maker_ = new met::METMaker("METMaker_H4l");
    met_maker_->initialize();
    met_maker_->msg().setLevel(MSG::INFO);

    // Muon Selection tool
    MS_ = new MuonSelection(m_debug);

    // Electron Selection tool
    ES_ = new ElectronSelection(m_debug);

    // Jet Selection tool
    // JS_ = new JetSelection(m_debug);
    JS_ = new JetSelection(false);
    cp_tools_ = new CPToolsHelper(); 

    event_info_filler_ = new EventInfoCreator();
    associate_tree_ = new TTree("associate", "associate");
    associate_tree_->SetDirectory(outFile);
    event_info_filler_->AttachMiniToTree(*associate_tree_);
    associate_tree_->Branch("nEventsProcessed", &total_evts_pro_, "nEventsProcessed/l");
    associate_tree_->Branch("nSumEventWeights", &sum_of_evt_w_, "nSumEventWeights/D");
    associate_tree_->Branch("nSumEventWeightsSquared", &sum_of_evt_w_sq_, 
            "nSumEventWeightsSquared/D");

    event_info_filler_->AttachBranchToTree(*(mt->m_tree_incl_all));

    trig_el_ = new vector<string>();
    trig_el_->push_back("HLT_e24_lhmedium_iloose_L1EM18VH");
    trig_el_->push_back("HLT_e60_lhmedium");
    trig_el_->push_back("HLT_e120_lhloose");
    trig_el_->push_back("HLT_e140_lhloose");
    trig_el_->push_back("HLT_2e12_lhloose_L12EM10VH");
    trig_el_->push_back("HLT_e17_lhloose_2e9_lhloose");

    trig_mu_ = new vector<string>();
    trig_mu_->push_back("HLT_mu20_iloose_L1MU15");
    trig_mu_->push_back("HLT_mu50");
    trig_mu_->push_back("HLT_mu60_0eta105_msonly");
    trig_mu_->push_back("HLT_2mu10");
    trig_mu_->push_back("HLT_mu18_mu8noL1");
    trig_mu_->push_back("HLT_3mu6");
    trig_mu_->push_back("HLT_3mu6_msonly");

    trig_mu_el_ = new vector<string>();
    trig_mu_el_->push_back("HLT_e17_lhloose_mu14");
    trig_mu_el_->push_back("HLT_2e12_lhloose_mu10");
    trig_mu_el_->push_back("HLT_e12_lhloose_2mu10");
    
    // electron fail codes
    m_es_failCodes.insert(std::make_pair<int, std::string>(0, "Pass"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(1, "Vertex"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(2, "Trigger"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(3, "Author"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(4, "Loose Likelihood"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(5, "eta"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(6, "Et"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(7, "OQ"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(8, "z0"));
    m_es_failCodes.insert(std::make_pair<int, std::string>(9, "Overlap"));

    // muon fail codes
    m_ms_failCodes.insert(std::make_pair<int, std::string>(0, "Pass"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(1, "Vertex"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(2, "Trigger"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(3, "Author"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(4, "pt"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(5, "eta"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(6, "SA MS"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(7, "Pix"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(8, "SCT"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(9, "Holes"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(10, "TRT"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(11, "d0/z0"));
    m_ms_failCodes.insert(std::make_pair<int, std::string>(12, "Overlap"));

    // jet fail codes
    m_js_failCodes.insert(std::make_pair<int, std::string>(0, "Pass"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(1, "Vertex"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(2, "Trigger"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(3, "pt"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(4, "eta"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(5, "Pile Up Removal"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(6, "Cleaning"));
    m_js_failCodes.insert(std::make_pair<int, std::string>(7, "Overlap"));

    // quadruplet fail codes
    m_qs_failCodes.insert(std::make_pair<int, std::string>(0, "Pass"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(1, "GRL"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(2, "LAr Error"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(3, "Vertex"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(4, "Trigger"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(5, "Selected Leptons"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(6, "SFOS"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(7, "Kinematics"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(8, "Trigger Match"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(9, "Z1 Mass"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(10, "Z2 Mass"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(11, "DeltaR / J/Psi"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(12, "Track Isolation"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(13, "Calorimeter Isolation"));
    m_qs_failCodes.insert(std::make_pair<int, std::string>(14, "d0 Significance"));
   
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: execute ()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.
    const char* APP_NAME = "H4lAnalysis";
    mt->cleanTree();
    event_info_filler_->ClearBranch();
    if(m_debug) cout << "In execute()" << endl;
    
    // if ((m_eventCounter % 1000) == 0) Info("execute()", "Event number = %i", m_eventCounter);
    m_eventCounter++;

    // Retrieve containers
    const xAOD::EventInfo *event_info = 0;
    if (!m_event->retrieve(event_info, "EventInfo").isSuccess()) {
        Error("execute()", "Failed to retrieve EventInfo. Exiting.");
        return EL::StatusCode::FAILURE;
    }
    event_info_filler_->Fill(*event_info);
    if(m_eventCounter == 1) associate_tree_->Fill();
    m_runNumber = event_info->runNumber();
    m_eventNumber = event_info->eventNumber();

    // pre-selections (event level)
    if(! cp_tools_->PassGRL(*event_info)) return EL::StatusCode::SUCCESS;
    if(! cp_tools_->PassEventCleaning(*event_info)) return EL::StatusCode::SUCCESS;

    //primary vertex
    const xAOD::VertexContainer* vertice = 0;
    CHECK( m_event->retrieve(vertice, "PrimaryVertices") );
    if(! cp_tools_->HasPrimaryVertex(*vertice, 1)) return EL::StatusCode::SUCCESS;
    if(m_debug) cout << "After primary vertex" << endl;
     
    // Trigger
    bool pass_trig_el = false;
    bool pass_trig_mu = false;
    bool pass_trig_mu_el = false;
    for(auto& trig_name : *trig_el_){
        if(cp_tools_->PassTrigger(trig_name)) {
            pass_trig_el = true;
            break;
        }
    }
    for(auto& trig_name : *trig_mu_){
        if(cp_tools_->PassTrigger(trig_name)){
            pass_trig_mu = true; break;
        }
    }
    for(auto& trig_name : *trig_mu_el_){
        if(cp_tools_->PassTrigger(trig_name)){
            pass_trig_mu_el = true; break;
        }
    }
    bool pass_trigger = pass_trig_el || pass_trig_mu || pass_trig_mu_el;

    const xAOD::ElectronContainer *electrons = 0;
    if (!m_event->retrieve(electrons, "Electrons").isSuccess()) {
        Error("execute()", "Failed to retrieve ElectronCollection. Exiting.");
        return EL::StatusCode::FAILURE;
    }
    if (m_debug) Info("execute()", "Number of electrons = %lu", electrons->size());
    
    const xAOD::PhotonContainer *photons = 0;
    if (!m_event->retrieve(photons, "Photons").isSuccess()) {
        Error("execute()", "Failed to retrieve PhotonCollection. Exiting.");
        return EL::StatusCode::FAILURE;
    }
    if (m_debug) Info("execute()", "Number of photons = %lu", photons->size());

    const xAOD::MuonContainer *muons = 0;
    if (!m_event->retrieve(muons, "Muons").isSuccess()) {
        Error("execute()", "Failed to retrieve Muons. Exiting.");
        return EL::StatusCode::FAILURE;
    }
    if (m_debug) Info("execute()", "Number of muons = %lu", muons->size());


    const xAOD::VertexContainer *vertices = 0;
    if (!m_event->retrieve(vertices, "PrimaryVertices").isSuccess()) {
        Error("execute()", "Failed to retrieve PrimaryVertices. Exiting.");
        return EL::StatusCode::FAILURE;
    }
    const xAOD::Vertex *pvx = vertices->at(0);

    const xAOD::JetContainer *jets = 0;
    if (!m_event->retrieve(jets, jet_con_name_).isSuccess()) 
    {
        Error("execute()", "Failed to retrieve %s. Exiting.", jet_con_name_.c_str());
        return EL::StatusCode::FAILURE;
    }
    if (m_debug) Info("execute()", "Number of jets = %lu", jets->size());



    // Electron selection

    std::pair<xAOD::ElectronContainer*,xAOD::ShallowAuxContainer*> el_shallowcopy = xAOD::shallowCopyContainer(*electrons);
    xAOD::ElectronContainer* el_copy = el_shallowcopy.first;
    xAOD::ShallowAuxContainer* el_copyaux = el_shallowcopy.second;
    xAOD::setOriginalObjectLink(*electrons, *el_copy);

    int el_count = 0;
    for (xAOD::ElectronContainer::const_iterator el_itr = el_copy->begin(); 
            el_itr != el_copy->end(); ++el_itr) {
        if (m_smearing && abs((*el_itr)->caloCluster()->eta()) <= 2.47 ) 
        {
            m_egammaCorrectionTool->setRandomSeed(event_info->eventNumber() + 100 * el_count);
            m_egammaCorrectionTool->applyCorrection(const_cast<xAOD::Electron&>(**el_itr));
        }
        int code = ES_->applyCuts(*el_itr, pvx);
        int binMax = (code == 0) ? h_es_fail->GetNbinsX() : code;
        if(!pass_trigger) binMax = 2;
        for (int ibin = 0; ibin < binMax; ibin++) {
            h_es_fail->Fill(ibin);
        }
        if (ES_->pass()) {
            dec_signal(**el_itr) = true;
        }else{
            dec_signal(**el_itr) = false;
        }
        el_count++;
    }

    // Photon smearing (for FSR correction)
    xAOD::PhotonContainer *ph_calib = new xAOD::PhotonContainer();
    xAOD::PhotonAuxContainer *ph_calib_aux = new xAOD::PhotonAuxContainer();
    ph_calib->setStore(ph_calib_aux);
    m_store->record(ph_calib, "photonCalib");
    m_store->record(ph_calib_aux, "photonCalibAux");
    int ph_count = 0;
    for (xAOD::PhotonContainer::const_iterator ph_itr = photons->begin(); ph_itr != photons->end(); ++ph_itr) {
        // Author cuts needed according to https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Photon_authors
        if ( !((*ph_itr)->author() & (xAOD::EgammaParameters::AuthorPhoton + xAOD::EgammaParameters::AuthorAmbiguous)) ) continue;
            
        xAOD::Photon *in_ph = new xAOD::Photon();
        in_ph->makePrivateStore(*ph_itr);
        if (m_smearing) {
            m_egammaCorrectionTool->setRandomSeed(event_info->eventNumber() + 100 * ph_count);
            m_egammaCorrectionTool->applyCorrection(*in_ph);
        }
        ph_calib->push_back(in_ph);
        ph_count++;
    }

    // Muon selection
    std::pair<xAOD::MuonContainer*,xAOD::ShallowAuxContainer*> mu_shallowcopy = xAOD::shallowCopyContainer(*muons);
    xAOD::MuonContainer* mu_copy = mu_shallowcopy.first;
    xAOD::ShallowAuxContainer* mu_copyaux = mu_shallowcopy.second;
    xAOD::setOriginalObjectLink(*muons, *mu_copy);

    for (xAOD::MuonContainer::const_iterator mu_itr = mu_copy->begin(); mu_itr != mu_copy->end(); ++mu_itr) {
        if (m_smearing){
            m_muonCorrectionTool->applyCorrection(const_cast<xAOD::Muon&>(**mu_itr));
        }

        int code = MS_->applyCuts(*mu_itr, pvx);
        int binMax = (code == 0) ? h_ms_fail->GetNbinsX() : code;
        if(!pass_trigger) binMax = 2;
        for (int ibin = 0; ibin < binMax; ibin++) {
            h_ms_fail->Fill(ibin);
        }
        if (MS_->pass()){ 
            dec_signal(**mu_itr) = true;
        }else{
            dec_signal(**mu_itr) = false;
        }
    }

    // Jet selection
    std::pair<xAOD::JetContainer*,xAOD::ShallowAuxContainer*> jets_shallowcopy = xAOD::shallowCopyContainer(*jets);
    xAOD::JetContainer* jet_copy = jets_shallowcopy.first;
    xAOD::ShallowAuxContainer* jet_copyaux = jets_shallowcopy.second;
    xAOD::setOriginalObjectLink(*jets, *jet_copy);
    int n_good_jets = 0;
    for (xAOD::JetContainer::const_iterator j_itr = jet_copy->begin(); 
            j_itr != jet_copy->end(); ++j_itr) 
    {
        if (m_smearing){ 
            m_jetCorrectionTool->applyCalibration(const_cast<xAOD::Jet&>(**j_itr));
        }

        int code = JS_->applyCuts(*j_itr, m_runNumber);
        int binMax = (code == 0) ? h_js_fail->GetNbinsX() : code;
        for (int ibin = 0; ibin < binMax; ibin++) {
            h_js_fail->Fill(ibin);
        }
        if (JS_->pass()){ 
            n_good_jets ++;
            dec_signal(**j_itr) = true;
        }else{
            dec_signal(**j_itr) = false;
        }
    }
    mt->n_good_jets = n_good_jets;
    mt->n_jets = (int)jet_copy->size();

    // recalculate the MET
    auto* metCon = new xAOD::MissingETContainer;
    auto* metAux = new xAOD::MissingETAuxContainer;
    metCon->setStore(metAux);
    m_store->record(metCon, "MyMETCon");
    m_store->record(metAux, "MyMETConAux");

    GetMET(*metCon, jet_copy, el_copy, mu_copy, 
            nullptr,  // gamma
            nullptr,  // taujet
            true // doTST
            ); 
    xAOD::MissingETContainer::const_iterator met_it = metCon->find("Final");
    double mpx = (*met_it)->mpx();
    double mpy = (*met_it)->mpy();
    double met_final = sqrt(mpx*mpx + mpy*mpy);
    float met_phi = (float) (*met_it)->phi();

    // Overlap removal
    if (m_debug) Info("execute()", "Going to overlap removal.");
    OverlapRemoval::applyOverlapRemoval(el_copy, mu_copy, jet_copy);
    
    //Get selected leptons
    auto* el_pass = new xAOD::ElectronContainer();
    auto* el_pass_aux = new xAOD::ElectronAuxContainer();
    el_pass->setStore(el_pass_aux);
    m_store->record(el_pass, "ElectronPass");
    m_store->record(el_pass_aux, "ElectronPassAux");
    for(auto el_iter = el_copy->begin(); el_iter != el_copy->end(); el_iter++)
    {
        if(dec_passOR(**el_iter) && dec_signal(**el_iter)){
            xAOD::Electron* el_new = new xAOD::Electron();
            el_new->makePrivateStore(**el_iter);
            el_pass->push_back(el_new);
        }
    }

    auto* mu_pass = new xAOD::MuonContainer();
    auto* mu_pass_aux = new xAOD::MuonAuxContainer();
    mu_pass->setStore(mu_pass_aux);
    m_store->record(mu_pass, "MuonPass");
    m_store->record(mu_pass_aux, "MuonPassAux");
    for(auto mu_iter = mu_copy->begin(); mu_iter != mu_copy->end(); mu_iter++)
    {
        if(dec_passOR(**mu_iter) && dec_signal(**mu_iter)){
            auto* muon_new = new xAOD::Muon();
            muon_new->makePrivateStore(**mu_iter);
            mu_pass->push_back(muon_new);
        }
    }
    if(m_debug){
        for(const auto& el : *el_pass){
            cout <<"Electron: " << el->p4().Pt() << endl;
        }
        for(const auto& mu : *mu_pass){
            cout <<"Muon: " << mu->p4().Pt() << endl;
        }
    }
    if(m_debug){
        Info("execute()", "Total good electrons: %ld", el_pass->size());
        Info("execute()", "Total good muons: %ld", mu_pass->size());
    }
    // Quadruplet selection
    if (m_debug) Info("execute()", "Going to quadruplet selection.");
    QuadrupletSelection *QS = new QuadrupletSelection(el_pass, mu_pass, m_debug);
    QS->SetCPTool(cp_tools_);
    QS->SetEventInfo(event_info);

    if(m_doCR) QS->doCR();
    std::vector<int> code;
    if (m_inclusive) code = QS->applyCuts_Inclusive();
    else code = QS->applyCuts_Type();
    int binMax_4mu = (code[_4mu] == 0) ? h_qs_fail_4mu->GetNbinsX() : code[_4mu];
    if(!pass_trigger) binMax_4mu = 4;
    for (int ibin = 0; ibin < binMax_4mu; ibin++) {
        h_qs_fail_4mu->Fill(ibin);
    }
    int binMax_4e = (code[_4e] == 0) ? h_qs_fail_4e->GetNbinsX() : code[_4e];
    if(!pass_trigger)  binMax_4e = 4;
    for (int ibin = 0; ibin < binMax_4e; ibin++) {
        h_qs_fail_4e->Fill(ibin);
    }
    int binMax_2e2mu = (code[_2e2mu] == 0) ? h_qs_fail_2e2mu->GetNbinsX() : code[_2e2mu];
    if(!pass_trigger)  binMax_2e2mu = 4;
    for (int ibin = 0; ibin < binMax_2e2mu; ibin++) {
        h_qs_fail_2e2mu->Fill(ibin);
    }
    int binMax_2mu2e = (code[_2mu2e] == 0) ? h_qs_fail_2mu2e->GetNbinsX() : code[_2mu2e];
    if(!pass_trigger)  binMax_2mu2e = 4;
    for (int ibin = 0; ibin < binMax_2mu2e; ibin++) {
        h_qs_fail_2mu2e->Fill(ibin);
    }
    
    if (!QS->pass()) {
        if (m_debug) Info("execute()", "No quadruplet passes full selection.");
    }
    else {
        quadruplet *quad = QS->getQuadruplet();

        // quad->met = (float) met->met();
        quad->met =  -9999.0;
        quad->met_final = (float) met_final;
        quad->met_phi = met_phi;
        quad->mpx = mpx;
        quad->mpy = mpy;

        // FSR correction
        if (m_debug) Info("execute()", "Applying FSR correction.");
        FSRCorrection *fsr = new FSRCorrection(quad, m_debug);
        fsr->applyFSR(m_fsrTool, el_copy, ph_calib);
        quadruplet *quad_fsr = fsr->getQuadruplet();

        // Fill minitree
        if (m_debug) Info("execute()", "Filling minitree.");
        mt->fillTree((int) m_runNumber, (int) m_eventNumber, quad, quad_fsr);

        // if(quad_fsr) delete quad_fsr;
        if(fsr) delete fsr;
        // if(quad) delete quad;
    }
    // delete QS;
    
    
    m_store->clear();

    delete el_copy;
    delete el_copyaux;
    delete mu_copy;
    delete mu_copyaux;
    delete jet_copy;
    delete jet_copyaux;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: postExecute ()
{
    // Here you do everything that needs to be done after the main event
    // processing.  This is typically very rare, particularly in user
    // code.  It is mainly used in implementing the NTupleSvc.
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: finalize ()
{
    // This method is the mirror image of initialize(), meaning it gets
    // called after the last event has been processed on the worker node
    // and allows you to finish up any objects you created in
    // initialize() before they are written to disk.  This is actually
    // fairly rare, since this happens separately for each worker node.
    // Most of the time you want to do your post-processing on the
    // submission node after all your histogram outputs have been
    // merged.  This is different from histFinalize() in that it only
    // gets called on worker nodes that processed input events.
    
    // delete mt;
    delete m_fsrTool;
    delete m_muonCorrectionTool;
    delete m_egammaCorrectionTool;
    delete m_jetCorrectionTool;
    delete met_maker_; 
    delete MS_;
    delete ES_;
    delete JS_;
    
    delete trig_el_;
    delete trig_mu_;
    delete trig_mu_el_;

    delete cp_tools_;
    delete event_info_filler_;

    return EL::StatusCode::SUCCESS;
}



EL::StatusCode H4lAnalysis :: histFinalize ()
{
    // This method is the mirror image of histInitialize(), meaning it
    // gets called after the last event has been processed on the worker
    // node and allows you to finish up any objects you created in
    // histInitialize() before they are written to disk.  This is
    // actually fairly rare, since this happens separately for each
    // worker node.  Most of the time you want to do your
    // post-processing on the submission node after all your histogram
    // outputs have been merged.  This is different from finalize() in
    // that it gets called on all worker nodes regardless of whether
    // they processed input events.
    //
    
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode H4lAnalysis::GetMET(xAOD::MissingETContainer &met,
				     const xAOD::JetContainer* jet,
				     const xAOD::ElectronContainer* elec,
				     const xAOD::MuonContainer* muon,
				     const xAOD::PhotonContainer* gamma,
				     const xAOD::TauJetContainer* taujet,
				     bool doTST) 
{
   if(!jet) {
        Warning("GetMET", "Invalid jet container specified for MET rebuilding!");
        return StatusCode::SUCCESS;
    }
    string m_eleTerm = "RefEle";
    string m_gammaTerm = "RefGamma";
    string m_tauTerm = "RefTau";
    string m_jetTerm = "RefJet";
    string m_muonTerm = "Muons";
    string m_outMETTerm = "Final";

    const xAOD::MissingETContainer* metcore(0);
    if ( m_event->retrieve( metcore, met_core_name_).isFailure() ) {
        Warning("GetMET", "Unable to retrieve MET core container: %s", met_core_name_.c_str());
        return StatusCode::FAILURE;
    }
    const xAOD::MissingETAssociationMap* metMap(0);
    if ( m_event->retrieve(metMap, met_map_name_).isFailure() ) {
        Warning("Unable to retrieve MissingETAssociationMap: %s", met_map_name_.c_str());
        return StatusCode::FAILURE;
    }
    metMap->resetObjSelectionFlags();
		
    // bool doTracks = false;
    bool doJVTCut = false;
    std::string softTerm = "SoftClus";

    if(doTST) { 
        // use Track soft term
        doJVTCut = true;
        softTerm = "PVSoftTrk";
    }
		
    if(elec) {
        ConstDataVector<xAOD::ElectronContainer> metelectron(SG::VIEW_ELEMENTS);
        for (const auto& el : *elec) {
            if(! dec_signal(*el)) continue;
            metelectron.push_back(el);
        }
        met_maker_->rebuildMET(m_eleTerm, xAOD::Type::Electron, &met, metelectron.asDataVector(), metMap);
    } 

    if(gamma) {
        ConstDataVector<xAOD::PhotonContainer> metgamma(SG::VIEW_ELEMENTS);
        for (const auto& ph : *gamma) {
            metgamma.push_back(ph);
        }
        met_maker_->rebuildMET(m_gammaTerm, xAOD::Type::Photon, &met, metgamma.asDataVector(), metMap);
    } 

    if(taujet) {
        ConstDataVector<xAOD::TauJetContainer> mettau(SG::VIEW_ELEMENTS);
        for (const auto& tau : *taujet) {
            mettau.push_back(tau);
        }
        met_maker_->rebuildMET(m_tauTerm, xAOD::Type::Tau, &met, mettau.asDataVector(), metMap);
    } 

    if(muon) {
        ConstDataVector<xAOD::MuonContainer> metmuon(SG::VIEW_ELEMENTS);
        for(const auto& mu : *muon) {
            if(! dec_signal(*mu)) continue;
            metmuon.push_back(mu);
        }
        met_maker_->rebuildMET(m_muonTerm, xAOD::Type::Muon, &met, metmuon.asDataVector(), metMap);
    } 
    met_maker_->rebuildJetMET(m_jetTerm, softTerm, &met, jet, metcore, metMap, doJVTCut);

    // if( m_metSystTool->applyCorrection(*met[softTerm]) != CP::CorrectionCode::Ok ) {
    //     Warning("GetMET", "Failed to apply MET soft term systematics.");
    // }
    met_maker_->buildMETSum(m_outMETTerm, &met, met[softTerm]->source());
    return EL::StatusCode::SUCCESS;
}
