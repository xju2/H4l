#ifndef H4l_H4lAnalysis_H
#define H4l_H4lAnalysis_H

#include <EventLoop/Algorithm.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "METUtilities/METMaker.h"
#include "METUtilities/METSystematicsTool.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/PhotonAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODTau/TauJetContainer.h"

#include "TH1.h"
#include "H4l/MuonSelection.h"
#include "H4l/ElectronSelection.h"
#include "H4l/JetSelection.h"
#include "MyXAODTools/CPToolsHelper.h"
#include "MyXAODTools/EventInfoCreator.h"

#include <string>
#include <vector>

#ifndef __CINT__
#include "FsrUtils/FsrPhotonTool.h"
#endif
using namespace std;
namespace CP {
    class MuonCalibrationAndSmearingTool;
    class EgammaCalibrationAndSmearingTool;
    class MuonSelectionTool;
}

class JetCalibrationTool;

class AsgElectronLikelihoodTool;
class Minitree;

class H4lAnalysis : public EL::Algorithm {
    // put your configuration variables here as public variables.
    // that way they can be set directly from CINT and python.

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)

    public:
        H4lAnalysis(bool inclusive = false, bool smearing = false, 
                bool weights = false, bool doCR = false, 
                bool is_data = false, bool atlfast = false,
                bool debug = false);
        virtual ~H4lAnalysis();

        // these are the functions inherited from Algorithm
        virtual EL::StatusCode setupJob (EL::Job& job);
        virtual EL::StatusCode fileExecute ();
        virtual EL::StatusCode histInitialize ();
        virtual EL::StatusCode changeInput (bool firstFile);
        virtual EL::StatusCode initialize ();
        virtual EL::StatusCode execute ();
        virtual EL::StatusCode postExecute ();
        virtual EL::StatusCode finalize ();
        virtual EL::StatusCode histFinalize ();

        /////
        EL::StatusCode GetMET(xAOD::MissingETContainer &met,
				     const xAOD::JetContainer* jet,
				     const xAOD::ElectronContainer* elec,
				     const xAOD::MuonContainer* muon,
				     const xAOD::PhotonContainer* gamma,
				     const xAOD::TauJetContainer* taujet,
				     bool doTST) ;
         

        // this is needed to distribute the algorithm to the workers
        ClassDef(H4lAnalysis, 1);

        bool m_inclusive; 
        bool m_smearing;
        bool m_weights;
        bool m_doCR;
        bool m_is_data;
        bool m_debug;
        bool m_atlfast;

        uint64_t total_evts_pro_ = 0;
        double sum_of_evt_w_ = 0;
        double sum_of_evt_w_sq_ = 0;

        static const char* APP_NAME;
protected:
        xAOD::JetInput::Type m_jetInputType; // permit switching between LC, PFlow, EM jets
        string jet_base_name_;
        string jet_con_name_;
        string met_core_name_;
        string met_map_name_;


        xAOD::TEvent *m_event; //!
        xAOD::TStore *m_store; //!
        int m_eventCounter; //!
        uint32_t m_runNumber; //!
        uint32_t m_eventNumber; //!
        std::map<int, std::string> m_es_failCodes; //!
        std::map<int, std::string> m_ms_failCodes; //!
        std::map<int, std::string> m_js_failCodes; //!  
        std::map<int, std::string> m_qs_failCodes; //!
        TH1 *h_es_fail; //!
        TH1 *h_ms_fail; //!
        TH1 *h_js_fail; //!
        TH1 *h_qs_fail_4mu; //!
        TH1 *h_qs_fail_4e; //!
        TH1 *h_qs_fail_2mu2e; //!
        TH1 *h_qs_fail_2e2mu; //!
        std::map<int, TH1*> h_qs_fail_map; //!
        vector<string>* trig_el_; //!
        vector<string>* trig_mu_; //!
        vector<string>* trig_mu_el_; //!

        
#ifndef __CINT__
        CP::MuonCalibrationAndSmearingTool *m_muonCorrectionTool; //!
        CP::EgammaCalibrationAndSmearingTool *m_egammaCorrectionTool; //!
        JetCalibrationTool *m_jetCorrectionTool; //!
        FSR::FsrPhotonTool *m_fsrTool; //!
        Minitree *mt; //!
        met::METMaker* met_maker_; //!
        MuonSelection *MS_ ; //!
        ElectronSelection* ES_; //!
        JetSelection* JS_; //!
        CPToolsHelper* cp_tools_;  //!   
        EventInfoCreator* event_info_filler_; //!
        TTree* associate_tree_; //!
#endif
};

#endif
