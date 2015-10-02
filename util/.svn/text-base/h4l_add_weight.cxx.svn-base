#include <stdlib.h>
#include <string>
#include <iostream>

#include <TH1F.h>
#include <TFile.h>

#include "H4l/physics.h"
#include "MyXAODTools/EventCounter.h"

using namespace std;
int main( int argc, char** argv)
{
    if(argc < 3) {
        cout << argv[0] << " toberun.txt/test.root output.root" << endl;
        return 1;
    }
    // everything is MC
    string input_name(argv[1]) ;
    string output_name(argv[2]);
    TFile* file_map = TFile::Open("/afs/cern.ch/user/x/xju/work/monoH4l/shape/map_weight.root","read");
    TH1F* h_signal = nullptr;
    if(file_map->IsOpen() && !file_map->IsZombie()){
        h_signal = (TH1F*) file_map->Get("h_sob");
        h_signal->SetDirectory(0);
        file_map->Close();
    }

    auto* event_counter = new EventCounter();
    event_counter->GetTotalEventsDic(input_name.c_str(), "associate");
    event_counter->printTotalEvents();

    auto* fc = EventCounter::getChain(input_name.c_str(), "tree_incl_all");
    Long64_t nentries = fc->GetEntries();
    cout << "total entries: " << nentries << endl;

    physics* ntuple = new physics(fc);
    TFile* output = new TFile(output_name.c_str(), "recreate");
    TTree* out_tree = fc->CloneTree(0);
    float s_over_b;
    out_tree->Branch("soverb", &s_over_b, "soverb/F");

    // Loop over the tree
    for(Long64_t ientry = 0; ientry < nentries; ientry ++){
        ntuple->LoadTree(ientry);
        ntuple->GetEntry(ientry);
        int mc_id = ntuple->mc_channel_number;
        double weight = event_counter->getCrossSection(mc_id)/event_counter->getTotalEvents(mc_id)*ntuple->MCWeight;
        ntuple->weight = weight;
        /***
        float met_over_hpt = ntuple->MET_Final/ntuple->Hpt;
        if(h_signal != nullptr){
            int bin = h_signal->FindBin(met_over_hpt);
            s_over_b = h_signal->GetBinContent(bin);
        }
        ***/
        out_tree->Fill();
    }
    out_tree->AutoSave();
    output->Close();
    delete event_counter;
    delete ntuple;
    if(h_signal != nullptr){
        delete h_signal;
    }
    return 0;
}
