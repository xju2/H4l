#include <stdlib.h>
#include "/afs/cern.ch/user/x/xju/tool/loader.c"
#include <map>
#include <string>
#include <stdlib.h>

using namespace std;

float in_m4l_unconstrained, in_weight; 
int in_run, in_event;

float out_m4l_truth, out_m4l, out_weight;
int out_category, out_event, out_run;
map<int, float> channel_truth_map;

TTree* tree;
bool is13TeV = false;
void GetTruthMap(const char* input_name){
    fstream input(input_name, fstream::in);
    string file_name;
    while(input >> file_name){
        size_t pos_mc12 = file_name.find("mc12_8TeV");
        int cut_down = 167978;
        int nskip = 10;
        if(is13TeV){ 
            cut_down = 206357;
            nskip = 11;
        }
        string channel_str = file_name.substr(pos_mc12+nskip, 6);
        int channel = atoi(channel_str.c_str());
        size_t pos_ggh = file_name.find("ggH");
        size_t length = 3;
        if(channel >= cut_down) length = 4;
        string m4ltruth_str = file_name.substr(pos_ggh+3,length);
        float m4l_truth = atof(m4ltruth_str.c_str());
        channel_truth_map[channel] = m4l_truth;
    }
    input.close();
}

void FillNewTree(const char* input_name, const char* chain_name, int category)
{
    TChain* input_4mu = loader(input_name, chain_name);
    input_4mu->SetBranchAddress("run", &in_run);
    input_4mu->SetBranchAddress("event", &in_event);
    input_4mu->SetBranchAddress("m4l_unconstrained", &in_m4l_unconstrained);
    input_4mu->SetBranchAddress("weight", &in_weight);

    int nentries = input_4mu ->GetEntries();
    for(int ientry = 0; ientry < nentries; ientry ++)
    {
        input_4mu->GetEntry(ientry);
        out_run = in_run;
        out_event = in_event;
        out_m4l = in_m4l_unconstrained;
        out_m4l_truth = channel_truth_map[in_run];
        out_weight = in_weight;
        out_category = category;
        tree->Fill();
    }
    delete input_4mu;
}

void merge_mini_tree(const char* input_name, const char* out_name){
    // const char* input_name = "/afs/cern.ch/user/x/xju/work/h4l/workspace/parametrisation/mc12_ggH_mini_tree.txt";
    // const char* out_name = "hzz_all.root"; 
    GetTruthMap(input_name);
    TFile* out_file = TFile::Open(out_name, "recreate");
    //output tree
    tree = new TTree("physics", "physics");
    tree->Branch("run", &out_run, "run/I");
    tree->Branch("event", &out_event, "event/I");
    tree->Branch("m4l", &out_m4l, "m4l/F");
    tree->Branch("m4l_truth", &out_m4l_truth, "m4l_truth/F");
    tree->Branch("weight", &out_weight, "weight/F");
    tree->Branch("category", &out_category, "category/I");
    
    FillNewTree(input_name, "tree_incl_4mu", 1);
    FillNewTree(input_name, "tree_incl_2mu2e", 2);
    FillNewTree(input_name, "tree_incl_2e2mu", 3);
    FillNewTree(input_name, "tree_incl_4e", 4);

    tree->Write();
    out_file->Close();
    cout << out_name << " is created." << endl;
}
