#include <stdlib.h>
#include <string>
#include <iostream>
#include <algorithm>

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
    string input_name(argv[1]) ;
    string output_name(argv[2]);

    auto* fc = EventCounter::getChain(input_name.c_str(), "physics");
    vector<double>* mc_weights = new vector<double>();
    double met;
    fc->SetBranchAddress("MCWeights", &mc_weights);
    fc->SetBranchAddress("MET_et", &met);

    Long64_t nentries = fc->GetEntries();
    cout << "total entries: " << nentries << endl;

    TFile* output = new TFile(output_name.c_str(), "recreate");
    TTree* out_tree = fc->CloneTree(0);
    double scale_low, scale_up;
    double pdf_norm, pdf_low, pdf_up;
    double norm_weight;

    out_tree->Branch("weight", &norm_weight, "weight/D");
    out_tree->Branch("scaleUp", &scale_up, "scaleUp/D");
    out_tree->Branch("scaleDown", &scale_low, "scaleDown/D");

    out_tree->Branch("pdfNorm", &pdf_norm, "pdfNorm/D");
    out_tree->Branch("pdfUp", &pdf_up, "pdfUp/D");
    out_tree->Branch("pdfDown", &pdf_low, "pdfDown/D");

    int n_up = 0, n_low = 0, n_eq = 0;
    out_tree->Branch("nUp", &n_up, "nUp/I");
    out_tree->Branch("nDown", &n_low, "nDown/I");
    out_tree->Branch("nEq", &n_eq, "nEq/I");

    vector<int> qcd_index_l1;  // NNPDF30lo_as_0118
    qcd_index_l1.push_back(25);
    qcd_index_l1.push_back(36);
    qcd_index_l1.push_back(80);
    qcd_index_l1.push_back(91);
    qcd_index_l1.push_back(102);

    vector<int> qcd_index_l2; // MMHT2014lo68cl
    qcd_index_l2.push_back(11);
    qcd_index_l2.push_back(22);
    qcd_index_l2.push_back(60);
    qcd_index_l2.push_back(61);
    qcd_index_l2.push_back(62);

    vector<int> qcd_index_l3; // CT14lo
    qcd_index_l2.push_back(5);
    qcd_index_l2.push_back(6);
    qcd_index_l2.push_back(10);
    qcd_index_l2.push_back(11);
    qcd_index_l2.push_back(12);

    vector<int> qcd_index_l4; // CT10
    qcd_index_l4.push_back(11);
    qcd_index_l4.push_back(22);
    qcd_index_l4.push_back(62);
    qcd_index_l4.push_back(63);
    qcd_index_l4.push_back(64);
    

    vector<int>* qcd_exclusion;
    fc->LoadTree(0);
    fc->GetEntry(0);
    bool is_pdf_ok = true;
    if(mc_weights->size() == 113) {
        qcd_exclusion = &qcd_index_l1;
    } else if(mc_weights->size() == 63) {
        qcd_exclusion = &qcd_index_l2;
    } else if(mc_weights->size() == 13) {
        qcd_exclusion = &qcd_index_l3;
    } else if(mc_weights->size() == 65) {
        qcd_exclusion = &qcd_index_l4;
    } else { 
        cout << "unknown pdfset!" << endl;
        is_pdf_ok = false; 
    }
    if(is_pdf_ok){
    for(Long64_t ientry = 0; ientry < nentries; ientry ++){
        fc->LoadTree(ientry);
        fc->GetEntry(ientry);
        
        norm_weight = mc_weights->at(1);
        double down_scale = mc_weights->at(0);
        double up_scale = mc_weights->at(qcd_exclusion->at(4));
        if(down_scale >= norm_weight && up_scale <= norm_weight){
            scale_low = up_scale;
            scale_up = down_scale;
        } else {
            scale_low = down_scale;
            scale_up = up_scale;
        }

        pdf_norm = 0;
        int n_varies = 0;
        for(int i = 2; i < (int)mc_weights->size(); i++){
            if(find(qcd_exclusion->begin(), qcd_exclusion->end(), i) != qcd_exclusion->end()) continue;
            pdf_norm += mc_weights->at(i); 
            n_varies += 1;
        }
        pdf_norm /= n_varies;

        if(mc_weights->size() == 113) { //NNPDF
            double sum_diff = 0;
            n_varies = 0;
            for(int i = 2; i < (int)mc_weights->size(); i++){
                if(find(qcd_exclusion->begin(), qcd_exclusion->end(), i) != qcd_exclusion->end()) continue;
                double diff = mc_weights->at(i) - norm_weight;
                sum_diff += diff*diff;
                n_varies ++;
            }
            pdf_up = norm_weight + sqrt(sum_diff)/(n_varies-1);
            pdf_low = norm_weight - sqrt(sum_diff)/(n_varies-1);
        } else { // CT10 and MMHT2014
            pdf_up = pdf_norm;
            pdf_low = pdf_norm;
            double pdf_up_sum = 0, pdf_low_sum = 0; 
            n_up = 0, n_low = 0, n_eq = 0;
            for(int i = 2; i < (int)mc_weights->size(); i++){
                if(find(qcd_exclusion->begin(), qcd_exclusion->end(), i) != qcd_exclusion->end()) continue;
                double diff = mc_weights->at(i) - norm_weight;
                if(diff > 0) {
                    pdf_up_sum += diff*diff;
                    n_up ++;
                }
                if(diff < 0) {
                    pdf_low_sum += diff*diff;
                    n_low ++;
                }
                if(diff == 0) n_eq ++;
            }
            pdf_up = norm_weight + sqrt(pdf_up_sum);
            pdf_low = norm_weight - sqrt(pdf_low_sum);
        }
        out_tree->Fill();
    }
    }
    out_tree->AutoSave();
    output->Close();
    delete fc;

    return 0;
}

