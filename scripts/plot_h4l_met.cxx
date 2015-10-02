#include "/afs/cern.ch/user/x/xju/tool/AtlasStyle.C"
#include "/afs/cern.ch/user/x/xju/tool/AtlasUtils.C"
#include "/afs/cern.ch/user/x/xju/tool/loader.c"

#include <TString.h>
#include <TF1.h>
#include <RooTFnPdfBinding.h>

#include <iostream>
#include <stdio.h>
#include <stack>
#include <utility>
#include <TColor.h>

using namespace std;
using namespace RooFit;

const char* ch_name = "tree_incl_all";
//////////////////////////////
// addtional selections
//////////////////////////////
// TCut cut("weight");
TCut cut("weight * (m4l_unconstrained > 110 && m4l_unconstrained < 140)");
// TCut cut("weight * (m4l_unconstrained > 110 && m4l_unconstrained < 140 && MET_Final/Hpt <  3 && MET_Final/Hpt > 0.6)");
// TCut cut("weight * (m4l_unconstrained > 110 && m4l_unconstrained < 140 && MET_Final > 80 && MET_Final/Hpt <  3 && MET_Final/Hpt > 0.6)");
// TCut cut("weight * (m4l_unconstrained > 110 && m4l_unconstrained < 140 && TMath::ATan(MET_Final/Hpt) > 0.5 && TMath::ATan(MET_Final/Hpt) < 1.1)");

//////////////////////////////
// analysis inputs
//////////////////////////////
const char* f_sig_name = "combined_zpxx_mzp200_mx1.root";
const char* f_bkg1_name = "combined_ggH125.root";
const char* f_bkg2_name = "combined_VBFH125.root";
const char* f_bkg3_name = "combined_zz.root";
const char* f_bkg4_name = "combined_WH125.root";
const char* f_bkg5_name = "combined_zee.root";
const char* f_bkg6_name = "combined_zmumu.root";
const char* f_bkg7_name = "combined_ZH125_llvv.root";
const char* f_bkg8_name = "combined_ZH125.root";
const char* f_bkg9_name = "combined_ZH125_lvlv.root";

//////////////////////////////
// which discriminate to use in the stack plot
//////////////////////////////
// #define _DPHI_
// #define _FIVE_
// #define _TAN_
// #define _RATIO_

#if defined(_DPHI_)
TH1F* h_template = new TH1F("h_template", "template", 50, 0, 3.2);
string var_name("dphi_met_higgs");
string x_title("#Delta#phi(E_{T}^{miss},h)");
#elif defined(_RATIO_)
TH1F* h_template = new TH1F("h_template", "template", 50, 0, 5);
string var_name("MET_Final/Hpt");
string x_title("E_{T}^{miss}/p_{T}^{h} [GeV]");
#elif defined(_TAN_)
TH1F* h_template = new TH1F("h_template", "template", 15, 0, 1.57);
string var_name("TMath::ATan(MET_Final/Hpt)");
string x_title("tan^{-1}(E_{T}^{miss}/p_{T}^{h})");
#else
TH1F* h_template = new TH1F("h_template", "template", 30, 0, 300);
string var_name("MET_Final");
string x_title("E_{T}^{miss} [GeV]");
#endif



Int_t font=42; // Helvetica
float t_size=0.04;

////////////
// define the color
//////////////
int color_signal = kRed;
int color_s2 = kRed+3;
int color_ZZ = kAzure+2;
int color_ggH = kOrange-2;
int color_VBF = kGreen+1;
int color_WH = kOrange+6;
int color_Zee = kAzure+5;
int color_Zmumu = kAzure+6;
// int color_ZH_llvv = kSpring-5;
int color_ZH_llvv = kViolet-4;
int color_ZH = kViolet-3;
int color_ZH_lvlv = kViolet-9;
// kAzure+4/5/6, kSpring+9, kYellow-9, kCyan-10, kVoilet-3, kOrange
// 206, 64, 95, 28, 209


// utilized functions
TH1F* add_sample(const char* file_name, const char* hist_name, int color, TCut cut_ );
void compare_full_atlfast(const char* var1, const char* x_title_, 
        int nbins_var1, float min_var1, float max_var1);
void shape_var(const char* var, const char* x_title_,
        int var_nbins, float var_min, float var_max);
TH1F* create_met(const char* fname, const TCut& cut, 
        const char* hist_name, int color);
void print_after_cut(const string& name, TH1F* h1, int cutbin);

// main functions
/* stacked backgrounds compare with signal
 * draw ROC to find optimized cut on the specific variable
 * */
void compare_h4l_met();
/* comapre full-sim sample with atl-fast-sim samples
 * 25 ns vs 50 ns
 * */
void compare_sim();
/* 
 * distributions are normalized to 1.
 * */
void shape();
/* compare MET distribution for different z' mass
 * */
void compare_zp_met();
void study_zz_met();
/* Fit the lineshape for ZZ and Z+jets
 * */
void fit_lineshape();

double get_sig_from_2d(TH2F* h_signal, TH2F* h_bkg, int x_low_bin, int y_low_bin, double& y_max_opt, double& cut_eff);
void optimize_met_hpt();

void compare_h4l_met()
{
    SetAtlasStyle();
    h_template->UseCurrentStyle();
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->SetLogy();
    

    TH1F* h_signal = add_sample(f_sig_name, "h_signal", color_signal, cut);
    TH1F* h_bkg1 = add_sample(f_bkg1_name, "h_bkg1", color_ggH, cut);
    TH1F* h_bkg2 = add_sample(f_bkg2_name, "h_bkg2", color_VBF, cut);
    TH1F* h_bkg3_raw = add_sample(f_bkg3_name, "h_bkg3_raw", color_ZZ, cut);
    TH1F* h_bkg4 = add_sample(f_bkg4_name, "h_bkg4", color_WH, cut);
    TH1F* h_bkg5_raw = add_sample(f_bkg5_name, "h_bkg5_raw", color_Zee, cut);
    TH1F* h_bkg6_raw = add_sample(f_bkg6_name, "h_bkg6_raw", color_Zmumu, cut);
    TH1F* h_bkg7 = add_sample(f_bkg7_name, "h_bkg7", color_ZH_llvv, cut);
    TH1F* h_bkg8 = add_sample(f_bkg8_name, "h_bkg8", color_ZH, cut);
    TH1F* h_bkg9 = add_sample(f_bkg9_name, "h_bkg9", color_ZH_lvlv, cut);

    // sampling low-stats samples
#ifndef _FIT_
    TH1F* h_bkg3 = generate_th(h_bkg3_raw, "h_bkg3", 10);
    TH1F* h_bkg5 = generate_th(h_bkg5_raw, "h_bkg5", 30);
    TH1F* h_bkg6 = generate_th(h_bkg6_raw, "h_bkg6", 30);
    canvas->cd();
#else

    TH1F* h_bkg3 = h_bkg3_raw;
    TH1F* h_bkg5 = h_bkg5_raw;
    TH1F* h_bkg6 = h_bkg6_raw;
#endif

    h_signal->SetFillColor(0);
    float lumi = 10e3; // 5/fb
#ifdef _FIVE_
    h_signal->Scale(lumi*1.25e-4); // lumi times BR
#else
    h_signal->Scale(lumi*1.25e-4); // lumi times BR
#endif
    h_bkg1->Scale(lumi);
    h_bkg2->Scale(lumi);
    h_bkg3->Scale(lumi);
    h_bkg4->Scale(lumi);
    h_bkg5->Scale(lumi);
    h_bkg6->Scale(lumi);
    h_bkg7->Scale(lumi);
    h_bkg8->Scale(lumi);
    h_bkg9->Scale(lumi);

    // summary of events
    printf("signal: %.2f\n", h_signal->Integral());
    printf("ggF: %.2f\n", h_bkg1->Integral());
    printf("VBF: %.2f\n", h_bkg2->Integral());
    printf("WH: %.2f\n", h_bkg4->Integral());
    printf("ZH: %.2f\n", h_bkg8->Integral());
    printf("ZH(llvv): %.2f\n", h_bkg7->Integral());
    printf("ZH(lvlv): %.2f\n", h_bkg9->Integral());
    printf("ZZ: %.2f\n", h_bkg3->Integral());
    printf("Zee: %.2f\n",   h_bkg5->Integral());
    printf("Zmumu: %.2f\n", h_bkg6->Integral());
    printf("---------------------------\n");
    printf("signal: %.f\n", h_signal->GetEntries());
    printf("ggF: %.f\n", h_bkg1->GetEntries());
    printf("VBF: %.f\n", h_bkg2->GetEntries());
    printf("WH: %.f\n", h_bkg4->GetEntries());
    printf("ZH: %.f\n", h_bkg8->GetEntries());
    printf("ZH(llvv): %.f\n", h_bkg7->GetEntries());
    printf("ZH(lvlv): %.f\n", h_bkg9->GetEntries());
    printf("ZZ: %.f\n",    h_bkg3->GetEntries());
    printf("Zee: %.f\n",   h_bkg5->GetEntries());
    printf("Zmumu: %.f\n", h_bkg6->GetEntries());
    printf("---------------------------\n");

    TList* objarray = new TList();
    for(int i = 1; i < 10; i ++){
        add_hist(objarray, Form("h_bkg%d",i));   
    }

    THStack* hs = new THStack("hs", "");
    for(int i=0; i < objarray->GetEntries(); i++){
        TH1F* h1 = (TH1F*)objarray->At(i);
        hs->Add(h1);
    }

#ifdef _FIVE_
    hs->Add(h_signal);
#endif

    float max_y = h_signal->GetMaximum()>hs->GetMaximum()?h_signal->GetMaximum():hs->GetMaximum();

    h_signal->SetXTitle(x_title.c_str());
    /**
    h_signal->GetXaxis()->SetTitleOffset(1.4);
    h_signal->GetXaxis()->SetTitleFont(font);
    h_signal->GetXaxis()->SetTitleSize(0.05);
    **/
    h_signal->GetYaxis()->SetRangeUser(1e-5, max_y*1e3);
    // h_signal->GetYaxis()->SetRangeUser(1e-6, max_y*1.1);
    h_signal->Draw("AXIS");
    hs->Draw("HISTsame");
#ifndef _FIVE_
    h_signal->Draw("AXISsame");
    h_signal->Draw("same");
#endif

    float x_offset = 0.75;
    float y_offset = 0.85; 
    float x_off_title = 0.2;
    if(h_signal->GetMaximumBin() > h_signal->GetNbinsX()/2.0) {
        x_offset = 0.30;
        x_off_title = 0.47;
    }
#ifdef _FIVE_
    myLineText(x_offset, y_offset, *h_signal, "signal #times 5");
#else
    myLineText(x_offset, y_offset, *h_signal, "signal");
#endif
    int item = 1;
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg1, "ggF");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg2, "VBF");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg4, "WH");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg8, "ZH");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg7, "ZH(llvv)");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg9, "ZH(lvlv)");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg3, "ZZ");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg5, "Zee");
    myLineText(x_offset, y_offset-t_size*item++, *h_bkg6, "Z#mu#mu");

    myText(x_off_title, 0.85, 1, "ATLAS Internal");
    myText(x_off_title, 0.80, 1, "m_{zp} = 200 GeV, m_{#chi} = 1 GeV");

    // do a roc scan
    TH1F* h_bkg = (TH1F*) h_bkg1->Clone("h_bkg");
    h_bkg->Add(h_bkg2);
    h_bkg->Add(h_bkg3);
    int cutbin = get_roc(h_signal, h_bkg, false);
    // summary with best cuts
    int tbins = h_signal->GetNbinsX();
    cout << "After cut of " << h_signal->GetBinLowEdge(cutbin) << endl;
    print_after_cut("signal", h_signal, cutbin);
    print_after_cut("ggF", h_bkg1, cutbin);
    print_after_cut("VBF", h_bkg2, cutbin);
    print_after_cut("WH", h_bkg4, cutbin);
    print_after_cut("ZH", h_bkg8, cutbin);
    print_after_cut("ZH(llvv)", h_bkg7, cutbin);
    print_after_cut("ZH(lvlv)", h_bkg9, cutbin);
    print_after_cut("ZZ", h_bkg3, cutbin);
    print_after_cut("Zee", h_bkg5, cutbin);
    print_after_cut("Zmm", h_bkg6, cutbin);
}

void compare_full_atlfast(const char* var1, const char* x_title_, 
        int nbins_var1, float min_var1, float max_var1, bool do_ratio)
{
#define _SIM_
    SetAtlasStyle();
    TCanvas* canvas = new TCanvas(Form("canvas_%s", var1), "canvas", 600, 600);

    const char* input_base_dir = "/afs/cern.ch/user/x/xju/work/h4l/minitrees/process";
    const char* f_sig1 = Form("%s/mini_fullsim_25ns_new.root", input_base_dir);
    const char* f_sig2 = Form("%s/mini_atlfast_25ns_new.root", input_base_dir);
    const char* f_sig3 = Form("%s/mini_fullsim_50ns_new.root", input_base_dir);
    const char* f_sig4 = Form("%s/mini_atlfast_50ns_new.root", input_base_dir);

    TH1F* h_s1 = create_hist(f_sig1, ch_name, var1, "",
            "h_s1", nbins_var1, min_var1, max_var1, 2);
    TH1F* h_s2 = create_hist(f_sig2, ch_name, var1, "", 
            "h_s2", nbins_var1, min_var1, max_var1, 4);
    TH1F* h_s3 = create_hist(f_sig3, ch_name, var1, "", 
            "h_s3", nbins_var1, min_var1, max_var1, 4);
    TH1F* h_s4 = create_hist(f_sig4, ch_name, var1, "", 
            "h_s4", nbins_var1, min_var1, max_var1, 4);

    
    h_s1->Sumw2();
    h_s2->Sumw2();
    h_s3->Sumw2();
    h_s4->Sumw2();

    h_s1->Scale(1./h_s1->Integral());
    h_s2->Scale(1./h_s2->Integral());
    h_s3->Scale(1./h_s3->Integral());
    h_s4->Scale(1./h_s4->Integral());

    float max_y = h_s1->GetMaximum()>h_s2->GetMaximum()?h_s1->GetMaximum():h_s2->GetMaximum();

    h_s1->SetXTitle(x_title_);
    h_s1->GetYaxis()->SetRangeUser(1e-6, max_y*1.4);
    
    if(do_ratio){
        h_s1->GetXaxis()->SetLabelSize(0);

        TH1F* h_copy = (TH1F*) h_s1->Clone("ratio");
#ifdef _SIM_
        add_ratio_pad(h_s1, h_s2);
#else
        add_ratio_pad(h_s1, h_s3);
#endif

    }
    //h_s3->SetLineStyle(2);
    h_s4->SetLineStyle(2);
    h_s1->Draw("HISTE");
    //h_s2->Draw("sameHIST");
#ifdef _SIM_
    h_s2->Draw("sameHISTE");
#else
    h_s3->Draw("sameHISTE");
#endif
    //h_s4->Draw("sameHIST");
    const char* tagname = "bs";

    float x_offset = 0.70;
    float y_offset = 0.85;
    float x_off_title = 0.2;
    if(h_s1->GetMaximumBin() > nbins_var1/2.0){ 
        x_offset = 0.27;
        x_off_title = 0.52;
    }
    myLineText(x_offset, y_offset, *h_s1, "Full Sim, 25 ns");
    // myLineText(x_offset, y_offset-t_size,  *h_s2, "AtlFast, 25 ns");
#ifdef _SIM_
    myLineText(x_offset, y_offset-t_size*2, *h_s2, "AtlFast, 25 ns");
    tagname = "sim";
#else
    myLineText(x_offset, y_offset-t_size*2, *h_s3, "Full Sim, 50 ns");
#endif
    //myLineText(x_offset, y_offset-t_size*3, *h_s4, "AtlFast, 50 ns");

    myText(x_off_title, 0.85, 1, "ATLAS Simulation, 13 TeV");
    myText(x_off_title, 0.80, 1, "m_{zp} = 200 GeV, m_{#chi} = 1 GeV");
        

    if(do_ratio){
        canvas->SaveAs(Form("eps/%s_%s_ratio.eps", var1, tagname));
    } else {
        canvas->SaveAs(Form("eps/%s_%s.eps", var1, tagname));
    }
    // add ratio
}

void compare_sim()
{
    // compare_full_atlfast("averageIPC", "<#mu>", 30, 0, 55, do_ratio); 
    bool do_ratio = true;
    compare_full_atlfast("MET_Final", "E_{T}^{miss} [GeV]", 50, 0, 600, do_ratio); 
    compare_full_atlfast("m4l_unconstrained", "m_{4l}", 30, 110, 140, do_ratio); 
    compare_full_atlfast("n_good_jets", "N_{jets}", 11, -0.5, 10.5, do_ratio); 
    compare_full_atlfast("Hpt", "p_{T}^{h}", 50, 0, 600, do_ratio); 
    // compare_full_atlfast("n_jets", "N_{jets}", 31, -0.5, 30.5, do_ratio); 

    int nbins_pt = 50;
    float low_pt = 0, high_pt = 250;
    compare_full_atlfast("Z1_lepplus_pt", "p_{T}^{Z_{1}^{+}} [GeV]", nbins_pt, low_pt, high_pt, do_ratio); 
    compare_full_atlfast("Z1_lepminus_pt", "p_{T}^{Z_{1}^{-} [GeV]", nbins_pt, low_pt, high_pt, do_ratio); 
    compare_full_atlfast("Z2_lepplus_pt", "p_{T}^{Z_{2}^{+}} [GeV]", nbins_pt, low_pt, high_pt, do_ratio); 
    compare_full_atlfast("Z2_lepminus_pt", "p_{T}^{Z_{2}^{-} [GeV]", nbins_pt, low_pt, high_pt, do_ratio); 

    compare_full_atlfast("mZ1_unconstrained", "m_{Z_{1}}", 70, 40, 110, do_ratio); 
    compare_full_atlfast("mZ2_unconstrained", "m_{Z_{2}}", 70, 10, 80, do_ratio); 
/***
    // eta
    int nbins_eta = 50;
    float low_eta = -3, high_eta = 3;
    compare_full_atlfast("Z1_lepplus_eta", "#eta Z1 +", nbins_eta, low_eta, high_eta, do_ratio); 
    compare_full_atlfast("Z1_lepminus_eta", "#eta Z1 -", nbins_eta, low_eta, high_eta, do_ratio); 
    compare_full_atlfast("Z2_lepplus_eta", "#eta Z2 +", nbins_eta, low_eta, high_eta, do_ratio); 
    compare_full_atlfast("Z2_lepminus_eta", "#eta Z2 -", nbins_eta, low_eta, high_eta, do_ratio); 
    
    // phi
    int nbins_phi = 50;
    float low_phi = -3, high_phi = 3;
    compare_full_atlfast("Z1_lepplus_phi", "#phi Z1 +", nbins_phi, low_phi, high_phi, do_ratio); 
    compare_full_atlfast("Z1_lepminus_phi", "#phi Z1 -", nbins_phi, low_phi, high_phi, do_ratio); 
    compare_full_atlfast("Z2_lepplus_phi", "#phi Z2 +", nbins_phi, low_phi, high_phi, do_ratio); 
    compare_full_atlfast("Z2_lepminus_phi", "#phi Z2 -", nbins_phi, low_phi, high_phi, do_ratio); 
**/
}

void shape_var(const char* var, const char* x_title_,
        int var_nbins, float var_min, float var_max)
{
    SetAtlasStyle();
    TCanvas* canvas = new TCanvas(Form("canvas_%s",var), "canvas", 600, 600);

    TH1F* h_signal = create_hist(f_sig_name, ch_name, var, cut,
            "h_signal", var_nbins, var_min, var_max, color_signal);
    TH1F* h_bkg1 = create_hist(f_bkg1_name, ch_name, var, cut,
            "h_bkg1", var_nbins, var_min, var_max, color_ggH);
    TH1F* h_bkg2 = create_hist(f_bkg2_name, ch_name, var, cut,
            "h_bkg2", var_nbins, var_min, var_max, color_VBF);
    TH1F* h_bkg3 = create_hist(f_bkg3_name, ch_name, var, cut,
            "h_bkg3", var_nbins, var_min, var_max, color_ZZ);
    TH1F* h_bkg4 = create_hist(f_bkg4_name, ch_name, var, cut,
            "h_bkg4", var_nbins, var_min, var_max, color_WH);
    TH1F* h_bkg5 = create_hist(f_bkg5_name, ch_name, var, cut,
            "h_bkg5", var_nbins, var_min, var_max, color_Zee);
    TH1F* h_bkg6 = create_hist(f_bkg6_name, ch_name, var, cut,
            "h_bkg6", var_nbins, var_min, var_max, color_Zmumu);
    TH1F* h_bkg7 = create_hist(f_bkg7_name, ch_name, var, cut,
            "h_bkg7", var_nbins, var_min, var_max, color_ZH_llvv);
    TH1F* h_bkg8 = create_hist(f_bkg8_name, ch_name, var, cut,
            "h_bkg8", var_nbins, var_min, var_max, color_ZH);

    TH1F* h_bkg = (TH1F*) h_bkg1->Clone("h_bkg");
    h_bkg->Add(h_bkg2);
    h_bkg->Add(h_bkg3);
    h_bkg->Add(h_bkg4);
    h_bkg->Add(h_bkg5);
    h_bkg->Add(h_bkg6);
    h_bkg->Add(h_bkg7);
    h_bkg->Add(h_bkg8);

    h_signal->Scale(1./h_signal->Integral());
    h_bkg1->Scale(1./h_bkg1->Integral());
    h_bkg2->Scale(1./h_bkg2->Integral());
    h_bkg3->Scale(1./h_bkg3->Integral());
    h_bkg4->Scale(1./h_bkg4->Integral());
    h_bkg5->Scale(1./h_bkg5->Integral());
    h_bkg6->Scale(1./h_bkg6->Integral());
    h_bkg7->Scale(1./h_bkg7->Integral());
    h_bkg8->Scale(1./h_bkg8->Integral());

    THStack* hs = new THStack("hs", "");
    hs->Add(h_bkg4);
    hs->Add(h_bkg2);
    hs->Add(h_bkg5);
    hs->Add(h_bkg6);
    hs->Add(h_bkg1);
    hs->Add(h_bkg3);
    hs->Add(h_bkg7);
    hs->Add(h_bkg8);
    hs->Add(h_signal);

    hs->Draw("nostackHIST");
    hs->GetXaxis()->SetTitle(x_title_);

    float x_offset = 0.75;
    float y_offset = 0.85;
    myLineText(x_offset, y_offset, *h_signal, "signal");
    myLineText(x_offset, y_offset-t_size, *h_bkg1, "ggF");
    myLineText(x_offset, y_offset-t_size*2, *h_bkg2, "VBF");
    myLineText(x_offset, y_offset-t_size*3, *h_bkg4, "WH");
    myLineText(x_offset, y_offset-t_size*4, *h_bkg8, "ZH");
    myLineText(x_offset, y_offset-t_size*5, *h_bkg7, "ZH(llvv)");
    myLineText(x_offset, y_offset-t_size*6, *h_bkg3, "ZZ");
    myLineText(x_offset, y_offset-t_size*7, *h_bkg5, "Zee");
    myLineText(x_offset, y_offset-t_size*8, *h_bkg6, "Z#mu#mu");

    myText(0.2, 0.85, 1, "ATLAS Internal");
    myText(0.2, 0.80, 1, "m_{zp} = 200 GeV, m_{#chi} = 1 GeV");
    
    // save to a histogram
    TFile* out_file = TFile::Open("map_weight.root","recreate");
    h_signal->Write();
    h_bkg->Scale(1./h_bkg->Integral());
    h_bkg->Write();
    TH1F* h_sob = (TH1F*) h_signal->Clone("h_sob");
    h_sob->Divide(h_bkg);
    h_sob->Scale(1./h_sob->Integral());
    h_sob->Write();
    out_file->Close();
}

void shape()
{
   // shape_var("MET_Final/Hpt", "MET/h_{pt}", 50, 0, 4);
   // shape_var("dphi_met_higgs", "#Delta#phi(MET,h)", 50, 0, 3.2);
   // shape_var("MET_Final", "MET [GeV]", 30, 0, 100);
   shape_var("TMath::ATan(MET_Final/Hpt)", "tan^{-1}(E_{T}^{miss}/p_{T}^h)", 30, 0, 1.57);
}

TH1F* create_met(const char* fname, const TCut& cut, const char* hist_name, int color)
{
    int bins_met = 100;
    float low_met = 0, high_met = 600;
    const char* met_name = "MET_Final";
    TH1F* h_temp = create_hist(fname, ch_name, met_name,
            cut, hist_name, 
            bins_met, low_met, high_met, color);
    return h_temp;
}

void compare_zp_met(){
    const char* f_zp_all = "combined_zp_all.root";
    SetAtlasStyle();
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->SetLogy();

    TH1F* h_med10_x1 = create_met(f_zp_all, "weight*(mc_channel_number == 341744)", 
            "h_med10_x1", color_signal);
    TH1F* h_med200_x1 = create_met(f_zp_all, "weight*(mc_channel_number == 341748)", 
            "h_med200_x1", color_ZZ);
    TH1F* h_med300_x1 = create_met(f_zp_all, "weight*(mc_channel_number == 341749)", 
            "h_med300_x1", color_ggH);
    TH1F* h_med2000_x1 = create_met(f_zp_all, "weight*(mc_channel_number == 341752)", 
            "h_med2000_x1", color_VBF);
    TH1F* h_med5000_x1 = create_met(f_zp_all, "weight*(mc_channel_number == 341753)", 
            "h_med5000_x1", color_WH);

    THStack* hs = new THStack("hs", "");
    norm_hist(h_med10_x1);
    norm_hist(h_med200_x1);
    norm_hist(h_med300_x1);
    norm_hist(h_med2000_x1);
    norm_hist(h_med5000_x1);

    hs->Add(h_med10_x1);
    hs->Add(h_med200_x1);
    hs->Add(h_med300_x1);
    hs->Add(h_med2000_x1);
    hs->Add(h_med5000_x1);

    float x_offset = 0.68;
    float y_offset = 0.85; 
    float x_off_title = 0.3;
   
    float max_y = h_med10_x1->GetMaximum();
    h_med10_x1->GetYaxis()->SetRangeUser(1e-5, max_y*1e2);
    h_med10_x1->GetXaxis()->SetTitleOffset(1.4);
    h_med10_x1->GetXaxis()->SetTitleFont(font);
    h_med10_x1->GetXaxis()->SetTitleSize(0.05);
    h_med10_x1->GetXaxis()->SetNdivisions(8);
    h_med10_x1->SetXTitle("MET [GeV]");
    h_med10_x1->Draw("AXIS");
    hs->Draw("nostackHISTsame");
    h_med10_x1->Draw("AXISsame");

    myLineText(x_offset, y_offset, *h_med10_x1, "m_{zp} = 10 GeV");
    myLineText(x_offset, y_offset-t_size, *h_med200_x1, "m_{zp} = 200 GeV");
    myLineText(x_offset, y_offset-t_size*2, *h_med300_x1, "m_{zp} = 300 GeV");
    myLineText(x_offset, y_offset-t_size*3, *h_med2000_x1, "m_{zp} = 2 TeV");
    myLineText(x_offset, y_offset-t_size*4, *h_med5000_x1, "m_{zp} = 5 TeV");

    myText(x_off_title, 0.85, 1, "ATLAS Internal");
    myText(x_off_title, 0.80, 1, "m_{#chi} = 1 GeV");
}

void study_zz_met()
{
    SetAtlasStyle();
    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    canvas->SetLogy();

    TH1F* h_s1 = add_sample(f_bkg3_name, "h_s1", color_signal, "weight*(m4l_unconstrained < 140 && m4l_unconstrained > 110)");
    TH1F* h_s2 = add_sample(f_bkg3_name, "h_s2", color_ZZ, "weight*(m4l_unconstrained >= 140 && m4l_unconstrained <= 300)");
    TH1F* h_s3 = add_sample(f_bkg3_name, "h_s3", color_WH, "weight*(m4l_unconstrained > 300)");
    printf("h_s1: %.f\n", h_s1->GetEntries());
    printf("h_s2: %.f\n", h_s2->GetEntries());
    printf("h_s3: %.f\n", h_s3->GetEntries());

    norm_hist(h_s1);
    norm_hist(h_s2);
    norm_hist(h_s3);

    h_s1->SetFillColor(0);
    h_s2->SetFillColor(0);
    h_s3->SetFillColor(0);

    THStack* hs = new THStack("hs", "");
    hs->Add(h_s1);
    hs->Add(h_s2);
    hs->Add(h_s3);
    
    h_s1->GetYaxis()->SetRangeUser(1e-6, 1);
    h_s1->GetXaxis()->SetTitle(x_title.c_str());
    h_s1->Draw("AXIS");
    
    cout <<"KolmogorovTest: " << h_s1->KolmogorovTest(h_s2) << endl;
    cout <<"Chi2Test: " << h_s1->Chi2Test(h_s2, "WW") << endl;
    hs->Draw("nostackHISTEsame");

    float x_offset = 0.6;
    float y_offset = 0.85; 
    myLineText(x_offset, y_offset, *h_s1, "110 < m4l < 140");
    myLineText(x_offset, y_offset-t_size*1, *h_s2, "140 <= m4l <= 300");
    myLineText(x_offset, y_offset-t_size*2, *h_s3, "m4l > 300");
    float x_off_title = 0.25;
    myText(x_off_title, 0.85, 1, "ATLAS Internal");
    myText(x_off_title, 0.80, 1, "qq#rightarrowZZ");

    TCanvas* canvas_new = new TCanvas("canvas_new", "canvas", 600, 600);
    canvas_new ->SetLogy();
    h_s1->GetYaxis()->SetRangeUser(1e-6, 10);
    h_s1->GetXaxis()->SetRangeUser(10, 150);
    h_s1->Draw();
    TH1F* h1 = generate_th(h_s1, "h1", 10);
    h1->SetLineStyle(2);
    h1->SetLineColor(4);
    h1->Draw("same");
    myLineText(x_offset, y_offset, *h_s1, "raw ZZ");
    myLineText(x_offset, y_offset-t_size*1, *h1, "fitted ZZ");
    myText(x_off_title, 0.85, 1, "ATLAS Internal");
    myText(x_off_title, 0.80, 1, "qq#rightarrowZZ");
}

TH1F* add_sample(const char* file_name, const char* hist_name, int color, TCut cut_)
{
    TH1F* h_tmp = (TH1F*) h_template->Clone(hist_name);
    h_tmp->Sumw2();
    TChain* chain = loader(file_name, ch_name);
    chain->Draw(Form("%s >> %s", var_name.c_str(), hist_name), cut_);
    // h_tmp->SetDirectory(0);
    delete chain;
    h_tmp->SetLineColor(color);
    h_tmp->SetFillColor(color);
    h_tmp->SetLineWidth(2);
    return h_tmp;
}


void fit_lineshape(const char* input, float low_value, const char* name)
{
    SetAtlasStyle();
    TCanvas* canvas = new TCanvas(Form("canvas_%s",name), "canvas", 600, 600);
    canvas->SetLogy();
    
    TH1F* h_s1 = add_sample(input, Form("h_%s",name), color_signal, cut);
    norm_hist(h_s1);
    string str_name(name);

    h_s1->GetYaxis()->SetRangeUser(1e-6, 10);
    h_s1->GetXaxis()->SetRangeUser(0, 200);
    h_s1->SetFillColor(0);
    h_s1->GetXaxis()->SetTitleOffset(1.4);
    h_s1->GetXaxis()->SetTitleFont(font);
    h_s1->GetXaxis()->SetTitleSize(0.05);
    h_s1->SetXTitle(x_title.c_str());
    h_s1->Draw();
    TH1F* h1 = generate_th(h_s1, Form("hfit_%s",name), low_value);
    h1->SetLineStyle(2);
    h1->SetLineColor(4);
    h1->Draw("same");
    if(false){
        const char* input_cr = nullptr;
        if(str_name == "zee"){
            input_cr = "combined_zee_CR.root";
        } else if(str_name == "zmm") {
            input_cr = "combined_zmumu_CR.root";
        }else {
        }
        TH1F* h_s1_cr = add_sample(input_cr, Form("h_%s_cr",name), color_VBF, cut);
        int bin60 = h_s1->FindBin(60);
        h_s1_cr->Scale(h_s1->Integral(1, bin60)/h_s1_cr->Integral(1,bin60));
        h_s1_cr->SetLineWidth(2);
        h_s1->Draw();
        h1->Draw("same");
        h_s1_cr->Draw("same");
    }

    float x_offset = 0.6;
    float y_offset = 0.85; 
    float x_off_title = 0.25;
    float t_size = 0.05;
    myLineText(x_offset, y_offset, *h_s1, "raw", 0.05);
    myLineText(x_offset, y_offset-t_size*1, *h1, "fitted", 0.05);
    myText(x_off_title, 0.85, 1, "ATLAS Internal");
    myText(x_off_title, 0.80, 1, name);
    canvas->SaveAs(Form("eps/fitted_%s.eps", name));
    delete h1;
    delete h_s1;
    delete canvas;
}

void fit(){
    // fit_lineshape(f_bkg3_name, 10, "ZZ");
    // fit_lineshape(f_bkg5_name, 30, "zee");
    // fit_lineshape(f_bkg6_name, 30, "zmm");
    fit_lineshape(f_bkg5_name, 30, "zee");
    fit_lineshape("combined_zee_CR.root", 30, "zee_cr");
    fit_lineshape(f_bkg6_name, 30, "zmm");
    fit_lineshape("combined_zmumu_CR.root", 30, "zmm_cr");
}

TH2F* get_2d_met(const string& file_name, const string& hist_name, TCut cut_)
{
    TH2F* h_tmp = new TH2F(hist_name.c_str(),"h_tmp;E_{T}^{miss};E_{T}^{miss}/p_{T}^{h}", 150, 0, 300, 100, 0, 5);
    TChain* chain = loader(file_name.c_str(), ch_name); 
    chain->Draw(Form("%s:%s >> %s","MET_Final/Hpt","MET_Final", h_tmp->GetName()), cut_);
    delete chain;
    return h_tmp;
}

double ratio_met_hpt_up = 4.0;
double get_sig_from_2d(TH2F* h_signal, TH2F* h_bkg, int x_low_bin, int y_low_bin, double& y_max_opt, double& cut_eff)
{
    double stotal = h_signal->Integral();
    // int x_low_bin = h_signal->GetXaxis()->FindBin(x_cut);
    int x_up_bin = h_signal->GetXaxis()->GetNbins();
    // int y_low_bin = h_signal->GetYaxis()->FindBin(y_cut);
    int ybins = h_signal->GetYaxis()->GetNbins();
    /***
    double max_sig = -99;
    int y_opt_bin = -1;
    for(int ybin = y_low_bin+1; ybin <= ybins; ybin ++)
    {
        double s = h_signal->Integral(x_low_bin, x_up_bin, y_low_bin, ybin);
        double b = h_bkg->Integral(x_low_bin, x_up_bin, y_low_bin, ybin);
        double sig = get_significance(s, b);
        if(sig < 0) continue;
        if(sig > max_sig){
            max_sig = sig;
            y_opt_bin = ybin;
        }
    }
    y_max_opt = h_signal->GetYaxis()->GetBinLowEdge(y_opt_bin);
    return max_sig;
    */
    int ybin = h_signal->GetYaxis()->FindBin(ratio_met_hpt_up);
    y_max_opt = ratio_met_hpt_up;
    double s = h_signal->Integral(x_low_bin, x_up_bin, y_low_bin, ybin);
    double b = h_bkg->Integral(x_low_bin, x_up_bin, y_low_bin, ybin);
    cut_eff = s/stotal;
    return get_significance(s, b);
}

class OptInfo{
public:
    OptInfo(double sig_, double eff_, double xv_, double yv_):
        sig(sig_), eff(eff_), xv(xv_), yv(yv_)
    {
    }
    OptInfo(){}
    ~OptInfo(){}
    double sig;
    double eff;
    double xv;
    double yv;
};

void optimize_met_hpt()
{
    // SetAtlasStyle();
    h_signal = get_2d_met(f_sig_name, "h_signal",cut);
    h_bkg1 = get_2d_met(f_bkg1_name, "h_bkg1",cut);
    h_bkg2 = get_2d_met(f_bkg2_name, "h_bkg2",cut);
    h_bkg3 = get_2d_met(f_bkg3_name, "h_bkg3",cut);
    h_bkg4 = get_2d_met(f_bkg4_name, "h_bkg4",cut);
    h_bkg5 = get_2d_met(f_bkg5_name, "h_bkg5",cut);
    h_bkg6 = get_2d_met(f_bkg6_name, "h_bkg6",cut);
    h_bkg7 = get_2d_met(f_bkg7_name, "h_bkg7",cut);
    h_bkg8 = get_2d_met(f_bkg8_name, "h_bkg8",cut);
    h_bkg9 = get_2d_met(f_bkg9_name, "h_bkg9",cut);

    float lumi = 4e3; // 4 ifb
    h_signal->Scale(lumi*1.25e-4);
    h_bkg1->Scale(lumi);
    h_bkg2->Scale(lumi);
    h_bkg3->Scale(lumi);
    h_bkg4->Scale(lumi);
    h_bkg5->Scale(lumi);
    h_bkg6->Scale(lumi);
    h_bkg7->Scale(lumi);
    h_bkg8->Scale(lumi);
    h_bkg9->Scale(lumi);
    
    h_bkg_all = (TH2F*) h_bkg1->Clone("h_bkg_all");
    h_bkg_all->Add(h_bkg2);
    h_bkg_all->Add(h_bkg3);
    h_bkg_all->Add(h_bkg4);
    h_bkg_all->Add(h_bkg5);
    h_bkg_all->Add(h_bkg6);
    h_bkg_all->Add(h_bkg7);
    h_bkg_all->Add(h_bkg8);
    h_bkg_all->Add(h_bkg9);
    
    fout = TFile::Open(Form("opt_%.f.root", ratio_met_hpt_up), "RECREATE");
    h_sig = (TH2F*) h_signal->Clone("h_sig");
    h_max_y = (TH2F*) h_signal->Clone("h_max_y");
    h_eff = (TH2F*) h_signal->Clone("h_eff");
    stack<OptInfo> max_sigs;
    stack<OptInfo> temp;
    for(int idx = 1; idx < h_signal->GetNbinsX(); idx++)
    {
        for(int idy = 1; idy < h_signal->GetNbinsY(); idy++)
        {
            double y_max= 0; 
            double cut_eff = 0;
            double sig = get_sig_from_2d(h_signal, h_bkg_all, idx, idy, y_max, cut_eff);
            // printf("%d,%d: %.2f\n", idx, idy, sig); 
            if(sig < 0) sig = 0;
            h_sig->SetBinContent(idx, idy, sig);
            h_max_y->SetBinContent(idx, idy, y_max);
            h_eff->SetBinContent(idx, idy, cut_eff);

            double x_value = h_signal->GetXaxis()->GetBinLowEdge(idx);
            double y_value = h_signal->GetYaxis()->GetBinLowEdge(idy);
            OptInfo info(sig, cut_eff, x_value, y_value); 

            if(!max_sigs.empty() && (max_sigs.size() <= 10 || sig > max_sigs.top().sig)){
                while(!max_sigs.empty() && sig > max_sigs.top().sig)
                {
                    temp.push(max_sigs.top());
                    max_sigs.pop();
                }
                max_sigs.push(info);
                while(!temp.empty() && max_sigs.size() <= 10){
                    max_sigs.push(temp.top());
                    temp.pop();
                }
            } else {
                max_sigs.push(info);
            }
        }
    }
    canvas_sig = new TCanvas("canvas_sig", "significance", 600, 600);
    canvas_sig->cd();
    h_sig->Draw("colz");
    canvas_sig->SaveAs(Form("eps/opt_sig_%.f.eps", ratio_met_hpt_up));

    canvas_maxy = new TCanvas("canvas_maxy", "max Y", 600, 600);
    canvas_maxy->cd();
    h_max_y->Draw("colz");
    canvas_maxy->SaveAs(Form("eps/opt_maxy_%.f.eps", ratio_met_hpt_up));

    canvas_eff = new TCanvas("canvas_eff", "signal efficiency", 600, 600);
    canvas_eff->cd();
    h_eff->Draw("colz");
    canvas_eff->SaveAs(Form("eps/opt_eff_%.f.eps", ratio_met_hpt_up));
    fout->cd();
    h_sig->Write();
    h_eff->Write();
    h_max_y->Write();
    cout <<"Here is the top 10: " << endl;
    while(!max_sigs.empty()){
        OptInfo& info = max_sigs.top();
        printf("%.4f : %.4f : %.4f : %.4f \n", info.sig, info.eff, info.xv, info.yv);
        max_sigs.pop();
    }
}
