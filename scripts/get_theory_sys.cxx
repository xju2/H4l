#include <stdlib.h>
#include "/afs/cern.ch/user/x/xju/tool/loader.c"
#include "/afs/cern.ch/user/x/xju/tool/AtlasStyle.C"
#include "/afs/cern.ch/user/x/xju/tool/AtlasUtils.C"

#include <TColor.h>

// #define _RATIO_
using namespace std;

double met_xbins[] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 600};
int met_nbins = 15;
TH1F* h_temp = new TH1F("h_temp","MET", met_nbins, met_xbins);
// string var_name = "MET_et";  // MET_et, m4l, Hpt

void get_pdf_sys(const char* file_name, TH1* h_nominal, TH1* h_up, TH1* h_down, int color, const string& var_name = "MET_et");
TH1F* get_var_cmp(const string&  var_name, const string& hist_name);

void get_theory_sys(const char* input_name = "new_8_MMHT2014.root")
{
    // h_temp->SetLineWidth(2);
    // SetXTitle(h_temp, "E_{T}^{miss} [GeV]");

    SetAtlasStyle();
    h_temp->UseCurrentStyle();

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
    // canvas->SetLogy();
    TChain* chain = loader(input_name);

    TH1F* h_norm = (TH1F*) h_temp->Clone("h_norm");
    TH1F* h_scale_up = (TH1F*) h_temp->Clone("h_scale_up");
    TH1F* h_scale_down = (TH1F*) h_temp->Clone("h_scale_down");

    TH1F* h_pdf_up = (TH1F*) h_temp->Clone("h_pdf_up");
    TH1F* h_pdf_down = (TH1F*) h_temp->Clone("h_pdf_down");
    const string var_name = "MET_et";
    chain->Draw(Form("%s/1e3>>h_norm", var_name.c_str()), "weight");
    chain->Draw(Form("%s/1e3>>h_scale_up", var_name.c_str()), "scaleUp");
    chain->Draw(Form("%s/1e3>>h_scale_down", var_name.c_str()), "scaleDown");
    chain->Draw(Form("%s/1e3>>h_pdf_up", var_name.c_str()), "pdfUp");
    chain->Draw(Form("%s/1e3>>h_pdf_down", var_name.c_str()), "pdfDown");
    h_scale_up->SetLineColor(kAzure+5);
    h_scale_down->SetLineColor(kAzure+5);
    // h_scale_down->SetLineStyle(2);

    h_pdf_up->SetLineColor(kOrange+5);
    h_pdf_down->SetLineColor(kOrange+5);
    // h_pdf_down->SetLineStyle(2);

    TH1F* h_ct10_nom = (TH1F*) h_temp->Clone("h_ct10_nom");
    TH1F* h_ct10_up = (TH1F*) h_temp->Clone("h_ct10_up");
    TH1F* h_ct10_down = (TH1F*) h_temp->Clone("h_ct10_down");
    get_pdf_sys("new_9_CT10.root", h_ct10_nom, h_ct10_up, h_ct10_down, kGreen+2);

    TH1F* h_nnpdf_nom = (TH1F*) h_temp->Clone("h_nnpdf_nom");
    TH1F* h_nnpdf_up = (TH1F*) h_temp->Clone("h_nnpdf_up");
    TH1F* h_nnpdf_down = (TH1F*) h_temp->Clone("h_nnpdf_down");
    get_pdf_sys("new_6_NNPDF30_0130.root", h_nnpdf_nom, h_nnpdf_up, h_nnpdf_down, kViolet-4);

    THStack* hs = new THStack("hs", "");
#ifdef _RATIO_

    h_scale_up->Divide(h_norm);
    h_scale_down->Divide(h_norm);
    h_pdf_up->Divide(h_norm);
    h_pdf_down->Divide(h_norm);
#else
    hs->Add(h_norm);
    // hs->Add(h_ct10_nom);
    // hs->Add(h_nnpdf_nom);
#endif

    hs->Add(h_scale_up);
    hs->Add(h_scale_down);

    hs->Add(h_pdf_up);
    hs->Add(h_pdf_down);

    hs->Add(h_ct10_up);
    hs->Add(h_ct10_down);
    hs->Add(h_nnpdf_up);
    hs->Add(h_nnpdf_down);

    
    h_temp->GetYaxis()->SetRangeUser(0.,0.55);
    h_temp->SetXTitle("E_{T}^{miss} [GeV]");
    h_temp->Draw("AXIS");
    hs->Draw("nostackHISTsame");
    // AddLine(h_temp, 1);

    legend = myLegend(0.5, 0.7, 0.85, 0.9);
    // legend->UseCurrentStyle();
    legend->AddEntry(h_scale_up, "QCD Scale", "l");
    legend->AddEntry(h_pdf_up, "MMHT2014LO68CL", "l");
    legend->AddEntry(h_ct10_up, "CT10", "l");
    legend->AddEntry(h_nnpdf_up, "NNPDF30LO", "l");
    legend->Draw();

    float x_off_title = 0.2;
    ATLAS(x_off_title, 0.85);
    myText(x_off_title, 0.80, 1, "m_{zp} = 200 GeV");
    myText(x_off_title, 0.75, 1, "m_{#chi} = 1 GeV");
    canvas->SaveAs("th_sys_new.eps");

    TCanvas* canvas_new = new TCanvas("canvas_new", "canvas", 600, 600);
    TH1F* h_nom_all = get_var_cmp("weight", "h_nom_all");
    TH1F* h_nom_up = get_var_cmp("pdfUp", "h_nom_up");
    TH1F* h_nom_down = get_var_cmp("pdfDown", "h_nom_down");
    TH1F* h_weight_up = get_var_cmp("scaleUp", "h_weight_up");
    TH1F* h_weight_down = get_var_cmp("scaleDown", "h_weight_down");

    h_nom_up->SetLineColor(kGreen+2);
    h_nom_up->SetMarkerColor(kGreen+2);
    h_nom_down->SetLineColor(kGreen+2);
    h_nom_down->SetMarkerColor(kGreen+2);

    // h_nom_down->SetLineColor(kOrange+5);
    // h_nom_down->SetMarkerColor(kOrange+5);

    h_weight_up->SetLineColor(kAzure+5);
    h_weight_down->SetLineColor(kAzure+5);

    double max_y = h_nom_up->GetMaximum();
    double min_y = h_nom_down->GetMinimum();
    h_nom_up->GetYaxis()->SetRangeUser(min_y*0.9, max_y*1.25);

    h_nom_up->GetXaxis()->SetBinLabel(1, "NNPDF30LO");
    h_nom_up->GetXaxis()->SetBinLabel(2, "CT10");
    h_nom_up->GetXaxis()->SetBinLabel(3, "MMHT2014LO68CL");
    h_nom_up->LabelsOption("h","X");
    h_nom_up->GetYaxis()->SetTitle("Weight per Event");

    h_nom_up->Draw("");
    h_nom_down->Draw("same");
    h_nom_all->Draw("same");
    h_weight_up->Draw("same");
    h_weight_down->Draw("same");
    AddLine(h_nom_all, h_nom_up->GetBinContent(2));
    AddLine(h_nom_all, h_nom_down->GetBinContent(3));
    AddLine(h_nom_all, h_nom_all->GetBinContent(1));

    leg2 = myLegend(0.65, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_nom_all, "Nominal", "l");
    leg2->AddEntry(h_nom_down, "PDF sys", "l");
    leg2->AddEntry(h_weight_up, "Scale sys", "l");
    leg2->Draw();
    ATLAS(x_off_title, 0.85);
    myText(x_off_title, 0.80, 1, "m_{zp} = 200 GeV");
    myText(x_off_title, 0.75, 1, "m_{#chi} = 1 GeV");
    double norm_value = h_nom_all->GetBinContent(1);
    double scale_up_val = h_weight_up->GetBinContent(1);
    double scale_down_val = h_weight_down->GetBinContent(1);
    myText(x_off_title, 0.70, 1, Form("PDF: + %.1f%% - %.1f%%", 
                100*(max_y-norm_value)/norm_value, 100*(norm_value-min_y)/norm_value));
    myText(x_off_title, 0.65, 1, Form("Scale: + %.1f%% - %.1f%%", 
                100*(scale_up_val-norm_value)/norm_value, 100*(norm_value-scale_down_val)/norm_value));
    canvas_new->SaveAs("th_sys_cmp.eps");
}

void get_pdf_sys(const char* file_name, TH1* h_nominal, TH1* h_up, TH1* h_down, int color, const string& var_name)
{
    TChain* chain = loader(file_name);
    // string norm_name(Form("norm_%s", h_up->GetName()));
    // TH1F* h_norm = (TH1F*) h_temp->Clone(norm_name.c_str());
    chain->Draw(Form("%s/1e3>>%s", var_name.c_str(), h_nominal->GetName()), "weight");
    chain->Draw(Form("%s/1e3>>%s", var_name.c_str(), h_up->GetName()), "pdfUp");
    chain->Draw(Form("%s/1e3>>%s", var_name.c_str(), h_down->GetName()), "pdfDown");
    h_up->SetLineColor(color);
    h_down->SetLineColor(color);
    // h_up->Divide(h_norm);
    // h_down->Divide(h_norm);
    delete chain;
    return;
}

TH1F* get_var_cmp(const string&  var_name, const string& hist_name)
{
    TH1F* h_nom_all = new TH1F(hist_name.c_str(), "h all", 3, 0.5, 3.5);

    TH1F* h_weight_temp = new TH1F("h_weight_temp", "weight temp", 100, 0, 100);

    TH1F* h_nnpdf_weight_nom = (TH1F*) h_weight_temp->Clone("h_nnpdf_weight_nom");
    TH1F* h_nnpdf_weight_up = (TH1F*) h_weight_temp->Clone("h_nnpdf_weight_up");
    TH1F* h_nnpdf_weight_down = (TH1F*) h_weight_temp->Clone("h_nnpdf_weight_down");

    TH1F* h_ct10_weight_nom = (TH1F*) h_weight_temp->Clone("h_ct10_weight_nom");
    TH1F* h_ct10_weight_up = (TH1F*) h_weight_temp->Clone("h_ct10_weight_up");
    TH1F* h_ct10_weight_down = (TH1F*) h_weight_temp->Clone("h_ct10_weight_down");

    TH1F* h_mmht_weight_nom = (TH1F*) h_weight_temp->Clone("h_mmht_weight_nom");
    TH1F* h_mmht_weight_up = (TH1F*) h_weight_temp->Clone("h_mmht_weight_up");
    TH1F* h_mmht_weight_down = (TH1F*) h_weight_temp->Clone("h_mmht_weight_down");

    get_pdf_sys("new_8_MMHT2014.root", h_mmht_weight_nom, h_mmht_weight_up, h_mmht_weight_down, kOrange+5, var_name.c_str());
    get_pdf_sys("new_9_CT10.root", h_ct10_weight_nom, h_ct10_weight_up, h_ct10_weight_down, kGreen+2, var_name.c_str());
    get_pdf_sys("new_6_NNPDF30_0130.root", h_nnpdf_weight_nom, h_nnpdf_weight_up, h_nnpdf_weight_down, kViolet-4, var_name.c_str());

    h_nom_all->SetBinContent(1, h_nnpdf_weight_nom->GetMean());
    h_nom_all->SetBinContent(2, h_ct10_weight_nom->GetMean());
    h_nom_all->SetBinContent(3, h_mmht_weight_nom->GetMean());
    /*
    h_nom_all->SetBinError(1, h_nnpdf_weight_nom->GetRMS());
    h_nom_all->SetBinError(2, h_ct10_weight_nom->GetRMS());
    h_nom_all->SetBinError(3, h_mmht_weight_nom->GetRMS());
    **/
    h_nom_all->SetBinError(1, h_nnpdf_weight_nom->GetMeanError());
    h_nom_all->SetBinError(2, h_ct10_weight_nom->GetMeanError());
    h_nom_all->SetBinError(3, h_mmht_weight_nom->GetMeanError());
   
    delete h_nnpdf_weight_nom;
    delete h_nnpdf_weight_up;
    delete h_nnpdf_weight_down;

    delete h_ct10_weight_nom;
    delete h_ct10_weight_up;
    delete h_ct10_weight_down;

    delete h_mmht_weight_nom;
    delete h_mmht_weight_up;
    delete h_mmht_weight_down;

    delete h_weight_temp;

    return h_nom_all;
}
