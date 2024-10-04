#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TH1.h>
#include <TF1.h>
#include <TApplication.h>
#include <TChain.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TFitResult.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "../include/vector_utility_functions.h"

using namespace std;

double CalculateStdDev(const std::vector<double>& data) {
    // Check if the vector is empty to avoid division by zero
    if (data.empty()) return 0.0;

    // Step 1: Calculate the mean
    double sum = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    double mean = sum / data.size();

    // Step 2: Calculate the variance
    double variance = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        variance += (data[i] - mean) * (data[i] - mean);
    }
    variance /= data.size();  // For population variance; use (data.size() - 1) for sample variance

    // Step 3: Return the standard deviation (square root of variance)
    return std::sqrt(variance);
}

double CalculateMean(const std::vector<double>& data) {
    // Check if the vector is empty to avoid division by zero
    if (data.empty()) return 0.0;

    // Step 1: Calculate the mean
    double sum = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    double mean = sum / data.size();

    return mean;
}

TGraph *grInitial;
TGraph *DDEInitial;
TGraph *grRefined;
TGraph *grDDERefined;
TGraph *grNominal;

TF1 *tf, *tf_expo, *tf_surv, *tf_irr, *tf_ln_expo, *tf_customPol4, *tf_gaus_custom;

vector<double> known_DDE, exper_DDE, known_SDE, exper_SDE, SDE_bias;

TH1D *h_SDE_bias;

vector<double> dataOW, dataPL, dataAl, dataCu, dataPLCU, dataOWAl, dataPLAl, dataOWCu_PLCu_ratio;
int Ndata;

void loadData(const char* filename, vector<double>& req_DDE, vector<double>& calc_DDE, vector<double>& req_SDE, vector<double>& calc_SDE) {
    ifstream file(filename);
    string line;

    bool header = true;

    int rowCounter = 0;
    while (getline(file, line)) {
        if( header ){
            if( rowCounter == 0 ){
                rowCounter++;
                continue;
            }
        }

        stringstream ss(line);
        double k_dde, e_dde, k_sde, e_sde, sde_bias_calc;
        sde_bias_calc = 0.0;

        ss >> k_dde >> e_dde >> k_sde >> e_sde;
        // if( sde > 525 || sde < 500){
        //     continue;
        // }
        req_SDE.push_back(k_dde);
        known_SDE.push_back(k_sde);
        exper_SDE.push_back(e_sde);

        known_DDE.push_back(k_dde);
        exper_DDE.push_back(e_dde);

        sde_bias_calc = (e_sde - k_sde)/k_sde;

        SDE_bias.push_back(sde_bias_calc);
        h_SDE_bias->Fill(sde_bias_calc);

        rowCounter++;
    }
}

vector<double> req_SDE, req_DDE, calc_SDE, calc_DDE;
int n = 0;
TString datafile;
void analyze() {

    const char *filename = datafile.Data();

    loadData(filename, req_DDE, calc_DDE, req_SDE, calc_SDE);


    n = req_SDE.size();
    if (n == 0) {
        cerr << "No data loaded." << endl;
        return;
    }


}

vector<double> oldCoefficients(4, 0.0);

double sde_bias_max, sde_bias_min, sde_bias_std_dev, sde_bias_mean, sde_Ndata, sde_min, sde_max, sde_mean, sde_final_mean;

TString bias_plot_title = "";
double bias_plot_xmin = 0.0, bias_plot_xmax = 0.0;
int bias_nbins;

double lower_95, upper_95, sig_level;
double t_value;

void plot_bias() {

    //Options:
    //photons
    //mixed
    //mixed_BL
    //mixed_BH

    // TString data_type = "photons";
    // TString data_type = "mixed";
    // TString data_type = "mixed_BL";
    TString data_type = "mixed_BH";

    TString dir = ".\\bias_data\\";

    if( data_type == "photons" ){
        bias_plot_title = "Hp(0.07) {SDE} Bias for Pure Photon Sources";
        bias_plot_xmin = -0.8;
        bias_plot_xmax = 1.0;
        bias_nbins = 50;
        // datafile = "pure_photons_DDE_SDE_bias_data.csv";
        datafile = "pure_photons_DDE_SDE_bias_data_02_10_2024.csv";
        t_value = 2.25;
    }
    else if( data_type == "mixed" ){
        bias_plot_title = "Hp(0.07) {SDE} Bias for Mixed Sources";
        bias_plot_xmin = -0.8;
        bias_plot_xmax = 0.4;
        bias_nbins = 40;
        // datafile = "all_mixed_DDE_SDE_bias_data.csv";
        datafile = "all_mixed_DDE_SDE_bias_data_02_10_2024.csv";
        t_value = 2.242;
    }
    else if( data_type == "mixed_BL" ){
        bias_plot_title = "Hp(0.07) {SDE} Bias for Mixed BL Sources";
        bias_plot_xmin = -0.8;
        bias_plot_xmax = 0.4;
        bias_nbins = 40;
        // datafile = "mixed_BL_DDE_SDE_bias_data.csv";
        datafile = "mixed_BL_DDE_SDE_bias_data_02_10_2024.csv";
        t_value = 2.2424;
    }
    else if( data_type == "mixed_BH" ){
        bias_plot_title = "Hp(0.07) {SDE} Bias for Mixed BH Sources";
        bias_plot_xmin = -0.8;
        bias_plot_xmax = 0.4;
        bias_nbins = 40;
        // datafile = "mixed_BH_DDE_SDE_bias_data.csv";
        datafile = "mixed_BH_DDE_SDE_bias_data_02_10_2024.csv";
        t_value = 2.2424;
    }
    else{
        bias_plot_title = "Hp(0.07) {SDE} Bias for Pure Beta Sources";
        bias_plot_xmin = -0.8;
        bias_plot_xmax = 0.4;
        bias_nbins = 100;
    }
    
    datafile = Form("%s%s", dir.Data(), datafile.Data());

    h_SDE_bias = new TH1D("h_SDE_bias", "SDE Bias", bias_nbins, bias_plot_xmin, bias_plot_xmax);

    analyze();

    TCanvas *c = new TCanvas("c", "c", 600, 500);

    sde_bias_max = *std::max_element(SDE_bias.begin(), SDE_bias.end());
    sde_bias_min = *std::min_element(SDE_bias.begin(), SDE_bias.end());
    sde_bias_std_dev = CalculateStdDev(SDE_bias);
    sde_bias_mean = CalculateMean(SDE_bias);
    sde_Ndata = SDE_bias.size();

    sde_min = *std::min_element(req_SDE.begin(), req_SDE.end());
    sde_max = *std::max_element(req_SDE.begin(), req_SDE.end());
    sde_mean = CalculateMean(req_SDE);

    h_SDE_bias->SetFillStyle(0);
    h_SDE_bias->SetFillColor(kRed);
    h_SDE_bias->SetTitle(bias_plot_title.Data());
    h_SDE_bias->GetXaxis()->SetTitle("Bias");
    h_SDE_bias->GetYaxis()->SetTitle("Entries (N)");
    h_SDE_bias->GetYaxis()->SetTitleOffset(1.3);
    h_SDE_bias->SetStats(false);

    h_SDE_bias->Draw("B");

    TH1D *h_SDE_bias_copy = (TH1D*)h_SDE_bias->Clone("h_SDE_bias_copy");
    h_SDE_bias_copy->SetFillStyle(1);
    h_SDE_bias_copy->SetLineColor(kBlack);
    h_SDE_bias_copy->Draw("Bsame");

    TPaveText *tpt_sde = new TPaveText(0.60, 0.60, 0.85, 0.85, "NDC");
    tpt_sde->AddText(Form("Min Bias: %0.3f", sde_bias_min));
    tpt_sde->AddText(Form("Max Bias: %0.3f", sde_bias_max));
    tpt_sde->AddText(Form("Mean Bias: %0.3f", sde_bias_mean));
    tpt_sde->AddText(Form("Bias StdDev: %0.4f", sde_bias_std_dev));
    tpt_sde->AddText(Form("N Data: %0.0f", sde_Ndata));

    lower_95 = sde_bias_mean - t_value*sde_bias_std_dev;
    upper_95 = sde_bias_mean + t_value*sde_bias_std_dev;

    // lower_95 = sde_bias_mean - 2*sde_bias_std_dev;
    // upper_95 = sde_bias_mean + 2*sde_bias_std_dev;
    
    TPaveText *tpt_ci = new TPaveText(0.60, 0.40, 0.85, 0.58, "NDC");
    tpt_ci->AddText("Dashed line represents");
    tpt_ci->AddText("95\% confidence interval");
    tpt_ci->AddText(Form("Lower limit: %0.03f", lower_95));
    tpt_ci->AddText(Form("Upper limit: %0.03f", upper_95));
    tpt_ci->Draw("same");

    if( sde_min == sde_max ){
        sde_final_mean = sde_min;
        tpt_sde->AddText(Form("All doses: %0.0f mrem", sde_min));
    }
    else{
        sde_final_mean = sde_mean;
        tpt_sde->AddText(Form("Dose range: %0.0f <= SDE <= %0.0f mrem", sde_min, sde_max));
    }


    tpt_sde->Draw("same");



    TLine *tl_lower_95 = new TLine(lower_95, 0, lower_95, h_SDE_bias->GetMaximum());
    tl_lower_95->SetLineStyle(10);
    tl_lower_95->SetLineColor(1);
    tl_lower_95->Draw("same");

    TLine *tl_upper_95 = new TLine(upper_95, 0, upper_95, h_SDE_bias->GetMaximum());
    tl_upper_95->SetLineStyle(10);
    tl_upper_95->SetLineColor(1);
    tl_upper_95->Draw("same");

    cout << endl;
    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl << endl;
    cout << "Finished analyzing datafile: " << datafile.Data() << endl << endl;
    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl << endl;
}
