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

using namespace std;

TGraph *grInitial;
TGraph *grRefined;
TGraph *grNominal;

TF1 *tf, *tf_expo, *tf_surv, *tf_irr, *tf_ln_expo, *tf_customPol4, *tf_gaus_custom;

vector<double> dataSDE, element_sums;
vector<double> M30_OW, M30_PL, M30_Al, M30_Cu, M30_SDE;

vector<double> dataOW, dataPL, dataAl, dataCu, dataPLCU, dataOWAl, dataPLAl, dataOWCu_PLCu_ratio;
int Ndata;

TGraph *gr_OW_PLCu_pure, *gr_OW_PLCu_mixed, *gr_OW_PLCu_beta, *gr_OW_PLAl_pure, *gr_OW_PLAl_mixed, *gr_OW_PLAl_beta;
TGraph *gr_sum_pure, *gr_sum_mixed;
TGraph *gr_sum_sde;
TGraph *gr_OW_sde, *gr_PL_sde, *gr_Al_sde, *gr_Cu_sde;
TGraph *gr_OW_Al, *gr_OW_Al_pure, *gr_OW_Al_beta, *gr_OW_Al_mixed;
TGraph *gr_OW_PL_pure, *gr_OW_PL_beta, *gr_OW_PL_mixed;
TGraph *gr_Cu_PL_pure, *gr_Cu_PL_beta, *gr_Cu_PL_mixed;
TGraph *gr_OWCu_PLCu_ratio_pure, *gr_OWCu_PLCu_ratio_beta, *gr_OWCu_PLCu_ratio_mixed;
TGraph *gr_sum_sde_mixed, *gr_sum_sde_photons, *gr_sum_sde_betas;
TGraph *gr_PLCu_sde;
TGraph *gr_M30_mixed_combined_sde;

vector<double> initialPredictedSDE;

Double_t gausFnOffset(Double_t *x, Double_t *par){
    double gausFn = 0.0;
    // if( x[0] < par[1] ){
    //     double diff = fabs(x[0] - par[1]);
    //     x[0] = par[1] + diff;
    // }
    result = par[0]*exp((-0.5)*pow((x[0] - par[1])/par[2], 2) ) + par[3];
    return result;
}

Double_t logFN(Double_t *x, Double_t *par){
    double result = 0;
    result = par[0] + par[1]*log(x[0]);

    return result;
}

Double_t expoFN(Double_t *x, Double_t *par){
    double result = 0;
    // result = par[0]*(pow(par[1], x[0])) + par[2];
    result = par[0] + par[1]*exp(par[2]*x[0]);

    return result;
}

Double_t sdeFN(Double_t *x, Double_t *par){
    double sdeFit = 0;
    double OWval = dataOW[x[0]];
    double PLval = dataPL[x[0]];
    double Alval = dataAl[x[0]];
    double Cuval = dataCu[x[0]];

    sdeFit = par[0]*OWval + par[1]*PLval + par[2]*Alval + par[3]*Cuval;

    return sdeFit;

}

Double_t survFN(Double_t *x, Double_t *par){
    double surv = 0;

    surv = par[0] + par[1]*pow( par[2]*x[0] + par[3], 3 );

    return surv;
}

Double_t pol3_PoI(Double_t *x, Double_t *par){
    double pol3 = 0;

    pol3 = par[0] + par[1]*pow(x[0] - par[2], 2);

    return pol3;
}

Double_t pol3_LinFF(Double_t *x, Double_t *par){
    double pol3 = 0;

    pol3 = par[0]*(x[0] - par[1])*(x[0] - par[2])*(x[0] - par[3]);

    return pol3;
}

Double_t pol3_IrrFF(Double_t *x, Double_t *par){
    double pol3 = 0;

    pol3 = (x[0] - par[0])*(par[1]*x[0]*x[0] + par[2]*x[0] + par[3]);

    return pol3;
}

Double_t ln_expo_fn(Double_t *x, Double_t *par){
    double ln_expo = 0.0;

    ln_expo = par[0]*log(x[0]/par[1]) + par[2]*exp((-0.5)*((x[0] - par[3])/par[4]));

    return ln_expo;
}

Double_t ln_expo_fn_6par(Double_t *x, Double_t *par){
    double ln_expo = 0.0;

    ln_expo = par[0]*log((x[0]-par[1])/par[2]) + par[3]*exp((-0.5)*((x[0] - par[4])/par[5]));

    return ln_expo;
}

Double_t customPol4(Double_t *x, Double_t *par){
    double pol4 = 0.0;

    pol4 = par[0] + par[1]*(x[0]-par[3]) + par[2]*pow(x[0] - par[3], 2);

    return pol4;
}

void loadData(const char* filename, vector<double>& OW, vector<double>& PL, vector<double>& Al, vector<double>& Cu, vector<double>& SDE, int percentage = 100) {
    ifstream file(filename);
    string line;
    int rowCounter = 0;
    int skipFactor = (percentage == 100) ? 0 : 100 / (100 - percentage);
    while (getline(file, line)) {
        if (skipFactor != 0 && rowCounter % skipFactor == 0) {
            rowCounter++;
            continue;
        }
        stringstream ss(line);
        double ow, pl, al, cu, sde;
        ss >> ow >> pl >> al >> cu >> sde;
        // if( sde > 525 || sde < 500){
        //     continue;
        // }
        OW.push_back(ow);
        PL.push_back(pl);
        Al.push_back(al);
        Cu.push_back(cu);
        SDE.push_back(sde);

        dataOW.push_back(ow);
        dataPL.push_back(pl);
        dataAl.push_back(al);
        dataCu.push_back(cu);
        dataPLCU.push_back(pl/cu);
        dataPLAl.push_back(pl/al);
        dataOWAl.push_back(ow/al);
        dataOWCu_PLCu_ratio.push_back( (ow - cu)/(pl - cu));
        dataSDE.push_back(sde);
        rowCounter++;
    }
}



double compositeFunction(double* x, double* params) {
    double OW = x[0];
    double PL = x[1];
    double Al = x[2];
    double Cu = x[3];
    return params[0] * OW + params[1] * PL + params[2] * Cu;
}

void analyze(const vector<double>& oldCoefficients, double percent_data_to_scan = 100) {

    cout << endl;
    cout << "Scanning " << int(percent_data_to_scan) << "% of the data in .csv file" << endl << endl;

    vector<double> OW, PL, Al, Cu, SDE;
    // loadData("bl_dde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("blm_sde_data_4coeff.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bl_kr_dde_data_4coeff.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bl_kr_mixed_sde_data_4coeff.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bl_sde_data_4coeff.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("mixed_sde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_mixed_sde_data_4coeff.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("mixed_dde_data_select_top_mixed2.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_points_mixed_and_pure_dde.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("mixed_dde_data_select_mixed.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_mixed_dde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    //loadData("bh_sde_data_orig_copy.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bh_sde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_values_dde.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_types_betaSDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("all_types_beta_SDE_data_grouping.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bl_kr_sde_data_pure.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bh_sr_sde_data_pure.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("beta_sde_data_mixed.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bl_kr_betaSDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_photon_SDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("ph_pure_sde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_beta_indicator_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("PL_sde_data_purePH.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("BH_sde_data_purePH.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_beta_sde_data_mixed.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_mixed_beta_dde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_mixed_beta_sde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("linear_beta_test_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("linear_mixed_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("BmixPL_betaSDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("mixed_combined_SDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("kr_mixed_betaSDE_data_PL.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_mixed_combinedSDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_mixedOnly_DDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);

    // loadData("petg_combined_SDE_NON_M30_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_combined_SDE_M30_ONLY_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("petg_mixed_betaSDE_M30_ONLY_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    loadData("petg_Sr_combined_SDE_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);


    int n = OW.size();
    if (n == 0) {
        cerr << "No data loaded." << endl;
        return;
    }

    // Calculate initial predicted SDE values
    initialPredictedSDE.resize(n);
    for (int i = 0; i < n; ++i) {
        initialPredictedSDE[i] = oldCoefficients[0] * OW[i] + oldCoefficients[1] * PL[i] + oldCoefficients[2] * Al[i] + oldCoefficients[3] * Cu[i];
    }

    // Perform linear regression to refine coefficients
    TMatrixD X(n, 4);  // Matrix for measurements
    TVectorD y(n);     // Vector for SDE values

    for (int i = 0; i < n; ++i) {
        X[i][0] = OW[i];
        X[i][1] = PL[i];
        X[i][2] = Al[i];
        X[i][3] = Cu[i];
        y[i] = SDE[i];
    }

    TMatrixD Xt(TMatrixD::kTransposed, X);
    TMatrixD XtX = Xt * X;
    TMatrixD XtX_inv = XtX.Invert();
    TVectorD beta = XtX_inv * (Xt * y);

    // Print out the refined coefficients
    cout << "Refined coefficients: OW = " << beta[0] << ", PL = " << beta[1] << ", Al = " << beta[2]<< ", Cu = " << beta[3] << endl;

    cout << endl;

    for( int coeff = 0; coeff < 4; coeff++ ){
        cout << "c" << coeff +1 << " = " << Form("%.3f", beta[coeff]) << endl;
    }

    // Calculate refined predicted SDE values
    vector<double> refinedPredictedSDE(n);
    for (int i = 0; i < n; ++i) {
        refinedPredictedSDE[i] = (beta[0] * OW[i] + beta[1] * PL[i] + beta[2] * Al[i]+ beta[3] * Cu[i]);
    }

    // Calculate initial and refined mean squared errors
    double initialError = 0, refinedError = 0;
    for (int i = 0; i < n; ++i) {
        initialError += TMath::Power(initialPredictedSDE[i] - SDE[i], 2);
        refinedError += TMath::Power(refinedPredictedSDE[i] - SDE[i], 2);
    }
    initialError /= n;
    refinedError /= n;
    cout << "Initial mean squared error: " << initialError << endl;
    cout << "Refined mean squared error: " << refinedError << endl;

    // Plot the comparison
    TCanvas *c1 = new TCanvas("c1", "Comparison", 800, 600);
    grInitial = new TGraph(n);
    grRefined = new TGraph(n);
    grNominal = new TGraph(n);


    for (int i = 0; i < n; ++i) {
        grInitial->SetPoint(i, i, dataSDE[i]);
        grRefined->SetPoint(i, i, refinedPredictedSDE[i]);
        grNominal->SetPoint(i, i, dataSDE[i]);
    }

    double initial_min = TMath::MinElement(n, grInitial->GetY());
    double initial_max = TMath::MaxElement(n, grInitial->GetY());

    double refined_min = TMath::MinElement(n, grRefined->GetY());
    double refined_max = TMath::MaxElement(n, grRefined->GetY());

    double graph_min = min(initial_min, refined_min);
    double graph_max = max(initial_max, refined_max);

    // TLine *tl_SDE_nominal = new TLine(0.0, 500, n, 500);
    // tl_SDE_nominal->SetLineStyle(10);

    // grInitial->SetMarkerStyle(21);
    // grInitial->SetMarkerColor(kRed);
    // grInitial->SetTitle(Form("High Energy Beta: SDE vs Measurement Entry (%0.0f%%);Measurement Entry;High Energy Beta SDE", percent_data_to_scan));
    // grInitial->GetYaxis()->SetRangeUser(0.95*graph_min, 1.05*graph_max);
    // grInitial->Draw("AP");

    grRefined->SetMarkerStyle(22);
    grRefined->SetMarkerColor(kBlue);
    grRefined->SetTitle("Refined vs Nominal SDE;Nominal SDE;Refined Predicted SDE");
    grRefined->Draw("AP");

    grNominal->SetMarkerColor(kRed);
    grNominal->SetMarkerStyle(20);
    grNominal->SetMarkerSize(0.5);
    grNominal->Draw("P same");

    // tl_SDE_nominal->Draw("same");

    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->AddEntry(grNominal, "Nominal", "p");
    // leg->AddEntry(grInitial, "Initial", "p");
    leg->AddEntry(grRefined, "Refined", "p");
    leg->Draw();

    c1->Update();
}

vector<double> oldCoefficients(4, 0.0);

void BH_sde_calc() {
    
    double percent_data_to_scan = 100;

    // Example old coefficients (OW, PL, Cu)
    // oldCoefficients[0] = 1.4624;
    // oldCoefficients[1] = -1.4983;
    // oldCoefficients[2] = 0.3246;
    // oldCoefficients[3] = 0.6824;

    oldCoefficients[0] = 1.0;
    oldCoefficients[1] = -1.0;
    oldCoefficients[2] = 0.0;
    oldCoefficients[3] = 0.0;

    // Call analyze with old coefficients
    analyze(oldCoefficients, percent_data_to_scan);

    Ndata = dataOW.size();
    gr_PLCu_sde = new TGraph(Ndata);
    gr_PLCu_sde->SetMarkerColor(kBlue);
    gr_PLCu_sde->SetMarkerStyle(6);
    gr_PLCu_sde->GetYaxis()->SetTitle("SDE");
    gr_PLCu_sde->GetXaxis()->SetTitle("PL/Cu");
    gr_PLCu_sde->SetTitle("SDE vs PL/Cu");

    gr_OW_PLCu_pure = new TGraph(Ndata);
    gr_OW_PLCu_pure->SetMarkerColor(kBlue);
    gr_OW_PLCu_pure->SetMarkerStyle(6);
    gr_OW_PLCu_pure->GetYaxis()->SetTitle("OW");
    gr_OW_PLCu_pure->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PLCu_pure->SetTitle("PL vs Cu for Pure and Mixed Sources");

    gr_OW_PLCu_mixed = new TGraph(Ndata);
    gr_OW_PLCu_mixed->SetMarkerColor(6);
    gr_OW_PLCu_mixed->SetMarkerStyle(7);
    gr_OW_PLCu_mixed->GetYaxis()->SetTitle("OW");
    gr_OW_PLCu_mixed->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PLCu_mixed->SetTitle("PL vs Cu for Pure and Mixed Sources");

    gr_OWCu_PLCu_ratio_beta = new TGraph(Ndata);
    gr_OWCu_PLCu_ratio_beta->SetMarkerColor(kGreen);
    gr_OWCu_PLCu_ratio_beta->SetMarkerStyle(7);
    gr_OWCu_PLCu_ratio_beta->GetYaxis()->SetTitle("OW");
    gr_OWCu_PLCu_ratio_beta->GetXaxis()->SetTitle("PL/Cu");
    gr_OWCu_PLCu_ratio_beta->SetTitle("PL vs Cu for Pure and Mixed Sources"); 

    gr_OWCu_PLCu_ratio_pure = new TGraph(Ndata);
    gr_OWCu_PLCu_ratio_pure->SetMarkerColor(kBlue);
    gr_OWCu_PLCu_ratio_pure->SetMarkerStyle(6);
    gr_OWCu_PLCu_ratio_pure->GetYaxis()->SetTitle("OW");
    gr_OWCu_PLCu_ratio_pure->GetXaxis()->SetTitle("PL/Cu");
    gr_OWCu_PLCu_ratio_pure->SetTitle("PL vs Cu for Pure and Mixed Sources");

    gr_OWCu_PLCu_ratio_mixed = new TGraph(Ndata);
    gr_OWCu_PLCu_ratio_mixed->SetMarkerColor(6);
    gr_OWCu_PLCu_ratio_mixed->SetMarkerStyle(7);
    gr_OWCu_PLCu_ratio_mixed->GetYaxis()->SetTitle("OW");
    gr_OWCu_PLCu_ratio_mixed->GetXaxis()->SetTitle("PL/Cu");
    gr_OWCu_PLCu_ratio_mixed->SetTitle("OWCu_PLCu for Pure and Mixed Sources");

    gr_OW_PLCu_beta = new TGraph(Ndata);
    gr_OW_PLCu_beta->SetMarkerColor(kGreen);
    gr_OW_PLCu_beta->SetMarkerStyle(7);
    gr_OW_PLCu_beta->GetYaxis()->SetTitle("OW");
    gr_OW_PLCu_beta->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PLCu_beta->SetTitle("PL/Cu vs OW for Pure and Mixed Sources"); 

    gr_Cu_PL_pure = new TGraph(Ndata);
    gr_Cu_PL_pure->SetMarkerColor(kBlue);
    gr_Cu_PL_pure->SetMarkerStyle(6);
    gr_Cu_PL_pure->GetYaxis()->SetTitle("OW");
    gr_Cu_PL_pure->GetXaxis()->SetTitle("PL/Cu");
    gr_Cu_PL_pure->SetTitle("PL/Cu vs OW for Pure and Mixed Sources");

    gr_Cu_PL_mixed = new TGraph(Ndata);
    gr_Cu_PL_mixed->SetMarkerColor(6);
    gr_Cu_PL_mixed->SetMarkerStyle(7);
    gr_Cu_PL_mixed->GetYaxis()->SetTitle("OW");
    gr_Cu_PL_mixed->GetXaxis()->SetTitle("PL/Cu");
    gr_Cu_PL_mixed->SetTitle("OW vs PL/Cu for Pure and Mixed Sources");

    gr_Cu_PL_beta = new TGraph(Ndata);
    gr_Cu_PL_beta->SetMarkerColor(kGreen);
    gr_Cu_PL_beta->SetMarkerStyle(7);
    gr_Cu_PL_beta->GetYaxis()->SetTitle("OW");
    gr_Cu_PL_beta->GetXaxis()->SetTitle("PL/Cu");
    gr_Cu_PL_beta->SetTitle("OW vs PL/Cu for Pure and Mixed Sources"); 

    gr_OW_PL_pure = new TGraph(Ndata);
    gr_OW_PL_pure->SetMarkerColor(kBlue);
    gr_OW_PL_pure->SetMarkerStyle(6);
    gr_OW_PL_pure->GetYaxis()->SetTitle("OW");
    gr_OW_PL_pure->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PL_pure->SetTitle("PL vs OW for Pure and Mixed Sources");

    gr_OW_PL_mixed = new TGraph(Ndata);
    gr_OW_PL_mixed->SetMarkerColor(6);
    gr_OW_PL_mixed->SetMarkerStyle(7);
    gr_OW_PL_mixed->GetYaxis()->SetTitle("OW");
    gr_OW_PL_mixed->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PL_mixed->SetTitle("PL vs OW for Pure and Mixed Sources");

    gr_OW_PL_beta = new TGraph(Ndata);
    gr_OW_PL_beta->SetMarkerColor(kGreen);
    gr_OW_PL_beta->SetMarkerStyle(7);
    gr_OW_PL_beta->GetYaxis()->SetTitle("OW");
    gr_OW_PL_beta->GetXaxis()->SetTitle("PL/Cu");
    gr_OW_PL_beta->SetTitle("PL vs OW for Pure and Mixed Sources"); 

    gr_OW_PLAl_pure = new TGraph(Ndata);
    gr_OW_PLAl_pure->SetMarkerColor(kGreen);
    gr_OW_PLAl_pure->SetMarkerStyle(6);
    gr_OW_PLAl_pure->GetYaxis()->SetTitle("OW");
    gr_OW_PLAl_pure->GetXaxis()->SetTitle("PL/Al");
    gr_OW_PLAl_pure->SetTitle("PL/Al vs OW for Pure and Mixed Sources");

    gr_OW_PLAl_mixed = new TGraph(Ndata);
    gr_OW_PLAl_mixed->SetMarkerColor(9);
    gr_OW_PLAl_mixed->SetMarkerStyle(7);
    gr_OW_PLAl_mixed->GetYaxis()->SetTitle("OW");
    gr_OW_PLAl_mixed->GetXaxis()->SetTitle("PL/Al");
    gr_OW_PLAl_mixed->SetTitle("PL/Al vs OW for Pure and Mixed Sources");

    gr_OW_PLAl_beta = new TGraph(Ndata);
    gr_OW_PLAl_beta->SetMarkerColor(2);
    gr_OW_PLAl_beta->SetMarkerStyle(7);
    gr_OW_PLAl_beta->GetYaxis()->SetTitle("OW");
    gr_OW_PLAl_beta->GetXaxis()->SetTitle("PL/Al");
    gr_OW_PLAl_beta->SetTitle("OW vs PL/Al for Pure and Mixed Sources"); 

    gr_OW_Al_pure = new TGraph(Ndata);
    gr_OW_Al_pure->SetMarkerColor(kGreen);
    gr_OW_Al_pure->SetMarkerStyle(6);
    gr_OW_Al_pure->GetYaxis()->SetTitle("OW");
    gr_OW_Al_pure->GetXaxis()->SetTitle("Al");
    gr_OW_Al_pure->SetTitle("Al vs OW for Pure and Mixed Sources");

    gr_OW_Al_mixed = new TGraph(Ndata);
    gr_OW_Al_mixed->SetMarkerColor(9);
    gr_OW_Al_mixed->SetMarkerStyle(7);
    gr_OW_Al_mixed->GetYaxis()->SetTitle("OW");
    gr_OW_Al_mixed->GetXaxis()->SetTitle("Al");
    gr_OW_Al_mixed->SetTitle("Al vs OW for Pure and Mixed Sources");

    gr_OW_Al_beta = new TGraph(Ndata);
    gr_OW_Al_beta->SetMarkerColor(2);
    gr_OW_Al_beta->SetMarkerStyle(7);
    gr_OW_Al_beta->GetYaxis()->SetTitle("OW");
    gr_OW_Al_beta->GetXaxis()->SetTitle("Al");
    gr_OW_Al_beta->SetTitle("Al vs OW for Pure and Mixed Sources"); 

    gr_M30_mixed_combined_sde = new TGraph(Ndata);
    gr_M30_mixed_combined_sde->SetMarkerColor(2);
    gr_M30_mixed_combined_sde->SetMarkerStyle(7);
    gr_M30_mixed_combined_sde->GetYaxis()->SetTitle("Combined SDE");
    gr_M30_mixed_combined_sde->GetXaxis()->SetTitle("(OW + PL + Al + Cu)");
    gr_M30_mixed_combined_sde->SetTitle("Combined SDE vs Element Sums - Mixed M30 Only"); 

    for( int i = 0; i < Ndata; i++ ){
        gr_OW_PLCu_pure->SetPoint(i, dataOW[i], dataPLCU[i]);
        gr_OW_PLAl_pure->SetPoint(i, dataOW[i], dataPLAl[i]);
        gr_OW_PL_pure->SetPoint(i, dataOW[i], dataPL[i]);
        gr_OW_Al_pure->SetPoint(i, dataOW[i], dataAl[i]);
        gr_Cu_PL_pure->SetPoint(i, dataCu[i], dataPL[i]);
        double OWCu_PLCu_ratio = 0;
        double newPL = dataPL[i];
        if( dataPL[i] = dataCu[i] ){
            newPL = 1.000000001*dataPL[i];
        }
        OWCu_PLCu_ratio = (dataOW[i] - dataCu[i])/(newPL - dataCu[i]);
        gr_OWCu_PLCu_ratio_pure->SetPoint(i, (newPL - dataCu[i]), (dataOW[i] - dataCu[i]));
    }
    // for( int i = 244; i < Ndata; i++ ){
    //     gr_OW_PLCu_beta->SetPoint(i, dataOW[i], dataPLCU[i]);
    //     gr_OW_PLAl_beta->SetPoint(i, dataOW[i], dataPLAl[i]);
    //     gr_OW_PL_beta->SetPoint(i, dataOW[i], dataPL[i]);
    //     gr_OW_Al_beta->SetPoint(i, dataOW[i], dataAl[i]);
    //     gr_Cu_PL_beta->SetPoint(i, dataCu[i], dataPL[i]);

    //     double OWCu_PLCu_ratio = 0;
    //     double newPL = dataPL[i];
    //     if( dataPL[i] = dataCu[i] ){
    //         newPL = 1.000000001*dataPL[i];
    //     }
    //     OWCu_PLCu_ratio = (dataOW[i] - dataCu[i])/(newPL - dataCu[i]);
    //     gr_OWCu_PLCu_ratio_beta->SetPoint(i, (newPL - dataCu[i]), (dataOW[i] - dataCu[i]));
    // }
    // for( int i = 99; i < Ndata; i++ ){
    //     gr_OW_PLCu_mixed->SetPoint(i, dataOW[i], dataPLCU[i]);
    //     gr_OW_PLAl_mixed->SetPoint(i, dataOW[i], dataPLAl[i]);
    //     gr_OW_PL_mixed->SetPoint(i, dataOW[i], dataPL[i]);
    //     gr_OW_Al_mixed->SetPoint(i, dataOW[i], dataAl[i]);
    //     gr_Cu_PL_mixed->SetPoint(i, dataCu[i], dataPL[i]);

    //     double OWCu_PLCu_ratio = 0;
    //     double newPL = dataPL[i];
    //     if( dataPL[i] = dataCu[i] ){
    //         newPL = 1.000000001*dataPL[i];
    //     }
    //     OWCu_PLCu_ratio = (dataOW[i] - dataCu[i])/(newPL - dataCu[i]);
    //     gr_OWCu_PLCu_ratio_mixed->SetPoint(i, (newPL - dataCu[i]), (dataOW[i] - dataCu[i]));
    // }

    TCanvas *c = new TCanvas("c", "c", 600, 500);
    gr_OW_PLCu_pure->Draw("AP");
    gr_OW_PLCu_mixed->Draw("P+same");
    gr_OW_PLCu_beta->Draw("P+same");

    TLegend *tl_mixed = new TLegend(0.7, 0.3, 0.85, 0.45);
    tl_mixed->AddEntry(gr_OW_PLCu_pure, "Pure Photons");
    tl_mixed->AddEntry(gr_OW_PLCu_beta, "Pure Beta");
    tl_mixed->AddEntry(gr_OW_PLCu_mixed, "Mixed");
    tl_mixed->Draw("same");

    c->Modified();
    c->Update();

    tf = new TF1("tf", "pol3", 0, 7000);

    // gr_sum_pure = new TGraph(10);
    // gr_sum_pure->SetMarkerColor(kRed);

    // gr_sum_mixed = new TGraph(Ndata - 10);
    // gr_sum_mixed->SetMarkerColor(kMagenta);

    // for( int i = 0; i < 10; i++ ){
    //     double sum = dataOW[i] + dataPL[i] + dataAl[i] + dataCu[i];
    //     gr_sum_pure->SetPoint(i, i, sum);

    // }
    // for( int i = 0; i < Ndata - 10; i++ ){
    //     double sum = dataOW[i] + dataPL[i] + dataAl[i] + dataCu[i];
    //     gr_sum_mixed->SetPoint(i, i, sum);

    // }
    double element_sum = 0;  
    for( int i = 0; i <  Ndata; i++ ){
        element_sum = 0;
        element_sum = dataOW[i] + dataPL[i] + dataAl[i] + dataCu[i];
        element_sums.push_back(element_sum);
    }
    gr_sum_sde = new TGraph(Ndata);
    gr_sum_sde->SetMarkerColor(9);
    gr_sum_sde->SetMarkerStyle(7);
    gr_sum_sde->GetYaxis()->SetTitle("SDE");
    gr_sum_sde->GetXaxis()->SetTitle("(OW+PL+Al+Cu)^{-1}");
    gr_sum_sde->SetTitle("SDE vs (OW+PL+Al+Cu) for Low E Mixed Beta");  

    gr_sum_sde_photons = new TGraph(Ndata);
    gr_sum_sde_photons->SetMarkerColor(kGreen);
    gr_sum_sde_photons->SetMarkerStyle(7);
    gr_sum_sde_photons->GetYaxis()->SetTitle("SDE");
    gr_sum_sde_photons->GetXaxis()->SetTitle("(OW+PL+Al+Cu)^{-1}");
    gr_sum_sde_photons->SetTitle("SDE vs (OW+PL+Al+Cu)");  

    gr_sum_sde_betas = new TGraph(Ndata);
    gr_sum_sde_betas->SetMarkerColor(6);
    gr_sum_sde_betas->SetMarkerStyle(7);
    gr_sum_sde_betas->GetYaxis()->SetTitle("SDE");
    gr_sum_sde_betas->GetXaxis()->SetTitle("(OW+PL+Al+Cu)^{-1}");
    gr_sum_sde_betas->SetTitle("SDE vs (OW+PL+Al+Cu)");  

    gr_sum_sde_mixed = new TGraph(Ndata);
    gr_sum_sde_mixed->SetMarkerColor(9);
    gr_sum_sde_mixed->SetMarkerStyle(7);
    gr_sum_sde_mixed->GetYaxis()->SetTitle("SDE");
    gr_sum_sde_mixed->GetXaxis()->SetTitle("(OW+PL+Al+Cu)^{-1}");
    gr_sum_sde_mixed->SetTitle("SDE vs (OW+PL+Al+Cu)");  

    int photons_lower = 0;
    int photons_upper = 60;

    int betas_lower = 60;
    int betas_upper = 98;

    int mixed_lower = 98;
    int mixed_upper = Ndata;

    for( int i = 0; i <  Ndata; i++ ){
        gr_sum_sde_photons->SetPoint(i, element_sums[i], dataSDE[i]);
    }
    for( int i = betas_lower; i <  betas_upper; i++ ){
        gr_sum_sde_betas->SetPoint(i, element_sums[i], dataSDE[i]);
    }
    for( int i = mixed_lower; i <  mixed_upper; i++ ){
        gr_sum_sde_mixed->SetPoint(i, element_sums[i], dataSDE[i]);
    }

    for( int i = 0; i <  Ndata; i++ ){
        gr_sum_sde->SetPoint(i, element_sum, dataSDE[i]);
        gr_PLCu_sde->SetPoint(i, sqrt(pow(dataOWCu_PLCu_ratio[i],2)), dataSDE[i]);
    }
    TCanvas *c_PLCu_sde = new TCanvas("c_PLCu_sde", "c_PLCu_sde", 600, 500);
    gr_sum_sde_photons->Draw("AP");
    // gr_sum_sde_betas->Draw("Psame");
    // gr_sum_sde_mixed->Draw("Psame");

    TLegend *tleg_sde = new TLegend(0.15, 0.6, 0.30, 0.75);
    tleg_sde->AddEntry(gr_sum_sde_photons, "Pure Photons");
    tleg_sde->AddEntry(gr_sum_sde_betas, "Pure Betas");
    tleg_sde->AddEntry(gr_sum_sde_mixed, "Mixed");
    tleg_sde->Draw("same");

    tf_exp = new TF1("tf_expo", expoFN, 0.00005,0.004,3);
    gr_sum_sde->Fit("tf_expo", "RMSE");
    gr_sum_sde->GetYaxis()->SetRangeUser(0, 3500);
    c_PLCu_sde->Modified();
    c_PLCu_sde->Update();

    gr_PL_sde = new TGraph(Ndata);
    gr_OW_sde = new TGraph(Ndata);
    gr_Al_sde = new TGraph(Ndata);
    gr_Cu_sde = new TGraph(Ndata);

    for( int i = 0; i < Ndata; i++ ){
        gr_OW_sde->SetPoint(i, dataOW[i], dataSDE[i]);
        gr_PL_sde->SetPoint(i, dataPL[i], dataSDE[i]);
        gr_Al_sde->SetPoint(i, dataAl[i], dataSDE[i]);
        gr_Cu_sde->SetPoint(i, dataCu[i], dataSDE[i]);
    }
    gr_OW_sde->SetLineColor(kBlue);
    gr_OW_sde->SetMarkerColor(kBlue);

    gr_PL_sde->SetLineColor(kGreen);
    gr_PL_sde->SetMarkerColor(kGreen);

    gr_Al_sde->SetLineColor(kCyan);
    gr_Al_sde->SetMarkerColor(kCyan);

    gr_Cu_sde->SetLineColor(kMagenta);
    gr_Cu_sde->SetMarkerColor(kMagenta);

    tf_surv = new TF1("tf_surv", survFN, 0, 13000, 4);
    // tf_surv->SetParLimits(0, 550, 650); //611
    // tf_surv->SetParLimits(1, 1.0e-7, 2.0e-7); //1.4681e-07, 1.23998e-07
    // tf_surv->SetParLimits(2, 0.3, 0.4); //0.3483, 0.358
    // tf_surv->SetParLimits(3, -1350, -1250); //-1295, -1304

    tf_surv->SetParLimits(0, 610, 612); //611
    tf_surv->SetParLimits(1, 1.2e-07, 1.3e-07); //1.4681e-07, 1.23998e-07
    tf_surv->SetParLimits(2, 0.355, 0.363); //0.3483, 0.358
    tf_surv->SetParLimits(3, -1310, -1300); //-1295, -1304

    // gr_sum_sde->Fit("tf_surv", "RMSE+");

    tf_irr = new TF1("tf_irr", pol3_IrrFF, 0, 8000, 4);
    tf_gaus_custom = new TF1("tf_gaus_custom", gausFnOffset, -2, 3, 4);
    tf_gaus_custom->SetParLimits(0, 2000, 2400);
    tf_gaus_custom->SetParLimits(1, 0, 2);
    tf_gaus_custom->SetParLimits(2, 0, 0.01);
    tf_gaus_custom->SetParLimits(3, 450, 550);
    // gr_PLCu_sde->Fit("tf_gaus_custom", "RMSE+");
    // gr_PLCu_sde->Draw("AP");

    tf_ln_expo = new TF1("tf_ln_expo", ln_expo_fn, -100, 8000, 5);
    tf_ln_expo->SetParLimits(0, 23, 47);
    // tf_ln_expo->SetParameter(0, 0.000096); //29.2469

    tf_ln_expo->SetParLimits(1, 0.0000000001, 0.0001);
    // tf_ln_expo->SetParameter(1, 34); //5e-5
    tf_ln_expo->SetParLimits(2, 1.0, 1.5);
    // tf_ln_expo->SetParameter(2, 1.2); //2.05

    tf_ln_expo->SetParLimits(3, -2500, -1000);
    // tf_ln_expo->SetParameter(3, 1540); //-750

    tf_ln_expo->SetParLimits(4, -700, -500);
    // tf_ln_expo->SetParameter(4, -254);//-557.5

    // gr_sum_sde->Fit("tf_irr", "RMSE+");

    tf_ln_expo_6par = new TF1("tf_ln_expo_6par", ln_expo_fn_6par, -100, 8000, 6);
    tf_ln_expo_6par->SetParLimits(0, 42.5, 43.5);
    tf_ln_expo_6par->SetParLimits(1, 0, 0.05);
    tf_ln_expo_6par->SetParLimits(2, 0.004, 0.006);
    tf_ln_expo_6par->SetParLimits(3, 0.03, 0.07);
    tf_ln_expo_6par->SetParLimits(4, -3250, -3150);
    tf_ln_expo_6par->SetParLimits(5, -500, -475);

    tf_customPol4 = new TF1("customPol4", customPol4, 0, 600, 4);

}
