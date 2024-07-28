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

vector<double> initialPredictedSDE;

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
        OW.push_back(ow);
        PL.push_back(pl);
        Al.push_back(al);
        Cu.push_back(cu);
        SDE.push_back(sde);
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
    loadData("bh_sde_data.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);
    // loadData("bh_sde_data_mod.csv", OW, PL, Al, Cu, SDE, percent_data_to_scan);

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

    // Calculate refined predicted SDE values
    vector<double> refinedPredictedSDE(n);
    for (int i = 0; i < n; ++i) {
        refinedPredictedSDE[i] = beta[0] * OW[i] + beta[1] * PL[i] + beta[2] * Al[i]+ beta[3] * Cu[i];
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

    for (int i = 0; i < n; ++i) {
        grInitial->SetPoint(i, i, initialPredictedSDE[i]);
        grRefined->SetPoint(i, i, refinedPredictedSDE[i]);
    }

    double initial_min = TMath::MinElement(n, grInitial->GetY());
    double initial_max = TMath::MaxElement(n, grInitial->GetY());

    double refined_min = TMath::MinElement(n, grRefined->GetY());
    double refined_max = TMath::MaxElement(n, grRefined->GetY());

    double graph_min = min(initial_min, refined_min);
    double graph_max = max(initial_max, refined_max);

    TLine *tl_SDE_nominal = new TLine(0.0, 500, n, 500);
    tl_SDE_nominal->SetLineStyle(10);

    grInitial->SetMarkerStyle(21);
    grInitial->SetMarkerColor(kRed);
    grInitial->SetTitle(Form("High Energy Beta: SDE vs Measurement Entry (%0.0f%%);Measurement Entry;High Energy Beta SDE", percent_data_to_scan));
    grInitial->GetYaxis()->SetRangeUser(0.95*graph_min, 1.05*graph_max);
    grInitial->Draw("AP");

    grRefined->SetMarkerStyle(22);
    grRefined->SetMarkerColor(kBlue);
    grRefined->SetTitle("Refined vs Nominal SDE;Nominal SDE;Refined Predicted SDE");
    grRefined->Draw("P same");

    tl_SDE_nominal->Draw("same");

    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->AddEntry(grInitial, "Initial", "p");
    leg->AddEntry(grRefined, "Refined", "p");
    leg->Draw();

    c1->Update();
}

vector<double> oldCoefficients;

void BH_sde_calc() {
    
    double percent_data_to_scan = 80;

    // Example old coefficients (OW, PL, Cu)
    oldCoefficients = {1.4624, -1.4983, 0.3246, 0.6824};

    // Call analyze with old coefficients
    analyze(oldCoefficients, percent_data_to_scan);

}
