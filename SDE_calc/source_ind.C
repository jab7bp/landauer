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

using namespace std;

double arr_boundary_x[] = {42.91, 42.91, 42.91, 42.91, 42.91, 42.91, 42.91, 42.91, 42.91, 90.20, 212.59, 403.93, 444.05, 453.02, 453.02, 453.02, 453.02, 453.02, 453.02};
double arr_boundary_y[] = {5026.68, 5026.68, 5026.68, 5026.68, 5026.68, 5026.68, 5026.68, 5026.68, 5026.68, 2610.08, 1851.75, 1401.94, 703.81, 397.72, 397.72, 397.72, 397.72, 397.72, 397.72};
vector<double> boundary_x, boundary_y;
TGraph *gr_boundary;
int boundary_n;
TF1 *tf_boundary = new TF1("tf_boundary", "pol5", 0, 500);

vector<double> dataSDE;
vector<double> pure_ph_M30_OW, pure_ph_M30_PL, pure_ph_M30_Cu, pure_ph_M30_Al;
vector<double> pure_ph_L50_OW, pure_ph_L50_PL, pure_ph_L50_Cu, pure_ph_L50_Al;
vector<double> pure_ph_Cs_OW, pure_ph_Cs_PL, pure_ph_Cs_Cu, pure_ph_Cs_Al;

vector<double> pure_b_kr_OW, pure_b_kr_PL, pure_b_kr_Al, pure_b_kr_Cu;
vector<double> pure_b_sr_OW, pure_b_sr_PL, pure_b_sr_Al, pure_b_sr_Cu;

vector<double> mixed_OW, mixed_PL, mixed_Al, mixed_Cu;

vector<double> v_max_PL, v_max_Cu, v_min_PL, v_min_Cu;

double M30_PL_max, M30_PL_min, M30_Cu_max, M30_cu_min;
double L50_PL_max, L50_PL_min, L50_Cu_max, L50_Cu_min;
double Cs_PL_max, Cs_PL_min, Cs_Cu_max, Cs_Cu_min;

double kr_PL_max, kr_PL_min, kr_Cu_max, kr_Cu_min;
double sr_PL_max, sr_PL_min, sr_Cu_max, sr_Cu_min;

double mixed_PL_max, mixed_PL_min, mixed_Cu_max, mixed_Cu_min;
double PL_max, PL_min, Cu_max, Cu_min;

TGraph *gr_pure_ph_M30_PL_Cu, *gr_pure_ph_L50_PL_Cu, *gr_pure_ph_Cs_PL_Cu;
TGraph *gr_pure_b_kr_PL_Cu, *gr_pure_b_sr_PL_Cu;
TGraph *gr_mixed_PL_Cu;

TGraph *copy_gr_pure_ph_M30_PL_Cu, *copy_gr_pure_ph_L50_PL_Cu, *copy_gr_pure_ph_Cs_PL_Cu;
TGraph *copy_gr_pure_b_kr_PL_Cu, *copy_gr_pure_b_sr_PL_Cu;
TGraph *copy_gr_mixed_PL_Cu;


int Ndata_M30, Ndata_L50, Ndata_Cs;
int Ndata_kr, Ndata_sr;
int Ndata_mixed;

TF1 *tf_poly, *tf_beta_upper, *tf_beta_lower;

// Double_t source_ind_poly(Double_t *x, Double_t *par){
//     double source_ind_poly = 0.0;

//     source_ind_poly = par[0] + par[1]*(x[0]-par[3]) + par[2]*pow((x[0]-par[3]), 2);

//     return source_ind_poly;
// }

Double_t source_ind_poly(Double_t *x, Double_t *par){
    double source_ind_poly = 0.0;

    source_ind_poly = par[0] + par[1]*x[0] + par[2]*pow(x[0], 2) + par[3]*pow(x[0], 3) + par[4]*pow(x[0], 4) + par[5]*pow(x[0], 5);

    return source_ind_poly;
}

Double_t pure_beta_upper(Double_t *x, Double_t *par){
    double pure_beta_upper = 0.0;

    pure_beta_upper = x[0]*par[0] + par[1];

    return pure_beta_upper;
}

Double_t pure_beta_lower(Double_t *x, Double_t *par){
    double pure_beta_lower = 0.0;

    pure_beta_lower = x[0]*par[0] + par[1];

    return pure_beta_lower;
}

void loadData(const char* filename, vector<double>& OW, vector<double>& PL, vector<double>& Al, vector<double>& Cu, int percentage = 100) {
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

        rowCounter++;
    }
}

void source_ind() {

    cout << endl;
    vector<double> OW, PL, Al, Cu, SDE;

    //Load Data from various source types:

    //Pure Photons:
    //M30
    // loadData("petg_M30_source_ind.csv", pure_ph_M30_OW, pure_ph_M30_PL, pure_ph_M30_Al, pure_ph_M30_Cu);
    loadData("source_ind_PL.csv", pure_ph_M30_OW, pure_ph_M30_PL, pure_ph_M30_Al, pure_ph_M30_Cu);

    //L50
	// loadData("petg_L50_source_ind.csv", pure_ph_L50_OW, pure_ph_L50_PL, pure_ph_L50_Al, pure_ph_L50_Cu);
    loadData("source_ind_PM.csv", pure_ph_L50_OW, pure_ph_L50_PL, pure_ph_L50_Al, pure_ph_L50_Cu);

    //Cs
	// loadData("petg_Cs_source_ind.csv", pure_ph_Cs_OW, pure_ph_Cs_PL, pure_ph_Cs_Al, pure_ph_Cs_Cu);
    loadData("source_ind_PH.csv", pure_ph_Cs_OW, pure_ph_Cs_PL, pure_ph_Cs_Al, pure_ph_Cs_Cu);

    //Pure Betas
    //Kr85
	// loadData("petg_all_Kr85_source_ind.csv", pure_b_kr_OW, pure_b_kr_PL, pure_b_kr_Al, pure_b_kr_Cu);
    loadData("source_ind_BL.csv", pure_b_kr_OW, pure_b_kr_PL, pure_b_kr_Al, pure_b_kr_Cu);

    //Sr90
	// loadData("petg_all_Sr90_source_ind.csv", pure_b_sr_OW, pure_b_sr_PL, pure_b_sr_Al, pure_b_sr_Cu);
    loadData("source_ind_BH_with_new.csv", pure_b_sr_OW, pure_b_sr_PL, pure_b_sr_Al, pure_b_sr_Cu);

    //Mixed Sources
	// loadData("petg_all_mixed_data_source_ind.csv", mixed_OW, mixed_PL, mixed_Al, mixed_Cu);
    loadData("source_ind_mixed.csv", mixed_OW, mixed_PL, mixed_Al, mixed_Cu);


    Ndata_M30 = pure_ph_M30_PL.size();
    Ndata_L50 = pure_ph_L50_PL.size();
    Ndata_Cs = pure_ph_Cs_PL.size();

    Ndata_kr = pure_b_kr_PL.size();
    Ndata_sr = pure_b_sr_PL.size();

    Ndata_mixed = mixed_PL.size();

    if (Ndata_M30 == 0){
    	cerr << "No M30 data loaded." << endl;
    	return;
    }
    if (Ndata_L50 == 0){
    	cerr << "No L50 data loaded." << endl;
    	return;
    }
    if (Ndata_Cs == 0){
    	cerr << "No Cs data loaded." << endl;
    	return;
    }
    if (Ndata_kr == 0){
    	cerr << "No Kr85 data loaded." << endl;
    	return;
    }
    if (Ndata_sr == 0){
    	cerr << "No Sr90 data loaded." << endl;
    	return;
    }
    if (Ndata_mixed == 0){
    	cerr << "No Mixed Source data loaded." << endl;
    	return;
    }

///Get maximum/minimum values in the element reading vectors

    // M30_PL_max = *std::max_element(pure_ph_M30_PL.begin(), pure_ph_M30_PL.end());
    // M30_PL_min = *std::min_element(pure_ph_M30_PL.begin(), pure_ph_M30_PL.end());
    // v_max_PL.push_back(M30_PL_max);
    // v_min_PL.push_back(M30_PL_min);

    // L50_PL_max = *std::max_element(pure_ph_L50_PL.begin(), pure_ph_L50_PL.end());
    // L50_PL_min = *std::min_element(pure_ph_L50_PL.begin(), pure_ph_L50_PL.end());
    // v_max_PL.push_back(L50_PL_max);
    // v_min_PL.push_back(L50_PL_min);

    // Cs_PL_max = *std::max_element(pure_ph_Cs_PL.begin(), pure_ph_Cs_PL.end());
    // Cs_PL_min = *std::min_element(pure_ph_Cs_PL.begin(), pure_ph_Cs_PL.end());
    // v_max_PL.push_back(Cs_PL_max);
    // v_min_PL.push_back(Cs_PL_min);

    // kr_PL_max = *std::max_element(pure_b_kr_PL.begin(), pure_b_kr_PL.end());
    // kr_PL_min = *std::min_element(pure_b_kr_PL.begin(), pure_b_kr_PL.end());
    // v_max_PL.push_back(kr_PL_max);
    // v_min_PL.push_back(kr_PL_min);

    // sr_PL_max = *std::max_element(pure_b_sr_PL.begin(), pure_b_sr_PL.end());
    // sr_PL_min = *std::min_element(pure_b_sr_PL.begin(), pure_b_sr_PL.end()); 
    // v_max_PL.push_back(sr_PL_max);
    // v_min_PL.push_back(sr_PL_min);

    // mixed_PL_max = *std::max_element(mixed_PL.begin(), mixed_PL.end());
	// mixed_PL_min = *std::min_element(mixed_PL.begin(), mixed_PL.end());
    // v_max_PL.push_back(mixed_PL_max);
    // v_min_PL.push_back(mixed_PL_min);

    // PL_max = *std::max_element(v_max_PL.begin(), v_max_PL.end());
    // PL_min = *std::min_element(v_min_PL.begin(), v_min_PL.end());

    // Cu_max = *std::max_element(v_max_Cu.begin(), v_max_Cu.end());
    // Cu_min = *std::min_element(v_min_Cu.begin(), v_min_Cu.end());

	double marker_size = 0.5;
	int marker_style = 20;

    gr_pure_ph_M30_PL_Cu = new TGraph(Ndata_M30);
    copy_gr_pure_ph_M30_PL_Cu = new TGraph(Ndata_M30);
    gr_pure_ph_M30_PL_Cu->SetTitle("PL vs Cu for Pure M30 Photons; Cu; PL");
    gr_pure_ph_M30_PL_Cu->SetMarkerColor(kCyan);
    copy_gr_pure_ph_M30_PL_Cu->SetMarkerColor(kCyan+2);
    gr_pure_ph_M30_PL_Cu->SetMarkerStyle(marker_style);
    gr_pure_ph_M30_PL_Cu->SetMarkerSize(marker_size);
    gr_pure_ph_M30_PL_Cu->SetFillStyle(0);
    gr_pure_ph_M30_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_M30; i++){
    	gr_pure_ph_M30_PL_Cu->SetPoint(i, pure_ph_M30_Cu[i], pure_ph_M30_PL[i]);
    	copy_gr_pure_ph_M30_PL_Cu->SetPoint(i, pure_ph_M30_Cu[i], pure_ph_M30_PL[i]);
    }
    gr_pure_ph_M30_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_pure_ph_M30_PL_Cu->GetXaxis()->SetLimits(0, 600);

    gr_pure_ph_L50_PL_Cu = new TGraph(Ndata_L50);
    copy_gr_pure_ph_L50_PL_Cu = new TGraph(Ndata_L50);
    gr_pure_ph_L50_PL_Cu->SetTitle("PL vs Cu for Pure L50 Photons; Cu; PL");
    gr_pure_ph_L50_PL_Cu->SetMarkerColor(kBlue);
    copy_gr_pure_ph_L50_PL_Cu->SetMarkerColor(kBlue+2);
    gr_pure_ph_L50_PL_Cu->SetMarkerStyle(marker_style);
    gr_pure_ph_L50_PL_Cu->SetMarkerSize(marker_size);
    gr_pure_ph_L50_PL_Cu->SetFillStyle(0);
    gr_pure_ph_L50_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_L50; i++){
    	gr_pure_ph_L50_PL_Cu->SetPoint(i, pure_ph_L50_Cu[i], pure_ph_L50_PL[i]);
    	copy_gr_pure_ph_L50_PL_Cu->SetPoint(i, pure_ph_L50_Cu[i], pure_ph_L50_PL[i]);
    }
    gr_pure_ph_L50_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_pure_ph_L50_PL_Cu->GetXaxis()->SetLimits(0, 600);


    gr_pure_ph_Cs_PL_Cu = new TGraph(Ndata_Cs);
    copy_gr_pure_ph_Cs_PL_Cu = new TGraph(Ndata_Cs);
    gr_pure_ph_Cs_PL_Cu->SetTitle("PL vs Cu for Pure Cs Photons; Cu; PL");
    gr_pure_ph_Cs_PL_Cu->SetMarkerColor(kMagenta);
    copy_gr_pure_ph_Cs_PL_Cu->SetMarkerColor(kMagenta+2);
    gr_pure_ph_Cs_PL_Cu->SetMarkerStyle(marker_style);
    gr_pure_ph_Cs_PL_Cu->SetMarkerSize(marker_size);
    gr_pure_ph_Cs_PL_Cu->SetFillStyle(0);
    gr_pure_ph_Cs_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_Cs; i++){
    	gr_pure_ph_Cs_PL_Cu->SetPoint(i, pure_ph_Cs_Cu[i], pure_ph_Cs_PL[i]);
    	copy_gr_pure_ph_Cs_PL_Cu->SetPoint(i, pure_ph_Cs_Cu[i], pure_ph_Cs_PL[i]);
    }
    gr_pure_ph_Cs_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_pure_ph_Cs_PL_Cu->GetXaxis()->SetLimits(0, 600);


    gr_pure_b_kr_PL_Cu = new TGraph(Ndata_kr);
    copy_gr_pure_b_kr_PL_Cu = new TGraph(Ndata_kr);
    gr_pure_b_kr_PL_Cu->SetTitle("PL vs Cu for Pure Kr85 Betas; Cu; PL");
    gr_pure_b_kr_PL_Cu->SetMarkerColor(kRed);
    copy_gr_pure_b_kr_PL_Cu->SetMarkerColor(kRed+2);
    gr_pure_b_kr_PL_Cu->SetMarkerStyle(marker_style);
    gr_pure_b_kr_PL_Cu->SetMarkerSize(marker_size);
    gr_pure_b_kr_PL_Cu->SetFillStyle(0);
    gr_pure_b_kr_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_kr; i++){
    	gr_pure_b_kr_PL_Cu->SetPoint(i, pure_b_kr_Cu[i], pure_b_kr_PL[i]);
    	copy_gr_pure_b_kr_PL_Cu->SetPoint(i, pure_b_kr_Cu[i], pure_b_kr_PL[i]);
    }
    gr_pure_b_kr_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_pure_b_kr_PL_Cu->GetXaxis()->SetLimits(0, 600);


    gr_pure_b_sr_PL_Cu = new TGraph(Ndata_sr);
    copy_gr_pure_b_sr_PL_Cu = new TGraph(Ndata_sr);
    gr_pure_b_sr_PL_Cu->SetTitle("PL vs Cu for Pure Sr90 Betas; Cu; PL");
    gr_pure_b_sr_PL_Cu->SetMarkerColor(kOrange);
    copy_gr_pure_b_sr_PL_Cu->SetMarkerColor(kOrange+2);
    gr_pure_b_sr_PL_Cu->SetMarkerStyle(marker_style);
    gr_pure_b_sr_PL_Cu->SetMarkerSize(marker_size);
    gr_pure_b_sr_PL_Cu->SetFillStyle(0);
    gr_pure_b_sr_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_sr; i++){
    	gr_pure_b_sr_PL_Cu->SetPoint(i, pure_b_sr_Cu[i], pure_b_sr_PL[i]);
    	copy_gr_pure_b_sr_PL_Cu->SetPoint(i, pure_b_sr_Cu[i], pure_b_sr_PL[i]);
    }
    gr_pure_b_sr_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_pure_b_sr_PL_Cu->GetXaxis()->SetLimits(0, 600);


    gr_mixed_PL_Cu = new TGraph(Ndata_mixed);
    copy_gr_mixed_PL_Cu = new TGraph(Ndata_mixed);
    gr_mixed_PL_Cu->SetTitle("PL vs Cu for Mixed Sources; Cu; PL");
    gr_mixed_PL_Cu->SetMarkerColor(kGreen);
    copy_gr_mixed_PL_Cu->SetMarkerColor(kGreen+2);
    gr_mixed_PL_Cu->SetMarkerStyle(21);
    gr_mixed_PL_Cu->SetMarkerSize(0.125);
    gr_mixed_PL_Cu->SetFillStyle(0);
    gr_mixed_PL_Cu->SetLineColor(0);
    for( int i = 0; i < Ndata_mixed; i++){
    	gr_mixed_PL_Cu->SetPoint(i, mixed_Cu[i], mixed_PL[i]);
    	copy_gr_mixed_PL_Cu->SetPoint(i, mixed_Cu[i], mixed_PL[i]);
    }
    gr_mixed_PL_Cu->GetYaxis()->SetRangeUser(0, 5200);
    gr_mixed_PL_Cu->GetXaxis()->SetLimits(0, 600);


    TCanvas *c_PL_Cu = new TCanvas("c_PL_Cu", "PL vs Cu for Various Source Types", 600, 500);
    gr_mixed_PL_Cu->Draw("AP");

    gr_pure_ph_M30_PL_Cu->Draw("Psame");
    gr_pure_ph_L50_PL_Cu->Draw("Psame");
    gr_pure_ph_Cs_PL_Cu->Draw("Psame");

    gr_pure_b_kr_PL_Cu->Draw("Psame");
    gr_pure_b_sr_PL_Cu->Draw("Psame");

    copy_gr_mixed_PL_Cu->Draw("Psame");
    copy_gr_pure_ph_M30_PL_Cu->Draw("Psame");
    copy_gr_pure_ph_L50_PL_Cu->Draw("Psame");
    copy_gr_pure_ph_Cs_PL_Cu->Draw("Psame");

    copy_gr_pure_b_kr_PL_Cu->Draw("Psame");
    copy_gr_pure_b_sr_PL_Cu->Draw("Psame");

    TLegend *tleg_PL_Cu = new TLegend(0.7, 0.7, 0.85, 0.85);
    tleg_PL_Cu->AddEntry(gr_pure_b_kr_PL_Cu, "Pure Kr85");
    tleg_PL_Cu->AddEntry(gr_pure_b_sr_PL_Cu, "Pure Sr90");
    tleg_PL_Cu->AddEntry(gr_pure_ph_M30_PL_Cu, "Pure M30");
    tleg_PL_Cu->AddEntry(gr_pure_ph_L50_PL_Cu, "Pure L50");
    tleg_PL_Cu->AddEntry(gr_pure_ph_Cs_PL_Cu, "Pure Cs");
    tleg_PL_Cu->AddEntry(gr_mixed_PL_Cu, "Mixed Photons Betas");
    tleg_PL_Cu->SetFillStyle(0);
    tleg_PL_Cu->Draw("same");

    tf_poly = new TF1("tf_poly", source_ind_poly, 0, 600, 6);
    // tf_poly->FixParameter(0, 311.6);
    // tf_poly->FixParameter(1, -2.6);
    // tf_poly->FixParameter(2, 0.0183);
    // tf_poly->FixParameter(3, 534);
    // tf_poly->SetParameter(0, 11521.3);
    // tf_poly->SetParameter(1, -156);
    // tf_poly->SetParameter(2, 1.0012);
    // tf_poly->SetParameter(3, -0.003141);
    // tf_poly->SetParameter(4, 4.85e-06);
    // tf_poly->SetParameter(5, -2.9946e-09);
    // tf_poly->SetLineColor(kRed);
    // tf_poly->Draw("same");

    // tf_beta_upper = new TF1("tf_beta_upper", pure_beta_upper, 0, 280, 2);
    // tf_beta_upper->FixParameter(0, 9.2329);
    // tf_beta_upper->FixParameter(1, 17.5);
    // tf_beta_upper->SetLineColor(kBlack);
    // tf_beta_upper->SetLineStyle(9);
    // tf_beta_upper->SetLineWidth(1.5);
    // tf_beta_upper->Draw("same");

    // tf_beta_lower = new TF1("tf_beta_lower", pure_beta_lower, 0, 280, 2);
    // tf_beta_lower->FixParameter(0, 9.2329);
    // tf_beta_lower->FixParameter(1, -250);
    // tf_beta_lower->SetLineColor(kBlack);
    // tf_beta_lower->SetLineStyle(9);
    // tf_beta_lower->SetLineWidth(1.5);
    // tf_beta_lower->Draw("same");



    boundary_n = sizeof(arr_boundary_x)/sizeof(double);

    gr_boundary = new TGraph(boundary_n);

    for( int i = 0; i < boundary_n; i++ ){
        boundary_x.push_back(arr_boundary_x[i]);
        boundary_y.push_back(arr_boundary_y[i]);
        gr_boundary->SetPoint(i, arr_boundary_x[i], arr_boundary_y[i]);
    }

    tf_boundary->FixParameter(0, 10569);
    tf_boundary->FixParameter(1, -184.612);
    tf_boundary->FixParameter(2, 1.5375);
    tf_boundary->FixParameter(3, -0.00621976);
    tf_boundary->FixParameter(4, 1.21312e-005);
    tf_boundary->FixParameter(5, -9.1689e-009);
    gr_boundary->Fit("tf_boundary", "RMSE");
    tf_boundary->Draw("same");

    TLine *tl_beta_upper = new TLine(0, 338.81, 37.895, 338.81);
    tl_beta_upper->SetLineStyle(8);
    tl_beta_upper->SetLineColor(kMagenta);
    tl_beta_upper->Draw("same");

    TLine *tl_beta_right = new TLine(37.895, 0, 37.895, 338.81);
    tl_beta_right->SetLineStyle(8);
    tl_beta_right->SetLineColor(kMagenta);
    tl_beta_right->Draw("same");

}