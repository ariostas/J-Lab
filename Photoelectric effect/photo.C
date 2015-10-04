#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <TString.h>
#include "TMath.h"
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TGraphErrors* analyze(TString);

/*
 * MAIN FUNCTION
 */

void photo(){
    
    cout << "\n\nStarting process...\n\n";

    TGraphErrors* photo365 = analyze("photoData365-0_1.txt");
    TGraphErrors* photo404 = analyze("photoData404-7.txt");
    TGraphErrors* photo435 = analyze("photoData435-8.txt");
    TGraphErrors* photo546 = analyze("photoData546-1.txt");
    TGraphErrors* photo577 = analyze("photoData577-0.txt");

    TCanvas *can = new TCanvas("Photoelectric", "Photoelectric", 1920, 1080);
    can->SetGrid();

    photo365->SetMarkerColor(4);
    photo365->SetMarkerStyle(21);
    photo404->SetMarkerColor(4);
    photo404->SetMarkerStyle(20);
    photo435->SetMarkerColor(2);
    photo435->SetMarkerStyle(21);
    photo546->SetMarkerColor(1);
    photo546->SetMarkerStyle(21);
    photo577->SetMarkerColor(1);
    photo577->SetMarkerStyle(20);

    photo365->SetLineWidth(2);
    photo404->SetLineWidth(2);
    photo435->SetLineWidth(2);
    photo546->SetLineWidth(2);
    photo577->SetLineWidth(2);
    photo365->SetLineColor(kBlack);
    photo404->SetLineColor(kRed);
    photo435->SetLineColor(kBlue);
    photo546->SetLineColor(kGreen);
    photo577->SetLineColor(kMagenta+2);
    // photo365->Draw("ALP");
    // photo404->Draw("ALP");
    // photo435->Draw("ALP");
    // photo546->Draw("ALP");
    // photo577->Draw("same");

    TMultiGraph *photoGraphs = new TMultiGraph();
    photoGraphs->Add(photo365);
    photoGraphs->Add(photo404);
    photoGraphs->Add(photo435);
    photoGraphs->Add(photo546);
    photoGraphs->Add(photo577);
    photoGraphs->Draw("ALP");
    photoGraphs->GetXaxis()->SetTitle("Voltage (V)");
    photoGraphs->GetYaxis()->SetTitle("Normalized Current");
    photoGraphs->SetTitle("Photo electric effect");

    TLegend *leg = new TLegend(0.7,0.6,0.885,0.875);
    leg->SetTextFont(72);
    leg->SetTextSize(0.04);
    leg->SetHeader("Wavelengths");
    leg->AddEntry(photo365, "365.0 nm","lep");
    leg->AddEntry(photo404, "404.7 nm","lep");
    leg->AddEntry(photo435, "435.8 nm","lep");
    leg->AddEntry(photo546, "546.1 nm","lep");
    leg->AddEntry(photo577, "577.0 nm","lep");
    leg->Draw();

    can->SaveAs("photoelectric.png");

    TCanvas *can2 = new TCanvas("Photoelectric2", "Photoelectric2", 1920, 1080);

    Double_t cutVoltage[10], cutWavelength[10];

    cutWavelength[0] = 8.213e14; cutVoltage[0] = 0.978524;
    cutWavelength[1] = 7.408e14; cutVoltage[1] = 0.777271;
    cutWavelength[2] = 6.879e14; cutVoltage[2] = 0.632909;
    cutWavelength[3] = 5.49e14; cutVoltage[3] = 0.178967;
    cutWavelength[4] = 5.196e14; cutVoltage[4] = 0.541412;

    TGraph *cut = new TGraph(5, cutWavelength, cutVoltage);
    cut->SetLineWidth(2);
    cut->SetMarkerColor(4);
    cut->SetMarkerStyle(21);
    TF1 *f = new TF1("f", "[1] * x + [0]");
    cut->Fit(f);
    cut->Draw("ALP");
    cut->GetXaxis()->SetTitle("Frequency (Hz)");
    cut->GetYaxis()->SetTitle("Cutoff voltage (V)");
    cut->SetTitle("Photo electric effect");
    can2->SaveAs("fit.png");


}

/*
 * ANALYZE DATA
 */

TGraphErrors* analyze(TString inputFile = "photoData365-0_2.txt"){
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return NULL;}

    Double_t voltage[100], current[100], voltageErr[100], currentErr[100], tempVoltage, tempCurrent;
    Int_t x=0, x2=0;
    Double_t scaleFactor=0;
    
    while(ifs >> tempVoltage >> tempCurrent){

        if(x == 0){
            scaleFactor = 100./tempCurrent;
        }
        if(scaleFactor*tempCurrent > 30){
            x2++;
        }
        
        voltage[x] = tempVoltage;
        current[x] = scaleFactor*tempCurrent;
        voltageErr[x] = 0.02;
        currentErr[x] = 2.*scaleFactor*(tempCurrent > 20. ? 1. : (tempCurrent > 2. ? 0.1 : 0.01));
        
        x++;
    }

    TGraphErrors *graph = new TGraphErrors(x, voltage, current, voltageErr, currentErr);
    TGraph *graph2 = new TGraphErrors(x2, current, voltage);
    TF1 *f = new TF1("f", "[1] * x + [0]");
    graph2->Fit(f);

    return graph;
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histo->SetLineWidth(3);
    histo->Draw("E1");
    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->SetTitle(histName); // title on top

    can->SaveAs(name + ".jpg");
}
