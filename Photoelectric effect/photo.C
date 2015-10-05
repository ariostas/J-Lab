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

    Double_t cutVoltage[10], cutFrequency[10], cutVoltageErr[10], cutFrequencyErr[10];

    cutFrequency[0] = 8.213e14; cutVoltage[0] = 9.78680e-01;
    cutFrequency[1] = 7.408e14; cutVoltage[1] = 7.77430e-01;
    cutFrequency[2] = 6.879e14; cutVoltage[2] = 6.32989e-01;
    cutFrequency[3] = 5.490e14; cutVoltage[3] = 1.81230e-01;
    cutFrequency[4] = 5.196e14; cutVoltage[4] = 6.07524e-01;

    cutFrequencyErr[0] = 9.051e12/2.0; cutVoltageErr[0] = 2.*8.20815e-03;
    cutFrequencyErr[1] = 7.358e12/2.0; cutVoltageErr[1] = 2.*9.82655e-03;
    cutFrequencyErr[2] = 6.343e12/2.0; cutVoltageErr[2] = 2.*9.83548e-03;
    cutFrequencyErr[3] = 4.036e12/2.0; cutVoltageErr[3] = 2.*2.45016e-02;
    cutFrequencyErr[4] = 3.614e12/2.0; cutVoltageErr[4] = 2.*1.49341e-02;

    TGraphErrors *cut = new TGraphErrors(4, cutFrequency, cutVoltage, cutFrequencyErr, cutVoltageErr);
    cut->SetLineWidth(2);
    cut->SetMarkerColor(4);
    cut->SetMarkerStyle(21);
    TF1 *f = new TF1("f", "[1] * x + [0]");
    f->SetParameter(1, 4e-15);
    f->SetParameter(0, -1.2);
    cut->Fit(f, "m");
    cut->Draw("ALP");
    cut->GetXaxis()->SetTitle("Light Frequency (Hz)");
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
        voltageErr[x] = 0.01;
        currentErr[x] = 2*scaleFactor*(tempCurrent > 20. ? 1. : (tempCurrent > 2. ? 0.1 : 0.01));
        
        x++;
    }

    TGraphErrors *graph = new TGraphErrors(x, voltage, current, voltageErr, currentErr);
    TGraphErrors *graph2 = new TGraphErrors(x2, current, voltage, currentErr, voltageErr);
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
