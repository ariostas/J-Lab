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
#include <TPaveStats.h>
#include <TText.h>
#include <TLatex.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);

// Initialize histograms
TH1D *histPoisson = new TH1D("histCubeX", "histCubeX", 50, 50, 150);

/*
 * MAIN FUNCTION
 */

void poisson(TString dataSet = "100"){

    TString inputFile = "poissonData" + dataSet + "ps1s_2.txt";

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}

    Double_t count;
    
    while(ifs >> count){
        
        histPoisson->Fill(count);

    }
    
    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1920, 1080);

    gStyle->SetOptStat(2210);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1111);


    histogram(histPoisson, "", c1, "Number of decays", "Count", "poisson_" + dataSet);
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histo->SetLineWidth(3);
    histo->Draw("E1");

    TF1 *f1 = new TF1(histName,TString::Format("%4.3f*TMath::Gaus(x,%4.3f,%4.3f,1)",
        histo->Integral()*histo->GetBinWidth(1), histo->GetMean(), histo->GetStdDev()),
        histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax()); 
    f1->SetLineWidth(3);
    f1->SetLineColor(kRed);
    //f1->Draw("same");

    TPaveText *pt = new TPaveText(0.605,0.675,0.885,0.875, "brNDC");
    // pt->AddText(TString::Format("Entries = %3i", histo->Integral()));
    pt->AddText("Total entries = 103");
    pt->AddText(TString::Format("Mean = %4.1f #pm %4.1f_{stat} #pm 0.5_{syst} mm", histo->GetMean(), Max(histo->GetMeanError(),0.1)));
    pt->AddText(TString::Format("Std dev = %4.2f mm", histo->GetStdDev()));
    //pt->Draw();

    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetXaxis()->SetTitleSize(0.045);
    histo->GetXaxis()->SetTitleOffset(1.05);
    histo->GetXaxis()->SetLabelOffset(0.010);
    histo->GetXaxis()->SetLabelSize(0.05);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->GetYaxis()->SetTitleSize(0.045);
    histo->GetYaxis()->SetTitleOffset(0.70);
    histo->GetYaxis()->SetLabelSize(0.05);
    histo->SetTitle(histName); // title on top
    can->Update();

    can->SaveAs(name + ".png");
}
