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
void histogram(Int_t, TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(Int_t, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void analyze(Int_t, TH1D*, TH1D*);

// Initialize histograms
TH1D *histPoisson1 = new TH1D("histPoisson1", "histPoisson1", 7, -.5, 6.5);
TH1D *histPoisson5 = new TH1D("histPoisson5", "histPoisson5", 16, -.5, 15.5);
TH1D *histPoisson10 = new TH1D("histPoisson10", "histPoisson10", 22, -.5, 21.5);
TH1D *histPoisson100 = new TH1D("histPoisson100", "histPoisson100", 25, 50, 150);

TH1D *histMC1 = new TH1D("histMC1", "histMC1", 7, -.5, 6.5);
TH1D *histMC5 = new TH1D("histMC5", "histMC5", 16, -.5, 15.5);
TH1D *histMC10 = new TH1D("histMC10", "histMC10", 22, -.5, 21.5);
TH1D *histMC100 = new TH1D("histMC100", "histMC100", 25, 50, 150);

/*
 * MAIN FUNCTION
 */

 void poisson(TString dataSet = "10"){

    cout << "\n\nStarting process...\n\n";

    TH1::StatOverflows(kTRUE);
    gStyle->SetOptStat(2210);
    //gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1100);

    analyze(1, histPoisson1, histMC1);
    analyze(5, histPoisson5, histMC5);
    analyze(10, histPoisson10, histMC10);
    analyze(100, histPoisson100, histMC100);

    TCanvas *can = new TCanvas("Histogram", "Histogram", 1920, 1080);
    can->Divide(2,2);

    can->cd(1);
    // histogram(5, histMC5, "5 Hz simulation", can, "Number of decays", "Count", "poisson_mc");
    histogram(1, histPoisson1, histMC1, "1 Hz run", can, "Number of decays", "Count", "poisson_1");
    can->cd(2);
    histogram(5, histPoisson5, histMC5, "5 Hz run", can, "Number of decays", "Count", "poisson_5");
    can->cd(3);
    histogram(10, histPoisson10, histMC10, "10 Hz run", can, "Number of decays", "Count", "poisson_10");
    can->cd(4);
    histogram(100, histPoisson100, histMC100, "100 Hz run", can, "Number of decays", "Count", "poisson_100");

 }

void analyze(Int_t dataSet = 1, TH1D *histoData = histPoisson1, TH1D *histoMC = histMC1){

    TString inputFile;
    if(dataSet < 10) inputFile = TString::Format("poissonData%1ips1s.txt", dataSet);
    else inputFile = TString::Format("poissonData%1ips1s_2.txt", dataSet);
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}

    Double_t count;
    
    while(ifs >> count){
        
        histoData->Fill(count);

    }

    ifs.close();

    inputFile = TString::Format("poissonMC%1ips.txt", dataSet);

    ifs.open(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}
    
    while(ifs >> count){

        histoMC->Fill(count);

    }

    ifs.close();
    
    Double_t scale = (dataSet < 10 ? 1.0/10.0 : 1.0/5.0);
    histoMC->Sumw2();
    histoMC->Scale(scale);
    
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(Int_t dataSet, TH1D *histoData, TH1D *histoMC, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoData->SetLineWidth(2);
    histoData->SetMarkerStyle(20);
    histoData->SetMarkerSize(1.4);
    histoData->SetMarkerColor(kRed);
    histoData->SetLineColor(kBlack);
    //TF1 *f = new TF1(histName, "[0]*TMath::Poisson(x,[1])", 0, histoData->GetXaxis()->GetXmax());
    TF1 *f = new TF1(histName, "gaus", 0, histoData->GetXaxis()->GetXmax());
    f->SetParameters(0, histoData->Integral()*histoData->GetBinWidth(1));
    f->SetParameters(1, dataSet);
    f->SetParameters(2, TMath::Sqrt(dataSet));
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    histoData->Fit(f);

    histoMC->SetLineWidth(2);
    histoMC->SetMarkerStyle(22);
    histoMC->SetMarkerSize(1.8);
    histoMC->SetMarkerColor(kBlue);
    histoMC->SetFillColor(kOrange-9);

    TF1 *f1 = new TF1(histName,TString::Format("%4.3f*TMath::Poisson(x,%1i)",
        histoData->Integral()*histoData->GetBinWidth(1), dataSet),
        0, histoData->GetXaxis()->GetXmax()); 
    f1->SetLineWidth(4);
    f1->SetLineColor(kBlue+1);
    f1->SetLineStyle(7);

    Double_t max = histoData->GetMaximum();
    if(histoMC->GetMaximum() > max) max = histoMC->GetMaximum();
    max*=1.25;
    histoMC->SetMaximum(max);
    histoMC->SetMinimum(0.);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.65,0.6,0.885,0.875);
    //leg->SetTextFont(72);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data","lep");
    //leg->AddEntry(histoMC, "Monte Carlo","pf");
    //leg->AddEntry(f1, "#splitline{Theoretical}{distribution}","l");
    leg->AddEntry(f, "Gaussian fit","l");
    
    histoData->Draw("E1");
    //histoMC->Draw("E2 same");
    //f1->Draw("same");
    //histoData->Draw("E1 same");
    leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.87);
    histoData->GetXaxis()->SetLabelOffset(0.010);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(0.70);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
    can->SaveAs(name + ".pdf");
}

void histogram(Int_t dataSet, TH1D *histoMC, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoMC->SetLineWidth(2);
    histoMC->SetMarkerStyle(20);
    histoMC->SetMarkerSize(1.4);
    histoMC->SetMarkerColor(kRed);
    histoMC->SetLineColor(kBlack);

    TF1 *f1 = new TF1(histName,TString::Format("%4.3f*TMath::Poisson(x,%1i)",
        histoMC->Integral()*histoMC->GetBinWidth(1), dataSet),
        0, histoMC->GetXaxis()->GetXmax()); 
    f1->SetLineWidth(3);
    f1->SetLineColor(kBlue+1);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.65,0.5,0.885,0.875);
    //leg->SetTextFont(72);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoMC, "Monte Carlo","lep");
    leg->AddEntry(f1, "#splitline{Theoretical}{distribution}","l");
    
    histoMC->Draw("E1");
    f1->Draw("same");
    leg->Draw("same");

    // add axis labels
    histoMC->GetXaxis()->SetTitle(xTitle);
    histoMC->GetXaxis()->CenterTitle();
    histoMC->GetXaxis()->SetTitleSize(0.055);
    histoMC->GetXaxis()->SetTitleOffset(0.85);
    histoMC->GetXaxis()->SetLabelOffset(0.010);
    histoMC->GetXaxis()->SetLabelSize(0.05);
    histoMC->GetYaxis()->SetTitle(yTitle);
    histoMC->GetYaxis()->CenterTitle();
    histoMC->GetYaxis()->SetTitleSize(0.055);
    histoMC->GetYaxis()->SetTitleOffset(0.70);
    histoMC->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoMC->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
    can->SaveAs(name + ".pdf");
}
