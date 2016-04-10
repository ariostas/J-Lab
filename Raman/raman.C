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
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TSpectrum.h"
#include <algorithm>

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
TH1D* readData(TString, TString, Double_t, Double_t);
void findPeaks(TH1D*);
void fitPeaks(TH1D*, vector<Double_t>, Double_t);
void printPeaks(TH1D*);
void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);

TGraphErrors *graphLarge, *graphSmall, *allPeaks;

/*
 * MAIN FUNCTION
 */

void raman(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    // TH1D *oxygen1 = readData("oxygen_0psi_25cm_mar15", "Oxygen", 25, 180.2);
    // TH1D *oxygen2 = readData("oxygen_15psi_25cm_mar15", "Oxygen", 25, 218.33);
    // TH1D *oxygen3 = readData("oxygen_30psi_25cm_mar15", "Oxygen", 25, 0);
    // TH1D *oxygen4 = readData("oxygen_0psi_25cm_mar24", "Oxygen", 25, 0);

    // TH1D *nitrogen1 = readData("nitrogen_15psi_25cm_mar17", "Nitrogen", 25, 0);
    // TH1D *nitrogen2 = readData("nitrogen_30psi_5cm_mar17", "Nitrogen", 5, 0);
    TH1D *nitrogen3 = readData("test_117", "Nitrogen", 5, 190.04);

    // TH1D *carbondio = readData("co2_0psi_5cm_mar24", "CO2", 5, 0);

    // findPeaks(oxygen1);
    // findPeaks(oxygen2);
    // findPeaks(oxygen3);
    // findPeaks(oxygen4);

    // findPeaks(nitrogen1);
    // findPeaks(nitrogen2);
    // findPeaks(nitrogen3);

    // findPeaks(carbondio);

    // double peaksOxygen1Arr[] = {-130.41, -118.74, -106.87, -95.62, -84.16, -72.49, -60.83, -48.95, -37.28, -25.41, 0, 26.26, 38.13, 49.80, 61.67, 73.34, 85.22, 96.88, 108.55, 119.80, 131.47};
    // vector<Double_t> peaksOxygen1 (peaksOxygen1Arr, peaksOxygen1Arr + sizeof(peaksOxygen1Arr) / sizeof(double) );
    // fitPeaks(oxygen1, peaksOxygen1, 5.);

    double peaksNitrogenArr[] = {-173.87, -165.58, -157.29, -149, -141.41, -133.08, -125.12, -117.08, -108.96, -100.58, -92.83, -85.33, -77.12, -69.08, -60.75, -52.37, -44.58, -36.36, -28.25, -20, 0, 20.00, 28.21, 36.33, 44.29, 52.88, 60.79, 69.00, 77.08, 85.08, 93.04, 101.25, 109.00, 117.00, 125.38, 132.54, 141.38, 149.67, 157.96};
    vector<Double_t> peaksNitrogen(peaksNitrogenArr, peaksNitrogenArr + sizeof(peaksNitrogenArr) / sizeof(double) );
    fitPeaks(nitrogen3, peaksNitrogen, 3.5);
    printPeaks(nitrogen3);


    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    // histogram(oxygen1, "Oxygen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen1", "line");
    // histogram(oxygen2, "Oxygen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen2", "line");
    // histogram(oxygen3, "Oxygen 3 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen3", "line");
    // histogram(oxygen4, "Oxygen 1 atm (fixed slit)", can, "#Delta k (cm^{-1})", "Photon count", "oxygen4", "line");

    // histogram(nitrogen1, "Nitrogen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen1", "line");
    // histogram(nitrogen2, "Nitrogen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen2", "line");
    histogram(nitrogen3, "Nitrogen 1 atm (small slits)", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen3", "line");

    // histogram(carbondio, "CO2 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "carbondio", "line");


    graph(graphSmall, "Small peaks", can, "#Delta k (cm^{-1})", "Photon count", "PeaksSmall");
    graph(graphLarge, "Large peaks", can, "#Delta k (cm^{-1})", "Photon count", "PeaksLarge");
    graph(allPeaks, "All peaks", can, "#Delta k (cm^{-1})", "J", "PeaksAll");

}

TH1D* readData(TString dataset, TString name, Double_t speed, Double_t centralLine){

    TString fileName = "./Data/" + dataset + ".txt";

    Int_t nBins = 0;
    Double_t offset = 0;

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}

    Double_t entry, intensity;
    
    bool offsetDone = false;
    while(ifs >> entry >> intensity){

        if(!offsetDone){offset = entry; offsetDone = true;}
        nBins++;

    }

    ifs.close();

    Double_t minHist, maxHist;

    if(centralLine == 0){

        minHist = +0.5*speed/60./2.;
        maxHist = (nBins+0.5)*speed/60./2.;

    }
    else{

        minHist = -centralLine+0.5*speed/60./2.;
        maxHist = -centralLine+(nBins+0.5)*speed/60./2.;
    }

    TH1D *histo = new TH1D(dataset, dataset, nBins/5, minHist, maxHist);


    ifs.open(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}
    
    while(ifs >> entry >> intensity){
        
        for(Int_t x = 0; x< intensity*1828.5323/2.; x++){
            if(intensity*1828.5323/2. > 2000) intensity = 150./1828.5323/2.;
            histo->Fill(-centralLine+(nBins/2.-entry+offset)*speed/60.);
        }

    }

    ifs.close();

    return histo;

}

void findPeaks(TH1D* hist){

    TSpectrum *s = new TSpectrum(60);
    Int_t nfound = s->Search(hist,8,"",0.005);
    TH1 *hb = s->Background(hist,20,"same");
    // if (hb) c1->Update();

    TF1 *fline = new TF1("fline","pol1",0,1000);
    hist->Fit(fline,"qn");

    Int_t npeaks = 0;
    Float_t *xpeaks = s->GetPositionX();
    std::sort(xpeaks, xpeaks + 50);
    for (Int_t p=0;p<nfound;p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = hist->GetXaxis()->FindBin(xp);
        Double_t yp = hist->GetBinContent(bin);
        // if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
        npeaks++;
        printf("%s: peak %i at position %4.2f\n",hist->GetTitle(), p, xp);
    }
    cout << endl;

}

void fitPeaks(TH1D* histo, vector<Double_t> peaks, Double_t width){

    if(!histo){
        cout << "Histogram not found" << endl;
        return;
    }

    for(UInt_t x = 0; x<peaks.size(); x++){

        Double_t height = histo->GetBinContent(histo->FindBin(peaks.at(x)));
        TF1* f = new TF1(TString::Format("Peak_%i", x), "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
        if(height > 1500){
            f->SetParameter(0,  250);
            f->SetParameter(1,  0.004*height);
            f->SetParameter(2, peaks.at(x));
            f->SetParameter(3,  1);
            f->SetParameter(4, 1);
        }
        else if(height > 750){
            f->SetParameter(0,  180);
            f->SetParameter(1,  0.004*height);
            f->SetParameter(2, peaks.at(x));
            f->SetParameter(3,  1);
            f->SetParameter(4, 1);
        }
        else{
            f->SetParameter(0,  100);
            f->SetParameter(1,  0.0025*height);
            f->SetParameter(2, peaks.at(x));
            f->SetParameter(3,  1);
            f->SetParameter(4, 1);
        }
        if(Abs(peaks.at(x)) < 5){
            f->SetParameter(0,  50);
            f->SetParameter(1,  0.004*height);
            f->SetParameter(2, peaks.at(x));
            f->SetParameter(3,  1);
            f->SetParameter(4, 1);
        }
        
        // f->SetParameter(1, 0.009*histo->GetBinContent(histo->FindBin(peaks.at(x))));
        // f->SetParameter(2, peaks.at(x));
        // f->SetParameter(3, 1);
        // f->SetParameter(4, 1);

        histo->Fit(f, "ME+", "", peaks.at(x)-width, peaks.at(x)+width);

    }
}

void printPeaks(TH1D* histo){

    if(!histo){
        cout << "Histogram not found" << endl;
        return;
    }

    Double_t xLarge[20], xErrorsLarge[20], xSmall[20], xErrorsSmall[20];
    Double_t yLarge[20], yErrorsLarge[20], ySmall[20], yErrorsSmall[20];
    Int_t nxLarge = 0, nxSmall = 0;
    Double_t xPeak[40], xErrorsPeak[40], nPeak[40], nErrorsPeak[40];
    Int_t npeaks = 0, central=0;

    for(Int_t x = 0; x < 1000; x++){

        TF1* func = histo->GetFunction(TString::Format("Peak_%i", x));
        if(!func) break;

        xPeak[x] = func->GetParameter(2);
        nPeak[x] = x-20;
        xErrorsPeak[x] = func->GetParError(2);
        nErrorsPeak[x] = 0.1;
        npeaks++;

        if((x%2 == 0) && (Abs(func->GetParameter(2)) > 5.)){

            xLarge[x/2-central] = func->GetParameter(2);
            yLarge[x/2-central] = func->Eval(func->GetParameter(2));
            xErrorsLarge[x/2-central] = func->GetParError(2);
            yErrorsLarge[x/2-central] = Sqrt(func->Eval(func->GetParameter(2)));
            nxLarge++;

        }
        else if(x%2 == 1 && Abs(func->GetParameter(2)) > 5){

            xSmall[(x-1)/2] = func->GetParameter(2);
            ySmall[(x-1)/2] = func->Eval(func->GetParameter(2));
            xErrorsSmall[(x-1)/2] = func->GetParError(2);
            yErrorsSmall[(x-1)/2] = Sqrt(func->Eval(func->GetParameter(2)));
            nxSmall++;

        }
        if(Abs(func->GetParameter(2)) < 5.) central = 1;

        cout << "Peak " << x << " at " << func->GetParameter(2) << " +- " << func->GetParError(2) << " with height " << func->Eval(func->GetParameter(2)) << " +- " << Sqrt(func->Eval(func->GetParameter(2))) << endl;

    }

    graphLarge = new TGraphErrors(nxLarge, xLarge,  yLarge, xErrorsLarge, yErrorsLarge);
    graphSmall = new TGraphErrors(nxSmall, xSmall,  ySmall, xErrorsSmall, yErrorsSmall);
    allPeaks = new TGraphErrors(npeaks, xPeak, nPeak, xErrorsPeak, nErrorsPeak);

}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );
    gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("4.3f");

    if(!histoData) return;

    // TF1 *f = histoData->GetFunction("4");

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.6);
        histoData->SetMarkerColor(kRed);
        histoData->SetLineColor(kBlack);
        // f->SetLineColor(kBlue);
        // f->SetLineWidth(4);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        // f->SetLineColor(kRed);
        // f->SetLineWidth(2);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.755,0.675,0.985,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data", (type.Contains("dot") ? "lep" : "l"));
    // leg->AddEntry(f, "Voigtian fits","l");
    
    histoData->Draw((type.Contains("dot") ? "E1" : ""));
    // leg->Draw("same");

    // add axis labels
    histoData->GetXaxis()->SetTitle(xTitle);
    histoData->GetXaxis()->CenterTitle();
    histoData->GetXaxis()->SetTitleSize(0.055);
    histoData->GetXaxis()->SetTitleOffset(0.84);
    histoData->GetXaxis()->SetLabelOffset(-0.001);
    histoData->GetXaxis()->SetLabelSize(0.05);
    histoData->GetYaxis()->SetTitle(yTitle);
    histoData->GetYaxis()->CenterTitle();
    histoData->GetYaxis()->SetTitleSize(0.055);
    histoData->GetYaxis()->SetTitleOffset(1.0);
    histoData->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");
}

void histogram(TH1D *histoData1, TH1D *histoData2, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    histoData1->Sumw2();
    histoData1->Scale(10.0/histoData1->Integral());
    histoData1->SetLineWidth(2);
    histoData1->SetMarkerStyle(20);
    histoData1->SetMarkerSize(1.5);
    histoData1->SetMarkerColor(kRed);
    histoData1->SetLineColor(kBlack);
    // TF1 *f1 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f1->SetParameter(0, 200);
    // f1->SetParameter(1, histoData1->GetMaximum()/10.);
    // f1->SetParameter(2, histoData1->GetMean());
    // f1->SetParameter(3, .05);
    // f1->SetParameter(4, .05);
    // f1->SetParName(0, "Offset");
    // f1->SetParName(1, "Normalization");
    // f1->SetParName(2, "Mean");
    // f1->SetParName(3, "Sigma");
    // f1->SetParName(4, "gamma");
    // f1->SetLineColor(kBlue);
    // f1->SetLineWidth(4);
    // histoData1->Fit(f1, "ME");
    histoData1->Fit("gaus", "0");

    histoData2->Sumw2();
    histoData2->Scale(10.0/histoData2->Integral());
    histoData2->SetLineWidth(2);
    histoData2->SetMarkerStyle(21);
    histoData2->SetMarkerSize(1.5);
    histoData2->SetMarkerColor(kBlue);
    histoData2->SetLineColor(kBlack);
    // TF1 *f2 = new TF1(histName + "1", "[0]+[1]*TMath::Voigt(x-[2],[3],[4])");
    // f2->SetParameter(0, 200);
    // f2->SetParameter(1, histoData1->GetMaximum()/10.);
    // f2->SetParameter(2, histoData1->GetMean());
    // f2->SetParameter(3, .05);
    // f2->SetParameter(4, .05);
    // f2->SetParName(0, "Offset");
    // f2->SetParName(1, "Normalization");
    // f2->SetParName(2, "Mean");
    // f2->SetParName(3, "Sigma");
    // f2->SetParName(4, "gamma");
    // f2->SetLineColor(kBlue);
    // f2->SetLineWidth(4);
    // histoData2->Fit(f2, "ME");
    histoData2->Fit("gaus", "0");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.15,0.6,0.385,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData1, "Hydrogen","lep");
    leg->AddEntry(histoData2, "Deuterium","lep");
    // leg->AddEntry(f, "Voigt fit","l");

    Double_t max = histoData1->GetMaximum();
    if(histoData2->GetMaximum() > max) max = histoData2->GetMaximum();
    max*=1.2;
    histoData1->SetMaximum(max);
    histoData1->SetMinimum(0.);
    histoData2->SetMaximum(max);
    histoData2->SetMinimum(0.);
    
    histoData1->Draw("E1");
    histoData2->Draw("same E1");
    leg->Draw("same");

    // add axis labels
    histoData1->GetXaxis()->SetTitle(xTitle);
    histoData1->GetXaxis()->CenterTitle();
    histoData1->GetXaxis()->SetTitleSize(0.055);
    histoData1->GetXaxis()->SetTitleOffset(0.90);
    histoData1->GetXaxis()->SetLabelOffset(0.010);
    histoData1->GetXaxis()->SetLabelSize(0.05);
    histoData1->GetYaxis()->SetTitle(yTitle);
    histoData1->GetYaxis()->CenterTitle();
    histoData1->GetYaxis()->SetTitleSize(0.055);
    histoData1->GetYaxis()->SetTitleOffset(1.0);
    histoData1->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histoData1->SetTitle(histName);
    can->Update();

    can->SaveAs(name + ".png");


}

void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    // TF1 *f = new TF1(graphName, "[0]+[1]*TMath::Power(TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])),0.5)");
    TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    f->SetParameter(0, 250);
    f->SetParameter(1, 100);
    f->SetParameter(2, 0.395);
    f->SetParameter(3, 1);
    f->SetParameter(4, 10000);
    f->SetParameter(5, 100);
    f->SetParName(0, "Offset");
    f->SetParName(1, "Peak intensity");
    f->SetParName(2, "Aperture");
    f->SetParName(3, "Center value");
    f->SetParName(4, "Amp");
    f->SetParLimits(0, 100, 300);
    f->SetParLimits(1, 10, 100);
    f->SetParLimits(2, 0.1, 1);
    f->SetParLimits(3, -0.2, 2);
    f->SetParLimits(4, 50, 10000);
    f->SetParLimits(5, 1, 10000);
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    //leg->AddEntry(f, "Airy disk fit","l");
    
    graph->Draw("ap");
    //leg->Draw("same");

    // add axis labels
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetYaxis()->SetTitle(yTitle);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(0.9);
    graph->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    graph->SetTitle(graphName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}
