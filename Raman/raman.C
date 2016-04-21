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
#include <TMultiGraph.h>
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
void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);
void computeRotConstant(TCanvas*);
void computeTemperature(TCanvas*);

TGraphErrors *graphLarge, *graphSmall, *allPeaks;

/*
 * MAIN FUNCTION
 */

void raman(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    // TFile f("histos.root","recreate");

    // TH1D *oxygen1 = readData("oxygen_0psi_25cm_mar15", "Oxygen", 25, 180.2);
    // TH1D *oxygen2 = readData("oxygen_15psi_25cm_mar15", "Oxygen", 25, 218.33);
    // TH1D *oxygen3 = readData("oxygen_30psi_25cm_mar15", "Oxygen", 25, 0);
    // TH1D *oxygen4 = readData("oxygen_0psi_25cm_mar24", "Oxygen", 25, 209);

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
    // double peaksNitrogenArr[] = {-60.75};
    vector<Double_t> peaksNitrogen(peaksNitrogenArr, peaksNitrogenArr + sizeof(peaksNitrogenArr) / sizeof(double) );
    fitPeaks(nitrogen3, peaksNitrogen, 3.5);
    printPeaks(nitrogen3);


    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    // histogram(oxygen1, "Oxygen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen1", "line");
    // histogram(oxygen2, "Oxygen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen2", "line");
    // histogram(oxygen3, "Oxygen 3 atm", can, "#Delta k (cm^{-1})", "Photon count", "oxygen3", "line");
    // histogram(oxygen4, "Oxygen", can, "Wavenumber shift #Deltak (cm^{-1})", "Photon count", "oxygen4", "line");

    // histogram(nitrogen1, "Nitrogen 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen1", "line");
    // histogram(nitrogen2, "Nitrogen 2 atm", can, "#Delta k (cm^{-1})", "Photon count", "nitrogen2", "line");
    histogram(nitrogen3, "Nitrogen", can, "Wavenumber shift #Deltak (cm^{-1})", "Photon count", "nitrogen3", "line");

    // histogram(carbondio, "CO2 1 atm", can, "#Delta k (cm^{-1})", "Photon count", "carbondio", "line");

    computeRotConstant(can);
    computeTemperature(can);

    // nitrogen3->Write();
    // allPeaks->Write();
    // graphSmall->Write();
    // graphLarge->Write();

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
            if(intensity*1828.5323 > 2000) intensity = 150./1828.5323;
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
    Int_t npeaks = 0, central=0, missing=0;

    for(Int_t x = 0; x < 1000; x++){

        TF1* func = histo->GetFunction(TString::Format("Peak_%i", x));
        if(!func) break;

        if(Abs(func->GetParameter(2)) < 0.5) missing++;
        xPeak[x] = func->GetParameter(2);// + 0.255018;
        nPeak[x] = x-21+missing;
        xErrorsPeak[x] = func->GetParError(2) + 0.0230688;
        nErrorsPeak[x] = 0.2;
        npeaks++;
        if(missing==1)missing++;

        if((x%2 == 0) && (Abs(func->GetParameter(2)) > 5.)){

            xLarge[x/2-central] = x-21+missing;
            yLarge[x/2-central] = func->Eval(func->GetParameter(2));
            xErrorsLarge[x/2-central] = 0.2;
            yErrorsLarge[x/2-central] = Sqrt(func->Eval(func->GetParameter(2)));
            nxLarge++;

        }
        else if(x%2 == 1 && Abs(func->GetParameter(2)) > 5){

            xSmall[(x-1)/2] = x-21+missing;
            ySmall[(x-1)/2] = func->Eval(func->GetParameter(2));
            xErrorsSmall[(x-1)/2] = 0.2;
            yErrorsSmall[(x-1)/2] = Sqrt(func->Eval(func->GetParameter(2)));
            nxSmall++;

        }
        if(Abs(func->GetParameter(2)) < 5.) central = 1;

        cout << "Peak " << x << " at " << func->GetParameter(2) << " +- " << func->GetParError(2) << " with height " << func->Eval(func->GetParameter(2)) << " +- " << Sqrt(func->Eval(func->GetParameter(2))) << endl;

    }

    graphLarge = new TGraphErrors(nxLarge, xLarge,  yLarge, xErrorsLarge, yErrorsLarge);
    graphSmall = new TGraphErrors(nxSmall, xSmall,  ySmall, xErrorsSmall, yErrorsSmall);
    allPeaks = new TGraphErrors(npeaks, nPeak, xPeak, nErrorsPeak, xErrorsPeak);

}

void computeRotConstant(TCanvas *can){

    if(!allPeaks){

        cout << "All peaks graph not defined" << endl;
        return;

    }

    TF1 *linFit1 = new TF1("linFit1", "[0]+[1]*x");
    TF1 *linFit2 = new TF1("linFit2", "[0]+[1]*x");
    linFit1->SetParameter(0, -4);
    linFit1->SetParameter(1, 8);
    linFit2->SetParameter(0, 4);
    linFit2->SetParameter(1, 8);
    linFit1->SetLineColor(kBlue);
    linFit2->SetLineColor(kGreen+2);
    linFit2->SetLineStyle(7);

    linFit1->SetParName(0, "Offset");
    linFit1->SetParName(1, "Slope");
    linFit2->SetParName(0, "Offset");
    linFit2->SetParName(1, "Slope");

    // allPeaks->Fit(linFit1, "ME", "", -21.5, -1.5);
    allPeaks->Fit(linFit2, "ME+", "", 1.5, 19.5);

    graph(allPeaks, "Nitrogen", can, "Angular momentum #font[12]{l}", "Wavenumber shift #Delta k (cm^{-1})", "PeaksAll");

    Double_t avgRotConst = (linFit1->GetParameter(1)+linFit2->GetParameter(1))/2./4.;
    Double_t avgRotError = (linFit1->GetParError(1)+linFit2->GetParError(1))/2./4.;

    cout << "The rotational constant from fit is " << avgRotConst << " +- " << avgRotError << " 1/cm" << endl;
    cout << "Stokes: " << linFit1->GetParameter(1)/4. << " +- " << linFit1->GetParError(1)/4. << endl;
    cout << "Anti-Stokes: " << linFit2->GetParameter(1)/4. << " +- " << linFit2->GetParError(1)/4. << endl;

    Double_t *waveNumbers = allPeaks->GetY(), *wnErrors = allPeaks->GetEX();

    Double_t sepFirstPeak = (waveNumbers[21]-waveNumbers[19])/2./10.;
    Double_t sepFirstPeakError = (wnErrors[21]+wnErrors[19])/10.;

    cout << "The rotational constant from central line is " << sepFirstPeak << " +- " << sepFirstPeakError << " 1/cm" << endl;

    Double_t h = 6.62607004e-34, mass = 14*1.660539e-27, c = 299792458, pi = 3.14159265359;

    Double_t radius = Sqrt(h/(16*pi*pi*mass*avgRotConst*100*c));
    Double_t radiusError = avgRotError*avgRotError*h/(16*pi*pi*mass*100*c)*0.5*0.5/(avgRotConst*avgRotConst*avgRotConst);
    radiusError = Sqrt(radiusError);
    cout << "Radius is " << radius*1e12 << " +- " << radiusError*1e12 << " pm" << endl;

    Double_t radiusTemp = Sqrt(h/(16*pi*pi*mass*linFit2->GetParameter(1)/4.*100*c));
    Double_t radiusErrorTemp = linFit2->GetParError(1)/4.*linFit2->GetParError(1)/4.*h/(16*pi*pi*mass*100*c)*0.5*0.5/(linFit2->GetParameter(1)/4.*linFit2->GetParameter(1)/4.*linFit2->GetParameter(1)/4.);
    radiusErrorTemp = Sqrt(radiusErrorTemp);

    cout << "temp: " << radiusTemp*1e12 << " +- " << radiusErrorTemp*1e12 << endl;
}

void computeTemperature(TCanvas *can){

    if(!graphLarge || !graphSmall){

        cout << "Large or small graph not defined" << endl;
        return;

    }

    TF1 *fitLarge = new TF1("fitLarge", "[0]+[1]*(2*TMath::Abs(x)+1)*TMath::Exp(-[2]*TMath::Abs(x)*(TMath::Abs(x)+1))");
    fitLarge->SetParameter(0, 450);
    fitLarge->SetParameter(1, 250);
    fitLarge->SetParameter(2, 0.011);
    fitLarge->SetLineColor(kBlue);
    fitLarge->SetParName(0, "Offset");
    fitLarge->SetParName(1, "Normalization");
    fitLarge->SetParName(2, "Decay parameter");

    TF1 *fitSmall = new TF1("fitSmall", "[0]+[1]*(2*TMath::Abs(x)+1)*TMath::Exp(-[2]*TMath::Abs(x)*(TMath::Abs(x)+1))");
    fitSmall->SetParameter(0, 450);
    fitSmall->SetParameter(1, 130);
    fitSmall->SetParameter(2, 0.011);
    fitSmall->SetLineColor(kBlue);
    fitSmall->SetParName(0, "Offset");
    fitSmall->SetParName(1, "Normalization");
    fitSmall->SetParName(2, "Decay parameter");

    // graphSmall->Fit(fitSmall, "ME", "", -21.5, -0.5);
    // graphLarge->Fit(fitLarge, "ME", "", -21.5, -0.5);
    graphSmall->Fit(fitSmall, "ME", "", 0.5, 21.5);
    graphLarge->Fit(fitLarge, "ME", "", 0.5, 21.5);

    // graph(graphSmall, "Small peaks", can, "Wavenumber shift #Delta k (cm^{-1})", "Photon count", "PeaksSmall");
    // graph(graphLarge, "Large peaks", can, "Wavenumber shift #Delta k (cm^{-1})", "Photon count", "PeaksLarge");
    graph(graphSmall, graphLarge, "Nitrogen", can, "Wavenumber shift #Delta k (cm^{-1})", "Photon count", "Temperature");

    Double_t c = 299792458, KB = 1.38064852e-23, h = 6.62607004e-34, rotConst = 201.825, sigmaRot = 1.70133, sigmaB = fitSmall->GetParError(2);

    Double_t temperature1 = h*c*rotConst/(KB*fitSmall->GetParameter(2));
    Double_t temperature2 = h*c*rotConst/(KB*fitLarge->GetParameter(2));

    Double_t B1 = fitSmall->GetParameter(2);
    Double_t B2 = fitLarge->GetParameter(2);

    Double_t temp1Error = sigmaB*sigmaB*h*h*c*c*rotConst*rotConst/(KB*KB*B1*B1*B1*B1)+sigmaRot*sigmaRot*h*h*c*c/(KB*KB*B1*B1);
    Double_t temp2Error = sigmaB*sigmaB*h*h*c*c*rotConst*rotConst/(KB*KB*B2*B2*B2*B2)+sigmaRot*sigmaRot*h*h*c*c/(KB*KB*B2*B2);
    temp1Error = Sqrt(temp1Error);
    temp2Error = Sqrt(temp2Error);

    cout << "Temperature from small peaks is " << temperature1 << " +- " << temp1Error << " K and from the large peaks is " << temperature2 << " +- " << temp2Error << " K" << endl;

    Double_t nuclearSpin = fitLarge->GetParameter(1)/fitSmall->GetParameter(1);
    Double_t spinError = (fitLarge->GetParError(1))*(fitLarge->GetParError(1))/((fitSmall->GetParameter(1))*(fitSmall->GetParameter(1)))+(fitSmall->GetParError(1))*(fitSmall->GetParError(1))*(fitLarge->GetParameter(1))*(fitLarge->GetParameter(1))/((fitSmall->GetParameter(1))*(fitSmall->GetParameter(1))*(fitSmall->GetParameter(1))*(fitSmall->GetParameter(1)));
    spinError = Sqrt(spinError);
    cout << "The nuclear ratio is " << nuclearSpin << " +- " << spinError << endl;
    cout << "The nuclear spin is " << 1/(nuclearSpin-1.) << " +- " << spinError/((nuclearSpin-1.)*(nuclearSpin-1.)) << endl;

    Double_t nuclearSpinFixed = (2958.76-450)/(1731.87-450);
    Double_t spinErrorFixed = (2958.76-450)/(1731.87-450)/(1731.87-450)+(1731.87-450)*(2958.76-450)*(2958.76-450)/(1731.87-450)/(1731.87-450)/(1731.87-450)/(1731.87-450);
    spinErrorFixed = Sqrt(spinErrorFixed);
    cout << "The fixed nuclear ratio is " << nuclearSpinFixed << " +- " << spinErrorFixed << endl;
    cout << "The fixed nuclear spin is " << 1/(nuclearSpinFixed-1.) << " +- " << spinErrorFixed/((nuclearSpinFixed-1.)*(nuclearSpinFixed-1.)) << endl;

}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );
    // gStyle->SetOptFit(kFALSE);
    gStyle->SetOptFit(1100);
    gStyle->SetOptStat(kFALSE);
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("4.3f");

    if(!histoData) return;

    // TF1 *f = histoData->GetFunction("Peak_0");

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.6);
        histoData->SetMarkerColor(kRed);
        // histoData->SetMarkerColor(kBlue);
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
    // TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    // f->SetParameter(0, 250);
    // f->SetParameter(1, 100);
    // f->SetParameter(2, 0.395);
    // f->SetParameter(3, 1);
    // f->SetParameter(4, 10000);
    // f->SetParameter(5, 100);
    // f->SetParName(0, "Offset");
    // f->SetParName(1, "Peak intensity");
    // f->SetParName(2, "Aperture");
    // f->SetParName(3, "Center value");
    // f->SetParName(4, "Amp");
    // f->SetParLimits(0, 100, 300);
    // f->SetParLimits(1, 10, 100);
    // f->SetParLimits(2, 0.1, 1);
    // f->SetParLimits(3, -0.2, 2);
    // f->SetParLimits(4, 50, 10000);
    // f->SetParLimits(5, 1, 10000);
    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    // TF1 *f1 = graph->GetFunction("linFit1");
    // TF1 *f2 = graph->GetFunction("linFit2");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    // leg->AddEntry(f1, "Stokes linear fit","l");
    // leg->AddEntry(f2, "Anti-Stokes linear fit","l");
    
    graph->Draw("ap");
    // leg->Draw("same");

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

void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph1 || !graph2){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph1->SetLineWidth(2);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(1.5);
    graph1->SetMarkerColor(kRed);
    graph1->SetLineColor(kBlack);
    graph2->SetLineWidth(2);
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(1.5);
    graph2->SetMarkerColor(kRed);
    graph2->SetLineColor(kBlack);

    TF1 *f1 = graph1->GetFunction("fitSmall");
    TF1 *f2 = graph2->GetFunction("fitLarge");

    f1->SetLineColor(kBlue);
    f1->SetLineWidth(4);
    f2->SetLineColor(kGreen+1);
    f2->SetLineWidth(4);
    f2->SetLineStyle(9);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph1, "Data","lep");
    leg->AddEntry(f1, "Small Stokes fit","l");
    leg->AddEntry(f2, "Large Stokes fit","l");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Nitrogen");

    mg->Add(graph1);
    mg->Add(graph2);
    
    // graph1->Draw("ap");
    // graph2->Draw("same ap");
    mg->Draw("ap");
    leg->Draw("same");

    // add axis labels
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitle(yTitle);
    mg->GetYaxis()->CenterTitle();
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(0.9);
    mg->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    mg->SetTitle(graphName);
    // can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}
