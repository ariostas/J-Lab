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
#include <TF2.h>
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
#include "TTimer.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString, const TString);
void graph(TGraph*, const TString, TCanvas*, const TString, const TString, const TString);
void plot(TF1*, const TString, TCanvas*, const TString, const TString, const TString, const TString type = "");
TH1D* readFile(TString, TF1* cal = 0);
void fitHist(TH1D*, TString, Double_t, Double_t, Double_t, Double_t);
TF1* calibration();
Double_t calibErrorSquared(TF1*, TF1*, Double_t);

/*
 * MAIN FUNCTION
 */

void xrays(){

    TH1::StatOverflows(kTRUE);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(10);

    vector<TString> spectrum;
    //spectrum.push_back("nothing");
    spectrum.push_back("copper");
    spectrum.push_back("rubidium");
    spectrum.push_back("molybdenum");
    spectrum.push_back("silver");
    spectrum.push_back("barium");
    spectrum.push_back("terbium");

    vector<TString> attenuation;
    vector<Int_t> nData;
    attenuation.push_back("nothing"); nData.push_back(1);
    attenuation.push_back("al0_062"); nData.push_back(4);
    attenuation.push_back("cu"); nData.push_back(6);
    attenuation.push_back("alloy"); nData.push_back(12);

    vector<TString> brem;
    brem.push_back("sr_23feb_al_100_0");
    brem.push_back("sr_23feb_al_200_0");
    brem.push_back("sr_23feb_lead_20_0");
    brem.push_back("sr_23feb_lead_50_0");
    brem.push_back("sr_23feb_lead_100_0");
    brem.push_back("sr_23feb_nothing_20_0");
    brem.push_back("sr_23feb_nothing_50_0");
    brem.push_back("sr_23feb_thinlead_200_0");

    vector<TString> random;
    // random.push_back("barium_4feb");
    // random.push_back("lead_4feb");
    // random.push_back("nothing_4feb");
    // random.push_back("sodium_4feb");

    // random.push_back("copper_9feb_final");
    // random.push_back("rubidium_9feb_final");
    // random.push_back("molybdenum_9feb_final");
    // random.push_back("silver_9feb_final");
    // random.push_back("barium_9feb_final");
    // random.push_back("terbium_9feb_final");

    // random.push_back("silver_11feb_al0_062_1");
    // random.push_back("silver_11feb_al0_062_4");
    // random.push_back("silver_11feb_alloy_1");
    // random.push_back("silver_11feb_alloy_11");
    // random.push_back("silver_11feb_cu_1");
    // random.push_back("silver_11feb_cu_6");

    // random.push_back("terbium_18feb_al_1");
    // random.push_back("terbium_18feb_al_4");
    // random.push_back("terbium_18feb_alloy_1");
    // random.push_back("terbium_18feb_alloy_11");
    // random.push_back("terbium_18feb_cu_1");
    // random.push_back("terbium_18feb_cu_6");

    // random.push_back("sodium_23feb_calib");
    // random.push_back("sr_23feb_al_200_0");
    // random.push_back("sr_23feb_lead_20_0");
    // random.push_back("sr_23feb_lead_100_0");

    // random.push_back("sodium_25feb");
    // random.push_back("cs_25feb");

    // random.push_back("sodium_feb26_c_20_f_10");
    // random.push_back("ba_feb26_c_100_f_5");
    // random.push_back("cs_26feb_c_20_f_10");

    random.push_back("bremsstrahlung_feb29");

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);

    TF1 *calib = calibration();
    // cout << calib->Eval(670) << endl;
    // cout << calib->Eval(340) << endl;
    // cout << calib->GetParError(0) << "  " << calib->GetParError(1) << endl;
    TF1 *calibError = new TF1("calibError", "[0]+[1]*x");
    calibError->FixParameter(0, calib->GetParameter(0)+calib->GetParError(0));
    calibError->FixParameter(1, calib->GetParameter(1)+calib->GetParError(1));
    // cout << calib->Eval(500) - calibError->Eval(500) << endl;

    // for(Int_t x=0; x<spectrum.size(); x++){
    //     TH1D *tempHist = readFile(spectrum.at(x)+"_9feb_final");
    //     // fitHist(tempHist, "1", 600, 750, 1e4, 1.);
    //     fitHist(tempHist, "2", 790, 830, 1e4, 0.2);
    //     fitHist(tempHist, "3", 90, 140, 1e3, 0.2);
    //     histogram(tempHist, "Calibration", can, "MCA bin", "Counts", TString::Format("Spectrum_%i_%s", x, spectrum.at(x).Data()), "dot log");
    // }

    // for(Int_t x=0; x<brem.size(); x++){
    //     TH1D *tempHist = readFile(brem.at(x));
    //     histogram(tempHist, "Data", can, "Decay energy", "Counts", TString::Format("Bremsstrahlung_%i_%s", x, brem.at(x).Data()), "line log");
    // }

    // for(Int_t x=0; x<attenuation.size(); x++){
    //     Double_t n[nData.at(x)], nError[nData.at(x)], integ[nData.at(x)], integError[nData.at(x)];
    //     for(Int_t y=1; y<=nData.at(x); y++){
    //         TH1D *tempHist = readFile("silver_11feb_"+attenuation.at(x)+TString::Format("_%i",y), calib);
    //         histogram(tempHist, "Data", can, "Decay energy", "Counts", "Attenuation_"+attenuation.at(x)+TString::Format("_%i",y), "dot log");
    //         n[y-1] = y;
    //         nError[y-1] = 0;
    //         integ[y-1] = tempHist->Integral(850, 1150); 
    //         integError[y-1] = Sqrt(tempHist->Integral(850, 1150));
    //     }
    //     TGraphErrors *tempgraph = new TGraphErrors(nData.at(x), n, integ, nError, integError);
    //     graph(tempgraph, "Attenuation "+attenuation.at(x), can, "Number of sheets", "Counts", "Plot_"+attenuation.at(x));
    // }

    // for(Int_t x=0; x<random.size(); x++){
    //     TH1D *tempHist = readFile(random.at(x));
    //     histogram(tempHist, "Bremsstrahlung x-rays", can, "MCA bin", "Counts", TString::Format("Random_%i_%s", x, random.at(x).Data()), "dot log");
    // }

    vector<TH1D*> histos;
    histos.push_back(readFile(spectrum.at(0)+"_9feb_final", calib));
    // fitHist(histos.at(0), "1", 58.5, 60.5, 1e4, 0.2);
    // fitHist(histos.at(0), "3", 7, 11, 1e3, 0.2);
    fitHist(histos.at(0), "4", 7, 11, 1e3, 0.2);

    histos.push_back(readFile(spectrum.at(1)+"_9feb_final", calib));
    fitHist(histos.at(1), "1", 3, 5, 1e3, 0.1);
    fitHist(histos.at(1), "2", 5, 7, 1e3, 0.1);
    fitHist(histos.at(1), "3", 13, 15, 1e4, 0.2);
    fitHist(histos.at(1), "4", 15.1, 16.5, 1e4, 0.3);

    histos.push_back(readFile(spectrum.at(2)+"_9feb_final", calib));
    fitHist(histos.at(2), "1", 7, 9, 1e4, 0.2);
    fitHist(histos.at(2), "2", 9, 12, 1e3, 0.2);
    fitHist(histos.at(2), "3", 17, 19.1, 1e4, 0.2);
    fitHist(histos.at(2), "4", 19.3, 21.5, 1e4, 0.2);

    cout << "silver" << endl;
    histos.push_back(readFile(spectrum.at(3)+"_9feb_final", calib));
    fitHist(histos.at(3), "1", 12, 14, 1e4, 0.2);
    fitHist(histos.at(3), "2", 15, 17, 1e4, 0.2);
    fitHist(histos.at(3), "3", 21, 24.1, 1e4, 0.2);
    fitHist(histos.at(3), "4", 24.4, 27, 1e4, 0.2);

    cout << "barium" << endl;
    histos.push_back(readFile(spectrum.at(4)+"_9feb_final", calib));
    fitHist(histos.at(4), "1", 21, 24.5, 1e4, 0.2);
    fitHist(histos.at(4), "2", 24.5, 28, 1e3, 0.2);
    fitHist(histos.at(4), "3", 31, 34.3, 1e5, 0.2);
    fitHist(histos.at(4), "4", 35.4, 39, 1e4, 0.2);

    cout << "terbum" << endl;
    histos.push_back(readFile(spectrum.at(5)+"_9feb_final", calib));
    fitHist(histos.at(5), "1", 32, 36, 1e4, 0.2);
    fitHist(histos.at(5), "3", 43, 46.1, 1e5, 0.1);
    fitHist(histos.at(5), "4", 48, 51.3, 1e4, 0.2);

    for(UInt_t x = 0; x < histos.size(); x++){
        histogram(histos.at(x), spectrum.at(x), can, "Energy (keV)", "Counts", TString::Format("Spectrum_%i_%s", x, spectrum.at(x).Data()), "dot log");

    }

    can->Clear();

    TString moseleyNames[4] = {"Reflection 1", "Reflection 2", "K#alpha lines", "K#beta lines"};

    Double_t zNumbers[6] = {29., 37., 42., 47., 56., 65.};
    Double_t zErrors[6] = {0, 0, 0, 0, 0, 0};

    for(Int_t x = 1; x <= 4; x++){

        Double_t sqrtEnergies[6], zNum[6], eErrors[6];
        Int_t el = 0;

        for(Int_t i = 0; i < histos.size(); i++){
            TF1 *fit = histos.at(i)->GetFunction(TString::Format("%i", x));
            if(fit == NULL) continue;

            sqrtEnergies[el] = Sqrt(fit->GetParameter(2));
            zNum[el] = zNumbers[i];
            // eErrors[el] = Sqrt(0.0569749*0.0569749+fit->GetParError(2)*fit->GetParError(2)+fit->GetParameter(3)*fit->GetParameter(3))/Sqrt(2.*fit->GetParameter(2));
            // if(i==0) eErrors[el] = Sqrt(fit->GetParError(2)*fit->GetParError(2)+fit->GetParameter(3)*fit->GetParameter(3))/Sqrt(2.*fit->GetParameter(2));
            // eErrors[el] = (0.0569749+fit->GetParError(2)+fit->GetParameter(3))/Sqrt(2.*fit->GetParameter(2));
            // if(i==0) eErrors[el] = (fit->GetParError(2)+fit->GetParameter(3))/Sqrt(2.*fit->GetParameter(2));
            eErrors[el] = Sqrt(calibErrorSquared(calib, calibError, fit->GetParameter(2))+fit->GetParError(2)*fit->GetParError(2))/Sqrt(2.*fit->GetParameter(2));
            // if(i==0) eErrors[el] = Sqrt(fit->GetParError(2)*fit->GetParError(2))/Sqrt(2.*fit->GetParameter(2));


            cout << "Z: " << zNum[el] << ", line: " << x << ". E = " << fit->GetParameter(2) << " +- " << Sqrt(0.0569749*0.0569749+fit->GetParError(2)*fit->GetParError(2)+fit->GetParameter(3)*fit->GetParameter(3)) << endl;
            // cout << "Z: " << zNum[el] << ", line: " << x << ". E = " << fit->GetParameter(2) << " +- " << 0.0569749+fit->GetParError(2)+fit->GetParameter(3) << endl;


            el++;

        }

        TGraphErrors *moseley = new TGraphErrors(el, zNum, sqrtEnergies, zErrors, eErrors);
        graph(moseley, moseleyNames[x-1], can, "Atomic number Z", "#sqrt{E} (#sqrt{keV})", TString::Format("Moseley_%i",x));

    }

}

TH1D* readFile(TString inputfile, TF1 *cal){

    ifstream ifs("./Data/" + inputfile + ".txt"); if(!ifs.is_open()){cout << "Error. File " << inputfile << " not found. Exiting...\n"; return NULL;}

    TH1D *p;
    if(cal == 0) p = new TH1D(inputfile, inputfile, 2048, -0.5, 2047.5);
    else p = new TH1D(inputfile, inputfile, 1024, cal->Eval(-0.5), cal->Eval(1023.5));

    Int_t bin = 0, count = 0;
    
    while(ifs >> count){
        
        if(bin < 30){bin++; continue;}
        for(Int_t x = 0; x < count; x++){
            if(cal == 0) p->Fill(bin);
            else p->Fill(cal->Eval(bin));
        }
        bin++;
                
    }

    ifs.close();

    //TF1 *f = new TF1(inputfile, "[0]*TMath::Landau(x,[1],[2])+[3]");
    TF1 *f = new TF1(inputfile, "[0]*TMath::Gaus(x,[1],[2])+[3]");
    f->SetParameter(0, p->GetMaximum());
    f->SetParameter(1, p->GetMean());
    f->SetParameter(2, p->GetStdDev());
    f->SetParameter(3, p->GetMinimum());
    f->SetParLimits(0, 4, 30);
    f->SetParLimits(1, 300, 600);
    f->SetParLimits(2, 20, 50);
    f->SetParLimits(3, -1, 1);
    //p->Fit(f, "ME");

    return p;

}

void fitHist(TH1D *histo, TString name, Double_t min, Double_t max, Double_t height, Double_t width){

    TF1 *func = new TF1(name, "[0]+[1]*TMath::Voigt(x-[2],[3],[4])", min, max);
    func->SetParameter(0, height/1000.);
    func->SetParameter(1, height);
    func->SetParameter(2, (max + min)/2.);
    func->SetParameter(3, width);
    func->SetParameter(4, width);
    func->SetParName(0, "Offset");
    func->SetParName(1, "Normalization");
    func->SetParName(2, "Mean");
    func->SetParName(3, "Sigma");
    func->SetParName(4, "Gamma");
    func->SetLineColor(kBlue);
    histo->Fit(func, "ME+", "", min, max);
}

TF1* calibration(){

    Double_t bin[2], realE[2], binError[2], realEError[2];

    realE[0] = 59.541; bin[0] = 809.611;
    //realE[1] = 8.904; bin[1] = 111.332;
    realE[1] = 8.048; bin[1] = 111.332;
    realEError[0] = 0.001; binError[0] = 1;
    realEError[1] = 0.001; binError[1] = 1;

    TGraphErrors *cal = new TGraphErrors(2, bin, realE, binError, realEError);
    TF1 *f = new TF1("fit", "[0]+[1]*x");
    cal->Fit(f, "ME");

    return f;
}

Double_t calibErrorSquared(TF1 *calib, TF1 *calibError, Double_t energy){

    Double_t bin = calib->GetX(energy);
    Double_t error = fabs(calib->Eval(energy)-calibError->Eval(energy));
    cout << error << endl;
    return error*error;

}

void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    // gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);
    gPad->SetLogy(0);
    gStyle->SetStatFormat("6.2g");

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    TF1 *f = new TF1(graphName, "[1]*x-[0]");
    f->SetParameter(0, 0);
    f->SetParameter(1, 1);
    f->SetParName(0, "C#sigma");
    f->SetParName(1, "C");
    f->SetLineColor(kBlue);
    f->SetLineWidth(4);
    graph->Fit(f, "MEF");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","lep");
    leg->AddEntry(f, "Linear fit","l");
    
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

void plot(TF1 *function, const TString fName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    if(!function){
        cout << "Error: Graph \"" << fName << "\" not defined" << endl;
        return;
    }

    function->SetLineWidth(2);
    function->SetLineColor(kBlack);
    
    //can->Clear();
    function->Draw();

    // add axis labels
    function->GetXaxis()->SetTitle(xTitle);
    function->GetXaxis()->CenterTitle();
    function->GetXaxis()->SetTitleSize(0.055);
    function->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    function->GetXaxis()->SetLabelSize(0.05);
    function->GetYaxis()->SetTitle(yTitle);
    function->GetYaxis()->CenterTitle();
    function->GetYaxis()->SetTitleSize(0.055);
    function->GetYaxis()->SetTitleOffset(0.8);
    function->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    function->SetTitle(fName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histoData, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );
    gStyle->SetOptFit(kFALSE);
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("4.3G");

    if(!histoData) return;

    TF1 *f = histoData->GetFunction("4");

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
    leg->AddEntry(f, "Voigtian fits","l");
    
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

void histogram(TH1D *histoData, TH1D *histoMC, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name, const TString type){

    gPad->SetLogy( (type.Contains("log") ? 1 : 0) );

    if(!histoData || !histoMC) return;

    if(type.Contains("dot")){
        histoData->SetLineWidth(1);
        histoData->SetMarkerStyle(20);
        histoData->SetMarkerSize(0.9);
        histoData->SetMarkerColor(kRed);
        histoData->SetLineColor(kBlack);
        histoMC->SetLineWidth(1);
        histoMC->SetMarkerStyle(21);
        histoMC->SetMarkerSize(0.9);
        histoMC->SetMarkerColor(kBlue);
        histoMC->SetLineColor(kBlack);
    }

    else if(type.Contains("line")){
        histoData->SetLineWidth(1);
        histoData->SetLineColor(kBlue);
        histoMC->SetLineWidth(1);
        histoMC->SetLineColor(kRed);
    }
    else return;

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.355,0.675,0.585,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(histoData, "Data", (type == "dot" ? "lep" : "l"));
    leg->AddEntry(histoData, "Monte Carlo", (type == "dot" ? "lep" : "l"));

    if(histoData->GetMaximum() > histoMC->GetMaximum()){
        histoData->Draw((type.Contains("dot") ? "E1" : ""));
        histoMC->Draw((type.Contains("dot") ? "same E1" : "same"));

        // add axis labels
        histoData->GetXaxis()->SetTitle(xTitle);
        histoData->GetXaxis()->CenterTitle();
        histoData->GetXaxis()->SetTitleSize(0.055);
        histoData->GetXaxis()->SetTitleOffset(0.84);
        histoData->GetXaxis()->SetLabelOffset(0.010);
        histoData->GetXaxis()->SetLabelSize(0.05);
        histoData->GetYaxis()->SetTitle(yTitle);
        histoData->GetYaxis()->CenterTitle();
        histoData->GetYaxis()->SetTitleSize(0.055);
        histoData->GetYaxis()->SetTitleOffset(1.0);
        histoData->GetYaxis()->SetLabelSize(0.05);
        gStyle->SetTitleSize(0.08, "t");
        histoData->SetTitle("");
    }
    else{
        histoMC->Draw((type == "dot" ? "E1" : ""));
        histoData->Draw((type == "dot" ? "same E1" : "same"));

        // add axis labels
        histoMC->GetXaxis()->SetTitle(xTitle);
        histoMC->GetXaxis()->CenterTitle();
        histoMC->GetXaxis()->SetTitleSize(0.055);
        histoMC->GetXaxis()->SetTitleOffset(0.84);
        histoMC->GetXaxis()->SetLabelOffset(0.010);
        histoMC->GetXaxis()->SetLabelSize(0.05);
        histoMC->GetYaxis()->SetTitle(yTitle);
        histoMC->GetYaxis()->CenterTitle();
        histoMC->GetYaxis()->SetTitleSize(0.055);
        histoMC->GetYaxis()->SetTitleOffset(1.0);
        histoMC->GetYaxis()->SetLabelSize(0.05);
        gStyle->SetTitleSize(0.08, "t");
        histoMC->SetTitle("");
    }
    
    leg->Draw("same");

    can->Update();

    can->SaveAs(name + ".png");
}
