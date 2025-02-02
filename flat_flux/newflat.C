#include <TRandom.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TAxis.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <Riostream.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <TBranch.h>
#include <TBasket.h>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <iterator>
#include <limits>

#define PI 3.14159265

using namespace std;


double FindMaximum(TGraph* graph) {
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }


    int nPoints = graph->GetN();
    double x, y;
    double maxY = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x, y);
        if (y > maxY) {
            maxY = y;
        }
    }

    return maxY;
}

double FindMinimum(TGraph* graph) {
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }


    int nPoints = graph->GetN();
    double x, y;
    double minY = -999.;

    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x, y);
        if (y < minY) {
            minY = y;
        }
        if (minY <= 0.){
            minY = 1.e-20;
        }
    }

    return minY;
}

double FindMeanWeight(TGraph* graph, double alpha, double maxVal, double minVal){
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    int nPoints = graph->GetN();
    double x, y;
    //double meanWeight = -std::numeric_limits<double>::infinity();
    double meanWeight = 0.;

    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x, y);
        double weight = 1./(y+alpha*maxVal+minVal);
        meanWeight += weight/nPoints;
        }
        return meanWeight;
}

void MakeReweightedGraph(TGraph* graph, double alpha, double maxVal, double minVal){
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return;
    }
    int nPoints = graph->GetN();
    double x, y;
    //double meanWeight = -std::numeric_limits<double>::infinity();
    double meanWeight = 0.;

    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x, y);
        double weight = 1./(y+alpha*maxVal+minVal);
        graph->SetPoint(i, x, y * weight); 

        if (i%100==0){std::cout << setprecision(10) <<"weight: " << weight << "    numubar_spectrum: " << y
                                     << "    ALPHA: " << alpha << "    maxVal: " << maxVal << std::endl;}
        }
}

void MakeReweightedGraphwClip(TGraph* graph, double alpha, double maxVal, double minVal){
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return;
    }
    int nPoints = graph->GetN();
    double x, y;
    //double meanWeight = -std::numeric_limits<double>::infinity();
    double meanWeight = 0.;

    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x, y);
        double weight = 1./(y+alpha*maxVal+minVal);
        if(weight < 0.2){
            weight = 0.2;
        }
        if(weight > 5.){
            weight = 5.;
        }
        // if (weight > 0.2 && weight < 5){
        //     weight = 1./(maxVal+alpha*maxVal+minVal);
        // }
        graph->SetPoint(i, x, y * weight); 

        if (i%100==0){std::cout << setprecision(10) <<"weight: " << weight << "    numubar_spectrum: " << y
                                     << "    ALPHA: " << alpha << "    maxVal: " << maxVal << std::endl;}
        }
}

void NormalizeTGraph(TGraph* graph) {
    if (!graph) {
        std::cerr << "TGraph object is null!" << std::endl;
        return;
    }

    int nPoints = graph->GetN();
    if (nPoints < 2) {
        std::cerr << "Not enough points in TGraph to normalize!" << std::endl;
        return;
    }

    // Calculate the integral using the trapezoidal rule
    double x1, y1, x2, y2;
    double integral = 0.0;

    for (int i = 0; i < nPoints - 1; ++i) {
        graph->GetPoint(i, x1, y1);
        graph->GetPoint(i + 1, x2, y2);
        double dx = x2 - x1;
        double avgY = (y1 + y2) / 2.0;
        integral += dx * avgY;
    }

    if (integral == 0.0) {
        std::cerr << "Integral of TGraph is zero, cannot normalize!" << std::endl;
        return;
    }

    // Normalize the y-values
    for (int i = 0; i < nPoints; ++i) {
        graph->GetPoint(i, x1, y1);
        graph->SetPoint(i, x1, y1 / integral); //SetPoint(i, x1, y1 * weight);
    }

    std::cout << "TGraph successfully normalized to have an integral of unity." << std::endl;
}



void newflat()
{
    //Create normalization factor for proper units in all spectra
    const double factor[] = {1.00};
    //Create variables for different cross section printouts
    double CCxsec_nuebar, NCxsec_nuebar, CCxsec_nue, NCxsec_nue, CCxsec_numubar, NCxsec_numubar, CCxsec_numu, NCxsec_numu, CCxsec_nutaubar, NCxsec_nutaubar, CCxsec_nutau, NCxsec_nutau;
    //Create strings of names for the file folders and text file outputs
    string elementnames[] = {"Ar40"};
    double elementpercentages[] = {1.0}/*{0.2,0.6,0.2}*/;
    //Read in from the cross sections from the gspl2root converted file
    TFile *file = new TFile("xsec_graphs.root");
    TFile *OUTFILE = new TFile("checks.root","RECREATE");
    //Create a large canvas for drawing things
    //TCanvas *c1 = new TCanvas("c1","c1",3000,2000);
    //Make output streams for the inverse cross sections for all neutrino flavors
    //In principle we can break this down by totals, CC, and NC--if we want to be precise,
    //we need to do this process by process, technically...
    //These output streams will output high precision flux values of the inverse cross sections
    //in a two column format that can be made sense by GENIE gevgen app
    ofstream numuCC_inverse_xsec_file("numu_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream numuNC_inverse_xsec_file("numu_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream numu_inverse_xsec_file("numu_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ofstream numubarCC_inverse_xsec_file("numubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream numubarNC_inverse_xsec_file("numubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream numubar_inverse_xsec_file("numubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ofstream nueCC_inverse_xsec_file("nue_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nueNC_inverse_xsec_file("nue_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nue_inverse_xsec_file("nue_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ofstream nuebarCC_inverse_xsec_file("nuebar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nuebarNC_inverse_xsec_file("nuebar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nuebar_inverse_xsec_file("nuebar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ofstream nutauCC_inverse_xsec_file("nutau_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nutauNC_inverse_xsec_file("nutau_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nutau_inverse_xsec_file("nutau_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ofstream nutaubarCC_inverse_xsec_file("nutaubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nutaubarNC_inverse_xsec_file("nutaubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ofstream nutaubar_inverse_xsec_file("nutaubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC 
    //Create a bin width for the energy values
    const double bin_width = 0.001;//1MeV in GENIE units (which are in GeV)
    //Loop through all of the various elements

    //Grab the TGraphs from the gspl2root converted cross section file
    // std::vector<TGraph*> graph
    // std::vector<TString*> names = ["",];
    //TGraph* graph[] = (TGraph*)file->Get(std:vector<names>);
    TGraph* tot_cc_nuebar = (TGraph*)file->Get("nu_e_bar_Ar40/tot_cc");
    TGraph* tot_nc_nuebar = (TGraph*)file->Get("nu_e_bar_Ar40/tot_nc");
    TGraph* tot_cc_nue = (TGraph*)file->Get("nu_e_Ar40/tot_cc");
    TGraph* tot_nc_nue = (TGraph*)file->Get("nu_e_Ar40/tot_nc");
    TGraph* tot_cc_numubar = (TGraph*)file->Get("nu_mu_bar_Ar40/tot_cc");
    TGraph* tot_nc_numubar = (TGraph*)file->Get("nu_mu_bar_Ar40/tot_nc");
    TGraph* tot_cc_numu = (TGraph*)file->Get("nu_mu_Ar40/tot_cc");
    TGraph* tot_nc_numu = (TGraph*)file->Get("nu_mu_Ar40/tot_nc");
    TGraph* tot_cc_nutaubar = (TGraph*)file->Get("nu_tau_bar_Ar40/tot_cc");
    TGraph* tot_nc_nutaubar = (TGraph*)file->Get("nu_tau_bar_Ar40/tot_nc");
    TGraph* tot_cc_nutau = (TGraph*)file->Get("nu_tau_Ar40/tot_cc");
    TGraph* tot_nc_nutau = (TGraph*)file->Get("nu_tau_Ar40/tot_nc");
    //Run through various energies (m) on which they will be interpolated using ROOT's TSpline fits to TGraphs
    int i;
    for (double m = 0.1; m < 10.000; m+=bin_width)
        {
            //i++;
            //Evaluate the cross sections at various energy values (via interpolation)
            CCxsec_nuebar = tot_cc_nuebar->Eval(m); if(CCxsec_nuebar<0.){CCxsec_nuebar=0.;}
            NCxsec_nuebar = tot_nc_nuebar->Eval(m); if(NCxsec_nuebar<0.){NCxsec_nuebar=0.;}
            CCxsec_nue = tot_cc_nue->Eval(m); if(CCxsec_nue<0.){CCxsec_nue=0.;}
            NCxsec_nue = tot_nc_nue->Eval(m); if(NCxsec_nue<0.){NCxsec_nue=0.;}
            CCxsec_numubar = tot_cc_numubar->Eval(m); if(CCxsec_numubar<0.){CCxsec_numubar=0.;}
            NCxsec_numubar = tot_nc_numubar->Eval(m); if(NCxsec_numubar<0.){NCxsec_numubar=0.;}
            CCxsec_numu = tot_cc_numu->Eval(m); if(CCxsec_numu<0.){CCxsec_numu=0.;}
            NCxsec_numu = tot_nc_numu->Eval(m); if(NCxsec_numu<0.){NCxsec_numu=0.;}
            CCxsec_nutaubar = tot_cc_nutaubar->Eval(m); if(CCxsec_nutaubar<0.){CCxsec_nutaubar=0.;}
            NCxsec_nutaubar = tot_nc_nutaubar->Eval(m); if(NCxsec_nutaubar<0.){NCxsec_nutaubar=0.;}
            CCxsec_nutau = tot_cc_nutau->Eval(m); if(CCxsec_nutau<0.){CCxsec_nutau=0.;}
            NCxsec_nutau = tot_nc_nutau->Eval(m); if(NCxsec_nutau<0.){NCxsec_nutau=0.;}
            //Print out those INVERSE cross sections to act as the fluxes for GENIE
            nuebarCC_inverse_xsec_file   << setprecision(8) << m << "   " << 1./CCxsec_nuebar                     << endl; 
            nuebarNC_inverse_xsec_file   << setprecision(8) << m << "   " << 1./NCxsec_nuebar                     << endl;
            nuebar_inverse_xsec_file     << setprecision(8) << m << "   " << 1./(CCxsec_nuebar+NCxsec_nuebar)     << endl;
            nueCC_inverse_xsec_file      << setprecision(8) << m << "   " << 1./(CCxsec_nue)                      << endl;
            nueNC_inverse_xsec_file      << setprecision(8) << m << "   " << 1./(NCxsec_nue)                      << endl;
            nue_inverse_xsec_file        << setprecision(8) << m << "   " << 1./(CCxsec_nue+NCxsec_nue)           << endl;
            numubarCC_inverse_xsec_file  << setprecision(8) << m << "   " << 1./CCxsec_numubar                    << endl;
            numubarNC_inverse_xsec_file  << setprecision(8) << m << "   " << 1./NCxsec_numubar                    << endl;
            numubar_inverse_xsec_file    << setprecision(8) << m << "   " << 1./(CCxsec_numubar+NCxsec_numubar)   << endl;
            numuCC_inverse_xsec_file     << setprecision(8) << m << "   " << 1./CCxsec_numu                       << endl;
            numuNC_inverse_xsec_file     << setprecision(8) << m << "   " << 1./NCxsec_numu                       << endl;
            numu_inverse_xsec_file       << setprecision(8) << m << "   " << 1./(CCxsec_numu+NCxsec_numu)         << endl;
            nutaubarCC_inverse_xsec_file << setprecision(8) << m << "   " << 1./CCxsec_nutaubar                   << endl;
            nutaubarNC_inverse_xsec_file << setprecision(8) << m << "   " << 1./NCxsec_nutaubar                   << endl;
            nutaubar_inverse_xsec_file   << setprecision(8) << m << "   " << 1./(CCxsec_nutaubar+NCxsec_nutaubar) << endl;
            nutauCC_inverse_xsec_file    << setprecision(8) << m << "   " << 1./CCxsec_nutau                      << endl;
            nutauNC_inverse_xsec_file    << setprecision(8) << m << "   " << 1./NCxsec_nutau                      << endl;
            nutau_inverse_xsec_file      << setprecision(8) << m << "   " << 1./(CCxsec_nutau+NCxsec_nutau)       << endl;                   
        
            //Energy[i-1] = m;
        }
        
    nuebarCC_inverse_xsec_file.close();
    nuebarNC_inverse_xsec_file.close();
    nuebar_inverse_xsec_file.close();
    nueCC_inverse_xsec_file.close();
    nueNC_inverse_xsec_file.close();
    nue_inverse_xsec_file.close();
    numuCC_inverse_xsec_file.close();
    numuNC_inverse_xsec_file.close();
    numu_inverse_xsec_file.close();
    nutaubarCC_inverse_xsec_file.close();
    nutaubarNC_inverse_xsec_file.close();
    nutaubar_inverse_xsec_file.close();
    nutauCC_inverse_xsec_file.close();
    nutauNC_inverse_xsec_file.close();
    nutau_inverse_xsec_file.close();

    ifstream numuCC_inverse_xsec_input("numu_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream numuNC_inverse_xsec_input("numu_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream numu_inverse_xsec_input("numu_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ifstream numubarCC_inverse_xsec_input("numubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream numubarNC_inverse_xsec_input("numubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream numubar_inverse_xsec_input("numubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ifstream nueCC_inverse_xsec_input("nue_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nueNC_inverse_xsec_input("nue_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nue_inverse_xsec_input("nue_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ifstream nuebarCC_inverse_xsec_input("nuebar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nuebarNC_inverse_xsec_input("nuebar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nuebar_inverse_xsec_input("nuebar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ifstream nutauCC_inverse_xsec_input("nutau_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nutauNC_inverse_xsec_input("nutau_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nutau_inverse_xsec_input("nutau_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC
    ifstream nutaubarCC_inverse_xsec_input("nutaubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nutaubarNC_inverse_xsec_input("nutaubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data");
    ifstream nutaubar_inverse_xsec_input("nutaubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data");//Total cross section, CC+NC 

    const int inverse_count = 9901;
    double energy, inverse_xsec_flux, flux;
    double Energy[inverse_count], Inverse_XSec_Flux[inverse_count];
    for(int i = 0; i<inverse_count; i++)
        {
            //Make vectors for TGraph to be filled
            numuCC_inverse_xsec_input >> energy >> inverse_xsec_flux;
            Energy[i] = energy;
            Inverse_XSec_Flux[i] = inverse_xsec_flux;
        }
    double XSEC_NC, XSEC_CC, FLAT_SPECTRA[inverse_count];
    for (int i = 0; i < inverse_count; i++)
        {
            //XSEC_NC = tot_nc_numu->Eval(Energy[i]);
            XSEC_CC = tot_cc_numu->Eval(Energy[i]);
            FLAT_SPECTRA[i] = (XSEC_CC/*+XSEC_NC*/)*Inverse_XSec_Flux[i];
        }
    TGraph* flat_spectra_check = new TGraph(inverse_count,Energy,FLAT_SPECTRA);
    flat_spectra_check->SetName("inverse_check");
    flat_spectra_check->SetTitle("inverse_check");
    flat_spectra_check->SetLineColor(kBlack);
    flat_spectra_check->Write();

    //ifstream numubar_flux_file("DUNE_OnAxis.d");
    ifstream numubar_flux_file("LE_NuMI.d");
    

    //const int nnn_DOn = 74;
    const int nnn_DOn = 101; //LE_NuMI

    double numubar_energy[nnn_DOn], numubar_flux[nnn_DOn];
    for(int i = 0; i<nnn_DOn; i++)
        {
            //Make vectors for TGraph to be filled
            numubar_flux_file >> energy >> flux;
            numubar_energy[i] = energy;
            //cout << "energy: " << numubar_energy[i] << endl;
            numubar_flux[i] = flux;
        }
    TGraph* numubar_flux_DOn = new TGraph(nnn_DOn, numubar_energy, numubar_flux);
    numubar_flux_DOn->SetName("DOn_check");
    numubar_flux_DOn->SetTitle("DOn_check");
    numubar_flux_DOn->SetLineColor(kRed);
    numubar_flux_DOn->Write();

    double flux_numubar, numubar_spectrum, UNWEIGHTED_SPECTRA[4900], WEIGHTED_SPECTRA[4900], CCXSEC[4900], FLUX[4900];
    for (int i = 0; i < 4900; i++) 
        {
            //Evaluate the cross sections at various energy values (via interpolation)
            CCxsec_numubar = tot_cc_numubar->Eval(Energy[i]); if(CCxsec_numubar<0.){CCxsec_numubar=0.;}
            flux_numubar = numubar_flux_DOn->Eval(Energy[i]); if(flux_numubar<0.){flux_numubar=0.;}
            numubar_spectrum = flux_numubar*CCxsec_numubar*bin_width;
            UNWEIGHTED_SPECTRA[i] = numubar_spectrum;
            CCXSEC[i] = CCxsec_numubar;
            FLUX[i] = flux_numubar;
            //cout << "numubar_spectrum: " << numubar_spectrum << endl;
        }
    
    TGraph* numubar_un_flux_DOn = new TGraph(4900, Energy, FLUX);
    numubar_un_flux_DOn->SetName("numubar_un_flux_DOn");
    numubar_un_flux_DOn->SetTitle("numubar_un_flux_DOn");
    numubar_un_flux_DOn->SetLineColor(kGreen);
    numubar_un_flux_DOn->Write();

    TGraph* numubar_un_ccxsec_DOn = new TGraph(4900, Energy, CCXSEC);
    numubar_un_ccxsec_DOn->SetName("numubar_un_ccxsec_DOn");
    numubar_un_ccxsec_DOn->SetTitle("numubar_un_ccxsec_DOn");
    numubar_un_ccxsec_DOn->SetLineColor(kGreen);
    numubar_un_ccxsec_DOn->Write();

    TGraph* numubar_un_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    numubar_un_spectrum_DOn->SetName("DOn_un_spectra_check");
    numubar_un_spectrum_DOn->SetTitle("DOn_un_spectra_check");
    numubar_un_spectrum_DOn->SetLineColor(kGreen);
    numubar_un_spectrum_DOn->Write();

    TGraph* numubar_normalized_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    NormalizeTGraph(numubar_normalized_spectrum_DOn);
    numubar_normalized_spectrum_DOn->SetName("numubar_normalized_spectrum_DOn_check");
    numubar_normalized_spectrum_DOn->SetTitle("numubar_normalized_spectrum_DOn_check");
    numubar_normalized_spectrum_DOn->SetLineColor(kGreen);
    numubar_normalized_spectrum_DOn->Write();

    const double maxVal = FindMaximum(numubar_un_spectrum_DOn);
    std::cout << "Maximum unweigted y-value: " << maxVal << std::endl;

    const double minVal = FindMinimum(numubar_un_spectrum_DOn);
    std::cout << "Minimum unweigted y-value: " << minVal << std::endl;

    double ALPHA = 0.0;

    std::cout << "before numubar_unN_reweighted_spectrum_DOn"  << std::endl;

    TGraph* numubar_unN_reweighted_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    MakeReweightedGraph(numubar_unN_reweighted_spectrum_DOn, ALPHA, maxVal, minVal);
    numubar_unN_reweighted_spectrum_DOn->SetName("numubar_unN_reweighted_spectrum_DOn");
    numubar_unN_reweighted_spectrum_DOn->SetTitle("numubar_unN_reweighted_spectrum_DOn");
    numubar_unN_reweighted_spectrum_DOn->SetLineColor(kGreen);
    numubar_unN_reweighted_spectrum_DOn->Write();

    std::cout << "before numubar_N_reweighted_spectrum_DOn y-value: " << std::endl;

    TGraph* numubar_N_reweighted_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    NormalizeTGraph(numubar_N_reweighted_spectrum_DOn);
    MakeReweightedGraph(numubar_N_reweighted_spectrum_DOn, ALPHA, maxVal, minVal);
    numubar_N_reweighted_spectrum_DOn->SetName("numubar_N_reweighted_spectrum_DOn");
    numubar_N_reweighted_spectrum_DOn->SetTitle("numubar_N_reweighted_spectrum_DOn");
    numubar_N_reweighted_spectrum_DOn->SetLineColor(kGreen);
    numubar_N_reweighted_spectrum_DOn->Write();

    
    TGraph* numubar_N_reweighted_wClip_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    NormalizeTGraph(numubar_N_reweighted_wClip_spectrum_DOn);
    MakeReweightedGraphwClip(numubar_N_reweighted_wClip_spectrum_DOn, ALPHA, maxVal, minVal);
    numubar_N_reweighted_wClip_spectrum_DOn->SetName("numubar_N_reweighted_wClip_spectrum_DOn");
    numubar_N_reweighted_wClip_spectrum_DOn->SetTitle("numubar_N_reweighted_wClip_spectrum_DOn");
    numubar_N_reweighted_wClip_spectrum_DOn->SetLineColor(kGreen);
    numubar_N_reweighted_wClip_spectrum_DOn->Write();

    ALPHA = 0.1;
    TGraph* numubar_N_reweighted_spectrum_wALPHA_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    NormalizeTGraph(numubar_N_reweighted_spectrum_wALPHA_DOn);
    MakeReweightedGraph(numubar_N_reweighted_spectrum_wALPHA_DOn, ALPHA, maxVal, minVal);
    numubar_N_reweighted_spectrum_wALPHA_DOn->SetName("numubar_N_reweighted_spectrum_wALPHA_DOn");
    numubar_N_reweighted_spectrum_wALPHA_DOn->SetTitle("numubar_N_reweighted_spectrum_wALPHA_DOn");
    numubar_N_reweighted_spectrum_wALPHA_DOn->SetLineColor(kGreen);
    numubar_N_reweighted_spectrum_wALPHA_DOn->Write();

    
    TGraph* numubar_N_reweighted_wALPHA_wClip_spectrum_DOn = new TGraph(4900, Energy, UNWEIGHTED_SPECTRA);
    NormalizeTGraph(numubar_N_reweighted_wALPHA_wClip_spectrum_DOn);
    MakeReweightedGraphwClip(numubar_N_reweighted_wALPHA_wClip_spectrum_DOn, ALPHA, maxVal, minVal);
    numubar_N_reweighted_wALPHA_wClip_spectrum_DOn->SetName("numubar_N_reweighted_wALPHA_wClip_spectrum_DOn");
    numubar_N_reweighted_wALPHA_wClip_spectrum_DOn->SetTitle("numubar_N_reweighted_wALPHA_wClip_spectrum_DOn");
    numubar_N_reweighted_wALPHA_wClip_spectrum_DOn->SetLineColor(kGreen);
    numubar_N_reweighted_wALPHA_wClip_spectrum_DOn->Write();
    
    
    const double mean_weight = FindMeanWeight(numubar_un_spectrum_DOn, ALPHA, maxVal, minVal);
    std::cout << "mean_weight: " << mean_weight << std::endl;
    for (int i = 0; i < 4900; i++) 
            {
                
                //Evaluate the cross sections at various energy values (via interpolation)
                CCxsec_numubar = tot_cc_numubar->Eval(Energy[i]); if(CCxsec_numubar<0.){CCxsec_numubar=0.;}
                flux_numubar = numubar_flux_DOn->Eval(Energy[i]); if(flux_numubar<0.){flux_numubar=0.;}
                numubar_spectrum = flux_numubar*CCxsec_numubar*bin_width; //has order E-10
                double weight = (1./(numubar_spectrum+ALPHA*maxVal+minVal))/mean_weight; //weight has order E10, mean weight has order E10
                //double weight = 1./(y+alpha*maxVal+1.);
                // if (i%100==0){std::cout << setprecision(10) <<"weight: " << weight << "    numubar_spectrum: " << numubar_spectrum
                //                      << "    ALPHA: " << ALPHA << "    maxVal: " << maxVal << "    mean_weight: " << mean_weight 
                //                      <<std::endl;}
                //double weight = weight/mean_weight; //has order 1
                if ( weight < 0.2 ){weight = 0.2;}
                if ( weight > 5.0 ){weight = 5.0;}
                //if (i%100==0){std::cout << setprecision(10) <<"         clipped weight: " << weight  <<std::endl;}
                WEIGHTED_SPECTRA[i] = numubar_spectrum*weight;

            }
    TGraph* numubar_w_spectrum_DOn = new TGraph(4900, Energy, WEIGHTED_SPECTRA);
    numubar_w_spectrum_DOn->SetName("DOn_w_spectra_check");
    numubar_w_spectrum_DOn->SetTitle("DOn_w_spectra_check");
    numubar_w_spectrum_DOn->SetLineColor(kBlue);
    numubar_w_spectrum_DOn->Write();

}

/*
void process_inverse_xsec_files()
{
    // List of filenames and corresponding graph names
    vector<pair<string, string>> files = {
        {"numu_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numuCC_flat_spectra"},
        {"numu_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numuNC_flat_spectra"},
        {"numu_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numu_total_flat_spectra"},
        {"numubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numubarCC_flat_spectra"},
        {"numubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numubarNC_flat_spectra"},
        {"numubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "numubar_total_flat_spectra"},
        {"nue_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nueCC_flat_spectra"},
        {"nue_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nueNC_flat_spectra"},
        {"nue_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nue_total_flat_spectra"},
        {"nuebar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nuebarCC_flat_spectra"},
        {"nuebar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nuebarNC_flat_spectra"},
        {"nuebar_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nuebar_total_flat_spectra"},
        {"nutau_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutauCC_flat_spectra"},
        {"nutau_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutauNC_flat_spectra"},
        {"nutau_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutau_total_flat_spectra"},
        {"nutaubar_CC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutaubarCC_flat_spectra"},
        {"nutaubar_NC_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutaubarNC_flat_spectra"},
        {"nutaubar_AR23_0p1-10GeV_inverse_xsec_as_flux.data", "nutaubar_total_flat_spectra"}
    };

    const int inverse_count = 9901; // Total number of data points
    double energy, inverse_xsec_flux;
    double Energy[inverse_count], Inverse_XSec_Flux[inverse_count];
    double FLAT_SPECTRA[inverse_count];

    // Open output ROOT file
    TFile* outputFile = new TFile("output_flat_spectra.root", "RECREATE");

    for (const auto& file : files)
    {
        const string& filename = file.first;
        const string& graphName = file.second;

        ifstream inputFile(filename);
        if (!inputFile.is_open())
        {
            cerr << "Error opening file: " << filename << endl;
            continue;
        }

        // Read data from file into arrays
        for (int i = 0; i < inverse_count; i++)
        {
            inputFile >> energy >> inverse_xsec_flux;
            Energy[i] = energy;
            Inverse_XSec_Flux[i] = inverse_xsec_flux;
        }
        inputFile.close();

        // Compute FLAT_SPECTRA values
        for (int i = 0; i < inverse_count; i++)
        {
            FLAT_SPECTRA[i] = Inverse_XSec_Flux[i]; // Adjust as needed for computations
        }

        // Create a TGraph for the current dataset
        TGraph* graph = new TGraph(inverse_count, Energy, FLAT_SPECTRA);
        graph->SetName(graphName.c_str());
        graph->SetTitle(graphName.c_str());

        // Write the graph to the ROOT file
        graph->Write();
    }

    // Close the ROOT file
    outputFile->Close();
    cout << "All graphs have been written to output_flat_spectra.root" << endl;
}
*/