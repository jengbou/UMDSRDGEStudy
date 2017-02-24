#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TLine.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace TMath;
using namespace std;

//---------------------------
//cheap sigmoid:
//0 if <= 0, 1 if >= 1
Double_t sigmoid(Double_t x){
	return Max(0.,Min(x,1.));
}

//---------------------------------------------
//model:
//from Wick (NIM B61 (1991) 472-486)
//z0^2 = g/R
//L = a(D) - f(D)*2z0/d = a(D) - f(D)*2/d*sqrt(g/R)
//if(2z0>=d) L = a(d) - f(D)
//parameters:
//g, f0, a0
//0,  1,  2
Double_t wick(Double_t *x, Double_t *par){
	//keep parameters in appropriate ranges
	par[0] = Abs(par[0]);
	par[1] = Abs(par[1]);
	par[2] = Abs(par[2]);

	double L = 0; //light yield
	double R = x[0]; //dose rate
	double D = x[1]; //total dose
	double f0 = par[1]; //param for f(D)
	double f = sigmoid(1 - exp(-f0*D)); //fraction of fluor destroyed, dependent on dose (poisson model)
	//cout << D << " -> " << (int)D/2 << ", " << f << endl;
	double g = par[0]; //scaling factor for z0 vs dose rate
	double z0 = Sqrt(g/R); //diffusion thickness
	double d = 4; //total thickness in mm
	double a0 = par[2]; //param for a(D)
	double a = sigmoid(1-a0*D); //light yield asymptote
	
	//if(2*z0 >= d) L = a - f;
	//else L = a - f*2*z0/d;
	L = Max(a - f*2*z0/d, a - f);
	
	return L;
}

//-----------------------------------------------------
//1d version of above function (only depends on R)
//parameters:
//g, f0, a0, D
//0,  1,  2, 3
Double_t wick1d(Double_t *x, Double_t *par){
	//keep parameters in appropriate ranges
	par[0] = Abs(par[0]);
	par[1] = Abs(par[1]);
	par[2] = Abs(par[2]);
	par[3] = Abs(par[3]);

	double L = 0; //light yield
	double R = x[0]; //dose rate
	double D = par[3]; //total dose
	double f0 = par[1]; //param for f(D)
	double f = sigmoid(1 - exp(-f0*D)); //fraction of fluor destroyed, dependent on dose (poisson model)
	//cout << D << " -> " << (int)D/2 << ", " << f << endl;
	double g = par[0]; //scaling factor for z0 vs dose rate
	double z0 = Sqrt(g/R); //diffusion thickness
	double d = 4; //total thickness in mm
	double a0 = par[2]; //param for a(D)
	double a = sigmoid(1-a0*D); //light yield asymptote
	
	//if(2*z0 >= d) L = a - f;
	//else L = a - f*2*z0/d;
	L = Max(a - f*2*z0/d, a - f);
	
	return L;
}

//----------------------------------------------
//fit the above function to Biagtan data
void fitWick5s(){
	//fitting data from Biagtan (NIM B108 (1996) 125-128)
	//points normalized to max light yield 100% (from 110% in paper)
	Double_t R[] = {0.014, 0.04, 1.5, 0.014, 0.04, 1.5, 0.014, 0.04, 0.14, 1.5, 0.014, 0.04, 0.14, 0.0022, 0.0043, 0.0064, 0.014, 0.14, 1.5};
	Double_t D[] = {2, 2, 2, 4, 4, 4, 6, 6, 6, 6, 8, 8, 8, 10, 10, 10, 10, 10, 10};
	//Bicron-499-35 
	//Double_t L[] = {0.8531, 0.8747, 0.9258, 0.9867, 0.7057, 0.7941, 0.8432, 0.9572, 0.5779, 0.6821, 0.7823, 0.9278, 0.4875, 0.5111, 0.6801, 0.8550, 0.1514, 0.2005, 0.2654, 0.3794, 0.3971, 0.5661, 0.7961};
	//SCSN-81
	Double_t L[] = {0.7850, 0.8549, 0.9361, 0.7444, 0.7624, 0.8977, 0.6564, 0.6992, 0.7602, 0.8436, 0.5820, 0.6271, 0.6722, 0.3451, 0.3925, 0.4602, 0.4737, 0.6564, 0.7511};
	
	//make graph
	TGraph2D* g2d = new TGraph2D("g2d","",19,R,D,L);
	g2d->GetXaxis()->SetTitle("dose rate [Mrad/hr]");
	g2d->GetYaxis()->SetTitle("total dose [Mrad]");
	g2d->GetZaxis()->SetTitle("light yield [%]");
	
	TF2* gfit = new TF2("wick",wick,1e-3,10,2,10,3);
	gfit->SetParameter(0,0.0099); //wick value, converted to Mrad from Gy
	gfit->SetParameter(1,0.02); //based on plot from wick3
	gfit->SetParameter(2,0.02); //based on plot from wick3
	//formatting
	gfit->SetLineColor(kRed);
	gfit->SetMarkerColor(kRed);
	gfit->SetLineWidth(2);
	
	g2d->Fit(gfit,"NR");
	
	//draw 2d stuff
	//TCanvas* can = new TCanvas("can","can",900,700);
	//can->SetLogx();
	//g2d->Draw("P");
	//gfit->Draw("same c");
	
	//draw 1d stuff
	TCanvas* can = new TCanvas("fitWick5s","fitWick5s",700,500);
	can->cd();
	can->SetLogx();
	
	TLegend* leg = new TLegend(0.7,0.25,0.8,0.5);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.05);
	leg->SetTextFont(42);
	
	TH1F* haxis = new TH1F("axis","",100,1e-3,10);
	haxis->GetYaxis()->SetRangeUser(0,1);
	haxis->GetXaxis()->SetTitle("dose rate [Mrad/hr]");
	haxis->GetYaxis()->SetTitle("light yield [%]");
	haxis->Draw("hist");
	
	TF1* gfit1d[5];
	TGraph* b1d[5];
	Double_t Dr[] = {2,4,6,8,10};
	Int_t npt[] = {3,3,4,3,6};
	Int_t ind = 0;
	Int_t marker[] = {21,25,20,24,22};
	//Double_t aD[5] = {};
	//Double_t aDe[5] = {};
	//Double_t fD[5] = {};
	//Double_t fDe[5] = {};
	for(int j = 0; j < 5; j++){
		//make 1d graph
		stringstream gname;
		gname << "b" << j+1;
		b1d[j] = new TGraph(npt[j],R+ind,L+ind);
		b1d[j]->SetTitle("");
		b1d[j]->GetXaxis()->SetTitle("dose rate [Mrad/hr]");
		b1d[j]->GetYaxis()->SetTitle("light yield [%]");
		b1d[j]->GetYaxis()->SetRangeUser(0,1);
		b1d[j]->SetMarkerStyle(marker[j]);
		stringstream legname;
		legname << Dr[j] << " Mrad";
		leg->AddEntry(b1d[j],(legname.str()).c_str(),"p");
		
		//make 1d fn
		stringstream fname;
		fname << "wick" << j+1;
		gfit1d[j] = new TF1((fname.str()).c_str(),wick1d,1e-3,10,4);
		gfit1d[j]->SetParameter(0,gfit->GetParameter(0));
		gfit1d[j]->SetParameter(1,gfit->GetParameter(1));
		gfit1d[j]->SetParameter(2,gfit->GetParameter(2));
		gfit1d[j]->SetParameter(3,Dr[j]);
		gfit1d[j]->SetLineColor(kRed);
		gfit1d[j]->SetMarkerColor(kRed);
		gfit1d[j]->SetLineWidth(2);
		
		b1d[j]->Draw("PZ same");
		gfit1d[j]->Draw("same");
		
		//aD[j] = gfit->GetParameter(j+1+5);
		//aDe[j] = gfit->GetParError(j+1+5);
		//fD[j] = gfit->GetParameter(j+1);
		//fDe[j] = gfit->GetParError(j+1);
		
		ind+=npt[j];
	}
	leg->Draw("same");
/*	
	TCanvas* can2 = new TCanvas("fitWick3a","fitWick3a",700,500);
	can2->cd();
	
	TGraphErrors* g_a = new TGraphErrors(5,Dr,aD,0,aDe);
	g_a->SetTitle("");
	g_a->GetXaxis()->SetTitle("total dose [Mrad]");
	g_a->GetYaxis()->SetTitle("asymptotic light yield [%]");
	g_a->Draw("APZ");
	
	TCanvas* can3 = new TCanvas("fitWick3f","fitWick3f",700,500);
	can3->cd();
	
	TGraphErrors* g_f = new TGraphErrors(5,Dr,fD,0,aDe);
	g_f->SetTitle("");
	g_f->GetXaxis()->SetTitle("total dose [Mrad]");
	g_f->GetYaxis()->SetTitle("fraction of fluors destroyed");
	g_f->Draw("APZ");
*/	
}

/*
results: SCSN-81
****************************************
Minimizer is Minuit / Migrad
Chi2                      =    0.0378756
NDf                       =           16
Edm                       =  7.39252e-10
NCalls                    =           70
p0                        =    0.0239684   +/-   0.00682766  
p1                        =    0.0434038   +/-   0.00730744  
p2                        =    0.0284216   +/-   0.00311842  

*/
