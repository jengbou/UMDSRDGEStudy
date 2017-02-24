#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TImageDump.h>
#include <TFile.h>
#include <TTree.h>
#include "Math/Interpolator.h"
#include <HEDarkening.C>
#include <fitWick5s.C>

using TMath::Exp;
using TMath::Power;
using TMath::Log;

void processDoseLine(string line, int& ieta, double& dose, double& fluence)
{
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {stringstream temp(field); temp >> ieta;}
		else if (ctr==4) {stringstream temp(field); temp >> dose;}
		else if (ctr==5) {stringstream temp(field); temp >> fluence;}
		else {}
		
		ctr++;
	}
}

void readDoseFile(string inputFile, std::map<int, double>& doseMap, std::map<int, double>& fluenceMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				int ieta;
				double dose, fluence;
				processDoseLine(line, ieta, dose, fluence);
				doseMap.insert(std::pair<int, double>(ieta, dose));
				fluenceMap.insert(std::pair<int, double>(ieta, fluence / 1e12));
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

Double_t ExponentialDecay(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = Exp(par[1]*(xx) + par[0]);//par[3]*pow(xx, 3) + par[2]*pow(xx, 2) + par[1]*(xx) + par[0]);
	return f;
}

float biagtan( float dose, float doserate ) //in Mrad, Mrad/hr
{
	float aD = 87.3 - 1.15*dose; //intercept
	float bD = 3.54 + 1.11*dose; //slope
	float part2 = 1;
	if(doserate>0) part2 = (aD + bD*log(doserate))/100;
	return part2;
}

float wickModel( float dose, float doserate ) //in Mrad, Mrad/hr
{
	double wickPar[3] = {0.0239684, 0.0434038, 0.0284216};
	double x[2] = {dose, doserate};
	return wick(x, wickPar);
}

void processLaserDataLine(string line, double& lumi, double* ratios)
{
	// ratio[0] through ratio[12] corresponds to ieta 17 to 29.
	// Make sure "ratios" has 13 elements.
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {
			stringstream temp(field);
			temp >> lumi;
			//cout << lumi;
		}
		else if (true) {//ctr >= 1 && ctr <= 13) {
			stringstream temp(field); 
			double ratio;
			temp >> ratio;
			//cout << ratio << "\t";
			ratios[ctr - 1] = ratio;
		}
		else {
			cout << endl;
			break;
		}
		
		ctr++;
	}
}

void readLaserDataFile(string inputFile, vector<double>& lumis, std::map< int, vector<double> > & ratiosMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				lumis.push_back(lumi);
				for (int i=0; i < 13; i++) {
					ratiosMap[i+17].push_back(ratios[i]);
					//cout << "ratios[i] " <<ratios[i] << endl;
					//cout << "ratiosMap[i+17] " << (ratiosMap[i+17])[1] << endl;
				}
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void ratiosMap(int ieta = 29) //bool doPlot = false)
{
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

	TCanvas* c1 = new TCanvas("c1","canvas1",800,600);
	TGraph* graph;
	graph = new TGraph(lumis.size(), &lumis[0], &((ratiosMap[ieta])[0]));
	c1->cd();
	graph->Draw("AP");
}

void plotRochCurve1()
{
	TCanvas* c1;
	c1 = new TCanvas("c1", "degradation vs dose Roch 1", 1200, 600);
	c1->SetLogx();
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.40);
	//c1->SetLogy();
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	TLegend* leg = new TLegend(0.62, 0.5, 0.95, 0.95);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);

	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

	HEDarkening heModel;
	// ieta tile closesest to 5cm x 8cm Rochester tile 
	{
		int ieta = 29;
		int n = lumis.size();
		vector<float> xV, yV; 
		cout << "ieta\tdose\tratio" << endl;
		for (int i=0; i<n; i++) {
			float lumi = lumis[i];
			float dose = lumi*heModel.towerDose[ieta-16][1+1];
			float ratio = ratiosMap[ieta][i];
			xV.push_back(dose);
			yV.push_back(ratio);
			cout << ieta << "\t" << dose << "\t" << ratio << endl;
		}
		TGraph* graph1 = new TGraph(n, &xV[0], &yV[0]);
		graph1->SetTitle("");
		graph1->GetXaxis()->SetTitle("dose [Mrad]");
		graph1->GetYaxis()->SetTitle("degradation in light output");
		graph1->GetXaxis()->SetLimits(1e-5, 1e2);
		graph1->GetXaxis()->SetTitleOffset(1.0);
		graph1->GetYaxis()->SetTitleOffset(0.8);
		graph1->GetYaxis()->SetRangeUser(-0.2, 1.2);
		graph1->SetMarkerColor(colors[8]);
		graph1->SetMarkerSize(1.0);
		graph1->SetMarkerStyle(21);
		graph1->Draw("AP");
		TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e+2, 2);
		Double_t par[] = {1.0, 1.0};
		functionLaserData->SetParameters(par);
		graph1->Fit(functionLaserData, "", "", 1e-4, 1e+2);

		leg->AddEntry(graph1, "i#eta 29 (~6cm x 8cm)", "p");
	}

	// Graph measurements from CMS IN 2001-022.
	// Rochester IN-2001-022: y1: 5cm x 8cm, y2: 12cm x 8cm, y3: 20cm x 20cm (x 0.4cm for all)
	int n = 5;
	float Rochx[5] = {0.1, 1, 2, 5, 10};
	float Roch1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
	float Roch2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
	float Roch3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch1[i];
		}
		TGraph* graphRoch1 = new TGraph(n, Rochx, yArray1);
		graphRoch1->SetMarkerColor(kRed);
		graphRoch1->SetMarkerSize(1.1);
		graphRoch1->SetMarkerStyle(21);
		graphRoch1->Draw("p same");
		leg->AddEntry(graphRoch1, "Rochester data (5cm x 8cm)", "p");
	}
	// Rochester data with naive model. (Extrapolation from Biagtan data points.)
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch1[i];
			yArray1[i] *= biagtan(Rochx[i], 1e-5) / biagtan(Rochx[i], 1e-1);
			//cout << "dose rate factor = " << biagtan(Rochx[i], 1e-1) / biagtan(Rochx[i], 1e-5) << endl;
			//cout << "biagtan @0.1 Mrad/hr = " << biagtan(Rochx[i], 1e-1) << endl;
			//cout << "biagtan @0.00001 Mrad/hr = " << biagtan(Rochx[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch1 = new TGraph(n, Rochx, yArray1);
		graphRoch1->SetMarkerColor(colors[5]);
		graphRoch1->SetMarkerSize(1.0);
		graphRoch1->SetMarkerStyle(21);
		graphRoch1->Draw("p same");
		//TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e+2, 2);
		//Double_t par[] = {1.0, 1.0};
		//functionLaserData->SetParameters(par);
		//graphRoch1->Fit(functionLaserData, "", "", 1e-4, 1e+2);
		leg->AddEntry(graphRoch1, "Rochester data w/ naive model", "p");
	}
	// Rochester data with Wick model.
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch1[i];
			yArray1[i] *= wickModel(Rochx[i], 1e-5) / wickModel(Rochx[i], 1e-1);
			//cout << "dose rate factor = " << wickModel(Rochx[i], 1e-1) / wickModel(Rochx[i], 1e-5) << endl;
			//cout << "wickModel @0.1 Mrad/hr = " << wickModel(Rochx[i], 1e-1) << endl;
			//cout << "wickModel @0.00001 Mrad/hr = " << wickModel(Rochx[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch1 = new TGraph(n, Rochx, yArray1);
		graphRoch1->SetMarkerColor(colors[6]);
		graphRoch1->SetMarkerSize(0.9);
		graphRoch1->SetMarkerStyle(22);
		graphRoch1->Draw("p same");
		leg->AddEntry(graphRoch1, "Rochester data w/ Wick model", "p");
	}
	leg->Draw();
}

void plotRochCurve2()
{
	TCanvas* c1;
	c1 = new TCanvas("c1", "degradation vs dose Roch 2", 1200, 600);
	c1->SetLogx();
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.40);
	//c1->SetLogy();
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	TLegend* leg = new TLegend(0.62, 0.5, 0.95, 0.95);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);

	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

	HEDarkening heModel;
	// ieta tile closesest to 12cm x 8cm Rochester tile 
	{
		int ieta = 28;
		int n = lumis.size();
		vector<float> xV, yV; 
		cout << "ieta\tdose\tratio" << endl;
		for (int i=0; i<n; i++) {
			float lumi = lumis[i];
			float dose = lumi*heModel.towerDose[ieta-16][1+1];
			float ratio = ratiosMap[ieta][i];
			xV.push_back(dose);
			yV.push_back(ratio);
			cout << ieta << "\t" << dose << "\t" << ratio << endl;
		}
		TGraph* graph1 = new TGraph(n, &xV[0], &yV[0]);
		graph1->SetTitle("");
		graph1->GetXaxis()->SetTitle("dose [Mrad]");
		graph1->GetYaxis()->SetTitle("degradation in light output");
		graph1->GetXaxis()->SetLimits(1e-5, 1e2);
		graph1->GetXaxis()->SetTitleOffset(1.0);
		graph1->GetYaxis()->SetTitleOffset(0.8);
		graph1->GetYaxis()->SetRangeUser(-0.2, 1.2);
		graph1->SetMarkerColor(colors[8]);
		graph1->SetMarkerSize(1.0);
		graph1->SetMarkerStyle(21);
		graph1->Draw("AP");
		TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e+2, 2);
		Double_t par[] = {1.0, 1.0};
		functionLaserData->SetParameters(par);
		graph1->Fit(functionLaserData, "", "", 1e-4, 1e+2);

		leg->AddEntry(graph1, "i#eta 28 (~12cm x 9cm)", "p");
	}

	// Graph measurements from CMS IN 2001-022.
	// Rochester IN-2001-022: y1: 5cm x 8cm, y2: 12cm x 8cm, y3: 20cm x 20cm (x 0.4cm for all)
	int n = 5;
	float Rochx[5] = {0.1, 1, 2, 5, 10};
	float Roch1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
	float Roch2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
	float Roch3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch2[i];
		}
		TGraph* graphRoch2 = new TGraph(n, Rochx, yArray1);
		graphRoch2->SetMarkerColor(kRed);
		graphRoch2->SetMarkerSize(1.1);
		graphRoch2->SetMarkerStyle(21);
		graphRoch2->Draw("p same");
		leg->AddEntry(graphRoch2, "Rochester data (12cm x 8cm)", "p");
	}
	// Rochester data with naive model. (Extrapolation from Biagtan data points.)
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch2[i];
			yArray1[i] *= biagtan(Rochx[i], 1e-5) / biagtan(Rochx[i], 1e-1);
			//cout << "dose rate factor = " << biagtan(Rochx[i], 1e-1) / biagtan(Rochx[i], 1e-5) << endl;
			//cout << "biagtan @0.1 Mrad/hr = " << biagtan(Rochx[i], 1e-1) << endl;
			//cout << "biagtan @0.00001 Mrad/hr = " << biagtan(Rochx[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch2 = new TGraph(n, Rochx, yArray1);
		graphRoch2->SetMarkerColor(colors[5]);
		graphRoch2->SetMarkerSize(1.0);
		graphRoch2->SetMarkerStyle(21);
		graphRoch2->Draw("p same");
		//TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e+2, 2);
		//Double_t par[] = {1.0, 1.0};
		//functionLaserData->SetParameters(par);
		//graphRoch2->Fit(functionLaserData, "", "", 1e-4, 1e+2);
		leg->AddEntry(graphRoch2, "Rochester data w/ naive model", "p");
	}
	// Rochester data with Wick model.
	{
		float yArray1[5];
		for (int i=0; i<5; i++) {
			yArray1[i] = Roch2[i];
			yArray1[i] *= wickModel(Rochx[i], 1e-5) / wickModel(Rochx[i], 1e-1);
			//cout << "dose rate factor = " << wickModel(Rochx[i], 1e-1) / wickModel(Rochx[i], 1e-5) << endl;
			//cout << "wickModel @0.1 Mrad/hr = " << wickModel(Rochx[i], 1e-1) << endl;
			//cout << "wickModel @0.00001 Mrad/hr = " << wickModel(Rochx[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch2 = new TGraph(n, Rochx, yArray1);
		graphRoch2->SetMarkerColor(colors[6]);
		graphRoch2->SetMarkerSize(0.9);
		graphRoch2->SetMarkerStyle(22);
		graphRoch2->Draw("p same");
		leg->AddEntry(graphRoch2, "Rochester data w/ Wick model", "p");
	}
	leg->Draw();
}

void plotDamageVsDose()
{
	TCanvas* c1;
	c1 = new TCanvas("c1", "degradation vs dose", 1200, 600);
	c1->SetLogx();
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.40);
	//c1->SetLogy();
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	TLegend* leg = new TLegend(0.62, 0.3, 0.95, 0.95);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);

	HEDarkening heModel;
	// Graph using old dose map.
	{
		std::map<int, double> doseMapLayer1, fluenceMapLayer1;
		readDoseFile("../doseLayer1.txt", doseMapLayer1, fluenceMapLayer1);
		vector<float> xV, yV;	
		cout << "ieta\tdose\tdegradation" << endl;
		for (int ieta = 17; ieta < 30; ieta ++) {
			float degradation = heModel.degradation(22, ieta, 1+1); // 22fb^-1, ieta, layer 1
			float dose = 0.22/10*doseMapLayer1[ieta]; //22fb^-1. Doses converted from kGy.
			cout << ieta << "\t" << dose << "\t" << degradation << endl;
			xV.push_back(dose);
			yV.push_back(degradation);
		}
		TGraph* graphLaserData = new TGraph(xV.size(), &xV[0], &yV[0]);
		graphLaserData->SetTitle("");
		graphLaserData->GetXaxis()->SetTitle("dose [Mrad]");
		graphLaserData->GetYaxis()->SetTitle("degradation in light output");
		graphLaserData->GetXaxis()->SetLimits(1e-5, 1e2);
		graphLaserData->GetXaxis()->SetTitleOffset(1.0);
		graphLaserData->GetYaxis()->SetTitleOffset(0.8);
		graphLaserData->GetYaxis()->SetRangeUser(-0.2, 1.2);
		graphLaserData->SetMarkerColor(colors[3]);
		graphLaserData->SetMarkerSize(1.0);
		graphLaserData->SetMarkerStyle(25);

		graphLaserData->Draw("AP");
		leg->AddEntry(graphLaserData, "laser data (previous dose map)", "p");
	}
	// Graph curve extracted from laser data with MARS dose.
	{
		vector<float> xV, yV;	
		for (int ieta = 16; ieta < 30; ieta ++) {
			float degradation = heModel.degradation(22, ieta, 1+1); // 22fb^-1, ieta, layer 1
			float dose = 22*heModel.towerDose[ieta-16][1+1];
			xV.push_back(dose);
			yV.push_back(degradation);
		}
		//for (int ieta = 16; ieta < 30; ieta ++) {
		//	float degradation = heModel.degradation(22, ieta, 7+1); // 22fb^-1, ieta, layer 7
		//	float dose = 22*heModel.towerDose[ieta-16][7+1];
		//	xV.push_back(dose);
		//	yV.push_back(degradation);
		//}
		TGraph* graphLaserData1 = new TGraph(xV.size(), &xV[0], &yV[0]);
		graphLaserData1->SetMarkerColor(colors[0]);
		graphLaserData1->SetMarkerSize(1.0);
		graphLaserData1->SetMarkerStyle(24);

		graphLaserData1->Draw("P same");
		leg->AddEntry(graphLaserData1, "laser data (MARS dose)", "p");

		//TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e-2, 2);
		//Double_t par[] = {1.0, 1.0};
		//functionLaserData->SetParameters(par);
		//graphLaserData1->Fit(functionLaserData, "", "", 1e-5, 1e-2);
	}
	// Graph curve extracted from laser data with MARS dose and corrected for geometry.
	{
		// p1 paramter from run11 simulation fits. ieta 22 to 29.
		float p1Array[8] = {-39.1601, -38.3167, -35.9355, -34.5196, -32.7535, -27.2942, -27.0732, -18.5917};
		float lumiFactor[8];
		for (int i = 0; i<8; i++) {
			lumiFactor[i] = p1Array[i]/p1Array[5];
		}
		
		vector<float> xV, yV;	
		for (int ieta = 22; ieta < 30; ieta++) {
			float degradation = heModel.degradation(22 / lumiFactor[ieta-22], ieta, 1+1); // 22fb^-1 * lumiFactor, ieta, layer 1
			float dose = 22*heModel.towerDose[ieta-16][1+1];
			xV.push_back(dose);
			yV.push_back(degradation);
		}
		TGraph* graphLaserGeomCorr = new TGraph(xV.size(), &xV[0], &yV[0]);
		graphLaserGeomCorr->SetTitle("");
		graphLaserGeomCorr->GetXaxis()->SetTitle("dose [Mrad]");
		graphLaserGeomCorr->GetYaxis()->SetTitle("degradation in light output");
		graphLaserGeomCorr->GetXaxis()->SetLimits(1e-5, 1e2);
		graphLaserGeomCorr->GetYaxis()->SetRangeUser(0, 1);
		graphLaserGeomCorr->SetMarkerColor(colors[8]);
		graphLaserGeomCorr->SetMarkerSize(1.0);
		graphLaserGeomCorr->SetMarkerStyle(26);

		graphLaserGeomCorr->Draw("P same");
		leg->AddEntry(graphLaserGeomCorr, "#splitline{laser data (MARS dose)}{+ geometry correction}", "p");

		TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e2, 2);
		Double_t par[] = {1.0, 1.0};
		functionLaserData->SetParameters(par);
		graphLaserGeomCorr->Fit(functionLaserData, "", "", 1e-3, 1e0);
	}

	// Graph Vasken measurements.
	{
		int n = 7;
		float xArray[7] = { 0.780743079, 1.863624049, 3.447142869, 4.831022809, 6.376175476, 8.415529239, 11.79400311 };
		float yArray[7] = { 0.859726027, 0.700273973, 0.598356164, 0.45369863, 0.394520548, 0.276164384, 0.170958904 };
		TGraph* graphVasken = new TGraph(n, xArray, yArray);
		graphVasken->Draw("p same");
		graphVasken->SetMarkerColor(colors[1]);
		graphVasken->SetMarkerSize(1.0);
		graphVasken->SetMarkerStyle(21);
		leg->AddEntry(graphVasken, "Vasken data", "p");
	}

	// Graph measurements from CMS IN 2001-022.
	// Rochester IN-2001-022: y1: 5cm x 8cm, y2: 12cm x 8cm, y3: 20cm x 20cm (x 0.4cm for all)
	{
		int n = 5;
		float xArray[5] = {0.1, 1, 2, 5, 10};
		float yArray1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
		float yArray2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
		float yArray3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
		
		TGraph* graphRoch1 = new TGraph(n, xArray, yArray1);
		graphRoch1->SetMarkerColor(kRed);
		graphRoch1->SetMarkerSize(1.1);
		graphRoch1->SetMarkerStyle(21);
		graphRoch1->Draw("p same");
		//TGraph* graphRoch2 = new TGraph(n, xArray, yArray2);
		//graphRoch2->SetMarkerColor(colors[4]);
		//graphRoch2->SetMarkerSize(0.5);
		//graphRoch2->SetMarkerStyle(21);
		//graphRoch2->Draw("p same");
		//TGraph* graphRoch3 = new TGraph(n, xArray, yArray3);
		//graphRoch3->SetMarkerColor(colors[4]);
		//graphRoch3->SetMarkerSize(0.5);
		//graphRoch3->SetMarkerStyle(3);
		//graphRoch3->Draw("p same");
		leg->AddEntry(graphRoch1, "Rochester data (5cm x 8cm)", "p");
	}
	// Rochester data with naive model. (Extrapolation from Biagtan data points.)
	{
		int n = 5;
		float xArray[5] = {0.1, 1, 2, 5, 10};
		float yArray1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
		float yArray2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
		float yArray3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
		for (int i=0; i<5; i++) {
			yArray1[i] *= biagtan(xArray[i], 1e-5) / biagtan(xArray[i], 1e-1);
			//cout << "dose rate factor = " << biagtan(xArray[i], 1e-1) / biagtan(xArray[i], 1e-5) << endl;
			//cout << "biagtan @0.1 Mrad/hr = " << biagtan(xArray[i], 1e-1) << endl;
			//cout << "biagtan @0.00001 Mrad/hr = " << biagtan(xArray[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch1 = new TGraph(n, xArray, yArray1);
		graphRoch1->SetMarkerColor(colors[5]);
		graphRoch1->SetMarkerSize(1.0);
		graphRoch1->SetMarkerStyle(21);
		graphRoch1->Draw("p same");
		//TGraph* graphRoch2 = new TGraph(n, xArray, yArray2);
		//graphRoch2->SetMarkerColor(colors[5]);
		//graphRoch2->SetMarkerSize(0.5);
		//graphRoch2->SetMarkerStyle(21);
		//graphRoch2->Draw("p same");
		//TGraph* graphRoch3 = new TGraph(n, xArray, yArray3);
		//graphRoch3->SetMarkerColor(colors[5]);
		//graphRoch3->SetMarkerSize(0.5);
		//graphRoch3->SetMarkerStyle(3);
		//graphRoch3->Draw("p same");
		//TF1* functionLaserData = new TF1("functionLaserData", ExponentialDecay, 1e-5, 1e+2, 2);
		//Double_t par[] = {1.0, 1.0};
		//functionLaserData->SetParameters(par);
		//graphRoch1->Fit(functionLaserData, "", "", 1e-4, 1e+2);
		leg->AddEntry(graphRoch1, "Rochester data w/ naive model", "p");
	}
	// Rochester data with Wick model.
	{
		int n = 5;
		float xArray[5] = {0.1, 1, 2, 5, 10};
		float yArray1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
		float yArray2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
		float yArray3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
		for (int i=0; i<5; i++) {
			yArray1[i] *= wickModel(xArray[i], 1e-5) / wickModel(xArray[i], 1e-1);
			//cout << "dose rate factor = " << wickModel(xArray[i], 1e-1) / wickModel(xArray[i], 1e-5) << endl;
			//cout << "wickModel @0.1 Mrad/hr = " << wickModel(xArray[i], 1e-1) << endl;
			//cout << "wickModel @0.00001 Mrad/hr = " << wickModel(xArray[i], 1e-5) << endl;
		}
		
		TGraph* graphRoch1 = new TGraph(n, xArray, yArray1);
		graphRoch1->SetMarkerColor(colors[6]);
		graphRoch1->SetMarkerSize(0.9);
		graphRoch1->SetMarkerStyle(22);
		graphRoch1->Draw("p same");
		//TGraph* graphRoch2 = new TGraph(n, xArray, yArray2);
		//graphRoch2->SetMarkerColor(colors[5]);
		//graphRoch2->SetMarkerSize(0.5);
		//graphRoch2->SetMarkerStyle(21);
		//graphRoch2->Draw("p same");
		//TGraph* graphRoch3 = new TGraph(n, xArray, yArray3);
		//graphRoch3->SetMarkerColor(colors[5]);
		//graphRoch3->SetMarkerSize(0.5);
		//graphRoch3->SetMarkerStyle(3);
		//graphRoch3->Draw("p same");
		leg->AddEntry(graphRoch1, "Rochester data w/ Wick model", "p");
	}
	leg->Draw();
}

void fitRochCurves()
{
	TCanvas* c1;
	c1 = new TCanvas("c1", "degradation vs dose", 1200, 600);
	c1->SetLogx();
	c1->SetLeftMargin(0.10);
	c1->SetRightMargin(0.40);

	// Graph measurements from CMS IN 2001-022.
	// Rochester IN-2001-022: y1: 5cm x 8cm, y2: 12cm x 8cm, y3: 20cm x 20cm (x 0.4cm for all)
	{
		int n = 5;
		float xArray[5] = {0.1, 1, 2, 5, 10};
		float yArray1[5] = {1, 0.702997275, 0.591280654, 0.201634877, 0.046321526};
		float yArray2[5] = {1, 0.670299728, 0.525885559, 0.133514986, 0.038147139};
		float yArray3[5] = {0.85013624, 0.365122616, 0.201634877, 0.103542234, 0.008174387};
		
		TGraph* graphRoch1 = new TGraph(n, xArray, yArray1);
		graphRoch1->SetTitle("");
		graphRoch1->GetXaxis()->SetTitle("dose [Mrad]");
		graphRoch1->GetYaxis()->SetTitle("degradation in light output");
		graphRoch1->GetXaxis()->SetLimits(1e-3, 1e2);
		graphRoch1->GetXaxis()->SetTitleOffset(1.0);
		graphRoch1->GetYaxis()->SetTitleOffset(0.8);
		graphRoch1->GetYaxis()->SetRangeUser(-0.2, 1.2);
		graphRoch1->SetMarkerColor(kRed);
		graphRoch1->SetMarkerSize(1.1);
		graphRoch1->SetMarkerStyle(21);
		graphRoch1->Draw("ap");

		TGraph* graphRoch2 = new TGraph(n, xArray, yArray2);
		graphRoch2->SetMarkerColor(kBlue);
		graphRoch2->SetMarkerSize(1.1);
		graphRoch2->SetMarkerStyle(21);
		graphRoch2->Draw("p same");

		TGraph* graphRoch3 = new TGraph(n, xArray, yArray3);
		graphRoch3->SetMarkerColor(kGreen);
		graphRoch3->SetMarkerSize(1.1);
		graphRoch3->SetMarkerStyle(21);
		graphRoch3->Draw("p same");

		Double_t par[] = {1.0, 1.0};
		double D0;

		TF1* functionRoch1 = new TF1("functionRoch1", ExponentialDecay, 1e-5, 1e2, 2);
		functionRoch1->SetParameters(par);
		graphRoch1->Fit(functionRoch1, "", "", 1e-2, 1e1);
		D0 = - 1 / functionRoch1->GetParameter(1);
		cout << "D0: " << D0 << endl;


		TF1* functionRoch2 = new TF1("functionRoch2", ExponentialDecay, 1e-5, 1e2, 2);
		functionRoch2->SetParameters(par);
		graphRoch2->Fit(functionRoch2, "", "", 1e-2, 1e1);
		D0 = - 1 / functionRoch2->GetParameter(1);
		cout << "D0: " << D0 << endl;


		TF1* functionRoch3 = new TF1("functionRoch3", ExponentialDecay, 1e-5, 1e2, 2);
		functionRoch3->SetParameters(par);
		graphRoch3->Fit(functionRoch3, "", "", 1e-2, 1e1);
		D0 = - 1 / functionRoch3->GetParameter(1);
		cout << "D0: " << D0 << endl;
	}
}
