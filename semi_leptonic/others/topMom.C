#include <unistd.h>
#include <iostream>
#include <string>
#define MAXV 8

using namespace std;
//void asymmetry(string filename = "TTBarProcessorLeft.root", TCanvas * c1 = NULL)
void topMom()
{
		int token=0;
		string filename0 = "/home/ilc/yokugawa/run/root_merge/";
		string filename1;

		/*
		cout << "0 = New/Small" 	  << endl;
		cout << "1 = New/Large" 	  << endl;
		cout << "2 = New/Large/QQbar" << endl;
		cout << "3 = Old      " 	  << endl;
		cout << "4 = Old/yyxylv      "       << endl;
		cout << "Choose from 0-4: ";
		cin  >> token;
		cout << endl;

		switch(token){
				case 0 : filename1 = "new/small/leptonic_yyxyev_eLeR_new_small.root";
						 break;
				case 1 : filename1 = "new/large/leptonic_yyxyev_eLeR_new_large.root";
						 break;
				case 2 : filename1 = "new/large/leptonic_yyxyev_eLeR_new_large_QQbar.root";
						 break;
				case 3 : filename1 = "old/leptonic_yyxyev_eLeR_old_lcut.root" ;
						 break;
				case 4 : filename1 = "old/leptonic_yyxylv_eLeR_iso_lep_lcut.root" ;
						 break;
		}
		*/

		filename1 = "new/large/leptonic_yyxyev_eLeR_new_large_QQbar_withBtrack.root";

		string filename = filename0 + filename1;
		cout << "Processing : " << filename << " ..." << endl;

		TFile * file = TFile::Open(filename.c_str());

		int bin_e = 30;
		int max_e = 1;

		TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,500,500);

		TH1F * hbMom = new TH1F("hbMom","b/W momentum (|cos#theta| = 1); p (GeV); nEvents", 100, 0, 200);
		TH1F * hWMom   = new TH1F("hWMom","b/W momentum (|cos#theta| = 1); p (GeV); nEvents", 100, 0, 200);

		//TH1F * cosReco = new TH1F("cosReco", "E(Ntracks)", bin_e,-1.0,max_e);
		//cosReco->Sumw2();

		TTree * normaltree = (TTree*) file->Get( "Stats" ) ;
		TTree * GenTree = (TTree*) file->Get( "GenTree" ) ;

		// Selection lists
		TCut thru = "Thrust < 0.9";
		TCut hadM = "hadMass > 180 && hadMass < 420";
		TCut rcTW = "Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50";
		TCut pcut = "Top1bmomentum > 15 && Top2bmomentum > 15";
		TCut gcut = "(Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2";

		// Methods selection
		TCut methodAll = "methodTaken > 0";
		TCut method1 = "methodTaken == 1";
		TCut method2 = "methodTaken == 2";
		TCut method3 = "methodTaken == 3";
		TCut method4 = "methodTaken == 4";
		TCut method5 = "methodTaken == 5";
		TCut method6 = "methodTaken == 6";
		TCut method7 = "methodTaken == 7";

		// Total cut applied
		TCut cuts = rcTW + hadM + pcut + gcut + methodAll;

		TCut cosP = "qCostheta == 1 " ;
		TCut cosN = "qCostheta == -1" ;

		TCut cos  = cosP || cosN ;

		//TCut fcuts = "qCostheta > 0" + cuts;
		//TCut bcuts = "qCostheta < 0 && qCostheta > -1.0 " + cuts;
		
		// Fill histograms from tree
		//normaltree->Draw("W1momentum    >> hWMom");
		normaltree->Draw("Top1bmomentum >> hbMom");
		
		c1->cd();
		hbMom->Draw();
		hWMom->Draw();
		c1->Update();

}

