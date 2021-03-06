#include <iostream>
#include <string>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TCut.h>
#include <TEventList.h>
#include <TStyle.h>

using namespace std ;

void FillandAdd( int b1 = 0, int b2 = 0, float tcos1 = 0, float tcos2 = 0, float bcos1 = 0, float bcos2 = 0, double weight = 0, TH1F* top = 0, TH1F* bottom = 0){

		if( b1 < 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			top->Fill( tcos1, weight ) ;
			//top->Fill( -tcos2, weight ) ;
			bottom->Fill( bcos1, weight ) ;
			bottom->Fill( -bcos2, weight ) ;
		}else if( b1 > 0 ){  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			top->Fill( -tcos1, weight ) ;
			//top->Fill( tcos2, weight ) ;
			bottom->Fill( -bcos1, weight ) ;
			bottom->Fill( bcos2, weight ) ;
		}else if( b2 > 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			top->Fill( tcos1, weight ) ;
			//top->Fill( -tcos2, weight ) ;
			bottom->Fill( bcos1, weight ) ;
			bottom->Fill( -bcos2, weight ) ;
		}else{  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			top->Fill( -tcos1, weight ) ;
			//top->Fill( tcos2, weight ) ;
			bottom->Fill( -bcos1, weight ) ;
			bottom->Fill( bcos2, weight ) ;
		}

}


#ifdef __CINT__
int polar_angle(){
#else
int main( int argc, char **argv ){
	TApplication app( "app", &argc, argv ) ;
#endif


	//string rootfilename = "rootfile/hadronic_eLpR_yyuyyc_00001.root" ;
	//string rootfilename = "rootfile/hadronic_eLpR_yyuyyc.root" ;
	//string rootfilename = "rootfile/hadronic_eRpL.root" ;
	//string rootfilename = "rootfile/hadronic_eLpR.root" ;
	//string rootfilename = "aftervertexrestore_hadronic_eLpR_yyuyyc_00001.root" ;
	//string rootfilename = "aftervertexrestore_hadronic_eRpL_yyuyyc_00001.root" ;
	//string rootfilename = "rootfile/aftervertexrestoretest.root" ;
	string rootfilename = "/home/ilc/yokugawa/TTBarAnalysis/rootfile/semi-leptonic.root" ;
	string treename = "Stats" ;

	TFile file( rootfilename.c_str() ) ;
	if( file.IsZombie() ){
		cout << "cannot open the file" << endl ;
		return 1 ;
	}

	gStyle->SetOptStat(0) ;

	int Top1bcharge, Top2bcharge ;
    int Top1TotalKaonCharge, Top2TotalKaonCharge ;
	float Top1costheta, Top1bcostheta, Top2costheta, Top2bcostheta ;
	float MCTopcostheta, MCTopBarcostheta ;
	float Top1bmomentum, Top2bmomentum;
	float qMCBcostheta[2] ;

	TCut btag = " ( Top1btag > 0.80 ) && ( Top2btag > 0.30 ) " ;
	TCut chi2_1 = " chiTopMass1 + chiTopE1 + chiPbstar1 < 30 " ;
	TCut chi2_2 = " chiTopMass2 + chiTopE2 + chiPbstar2 < 30 " ;
	//TCut chi2 = chi2_1 + chi2_2 ;
	TCut chi2 = chi2_1 ;
	TCut kinematic = " ( Top1mass > 140 ) && ( Top1mass < 210 ) && ( Top2mass > 140 ) && ( Top2mass < 210 ) " ;
	//TCut samecharge = " Top1bcharge * Top2bcharge > 0 " ;
	//TCut both0charge = " ( Top1bcharge == 0 ) || ( Top2bcharge == 0 ) " ;

	TTree* tree = (TTree*)file.Get( treename.c_str() ) ;
	int eventnum = tree->GetEntries() ;
	cout << "eventnum            = " << eventnum << " (100%)" << endl ;
	int afterbtag = tree->GetEntries( btag ) ;
	cout << "after b-tag cut     = " << afterbtag << " (" << (float)100*afterbtag/eventnum << "%)" << endl ;
	//int afterkinematic = tree->GetEntries( btag && kinematic ) ;
	//cout << "atfer kinematic cut = " << afterkinematic << " (" << (float)100*afterkinematic/eventnum << "%)" << endl;
	//int afterchi2 = tree->GetEntries( btag && kinematic && chi2 ) ;
	int afterchi2 = tree->GetEntries( btag && chi2 ) ;
	cout << "after chi2 cut      = " << afterchi2 << " (" << (float)100*afterchi2/eventnum << "%)" << endl;

	//int samesignnum = tree->GetEntries( btag && kinematic && chi2 && samecharge ) ;
	//int both0num = tree->GetEntries( btag && kinematic && chi2 && both0charge ) ;
	//int usednum = afterchi2 - samesignnum - both0num ;
	int usednum = afterchi2;
	cout << endl ;
	cout << "used number = " << usednum << endl ;
	//cout << "same charge sign number = " << samesignnum << endl ;
	//cout << "both charge 0 number = " << both0num << endl ;

	TEventList *elist = new TEventList( "elist", "Reconstructed Event List" ) ;
	//tree->Draw( ">>elist", btag && kinematic && chi2) ;
	tree->Draw( ">>elist", btag && chi2) ;
	//usednum = elist->GetN() ;
	//cout << "used number = " << usednum << endl ;
	
	tree->SetBranchAddress( "Top1TotalKaonCharge", &Top1TotalKaonCharge ) ;
	tree->SetBranchAddress( "Top2TotalKaonCharge", &Top2TotalKaonCharge ) ;
	tree->SetBranchAddress( "Top1bcharge", &Top1bcharge ) ;
	tree->SetBranchAddress( "Top1costheta", &Top1costheta ) ;
	tree->SetBranchAddress( "Top1bcostheta", &Top1bcostheta ) ;
	tree->SetBranchAddress( "Top2bcharge", &Top2bcharge ) ;
	tree->SetBranchAddress( "Top2costheta", &Top2costheta ) ;
	tree->SetBranchAddress( "Top2bcostheta", &Top2bcostheta ) ;
	tree->SetBranchAddress( "Top1bmomentum", &Top1bmomentum ) ;
	tree->SetBranchAddress( "Top2bmomentum", &Top2bmomentum ) ;

	TCanvas* c1 = new TCanvas( "c1", "c1", 1280, 480 ) ;
	TH1F* htop = new TH1F( "htop", "top polar angle (reconstructed)", 40, -1, 1 ) ;
	htop->SetMinimum(0) ;
	htop->SetLineColor(2) ;
	TH1F* hbottom = new TH1F( "hbottom", "bottom polar angle (reconstructed)", 40, -1, 1 ) ;
	hbottom->SetMinimum(0) ;
	hbottom->SetLineColor(2) ;
	TH1F* htopgen = new TH1F( "htopgen", "top polar angle (generated)", 40, -1, 1 ) ;
	htopgen->SetMinimum(0) ;
	htopgen->SetMaximum(0.11) ;
	htopgen->SetLineColor(4) ;
	TH1F* hbottomgen = new TH1F( "hbottomgen", "bottom polar angle (generated)", 40, -1, 1 ) ;
	hbottomgen->SetMinimum(0) ;
	hbottomgen->SetMaximum(0.11) ;
	hbottomgen->SetLineColor(4) ;

	double fillweight = 1 ;
    
	// in-loop varibale declaration
	bool isbsign = false, iskaonsign = false;
	bool isb1k2sign = false;
	bool isb2k1sign = false;
	bool isbMomentum = false;
    
	for( int i=0 ; i<usednum ; i++ ){

		tree->GetEntry( elist->GetEntry(i) ) ;

		isbsign = false ;
		iskaonsign = false;
		isb1k2sign = false ;
		isb2k1sign = false ;
		isbMomentum = false ;

		if( Top1bcharge*Top2bcharge < 0 ) isbsign = true;
		if( Top1TotalKaonCharge*Top2TotalKaonCharge < 0 ) iskaonsign = true;
		if( Top1bcharge*Top2TotalKaonCharge < 0 ) isb1k2sign = true;
		if( Top2bcharge*Top1TotalKaonCharge < 0 ) isb2k1sign = true;
		if( Top1bmomentum > 30 && Top2bmomentum > 30 ) isbMomentum = true;
	
		if( isbsign ){
				FillandAdd(Top1bcharge,Top2bcharge,Top1costheta,Top2costheta,Top1bcostheta,Top2bcostheta,fillweight,htop,hbottom);
		}else if( iskaonsign ){
				FillandAdd(Top1TotalKaonCharge,Top2TotalKaonCharge,Top1costheta,Top2costheta,Top1bcostheta,Top2bcostheta,fillweight,htop,hbottom);
		}else if( isb1k2sign ){
				FillandAdd(Top1bcharge,Top2TotalKaonCharge,Top1costheta,Top2costheta,Top1bcostheta,Top2bcostheta,fillweight,htop,hbottom);
		}else if( isb2k1sign ){
				FillandAdd(Top1TotalKaonCharge,Top2bcharge,Top1costheta,Top2costheta,Top1bcostheta,Top2bcostheta,fillweight,htop,hbottom);
		}


	    /*
		if( Top1bcharge < 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			htop->Fill( Top1costheta, fillweight ) ;
			htop->Fill( -Top2costheta, fillweight ) ;
			hbottom->Fill( Top1bcostheta, fillweight ) ;
			hbottom->Fill( -Top2bcostheta, fillweight ) ;
		}else if( Top1bcharge > 0 ){  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			htop->Fill( -Top1costheta, fillweight ) ;
			htop->Fill( Top2costheta, fillweight ) ;
			hbottom->Fill( -Top1bcostheta, fillweight ) ;
			hbottom->Fill( Top2bcostheta, fillweight ) ;
		}else if( Top2bcharge > 0 ){  //Top1 = top, Top2 = antitop, Top1b = bottom, Top2b = antibottom
			htop->Fill( Top1costheta, fillweight ) ;
			htop->Fill( -Top2costheta, fillweight ) ;
			hbottom->Fill( Top1bcostheta, fillweight ) ;
			hbottom->Fill( -Top2bcostheta, fillweight ) ;
		}else{  //Top1 = antitop, Top2 = top, Top1b = antibottom, Top2b = bottom
			htop->Fill( -Top1costheta, fillweight ) ;
			htop->Fill( Top2costheta, fillweight ) ;
			hbottom->Fill( -Top1bcostheta, fillweight ) ;
			hbottom->Fill( Top2bcostheta, fillweight ) ;
		}*/

	}

	tree->SetBranchStatus( "Top*", 0 ) ;
	tree->SetBranchAddress( "MCTopcostheta", &MCTopcostheta ) ;
	tree->SetBranchAddress( "MCTopBarcostheta", &MCTopBarcostheta ) ;
	tree->SetBranchAddress( "qMCBcostheta", &qMCBcostheta ) ;

	TEventList* elistgen = new TEventList( "elistgen", "Generated Event List" ) ;
	tree->Draw( ">>elistgen", " ( MCTopcostheta != -2 ) && ( MCTopBarcostheta != -2 ) && ( qMCBcostheta[0] != -2 ) && ( qMCBcostheta[1] != -2 ) " ) ;
	int usedgennum = elistgen->GetN() ;
	cout << endl ;
	cout << "used genarated number = " << usedgennum << endl ;

	fillweight = (double)1/(2*usedgennum) ;

	for( int i=0 ; i<usedgennum ; i++ ){

		tree->GetEntry( elistgen->GetEntry(i) ) ;

		htopgen->Fill( MCTopcostheta, fillweight ) ;
		htopgen->Fill( -MCTopBarcostheta, fillweight ) ;
		hbottomgen->Fill( qMCBcostheta[0], fillweight ) ;
		hbottomgen->Fill( -qMCBcostheta[1], fillweight ) ;

	}

	double norm = htop->GetEntries() ;
	htop->Scale(1/norm);

	c1->Divide(2,1) ;
	c1->cd(1) ;
	htopgen->SetTitle( "top polar angle (Rec:Red Gen:Blue)" ) ;
	htopgen->Draw() ;
	htop->Draw( "SAME" ) ;

	norm = hbottom->GetEntries()/2 ;
	hbottom->Scale(1/norm);

	c1->cd(2) ;
	hbottomgen->SetTitle( "bottom polar angle (Rec:Red Gen:Blue)" ) ;
	hbottomgen->Draw() ;
	hbottom->Draw( "SAME" ) ;

	c1->WaitPrimitive() ;

	c1->Print( "picture/newpicture.png" ) ;
	c1->Close() ;

	file.Close() ;

	return 0;

}




