#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include "../../style/Style.C"
#include "../../style/Labels.C"

//#include "rootlogon.C"
//#include "bilo_sty.C"
#define MAXV 8

using namespace std;

double cosb_topframe(float Bfourvect[4], float Topfourvect[4]){

  /*float bboosted[3];
    float betat[3];
    for(unsigned int i = 0; i<3; i++) betat[i] = Topfourvect[i]/Topfourvect[3];
    float betat2 = betat[0]*betat[0]+betat[1]*betat[1]+betat[2]*betat[2];
    float gammat = 1/sqrt(1-betat2);
    for(unsigned int i = 0; i<3; i++) bboosted[i] = Bfourvect[i] + betat[i] * (- gammat*Bfourvect[3] + (gammat-1)*(betat[0]/betat2)*Bfourvect[0] + (gammat-1)*(betat[1]/betat2)*Bfourvect[1] + (gammat-1)*(betat[2]/betat2)*Bfourvect[2]);
    double dot = 0;
    for(unsigned int i = 0; i<3; i++) dot += bboosted[i]*Topfourvect[i];
    double mag2b = 0;
    for(unsigned int i = 0; i<3; i++) mag2b += bboosted[i]*bboosted[i];
    double mag2t = 0;
    for(unsigned int i = 0; i<3; i++) mag2t += Topfourvect[i]*Topfourvect[i];
    double bboostedcos = dot/sqrt(mag2b*mag2t);*/

  TLorentzVector bboosted(Bfourvect);
  TLorentzVector cm(0.0,0.0,0.0,500.0);
  TLorentzVector top(Topfourvect);
  TVector3 betat = top.BoostVector();
  bboosted.Boost( - betat );
  cm.Boost( - betat );
  double bboostedcos = cos( (-cm).Angle(bboosted.Vect()) );

  return bboostedcos;

}

double theta_topframe(float Bfourvect[4], float Topfourvect[4]){


  TLorentzVector bboosted(Bfourvect);
  TLorentzVector cm(0.0,0.0,0.0,500.0);
  TLorentzVector top(Topfourvect);
  TVector3 betat = top.BoostVector();
  bboosted.Boost( - betat );
  cm.Boost( - betat );
  double theta = (-cm).Angle(bboosted.Vect());

  return theta;

}

double computeGamma(float fourvect[4]){


  TLorentzVector v(fourvect);

  return v.Gamma();

}

double computeMag(float fourvect[4]){

  TLorentzVector v(fourvect);
  return v.Vect().Mag();

}

void asymmetry_sl_btopframe()
{

  // initialize variables
  int styl = 0;
  int cx   = 500;
  double Legx1 = 0.20;
  double Legx2 = 0.6;

  int token_l5=0;
  int token_s5=0;

  // set plot style
  SetQQbarStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);  
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetMarkerSize(1);
  gStyle->SetTitleX(0.2); 
  gStyle->SetTitleY(0.9); 

	
  FileSelector fs;
  std::vector<FileSelector> rootfiles;
  std::ifstream in( "/home/ilc/vlohezic/working/VertexMacro/semi_leptonic/input/record.txt" );

  while( fs.input(in) ){
    rootfiles.push_back(fs);
  }

  int nrootfiles = 0;
  nrootfiles = rootfiles.size();

  std::cout << "Choose a file from below:" << std::endl;
  for( int i=0; i < nrootfiles; i++){
    std::cout << i << ": " << rootfiles[i].info() << endl;
  }

  std::cout << "Enter code for l5: ";
  std::cin >> token_l5;
  std::cout << "Enter code for s5: ";
  std::cin >> token_s5;

  std::cout << std::endl;

  TCanvas * c1 = new TCanvas("c1", "Data-MC",0,0,cx,500);
  TCanvas * c2 = new TCanvas("c2", "Data-MC",0,0,cx,500);
  TCanvas * c3 = new TCanvas("c3", "Data-MC",0,0,cx,500);
  

  ////////////// Cuts //////////////

  float MCTopfourvect[2][4];
  float MCWfourvect[2][4];
  float qMCBfourvect[2][4];
  float Topfourvect[2][4];
  float Wfourvect[2][4];
  float qBfourvect[2][4];
  float Thrust;
  float hadMass;
  float Top1mass;
  float W1mass;
  float Top1bmomentum;
  float Top2bmomentum;
  float Top1gamma;
  float Top2gamma;
  int methodTaken[12];


  ////////////// Detector Model LARGE //////////////

  std::string filename_l5 = rootfiles[token_l5].filename();
  cout << "Processing (l5) : " << filename_l5 << " ..." << endl;

  TFile * file_l5 = TFile::Open(filename_l5.c_str());

  int bin_e = 20;
  int max_e = 1;

  TH1F * cosReco_l5 = new TH1F("cosReco_l5", "E(Ntracks)", bin_e,-1.0,max_e);
  cosReco_l5->Sumw2();
  TH1F * cosGen_l5 = new TH1F("cosGen_l5", ";cos#theta_{b};Entries", bin_e,-1.0,max_e);
  cosGen_l5->Sumw2();
  TH1F * cosRecobar_l5 = new TH1F("cosRecobar_l5", "E(Ntracks)", bin_e,-1.0,max_e);
  cosRecobar_l5->Sumw2();
  TH1F * cosGenbar_l5 = new TH1F("cosGenbar_l5", ";cos#theta_{#bar{b}};Entries", bin_e,-1.0,max_e);
  cosGenbar_l5->Sumw2();

  TGaxis::SetMaxDigits(3);

  TTree * normaltree_l5 = (TTree*) file_l5->Get( "Stats" ) ;
  TTree * GenTree_l5 = (TTree*) file_l5->Get( "GenTree" ) ;

  int neventGen_l5 = GenTree_l5->GetEntries();
  int neventReco_l5 = normaltree_l5->GetEntries();

  cout << "l5Reco Entry = " << neventReco_l5 << endl;
  cout << "l5Gen Entry = " << neventGen_l5 << endl;

  int ngen_l5 = 0;
  int ngenpos_l5 = 0;
  int ngenneg_l5 = 0;
  int nreco_l5 = 0;
  int nrecopos_l5 = 0;
  int nreconeg_l5 = 0;
  int ngenbar_l5 = 0;
  int nrecobar_l5 = 0;

  GenTree_l5->SetBranchAddress("MCTopfourvect",&MCTopfourvect);
  GenTree_l5->SetBranchAddress("MCWfourvect",&MCWfourvect);
  GenTree_l5->SetBranchAddress("qMCBfourvect",&qMCBfourvect);
	
  for(int ev = 0; ev<neventGen_l5; ev++){
    GenTree_l5->GetEntry(ev);

    double bboostedcos = cosb_topframe(qMCBfourvect[0],MCTopfourvect[0]);

    double Wtheta = theta_topframe(MCWfourvect[0],MCTopfourvect[0]);
    double cosWb = cos( Wtheta + acos(bboostedcos) );

    float betat[3];
    for(unsigned int i = 0; i<3; i++) betat[i] = Topfourvect[0][i]/Topfourvect[0][3];
    float betat2 = betat[0]*betat[0]+betat[1]*betat[1]+betat[2]*betat[2];
    float gammat = 1/sqrt(1-betat2);

    double bmomentum = sqrt( qMCBfourvect[0][0]*qMCBfourvect[0][0]+qMCBfourvect[0][1]*qMCBfourvect[0][1]+qMCBfourvect[0][2]*qMCBfourvect[0][2] );
    double bbarmomentum = sqrt( qMCBfourvect[1][0]*qMCBfourvect[1][0]+qMCBfourvect[1][1]*qMCBfourvect[1][1]+qMCBfourvect[1][2]*qMCBfourvect[1][2] );

    if(bmomentum>15 && bbarmomentum>15 && bboostedcos>-1 && bboostedcos<1){

      (bboostedcos > 0)? ngenpos_l5++ : ngenneg_l5++;

      ngen_l5++;
      cosGen_l5->Fill(bboostedcos);
    }


    double bbarboostedcos = cosb_topframe(qMCBfourvect[1],MCTopfourvect[1]);

    float betatbar[3];
    for(unsigned int i = 0; i<3; i++) betatbar[i] = MCTopfourvect[1][i]/MCTopfourvect[1][3];
    float betatbar2 = betatbar[0]*betatbar[0]+betatbar[1]*betatbar[1]+betatbar[2]*betatbar[2];
    float gammatbar = 1/sqrt(1-betatbar2);

    if(bbarboostedcos>-1 && bbarboostedcos<1){
      ngenbar_l5++;
      cosGenbar_l5->Fill(bbarboostedcos);

    }
  }

  normaltree_l5->SetBranchAddress("MCTopfourvect",&MCTopfourvect);
  normaltree_l5->SetBranchAddress("qMCBfourvect",&qMCBfourvect);
  normaltree_l5->SetBranchAddress("Topfourvect",&Topfourvect);
  normaltree_l5->SetBranchAddress("Wfourvect",&Wfourvect);
  normaltree_l5->SetBranchAddress("qBfourvect",&qBfourvect);
  normaltree_l5->SetBranchAddress("Thrust",&Thrust);
  normaltree_l5->SetBranchAddress("hadMass",&hadMass);
  normaltree_l5->SetBranchAddress("Top1mass",&Top1mass);
  normaltree_l5->SetBranchAddress("W1mass",&W1mass);
  normaltree_l5->SetBranchAddress("Top1bmomentum",&Top1bmomentum);
  normaltree_l5->SetBranchAddress("Top2bmomentum",&Top2bmomentum);
  normaltree_l5->SetBranchAddress("Top1gamma",&Top1gamma);
  normaltree_l5->SetBranchAddress("Top2gamma",&Top2gamma);
  normaltree_l5->SetBranchAddress("methodTaken",&methodTaken);
    
  for(int ev = 0; ev<neventReco_l5; ev++){
    normaltree_l5->GetEntry(ev);

    bool methodok = false;
    for(unsigned int i = 0; i<12; i++){
      if(methodTaken[i] == 5 || methodTaken[i] == 6 || methodTaken[i] == 7) methodok = true;
    }

    double bboostedcos = cosb_topframe(qBfourvect[0],Topfourvect[0]);
    
    double Wtheta = theta_topframe(Wfourvect[0],MCTopfourvect[0]);
    double cosWb = cos( Wtheta + acos(bboostedcos) );
    
    if( Thrust < 0.9 && hadMass > 180 && hadMass < 420 && Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50 && computeMag(qBfourvect[0])>15 && computeMag(qMCBfourvect[0])>15 && Top1bmomentum > 15 && Top2bmomentum > 15 && (Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2 && methodok && bboostedcos>-1 && bboostedcos<1 /*&& cosWb<-0.99 */){
  
      (bboostedcos > 0)? nrecopos_l5++ : nreconeg_l5++;
	  
      nreco_l5++;
      cosReco_l5->Fill(bboostedcos);

    }

    double bbarboostedcos = cosb_topframe(qBfourvect[1],Topfourvect[1]);
  
    if( Thrust < 0.9 && hadMass > 180 && hadMass < 420 && Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50 && Top1bmomentum > 15 && Top2bmomentum > 15 && (Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2 && methodok && bbarboostedcos>-1 && bbarboostedcos<1 ){
      nrecobar_l5++;
      cosRecobar_l5->Fill(bbarboostedcos);

    }

  }

  TH1F * cosRecosum_l5 = new TH1F("cosRecosum_l5", "E(Ntracks)", bin_e,-1.0,max_e);
  TH1F * cosGensum_l5 = new TH1F("cosGensum_l5", ";cos#theta_{b#bar{b}};Entries", bin_e,-1.0,max_e);
  cosRecosum_l5->Add(cosReco_l5,cosRecobar_l5);
  cosGensum_l5->Add(cosGen_l5,cosGenbar_l5);

  ////////////// Detector Model SMALL //////////////

  std::string filename_s5 = rootfiles[token_s5].filename();
  cout << "Processing (s5) : " << filename_s5 << " ..." << endl;

  TFile * file_s5 = TFile::Open(filename_s5.c_str());

  TH1F * cosReco_s5 = new TH1F("cosReco_s5", "E(Ntracks)", bin_e,-1.0,max_e);
  cosReco_s5->Sumw2();
  TH1F * cosGen_s5 = new TH1F("cosGen_s5", ";cos#theta_{b};Entries", bin_e,-1.0,max_e);
  cosGen_s5->Sumw2();
  TH1F * cosRecobar_s5 = new TH1F("cosRecobar_s5", "E(Ntracks)", bin_e,-1.0,max_e);
  cosRecobar_s5->Sumw2();
  TH1F * cosGenbar_s5 = new TH1F("cosGenbar_s5", ";cos#theta_{#bar{b}};Entries", bin_e,-1.0,max_e);
  cosGenbar_s5->Sumw2();

  TTree * normaltree_s5 = (TTree*) file_s5->Get( "Stats" ) ;
  TTree * GenTree_s5 = (TTree*) file_s5->Get( "GenTree" ) ;

  int neventGen_s5 = GenTree_s5->GetEntries();
  int neventReco_s5 = normaltree_s5->GetEntries();

  cout << "s5Reco Entry = " << neventReco_s5 << endl;
  cout << "s5Gen Entry = " << neventGen_s5 << endl;

  int ngen_s5 = 0;
  int nreco_s5 = 0;
  int ngenbar_s5 = 0;
  int nrecobar_s5 = 0;

  GenTree_s5->SetBranchAddress("MCTopfourvect",&MCTopfourvect);
  GenTree_s5->SetBranchAddress("MCWfourvect",&MCWfourvect);
  GenTree_s5->SetBranchAddress("qMCBfourvect",&qMCBfourvect);
	
  for(int ev = 0; ev<neventGen_s5; ev++){
    GenTree_s5->GetEntry(ev);

    double bboostedcos = cosb_topframe(qMCBfourvect[0],MCTopfourvect[0]);

    if(bboostedcos>-1 && bboostedcos<1){
      ngen_s5++;
      cosGen_s5->Fill(bboostedcos);
    }

    double bbarboostedcos = cosb_topframe(qMCBfourvect[1],MCTopfourvect[1]);

    if(bbarboostedcos>-1 && bbarboostedcos<1){
      ngenbar_s5++;
      cosGenbar_s5->Fill(bbarboostedcos);
    }
  }

  normaltree_s5->SetBranchAddress("MCTopfourvect",&MCTopfourvect);
  normaltree_s5->SetBranchAddress("qMCBfourvect",&qMCBfourvect);
  normaltree_s5->SetBranchAddress("Topfourvect",&Topfourvect);
  normaltree_s5->SetBranchAddress("Wfourvect",&Wfourvect);
  normaltree_s5->SetBranchAddress("qBfourvect",&qBfourvect);
  normaltree_s5->SetBranchAddress("Thrust",&Thrust);
  normaltree_s5->SetBranchAddress("hadMass",&hadMass);
  normaltree_s5->SetBranchAddress("Top1mass",&Top1mass);
  normaltree_s5->SetBranchAddress("W1mass",&W1mass);
  normaltree_s5->SetBranchAddress("Top1bmomentum",&Top1bmomentum);
  normaltree_s5->SetBranchAddress("Top2bmomentum",&Top2bmomentum);
  normaltree_s5->SetBranchAddress("Top1gamma",&Top1gamma);
  normaltree_s5->SetBranchAddress("Top2gamma",&Top2gamma);
  normaltree_s5->SetBranchAddress("methodTaken",&methodTaken);
    
  for(int ev = 0; ev<neventReco_s5; ev++){
    normaltree_s5->GetEntry(ev);

    bool methodok = false;
    for(unsigned int i = 0; i<12; i++){
      if(methodTaken[i] == 5 || methodTaken[i] == 6 || methodTaken[i] == 7) methodok = true;
    }

    double bboostedcos = cosb_topframe(qBfourvect[0],Topfourvect[0]);

    double Wtheta = theta_topframe(Wfourvect[0],MCTopfourvect[0]);
    double cosWb = cos( Wtheta + acos(bboostedcos) );

    if( Thrust < 0.9 && hadMass > 180 && hadMass < 420 && Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50 && /*Top1bmomentum > 15 && Top2bmomentum > 15 &&*/ (Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2 && methodok && bboostedcos>-1 && bboostedcos<1 && cosWb<-0.99){
      nreco_s5++;
      cosReco_s5->Fill(bboostedcos);
    }

    double bbarboostedcos = cosb_topframe(qBfourvect[1],Topfourvect[1]);

    if( Thrust < 0.9 && hadMass > 180 && hadMass < 420 && Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50 && Top1bmomentum > 15 && Top2bmomentum > 15 && (Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2 && methodok && bbarboostedcos>-1 && bbarboostedcos<1){
      nrecobar_s5++;
      cosRecobar_s5->Fill(bbarboostedcos);
    }
  }

  TH1F * cosRecosum_s5 = new TH1F("cosRecosum_s5", "E(Ntracks)", bin_e,-1.0,max_e);
  cosRecosum_s5->Add(cosReco_s5,cosRecobar_s5);

  /////////////// Scaling //////////////////
  
  /*double intCosReco_l5 = cosReco_l5->Integral(2,19);
    double intCosGen  = cosGen_l5->Integral(2,19);
    cosGen_l5->Scale(intCosReco_l5 / intCosGen);
    double intCosReco_s5 = cosReco_s5->Integral(2,19);
    cosReco_s5->Scale(intCosReco_l5 / intCosReco_s5);

    double intCosRecobar_l5 = cosRecobar_l5->Integral(2,19);
    double intCosGenbar  = cosGenbar_l5->Integral(2,19);
    cosGenbar_l5->Scale(intCosRecobar_l5 / intCosGenbar);
    double intCosRecobar_s5 = cosRecobar_s5->Integral(2,19);
    cosRecobar_s5->Scale(intCosRecobar_l5 / intCosRecobar_s5);

    double intCosRecosum_l5 = cosRecosum_l5->Integral(2,19);
    double intCosGensum  = cosGensum_l5->Integral(2,19);
    cosGensum_l5->Scale(intCosRecosum_l5 / intCosGensum);
    double intCosRecosum_s5 = cosRecosum_s5->Integral(2,19);
    cosRecosum_s5->Scale(intCosRecosum_l5 / intCosRecosum_s5);*/

  ////////////// Style Setting //////////////

  /////For b
  cosReco_l5->SetLineWidth(3);
  cosReco_l5->SetLineColor(kBlue);
  cosReco_l5->SetMarkerColor(kBlue);
  cosReco_l5->SetMarkerStyle(21);

  cosReco_s5->SetLineWidth(3);
  cosReco_s5->SetLineColor(kRed);
  cosReco_s5->SetMarkerColor(kRed);
  cosReco_s5->SetMarkerStyle(22);
  cosReco_s5->SetLineStyle(2);

  cosGen_l5->SetLineWidth(3);
  //cosGen_l5->SetLineStyle(2);
  //cosGen_l5->SetLineColor(kGreen+1);
  //cosGen_l5->SetFillColor(kGreen+1);
  cosGen_l5->SetLineColor(kGray+1);
  cosGen_l5->SetFillColor(kGray+1);

  cosGen_l5->SetFillStyle(3004);
  cosGen_l5->SetStats(0);
  cosGen_l5->SetMinimum(0);
  //cosGen_l5->SetMaximum(7000);

  //cosGen_l5->SetTitle("e_{L}^{+}e_{R}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGen_l5->SetTitle("e_{R}^{+}e_{L}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGen_l5->GetXaxis()->SetTitleOffset(1.1);
  cosGen_l5->GetXaxis()->SetTitleFont(42);
  cosGen_l5->GetXaxis()->SetTitleSize(0.05);
  cosGen_l5->GetXaxis()->SetLabelSize(0.05);
  cosGen_l5->GetXaxis()->SetLabelOffset(0.015);

  cosGen_l5->GetYaxis()->SetTitle("entries / 0.1 rad");
  cosGen_l5->GetYaxis()->SetTitleOffset(1.4);
  cosGen_l5->GetYaxis()->SetTitleFont(42);
  cosGen_l5->GetYaxis()->SetTitleSize(0.05);
  cosGen_l5->GetYaxis()->SetLabelSize(0.05);
  cosGen_l5->GetYaxis()->SetLabelOffset(0.015);

  /////For bbar
  cosRecobar_l5->SetLineWidth(3);
  cosRecobar_l5->SetLineColor(kBlue);
  cosRecobar_l5->SetMarkerColor(kBlue);
  cosRecobar_l5->SetMarkerStyle(21);

  cosRecobar_s5->SetLineWidth(3);
  cosRecobar_s5->SetLineColor(kRed);
  cosRecobar_s5->SetMarkerColor(kRed);
  cosRecobar_s5->SetMarkerStyle(22);
  cosRecobar_s5->SetLineStyle(2);

  cosGenbar_l5->SetLineWidth(3);
  //cosGenbar_l5->SetLineStyle(2);
  //cosGenbar_l5->SetLineColor(kGreen+1);
  //cosGenbar_l5->SetFillColor(kGreen+1);
  cosGenbar_l5->SetLineColor(kGray+1);
  cosGenbar_l5->SetFillColor(kGray+1);

  cosGenbar_l5->SetFillStyle(3004);
  cosGenbar_l5->SetStats(0);
  cosGenbar_l5->SetMinimum(0);
  //cosGenbar_l5->SetMaximum(7000);

  //cosGenbar_l5->SetTitle("e_{L}^{+}e_{R}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGenbar_l5->SetTitle("e_{R}^{+}e_{L}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGenbar_l5->GetXaxis()->SetTitleOffset(1.1);
  cosGenbar_l5->GetXaxis()->SetTitleFont(42);
  cosGenbar_l5->GetXaxis()->SetTitleSize(0.05);
  cosGenbar_l5->GetXaxis()->SetLabelSize(0.05);
  cosGenbar_l5->GetXaxis()->SetLabelOffset(0.015);

  cosGenbar_l5->GetYaxis()->SetTitle("entries / 0.1 rad");
  cosGenbar_l5->GetYaxis()->SetTitleOffset(1.4);
  cosGenbar_l5->GetYaxis()->SetTitleFont(42);
  cosGenbar_l5->GetYaxis()->SetTitleSize(0.05);
  cosGenbar_l5->GetYaxis()->SetLabelSize(0.05);
  cosGenbar_l5->GetYaxis()->SetLabelOffset(0.015);

  /////For sum
  cosRecosum_l5->SetLineWidth(3);
  cosRecosum_l5->SetLineColor(kBlue);
  cosRecosum_l5->SetMarkerColor(kBlue);
  cosRecosum_l5->SetMarkerStyle(21);

  cosRecosum_s5->SetLineWidth(3);
  cosRecosum_s5->SetLineColor(kRed);
  cosRecosum_s5->SetMarkerColor(kRed);
  cosRecosum_s5->SetMarkerStyle(22);
  cosRecosum_s5->SetLineStyle(2);

  cosGensum_l5->SetLineWidth(3);
  //cosGensum_l5->SetLineStyle(2);
  //cosGensum_l5->SetLineColor(kGreen+1);
  //cosGensum_l5->SetFillColor(kGreen+1);
  cosGensum_l5->SetLineColor(kGray+1);
  cosGensum_l5->SetFillColor(kGray+1);

  cosGensum_l5->SetFillStyle(3004);
  cosGensum_l5->SetStats(0);
  cosGensum_l5->SetMinimum(0);
  //cosGensum_l5->SetMaximum(7000);

  //cosGensum_l5->SetTitle("e_{L}^{+}e_{R}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGensum_l5->SetTitle("e_{R}^{+}e_{L}^{-}#rightarrow t#bar{t} @ 500GeV");
  cosGensum_l5->GetXaxis()->SetTitleOffset(1.1);
  cosGensum_l5->GetXaxis()->SetTitleFont(42);
  cosGensum_l5->GetXaxis()->SetTitleSize(0.05);
  cosGensum_l5->GetXaxis()->SetLabelSize(0.05);
  cosGensum_l5->GetXaxis()->SetLabelOffset(0.015);

  cosGensum_l5->GetYaxis()->SetTitle("entries / 0.1 rad");
  cosGensum_l5->GetYaxis()->SetTitleOffset(1.4);
  cosGensum_l5->GetYaxis()->SetTitleFont(42);
  cosGensum_l5->GetYaxis()->SetTitleSize(0.05);
  cosGensum_l5->GetYaxis()->SetLabelSize(0.05);
  cosGensum_l5->GetYaxis()->SetLabelOffset(0.015);


  ////////////// Fitting //////////////

  /////For b
  TF1 * fgen = new TF1("fgen","pol2",-1,1);
  TF1 * freco_l5 = new TF1("freco_l5","pol2",-0.9,0.9);
  TF1 * freco_s5 = new TF1("freco_s5","pol2",-0.9,0.9);
  //fgen->SetLineColor(kGreen);
  fgen->SetLineColor(kGray);
  fgen->SetLineStyle(3);
  freco_l5->SetLineStyle(3);
  freco_l5->SetLineColor(kBlue+1);
  freco_s5->SetLineStyle(3);
  freco_s5->SetLineColor(kRed+1);
	
  cosGen_l5->Fit("fgen","Q");
  cosReco_l5->Fit("freco_l5", "QR");
  cosReco_s5->Fit("freco_s5", "QR");

  /////For bbar
  TF1 * fgenbar = new TF1("fgenbar","pol2",-1,1);
  TF1 * frecobar_l5 = new TF1("frecobar_l5","pol2",-0.9,0.9);
  TF1 * frecobar_s5 = new TF1("frecobar_s5","pol2",-0.9,0.9);
  //fgenbar->SetLineColor(kGreen);
  fgenbar->SetLineColor(kGray);
  fgenbar->SetLineStyle(3);
  frecobar_l5->SetLineStyle(3);
  frecobar_l5->SetLineColor(kBlue+1);
  frecobar_s5->SetLineStyle(3);
  frecobar_s5->SetLineColor(kRed+1);
	
  cosGenbar_l5->Fit("fgenbar","Q");
  cosRecobar_l5->Fit("frecobar_l5", "QR");
  cosRecobar_s5->Fit("frecobar_s5", "QR");

  /////For sum
  TF1 * fgensum = new TF1("fgensum","pol2",-1,1);
  TF1 * frecosum_l5 = new TF1("frecosum_l5","pol2",-0.9,0.9);
  TF1 * frecosum_s5 = new TF1("frecosum_s5","pol2",-0.9,0.9);
  //fgensum->SetLineColor(kGreen);
  fgensum->SetLineColor(kGray);
  fgensum->SetLineStyle(3);
  frecosum_l5->SetLineStyle(3);
  frecosum_l5->SetLineColor(kBlue+1);
  frecosum_s5->SetLineStyle(3);
  frecosum_s5->SetLineColor(kRed+1);
	
  cosGensum_l5->Fit("fgensum","Q");
  cosRecosum_l5->Fit("frecosum_l5", "QR");
  cosRecosum_s5->Fit("frecosum_s5", "QR");


  ////////////// Drawing //////////////

  /////For b
  c1->cd();
  cosGen_l5->Draw("he");
  fgen->Draw("same");
  cosReco_l5->Draw("samee");
  cosReco_s5->Draw("samee");

  //TLegend *leg = new TLegend(0.2,0.65,0.6,0.75); //set here your x_0,y_0, x_1,y_1 options
  TLegend *leg = new TLegend(0.2,0.65,0.55,0.8);
  leg->SetTextFont(42);
  leg->AddEntry(cosGen_l5,"Parton level","l");
  leg->AddEntry(cosReco_l5,"IDR-L","lep");
  leg->AddEntry(cosReco_s5,"IDR-S","lep");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->Draw();

  QQBARLabel(0.8,0.2,NULL,1);

  c1->Update();

  /////For bbar
  c2->cd();
  cosGenbar_l5->Draw("he");
  fgenbar->Draw("same");
  cosRecobar_l5->Draw("samee");
  cosRecobar_s5->Draw("samee");

  leg->Draw();

  QQBARLabel(0.8,0.2,NULL,1);

  c2->Update();

  /////For sum
  c3->cd();
  cosGensum_l5->Draw("he");
  fgensum->Draw("same");
  cosRecosum_l5->Draw("samee");
  cosRecosum_s5->Draw("samee");

  leg->Draw();

  QQBARLabel(0.8,0.2,NULL,1);

  c3->Update();

  ////////////// Calculation //////////////
  float nominal = 30.8;

  float efficiency_l5 = ((float)nreco_l5/(float)ngen_l5) * 100;
  float efficiency_s5 = ((float)nreco_s5/(float)ngen_s5) * 100;
  cout << "--------------------------------------------------------------\n";
  cout << "--------------------------------------------------------------\n";
  cout << endl;
  cout << "LARGE Final efficiencies: " << efficiency_l5 << "% (+" << efficiency_l5 / nominal *100 -100 << "%)" << endl;
  cout << "Ratio generated : forward -> " << (float)ngenpos_l5*100/(float)ngen_l5 << "%   backward -> " << (float)ngenneg_l5*100/(float)ngen_l5 << "%" << endl; 
  cout << "Ratio reconstructed : forward -> " << (float)nrecopos_l5*100/(float)nreco_l5 << "%   backward -> " << (float)nreconeg_l5*100/(float)nreco_l5 << "%" << endl; 
  cout << endl;
  cout << "--------------------------------------------------------------\n";
  cout << "--------------------------------------------------------------\n";
  cout << endl;
  cout << "SMALL Final efficiencies: " << efficiency_s5 << "% (+" << efficiency_s5 / nominal *100 -100 << "%) " << endl;
  cout << endl;
  cout << "--------------------------------------------------------------\n";
  cout << "--------------------------------------------------------------\n";
  //file->Close();
  /*
   */
	
}

