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

void comparison()
{

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
  //std::cout << "Enter code for s5: ";
  //std::cin >> token_s5;

  std::cout << std::endl;

  TCanvas * c1 = new TCanvas("c1", "gammaerr");
  TCanvas * c2 = new TCanvas("c2", "tmomentum",0,0,1500,500);
  TCanvas * c3 = new TCanvas("c3", "bmomentum",0,0,1500,500);
  TCanvas * c4 = new TCanvas("c4", "Wmomentum",0,0,1500,500);
  TCanvas * c5 = new TCanvas("c5", "top12b");
  
  TH2F * gammaerrb = new TH2F("gammaerrb",NULL,50,1,3,90,1,5);

  TH1F * top12b = new TH1F("top12b", NULL, 200, 0, 200);
  TH1F * bbbar = new TH1F("bbbar", NULL, 200, 0, 200);
  TH1F * mcbbbar = new TH1F("mcbbbar", NULL, 200, 0, 200);

  vector< TH2F* > t_P2;
  vector< TH2F* > b_P2;
  vector< TH2F* > W_P2;
  
  for( int i=0; i<3; i++){
    
    t_P2.push_back( new TH2F(Form("t_p%d",i),NULL,1000,0,4000,1000,0,4000) );

    b_P2.push_back( new TH2F(Form("b_p%d",i),NULL,1000,0,4000,1000,0,4000) );

    W_P2.push_back( new TH2F(Form("W_p%d",i),NULL,1000,0,4000,1000,0,4000) );

  }

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

  TGaxis::SetMaxDigits(3);

  TTree * normaltree_l5 = (TTree*) file_l5->Get( "Stats" ) ;
  TTree * GenTree_l5 = (TTree*) file_l5->Get( "GenTree" ) ;

  int neventGen_l5 = GenTree_l5->GetEntries();
  int neventReco_l5 = normaltree_l5->GetEntries();

  cout << "l5Reco Entry = " << neventReco_l5 << endl;
  cout << "l5Gen Entry = " << neventGen_l5 << endl;
  cout << endl;

  normaltree_l5->SetBranchAddress("MCTopfourvect",&MCTopfourvect);
  normaltree_l5->SetBranchAddress("MCWfourvect",&MCWfourvect);
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
    
    if( /*Thrust < 0.9 && hadMass > 180 && hadMass < 420 && Top1mass < 270 && W1mass < 250 && Top1mass > 120 && W1mass > 50 && Top1bmomentum > 15 && Top2bmomentum > 15 && (Top1gamma + Top2gamma) > 2.4  && Top2gamma < 2 && methodok && bboostedcos>-1 && bboostedcos<1 /*&& cosWb < -0.99*/ true){

      

      gammaerrb->Fill( computeGamma(MCTopfourvect[0]) , computeGamma(Topfourvect[0]) );

      top12b->Fill(Top1bmomentum);
      top12b->Fill(Top2bmomentum);
      bbbar->Fill(computeMag(qBfourvect[0]));
      bbbar->Fill(computeMag(qBfourvect[1]));
      mcbbbar->Fill(computeMag(qMCBfourvect[0]));
      mcbbbar->Fill(computeMag(qMCBfourvect[1]));

      for( int i=0; i<3; i++){

	t_P2.at(i)->Fill(MCTopfourvect[0][i]*MCTopfourvect[0][i],Topfourvect[0][i]*Topfourvect[0][i]);

	b_P2.at(i)->Fill(qMCBfourvect[0][i]*qMCBfourvect[0][i],qBfourvect[0][i]*qBfourvect[0][i]);

	W_P2.at(i)->Fill(MCWfourvect[0][i]*MCWfourvect[0][i],Wfourvect[0][i]*Wfourvect[0][i]);
  
      }

    }

  }

  int n_top12b = top12b->GetEntries();
  int n_bbbar = bbbar->GetEntries();
  
  cout << "top12b entries : " << n_top12b << endl;
  cout << "bbbar entries : " << n_bbbar << endl;
  cout << endl;

  ////////////// Detector Model SMALL //////////////



  ////////////// Style Setting //////////////

 

  ////////////// Drawing //////////////

  c1->cd();
  c1->cd()->SetRightMargin(0.18);
  gammaerrb->GetXaxis()->SetTitle("gamma_{MC}");
  gammaerrb->GetYaxis()->SetTitle("gamma_{Reco}");
  gammaerrb->Draw("colz");

  c2->Divide(3,1);
  c3->Divide(3,1);
  c4->Divide(3,1);

  for(int i=0; i<3; i++){

    c2->cd(i+1);
    c2->cd(i+1)->SetRightMargin(0.18);
    t_P2.at(i)->SetMaximum(5);
    t_P2.at(i)->GetXaxis()->SetTitle("p^{2}_{MC}");
    t_P2.at(i)->GetYaxis()->SetTitle("p^{2}_{Reco}");
    t_P2.at(i)->Draw("colz");

    c3->cd(i+1);
    c3->cd(i+1)->SetRightMargin(0.18);
    b_P2.at(i)->GetXaxis()->SetTitle("p^{2}_{MC}");
    b_P2.at(i)->GetYaxis()->SetTitle("p^{2}_{Reco}");
    b_P2.at(i)->SetMaximum(5);
    b_P2.at(i)->Draw("colz");

    c4->cd(i+1);
    c4->cd(i+1)->SetRightMargin(0.18);
    W_P2.at(i)->GetXaxis()->SetTitle("p^{2}_{MC}");
    W_P2.at(i)->GetYaxis()->SetTitle("p^{2}_{Reco}");
    W_P2.at(i)->SetMaximum(5);
    W_P2.at(i)->Draw("colz");

  }

  c5->cd();
  top12b->SetLineColor(kRed);
  top12b->SetMaximum(30000);
  top12b->Draw();
  bbbar->SetLineColor(kBlue);
  bbbar->Draw("same");
  mcbbbar->SetLineColor(kGray+1);
  mcbbbar->Draw("same");

}

