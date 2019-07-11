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

void class_asymmetry_sl()
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
  gStyle->SetMarkerSize(0);
  gStyle->SetTitleX(0.2); 
  gStyle->SetTitleY(0.9); 
	
  FileSelector fs;
  std::vector<FileSelector> rootfiles;
  std::ifstream in ( "/home/ilc/vlohezic/working/VertexMacro/semi_leptonic/input/record.txt" ); 

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

  //////////////// Detector Model LARGE : Fit and analysis ////////////////

  std::string filename_l5 = rootfiles[token_l5].filename();
  std::cout << "Processing (l5) : " << filename_l5 << " ..." << endl;

  asymmetry *res_l5 = new asymmetry(filename_l5);

  std::cout << "---------------------LARGE Detector Model---------------------" << std::endl;

  res_l5->analyse();

  //////////////// Detector Model SMALL : Fit and analysis ///////////////
  
  std::string filename_s5 = rootfiles[token_s5].filename();
  std::cout << "Processing (s5) : " << filename_s5 << " ..." << endl;

  asymmetry *res_s5 = new asymmetry(filename_s5);

  std::cout << "---------------------SMALL Detector Model---------------------" << std::endl;

  res_s5->analyse();

  /////////////// Drawing ///////////////

  res_l5->show(res_s5);

}
	
