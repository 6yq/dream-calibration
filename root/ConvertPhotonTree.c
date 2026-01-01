#include "TFile.h"
#include "TTree.h"
#include <cstring>
#include <iostream>
#include <regex>
#include <string>

static double parse_Etrue_GeV(const std::string &path) {
  std::regex re("E([0-9]+(\\.[0-9]+)?)");
  std::smatch m;
  if (std::regex_search(path, m, re))
    return std::stod(m[1].str());
  std::cerr << "[ConvertPhotonTree] cannot parse E from filename: " << path
            << std::endl;
  return -1.0;
}

void ConvertPhotonTree(const char *inpath, const char *outpath,
                       const char *intree = "PhotonTree") {
  TFile *fin = TFile::Open(inpath, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input: " << inpath << "\n";
    return;
  }

  TTree *tin = nullptr;
  fin->GetObject(intree, tin);
  if (!tin) {
    std::cerr << "Cannot find tree '" << intree << "' in " << inpath << "\n";
    return;
  }

  Int_t ScintPhoton = 0;
  Int_t CerenkovPhoton = 0;
  tin->SetBranchStatus("*", 0);
  tin->SetBranchStatus("ScintPhoton", 1);
  tin->SetBranchStatus("CerenkovPhoton", 1);
  tin->SetBranchAddress("ScintPhoton", &ScintPhoton);
  tin->SetBranchAddress("CerenkovPhoton", &CerenkovPhoton);

  TFile *fout = TFile::Open(outpath, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Cannot create output: " << outpath << "\n";
    return;
  }

  TTree *tout = new TTree("t", "compat tree for analysis");

  double Etrue_GeV = parse_Etrue_GeV(std::string(inpath));
  Long64_t S[16];
  Long64_t C[16];

  tout->Branch("Etrue_GeV", &Etrue_GeV, "Etrue_GeV/D");
  tout->Branch("S", S, "S[16]/L");
  tout->Branch("C", C, "C[16]/L");

  const Long64_t n = tin->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);

    std::memset(S, 0, sizeof(S));
    std::memset(C, 0, sizeof(C));

    S[0] = (Long64_t)ScintPhoton;
    C[0] = (Long64_t)CerenkovPhoton;

    tout->Fill();
  }

  fout->Write();
  fout->Close();
  fin->Close();
}
