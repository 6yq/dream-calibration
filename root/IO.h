#pragma once
#include "types.h"

#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"

#include <algorithm>
#include <string>
#include <vector>

static inline void read_one_file(const std::string &path, OneEnergy &out) {
  TFile f(path.c_str(), "READ");
  if (f.IsZombie())
    return;

  TTree *t = (TTree *)f.Get("t");
  if (!t)
    return;

  double Etrue = 0.0;
  Long64_t S[16];
  Long64_t C[16];

  t->SetBranchAddress("Etrue_GeV", &Etrue);
  t->SetBranchAddress("S", S);
  t->SetBranchAddress("C", C);

  const Long64_t n = t->GetEntries();
  out.Ns.reserve((size_t)n);
  out.Nc.reserve((size_t)n);

  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    long long ns = 0, nc = 0;
    for (int k = 0; k < 16; ++k) {
      ns += (long long)S[k];
      nc += (long long)C[k];
    }
    out.Ns.push_back((double)ns);
    out.Nc.push_back((double)nc);
  }
  out.E = Etrue;
}

static inline std::vector<std::pair<double, std::string>>
collect_input_paths(const char *simDir) {
  std::vector<std::pair<double, std::string>> paths;

  TSystemDirectory dir("simDir", simDir);
  TList *files = dir.GetListOfFiles();
  if (!files)
    return paths;

  TIter it(files);
  while (TSystemFile *f = (TSystemFile *)it()) {
    TString name = f->GetName();
    if (f->IsDirectory())
      continue;
    if (!name.EndsWith(".root"))
      continue;
    if (!name.BeginsWith("E"))
      continue;

    // E20.root -> "20"
    TString Est = name(1, name.Length() - 1 - 5);
    double E = Est.Atof();
    std::string full = std::string(simDir) + "/" + name.Data();
    paths.push_back({E, full});
  }

  std::sort(paths.begin(), paths.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });
  return paths;
}
