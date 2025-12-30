#include "AnalysisManager.hh"
#include "TFile.h"
#include "TTree.h"
#include <stdexcept>

AnalysisManager &AnalysisManager::Instance() {
  static AnalysisManager inst;
  return inst;
}

void AnalysisManager::SetOutput(const std::string &path) { m_outPath = path; }

void AnalysisManager::BeginRun() {
  if (m_file)
    return;
  m_file = TFile::Open(m_outPath.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
    throw std::runtime_error("Failed to open ROOT output");
  m_tree = new TTree("t", "DREAM tower readout");
  m_tree->Branch("Etrue_GeV", &m_Etrue, "Etrue_GeV/D");
  m_tree->Branch("S", m_S.data(), "S[16]/L");
  m_tree->Branch("C", m_C.data(), "C[16]/L");
}

void AnalysisManager::EndRun() {
  if (!m_file)
    return;
  m_file->cd();
  m_tree->Write();
  m_file->Close();
  delete m_file;
  m_file = nullptr;
  m_tree = nullptr;
}

void AnalysisManager::BeginEvent(double Etrue_GeV) {
  m_Etrue = Etrue_GeV;
  m_S.fill(0);
  m_C.fill(0);
}

void AnalysisManager::AddS(int towerId, long n) {
  if (towerId >= 0 && towerId < 16)
    m_S[towerId] += n;
}

void AnalysisManager::AddC(int towerId, long n) {
  if (towerId >= 0 && towerId < 16)
    m_C[towerId] += n;
}

void AnalysisManager::EndEvent() {
  if (m_tree)
    m_tree->Fill();
}
