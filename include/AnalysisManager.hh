#pragma once
#include <array>
#include <string>

class TFile;
class TTree;

class AnalysisManager {
public:
  static AnalysisManager &Instance();

  void SetOutput(const std::string &path);
  void BeginRun();
  void EndRun();

  void BeginEvent(double Etrue_GeV);
  void AddS(int towerId, long n);
  void AddC(int towerId, long n);
  void EndEvent();

private:
  AnalysisManager() = default;
  ~AnalysisManager() = default;

  std::string m_outPath = "sim/out.root";
  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;

  double m_Etrue = 0.0;
  std::array<long, 16> m_S{};
  std::array<long, 16> m_C{};
};
