#pragma once
#include "types.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

static inline void write_all_results_txt(
    const std::string &path, const CalibParams &calS, const CalibParams &calC,
    const CalibParams &calComb, const ResoParams &resS, const ResoParams &resC,
    const ResoParams &resComb, const std::vector<double> &E,
    const std::vector<double> &muS, const std::vector<double> &sigS,
    const std::vector<double> &muC, const std::vector<double> &sigC,
    const std::vector<double> &muComb, const std::vector<double> &muCombErr,
    const std::vector<double> &sigComb, const std::vector<double> &sigCombErr,
    const std::vector<double> &wStar, const std::vector<CombMoments> &combMM,
    const std::vector<long long> &NpComb, const std::vector<double> &sigSerr,
    const std::vector<double> &sigCerr) {

  std::ofstream ofs(path.c_str());
  ofs << std::setprecision(12);

  ofs << "# Calibration: E = a*Mean + b\n";
  ofs << "S_calib  a " << calS.a << "  aErr " << calS.aErr << "  b " << calS.b
      << "  bErr " << calS.bErr << "  chi2 " << calS.chi2 << "  ndf "
      << calS.ndf << "  p " << calS.prob << "\n";
  ofs << "C_calib  a " << calC.a << "  aErr " << calC.aErr << "  b " << calC.b
      << "  bErr " << calC.bErr << "  chi2 " << calC.chi2 << "  ndf "
      << calC.ndf << "  p " << calC.prob << "\n";
  ofs << "Comb_calib  a " << calComb.a << "  aErr " << calComb.aErr << "  b "
      << calComb.b << "  bErr " << calComb.bErr << "  chi2 " << calComb.chi2
      << "  ndf " << calComb.ndf << "  p " << calComb.prob << "\n";

  ofs << "\n# Resolution: sigma(E)/E = sqrt((alpha/sqrt(E))^2 + beta^2)\n";
  ofs << "S_res  alpha " << resS.alpha << "  alphaErr " << resS.alphaErr
      << "  beta " << resS.beta << "  betaErr " << resS.betaErr << "  chi2 "
      << resS.chi2 << "  ndf " << resS.ndf << "  p " << resS.prob << "\n";
  ofs << "C_res  alpha " << resC.alpha << "  alphaErr " << resC.alphaErr
      << "  beta " << resC.beta << "  betaErr " << resC.betaErr << "  chi2 "
      << resC.chi2 << "  ndf " << resC.ndf << "  p " << resC.prob << "\n";
  ofs << "Comb_res  alpha " << resComb.alpha << "  alphaErr "
      << resComb.alphaErr << "  beta " << resComb.beta << "  betaErr "
      << resComb.betaErr << "  chi2 " << resComb.chi2 << "  ndf " << resComb.ndf
      << "  p " << resComb.prob << "\n";

  ofs << "\n# Per-energy summary\n";
  ofs << "# E  w*  muS  sigS  muC  sigC  muComb  sigComb  NpComb\n";
  for (int i = 0; i < (int)E.size(); ++i) {
    ofs << E[i] << "  " << wStar[i] << "  " << muS[i] << "  " << sigS[i] << "  "
        << muC[i] << "  " << sigC[i] << "  " << muComb[i] << "  " << sigComb[i]
        << "  " << NpComb[i] << "\n";
  }

  ofs << "\n# Moments used in w*(E)\n";
  ofs << "# E  meanEs  meanEc  varEs  varEc  covSC  biasS  biasC  denom\n";
  for (int i = 0; i < (int)E.size(); ++i) {
    const auto &m = combMM[i];
    ofs << E[i] << "  " << m.mS << "  " << m.mC << "  " << m.vS << "  " << m.vC
        << "  " << m.covSC << "  " << m.biasS << "  " << m.biasC << "  "
        << m.denom << "\n";
  }

  ofs << "\n# Per-energy resolution (true E in denominator)\n";
  ofs << "# E  dE/E_S  errS  dE/E_C  errC  dE/E_Comb  errComb\n";

  auto rel_and_err = [&](double Etrue, double a, double aErr, double sig,
                         double sigErr, double &rel, double &relErr) {
    const double sigE = a * sig;
    const double sigEerr =
        std::sqrt(std::pow(sig * aErr, 2) + std::pow(a * sigErr, 2));
    rel = sigE / Etrue;
    relErr = sigEerr / Etrue;
  };

  for (int i = 0; i < (int)E.size(); ++i) {
    double rS = 0, rSe = 0, rC = 0, rCe = 0, rB = 0, rBe = 0;
    rel_and_err(E[i], calS.a, calS.aErr, sigS[i], sigSerr[i], rS, rSe);
    rel_and_err(E[i], calC.a, calC.aErr, sigC[i], sigCerr[i], rC, rCe);
    rel_and_err(E[i], calComb.a, calComb.aErr, sigComb[i], sigCombErr[i], rB,
                rBe);

    ofs << E[i] << "  " << rS << "  " << rSe << "  " << rC << "  " << rCe
        << "  " << rB << "  " << rBe << "\n";
  }
}
