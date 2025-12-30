#pragma once
#include <string>
#include <vector>

struct FitInfo {
  double A = 0, AErr = 0;
  double mu = 0, muErr = 0;
  double sigma = 0, sigmaErr = 0;
  double chi2 = 0;
  int ndf = 0;
  double prob = 0;
  int status = 999;
};

struct CalibParams {
  double a = 0, aErr = 0;
  double b = 0, bErr = 0;
  double chi2 = 0;
  int ndf = 0;
  double prob = 0;
};

struct ResoParams {
  double alpha = 0, alphaErr = 0;
  double beta = 0, betaErr = 0;
  double chi2 = 0;
  int ndf = 0;
  double prob = 0;
};

struct CombMoments {
  double mS = 0, mC = 0;
  double vS = 0, vC = 0;
  double covSC = 0;
  double biasS = 0, biasC = 0;
  double denom = 0;
};

struct TwoFits {
  FitInfo full;
  FitInfo peak;
  double peakLo = 0, peakHi = 0;
  long long Nprime = 0; // N'
};

struct OneEnergy {
  // Etrue_GeV from simulation files
  double E = 0;
  std::string path;
  std::vector<double> Ns, Nc;

  // peak-fit parameters
  double muS = 0, muSerr = 0, sigS = 0, sigSerr = 0;
  double muC = 0, muCerr = 0, sigC = 0, sigCerr = 0;
};
