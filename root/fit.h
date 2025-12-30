#pragma once
#include "bins.h"
#include "types.h"

#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>

static inline int nbins_in_range(TH1D *h, double lo, double hi) {
  if (!h)
    return 0;
  int b1 = h->GetXaxis()->FindBin(lo);
  int b2 = h->GetXaxis()->FindBin(hi);
  if (b1 < 1)
    b1 = 1;
  if (b2 > h->GetNbinsX())
    b2 = h->GetNbinsX();
  return std::max(0, b2 - b1 + 1);
}

static inline long long integral_in_range(TH1D *h, double lo, double hi) {
  if (!h)
    return 0;
  int b1 = h->GetXaxis()->FindBin(lo);
  int b2 = h->GetXaxis()->FindBin(hi);
  if (b1 < 1)
    b1 = 1;
  if (b2 > h->GetNbinsX())
    b2 = h->GetNbinsX();
  double s = h->Integral(b1, b2);
  if (s < 0)
    s = 0;
  return (long long)llround(s);
}

// L: likelihood, I: bin integral, R: range, S: save, Q: quiet
static inline FitInfo fit_hist_gaus_LI(TH1D *h, TF1 *f, double lo, double hi) {
  FitInfo info;
  if (!h || !f)
    return info;

  f->SetRange(lo, hi);
  TFitResultPtr r = h->Fit(f, "LIRSQ", "", lo, hi);
  info.status = (int)r;

  info.A = f->GetParameter(0);
  info.AErr = f->GetParError(0);
  info.mu = f->GetParameter(1);
  info.muErr = f->GetParError(1);
  info.sigma = f->GetParameter(2);
  info.sigmaErr = f->GetParError(2);

  info.chi2 = h->Chisquare(f, "IR");
  const int nBins = nbins_in_range(h, lo, hi);
  const int nPar = f->GetNpar();
  info.ndf = std::max(1, nBins - nPar);
  info.prob = TMath::Prob(info.chi2, info.ndf);

  return info;
}

// the full data are not precisely norm,
// so the full dataset MLE is only to provide a cut
static inline TwoFits do_full_then_peak(TH1D *h, double xmin, double xmax,
                                        double kWin, const char *tagFull,
                                        const char *tagPeak) {
  TwoFits out;

  TF1 fFull(tagFull, "gaus", xmin, xmax);
  fFull.SetNpx(500);
  fFull.SetParameters(h->GetMaximum(), h->GetMean(),
                      std::max(1e-9, h->GetRMS()));
  out.full = fit_hist_gaus_LI(h, &fFull, xmin, xmax);

  double lo = out.full.mu - kWin * std::fabs(out.full.sigma);
  double hi = out.full.mu + kWin * std::fabs(out.full.sigma);
  if (lo < xmin)
    lo = xmin;
  if (hi > xmax)
    hi = xmax;
  if (!(hi > lo)) {
    lo = xmin;
    hi = xmax;
  }

  TF1 fPeak(tagPeak, "gaus", lo, hi);
  fPeak.SetNpx(500);
  fPeak.SetParameters(h->GetMaximum(), out.full.mu,
                      std::max(1e-9, std::fabs(out.full.sigma)));
  out.peak = fit_hist_gaus_LI(h, &fPeak, lo, hi);

  out.peakLo = lo;
  out.peakHi = hi;
  out.Nprime = integral_in_range(h, lo, hi);

  return out;
}
