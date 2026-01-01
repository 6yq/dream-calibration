#pragma once
#include "bins.h"
#include "fit.h"
#include "plot_per_energy.h"
#include "types.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"

#include <algorithm>
#include <cmath>
#include <vector>

static inline void mean_var(const std::vector<double> &x, double &m,
                            double &v) {
  const int n = (int)x.size();
  if (n <= 0) {
    m = 0;
    v = 0;
    return;
  }
  double s = 0;
  for (double xi : x)
    s += xi;
  m = s / n;
  double ss = 0;
  for (double xi : x) {
    double d = xi - m;
    ss += d * d;
  }
  v = (n > 1 ? ss / (n - 1) : 0.0);
}

static inline double cov_xy(const std::vector<double> &x,
                            const std::vector<double> &y, double mx,
                            double my) {
  const int n = (int)std::min(x.size(), y.size());
  if (n <= 1)
    return 0.0;
  double s = 0;
  for (int i = 0; i < n; ++i)
    s += (x[i] - mx) * (y[i] - my);
  return s / (n - 1);
}

static inline double w_mse_star(const std::vector<double> &Es,
                                const std::vector<double> &Ec, double Etrue,
                                CombMoments *out = nullptr) {
  CombMoments mm;
  mean_var(Es, mm.mS, mm.vS);
  mean_var(Ec, mm.mC, mm.vC);
  mm.covSC = cov_xy(Es, Ec, mm.mS, mm.mC);
  mm.biasS = mm.mS - Etrue;
  mm.biasC = mm.mC - Etrue;

  const double Var_d = mm.vS + mm.vC - 2.0 * mm.covSC;
  const double Mean_d = (mm.biasS - mm.biasC);
  mm.denom = Var_d + Mean_d * Mean_d;

  const double numer = (mm.covSC - mm.vC) + mm.biasC * (mm.biasS - mm.biasC);

  double w = 0.5;
  if (std::fabs(mm.denom) > 1e-18)
    w = -numer / mm.denom;

  if (out)
    *out = mm;
  return w;
}

static inline void
build_Es_Ec_Ecomb(const OneEnergy &d, const CalibParams &calS,
                  const CalibParams &calC, double w, std::vector<double> &Es,
                  std::vector<double> &Ec, std::vector<double> &Ecomb) {
  const int n = (int)std::min(d.Ns.size(), d.Nc.size());
  Es.resize(n);
  Ec.resize(n);
  Ecomb.resize(n);
  for (int i = 0; i < n; ++i) {
    const double es = calS.a * d.Ns[i] + calS.b;
    const double ec = calC.a * d.Nc[i] + calC.b;
    Es[i] = es;
    Ec[i] = ec;
    Ecomb[i] = w * es + (1.0 - w) * ec;
  }
}

static inline void fit_vector_gaus_and_draw(
    TCanvas &c, const std::vector<double> &v, double Etrue,
    const std::string &pdfName, const char *xTitle, double &out_mu,
    double &out_muErr, double &out_sig, double &out_sigErr, double &out_w,
    long long &out_Nprime, const CombMoments &mm, bool draw = true) {
  double xmin = 0, xmax = 0;
  int nb = fd_nbins(v, xmin, xmax);
  double pad = 0.02 * (xmax - xmin);
  if (!(pad > 0))
    pad = 1.0;
  double x1 = std::max(0.0, xmin - pad);
  double x2 = xmax + pad;

  TH1D *h = new TH1D(Form("hComb_E%.0f", Etrue), "", nb, x1, x2);
  h->SetDirectory(nullptr);
  for (double x : v)
    h->Fill(x);
  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle("Events");

  const double kWin = 1.5;
  TwoFits F =
      do_full_then_peak(h, x1, x2, kWin, Form("fComb_full_E%.0f", Etrue),
                        Form("fComb_peak_E%.0f", Etrue));

  out_mu = F.peak.mu;
  out_muErr = F.peak.muErr;
  out_sig = F.peak.sigma;
  out_sigErr = F.peak.sigmaErr;
  out_Nprime = F.Nprime;

  if (draw) {
    TF1 fFull(Form("draw_fComb_full_E%.0f", Etrue), "gaus", x1, x2);
    fFull.SetNpx(800);
    fFull.SetParameters(F.full.A, F.full.mu, F.full.sigma);
    fFull.SetLineColor(kRed);
    fFull.SetLineStyle(1);
    fFull.SetLineWidth(3);

    TF1 fPeak(Form("draw_fComb_peak_E%.0f", Etrue), "gaus", F.peakLo, F.peakHi);
    fPeak.SetNpx(800);
    fPeak.SetParameters(F.peak.A, F.peak.mu, F.peak.sigma);
    fPeak.SetLineColor(kRed);
    fPeak.SetLineStyle(2);
    fPeak.SetLineWidth(3);

    c.Clear();
    c.cd();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.12);
    gPad->SetBottomMargin(0.12);

    h->SetTitle(Form("E = %.0f GeV : E_{comb}", Etrue));
    h->Draw("hist");
    fFull.Draw("same");
    fPeak.Draw("same");

    const double marginL = gPad->GetLeftMargin();
    const double marginT = gPad->GetTopMargin();
    const double w = 0.34;
    const double hh = 0.36;

    const double x1_leg = marginL + 0.03;
    const double x2_leg = x1_leg + w;
    const double y2_leg = 1.0 - marginT;
    const double y1_leg = y2_leg - hh;

    auto *leg = new TLegend(x1_leg, y1_leg, x2_leg, y2_leg);
    leg->SetHeader(Form("w*=%.4g", out_w), "");
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->SetTextFont(42);

    add_fit_legend(leg, &fFull, F.full, &fPeak, F.peak, F.Nprime);

    leg->Draw();
    c.Print(pdfName.c_str());
  }

  delete h;
}
