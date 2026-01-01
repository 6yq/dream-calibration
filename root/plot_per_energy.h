#pragma once
#include "bins.h"
#include "fit.h"
#include "types.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"

#include <algorithm>
#include <cmath>

static inline void add_fit_legend(TLegend *leg, TF1 *fFull, const FitInfo &FF,
                                  TF1 *fPeak, const FitInfo &FP,
                                  long long Nprime) {
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);

  leg->SetTextFont(42);
  const char *MU = "#font[122]{m}";

  leg->AddEntry(fFull, "Full data", "l");
  leg->AddEntry((TObject *)0, Form("%s=%.4g, #sigma=%.4g", MU, FF.mu, FF.sigma),
                "");
  leg->AddEntry((TObject *)0,
                Form("#chi^{2}/ndf=%.1f/%d, p=%.3g", FF.chi2, FF.ndf, FF.prob),
                "");

  leg->AddEntry(fPeak,
                Form("Peak (%s#pm1.5#sigma): N'=%lld", MU, (long long)Nprime),
                "l");
  leg->AddEntry((TObject *)0,
                Form("%s'=%.4g, #sigma'=%.4g", MU, FP.mu, FP.sigma), "");
  leg->AddEntry((TObject *)0,
                Form("#chi^{2}/ndf=%.1f/%d, p=%.3g", FP.chi2, FP.ndf, FP.prob),
                "");
}

static inline void draw_one_energy_page(TCanvas &c, const OneEnergy &d,
                                        const std::string &pdfName,
                                        double &out_muC, double &out_muCerr,
                                        double &out_sigC, double &out_sigCerr,
                                        double &out_muS, double &out_muSerr,
                                        double &out_sigS, double &out_sigSerr) {
  const int N = (int)d.Ns.size();
  if (N <= 0)
    return;

  const double kWin = 1.5;

  auto build_hist = [&](const std::vector<double> &v, const TString &hname,
                        const char *xTitle, double &xminPad,
                        double &xmaxPad) -> TH1D * {
    double xmin = 0, xmax = 0;
    int nb = fd_nbins(v, xmin, xmax);
    double pad = 0.02 * (xmax - xmin);
    if (!(pad > 0))
      pad = 1.0;

    xminPad = xmin - pad;
    if (xminPad < 0)
      xminPad = 0;
    xmaxPad = xmax + pad;

    TH1D *h = new TH1D(hname.Data(), "", nb, xminPad, xmaxPad);
    h->SetDirectory(nullptr);
    for (double x : v)
      h->Fill(x);

    h->SetLineColor(kBlue + 1);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle("Events");
    return h;
  };

  double xC1, xC2, xS1, xS2;
  TH1D *hC = build_hist(d.Nc, Form("hC_E%.0f", d.E), "N_{C}", xC1, xC2);
  TH1D *hS = build_hist(d.Ns, Form("hS_E%.0f", d.E), "N_{S}", xS1, xS2);

  TwoFits FC = do_full_then_peak(hC, xC1, xC2, kWin, Form("fC_full_E%.0f", d.E),
                                 Form("fC_peak_E%.0f", d.E));
  TwoFits FS = do_full_then_peak(hS, xS1, xS2, kWin, Form("fS_full_E%.0f", d.E),
                                 Form("fS_peak_E%.0f", d.E));

  out_muC = FC.peak.mu;
  out_muCerr = FC.peak.muErr;
  out_sigC = FC.peak.sigma;
  out_sigCerr = FC.peak.sigmaErr;
  out_muS = FS.peak.mu;
  out_muSerr = FS.peak.muErr;
  out_sigS = FS.peak.sigma;
  out_sigSerr = FS.peak.sigmaErr;

  TF1 fC_full(Form("draw_fC_full_E%.0f", d.E), "gaus", xC1, xC2);
  fC_full.SetNpx(800);
  fC_full.SetParameters(FC.full.A, FC.full.mu, FC.full.sigma);
  fC_full.SetLineColor(kRed);
  fC_full.SetLineStyle(1);
  fC_full.SetLineWidth(3);

  TF1 fC_peak(Form("draw_fC_peak_E%.0f", d.E), "gaus", FC.peakLo, FC.peakHi);
  fC_peak.SetNpx(800);
  fC_peak.SetParameters(FC.peak.A, FC.peak.mu, FC.peak.sigma);
  fC_peak.SetLineColor(kRed);
  fC_peak.SetLineStyle(2);
  fC_peak.SetLineWidth(3);

  TF1 fS_full(Form("draw_fS_full_E%.0f", d.E), "gaus", xS1, xS2);
  fS_full.SetNpx(800);
  fS_full.SetParameters(FS.full.A, FS.full.mu, FS.full.sigma);
  fS_full.SetLineColor(kRed);
  fS_full.SetLineStyle(1);
  fS_full.SetLineWidth(3);

  TF1 fS_peak(Form("draw_fS_peak_E%.0f", d.E), "gaus", FS.peakLo, FS.peakHi);
  fS_peak.SetNpx(800);
  fS_peak.SetParameters(FS.peak.A, FS.peak.mu, FS.peak.sigma);
  fS_peak.SetLineColor(kRed);
  fS_peak.SetLineStyle(2);
  fS_peak.SetLineWidth(3);

  c.Clear();
  c.Divide(2, 1);

  auto draw_pad = [&](int ipad, TH1D *h, TF1 &fFull, const FitInfo &FF,
                      TF1 &fPeak, const FitInfo &FP, long long Np,
                      const char *title) {
    c.cd(ipad);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.12);
    gPad->SetBottomMargin(0.12);

    h->SetTitle(title);
    h->Draw("hist");
    fFull.Draw("same");
    fPeak.Draw("same");

    double marginR = 0.10;
    double marginT = 0.12;
    double marginL = 0.15;
    double w = 0.36;
    double hh = 0.36;

    double x2 = 1.0 - marginR;
    double y2 = 1.0 - marginT;
    double x1 = x2 - w;
    double y1 = y2 - hh;

    if (ipad == 2) {
      x1 = marginL;
      x2 = x1 + w;
    } else {
      x2 = 1.0 - marginR;
      x1 = x2 - w;
    }

    auto *leg = new TLegend(x1, y1, x2, y2);

    leg->SetHeader(Form("N=%d", N), "");
    add_fit_legend(leg, &fFull, FF, &fPeak, FP, Np);
    leg->Draw();

    gPad->Modified();
    gPad->Update();
  };

  draw_pad(1, hC, fC_full, FC.full, fC_peak, FC.peak, FC.Nprime,
           Form("E = %.0f GeV : N_{C}", d.E));
  draw_pad(2, hS, fS_full, FS.full, fS_peak, FS.peak, FS.Nprime,
           Form("E = %.0f GeV : N_{S}", d.E));

  c.Print(pdfName.c_str());

  delete hC;
  delete hS;
}
