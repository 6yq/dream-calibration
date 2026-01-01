#pragma once
#include "types.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <vector>

static inline void make_calib_canvas(TCanvas &c, const char *title,
                                     const std::vector<double> &E,
                                     const std::vector<double> &mu,
                                     const std::vector<double> &muErr,
                                     CalibParams &out) {
  const int n = (int)E.size();
  auto *gr = new TGraphErrors(n);

  for (int i = 0; i < n; ++i) {
    gr->SetPoint(i, mu[i], E[i]);
    gr->SetPointError(i, 0.0, 1.0);
  }
  gr->SetTitle(title);
  gr->SetMarkerStyle(20);

  const double xmin = *std::min_element(mu.begin(), mu.end());
  const double xmax = *std::max_element(mu.begin(), mu.end());

  auto *fline = new TF1("fline_cal", "[0]*x+[1]", xmin, xmax);
  double a = (E.back() - E.front()) / std::max(1e-12, (mu.back() - mu.front()));
  double b = 0.0;
  fline->SetParameters(a, b);

  TFitResultPtr r;
  for (int it = 0; it < 3; ++it) {
    for (int i = 0; i < n; ++i) {
      double yerr = std::fabs(a) * muErr[i];
      if (!(yerr > 0))
        yerr = 1.0;
      gr->SetPointError(i, 0.0, yerr);
    }
    r = gr->Fit(fline, "SQ");
    a = fline->GetParameter(0);
    b = fline->GetParameter(1);
  }

  out.a = fline->GetParameter(0);
  out.b = fline->GetParameter(1);
  out.aErr = fline->GetParError(0);
  out.bErr = fline->GetParError(1);

  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf = r->Ndf();
    out.prob = (out.ndf > 0 ? TMath::Prob(out.chi2, out.ndf) : 0.0);
  }

  c.cd();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.10);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.12);

  gr->GetXaxis()->SetTitle("Mean(N)");
  gr->GetYaxis()->SetTitle("E [GeV]");
  gr->Draw("AP");

  fline->SetLineColor(kRed);
  fline->SetLineWidth(3);
  fline->Draw("same");

  const double marginL = gPad->GetLeftMargin();
  const double marginT = gPad->GetTopMargin();
  const double w = 0.32;
  const double hh = 0.24;

  const double x1 = marginL + 0.03;
  const double x2 = x1 + w;
  const double y2 = 1.0 - marginT;
  const double y1 = y2 - hh;

  auto *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  leg->SetTextFont(42);

  leg->AddEntry(gr, "data", "p");
  leg->AddEntry(fline, "fit: E = a#timesMean + b", "l");
  leg->AddEntry((TObject *)0, Form("a=%.3g #pm %.3g", out.a, out.aErr), "");
  leg->AddEntry((TObject *)0, Form("b=%.3g #pm %.3g", out.b, out.bErr), "");
  leg->AddEntry(
      (TObject *)0,
      Form("#chi^{2}/ndf=%.1f/%d, p=%.3g", out.chi2, out.ndf, out.prob), "");

  leg->Draw();

  gPad->Modified();
  gPad->Update();
}

static inline void make_res_canvas(TCanvas &c, const char *title,
                                   const std::vector<double> &E,
                                   const std::vector<double> &sigN,
                                   const std::vector<double> &sigNerr,
                                   const CalibParams &cal, ResoParams &out) {
  const int n = (int)E.size();
  auto *gr = new TGraphErrors(n);

  for (int i = 0; i < n; ++i) {
    const double Etrue = E[i];
    const double sigmaE = cal.a * sigN[i];

    const double term_a = (sigN[i] / Etrue) * cal.aErr;
    const double term_s = (cal.a / Etrue) * sigNerr[i];
    const double y = sigmaE / Etrue;
    const double yerr = std::sqrt(term_a * term_a + term_s * term_s);

    gr->SetPoint(i, Etrue, y);
    gr->SetPointError(i, 0.0, yerr);
  }

  gr->SetTitle(title);
  gr->SetMarkerStyle(20);

  auto *fres = new TF1("fres", "sqrt(([0]/sqrt(x))*([0]/sqrt(x)) + [1]*[1])",
                       *std::min_element(E.begin(), E.end()),
                       *std::max_element(E.begin(), E.end()));
  fres->SetParameters(0.3, 0.01);

  TFitResultPtr r = gr->Fit(fres, "SQ");

  out.alpha = fres->GetParameter(0);
  out.beta = fres->GetParameter(1);
  out.alphaErr = fres->GetParError(0);
  out.betaErr = fres->GetParError(1);

  if (r.Get()) {
    out.chi2 = r->Chi2();
    out.ndf = r->Ndf();
    out.prob = (out.ndf > 0 ? TMath::Prob(out.chi2, out.ndf) : 0.0);
  }

  c.cd();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.12);

  gr->GetXaxis()->SetTitle("E [GeV]");
  gr->GetYaxis()->SetTitle("#sigma(E)/E");
  gr->Draw("AP");

  fres->SetLineColor(kRed);
  fres->SetLineWidth(3);
  fres->Draw("same");

  const double marginR = gPad->GetRightMargin();
  const double marginT = gPad->GetTopMargin();
  const double w = 0.32;
  const double hh = 0.24;

  const double x2 = 1.0 - marginR;
  const double x1 = x2 - w;
  const double y2 = 1.0 - marginT;
  const double y1 = y2 - hh;

  auto *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->SetTextFont(42);

  leg->AddEntry(gr, "data", "p");
  leg->AddEntry(fres, "fit: #sqrt{(#alpha/#sqrt{E})^{2} + #beta^{2}}", "l");
  leg->AddEntry((TObject *)0,
                Form("#alpha=%.3g #pm %.3g", out.alpha, out.alphaErr), "");
  leg->AddEntry((TObject *)0,
                Form("#beta=%.3g #pm %.3g", out.beta, out.betaErr), "");
  leg->AddEntry(
      (TObject *)0,
      Form("#chi^{2}/ndf=%.1f/%d, p=%.3g", out.chi2, out.ndf, out.prob), "");

  leg->Draw();

  gPad->Modified();
  gPad->Update();
}
