#include "IO.h"
#include "comb.h"
#include "common.h"
#include "plot_per_energy.h"
#include "plot_summary.h"
#include "types.h"
#include "write.h"

void analysis(const char *simDir, const char *resDir) {
  style_global();
  gSystem->mkdir(resDir, true);

  // read in
  auto paths = collect_input_paths(simDir);
  if (paths.empty())
    return;

  std::vector<OneEnergy> all;
  all.reserve(paths.size());
  for (auto &p : paths) {
    OneEnergy d;
    d.path = p.second;
    read_one_file(d.path, d);
    if (!d.Ns.empty())
      all.push_back(std::move(d));
  }
  if (all.empty())
    return;

  // per-energy fits
  const std::string perPdf = std::string(resDir) + "/per_energy_fits.pdf";
  TCanvas cper("cper", "per-energy", 1200, 450);
  cper.Print((perPdf + "[").c_str());

  for (auto &d : all) {
    double muC, muCerr, sigC, sigCerr;
    double muS, muSerr, sigS, sigSerr;
    draw_one_energy_page(cper, d, perPdf, muC, muCerr, sigC, sigCerr, muS,
                         muSerr, sigS, sigSerr);
    d.muC = muC;
    d.muCerr = muCerr;
    d.sigC = sigC;
    d.sigCerr = sigCerr;
    d.muS = muS;
    d.muSerr = muSerr;
    d.sigS = sigS;
    d.sigSerr = sigSerr;
  }
  cper.Print((perPdf + "]").c_str());

  // results
  std::vector<double> E, muS, muSerr, sigS, sigSerr, muC, muCerr, sigC, sigCerr;
  for (auto &d : all) {
    E.push_back(d.E);
    muS.push_back(d.muS);
    muSerr.push_back(d.muSerr);
    sigS.push_back(d.sigS);
    sigSerr.push_back(d.sigSerr);
    muC.push_back(d.muC);
    muCerr.push_back(d.muCerr);
    sigC.push_back(d.sigC);
    sigCerr.push_back(d.sigCerr);
  }
  if (E.size() < 2)
    return;

  // summary
  const std::string sumPdf = std::string(resDir) + "/summary.pdf";
  TCanvas cs("cs", "summary", 900, 650);
  cs.Print((sumPdf + "[").c_str());

  CalibParams calS, calC;
  ResoParams resS, resC;

  cs.Clear();
  make_calib_canvas(cs, "Calibration (S): E vs Mean(N_{S})", E, muS, muSerr,
                    calS);
  cs.Print(sumPdf.c_str());
  cs.Clear();
  make_calib_canvas(cs, "Calibration (C): E vs Mean(N_{C})", E, muC, muCerr,
                    calC);
  cs.Print(sumPdf.c_str());
  cs.Clear();
  make_res_canvas(cs, "Resolution (S): #sigma(E)/E vs E", E, sigS, sigSerr,
                  calS, resS);
  cs.Print(sumPdf.c_str());
  cs.Clear();
  make_res_canvas(cs, "Resolution (C): #sigma(E)/E vs E", E, sigC, sigCerr,
                  calC, resC);
  cs.Print(sumPdf.c_str());

  cs.Print((sumPdf + "]").c_str());

  // per-energy comb
  const std::string combPerPdf =
      std::string(resDir) + "/per_energy_comb_fits.pdf";
  TCanvas ccomb("ccomb", "per-energy-comb", 900, 650);
  ccomb.Print((combPerPdf + "[").c_str());

  std::vector<double> wStar, muComb, muCombErr, sigComb, sigCombErr;
  std::vector<long long> NpComb;
  std::vector<CombMoments> combMM;

  wStar.reserve(all.size());
  muComb.reserve(all.size());
  muCombErr.reserve(all.size());
  sigComb.reserve(all.size());
  sigCombErr.reserve(all.size());
  NpComb.reserve(all.size());
  combMM.reserve(all.size());

  for (auto &d : all) {
    std::vector<double> Es, Ec, tmp;
    build_Es_Ec_Ecomb(d, calS, calC, 0.0, Es, Ec, tmp);

    CombMoments mm;
    const double w = w_mse_star(Es, Ec, d.E, &mm);

    std::vector<double> Ecomb;
    build_Es_Ec_Ecomb(d, calS, calC, w, Es, Ec, Ecomb);

    double mu = 0, muE = 0, sg = 0, sgE = 0;
    long long Np = 0;
    double wForPlot = w;

    fit_vector_gaus_and_draw(ccomb, Ecomb, d.E, combPerPdf, "E_{comb}  [GeV]",
                             mu, muE, sg, sgE, wForPlot, Np, mm, true);

    wStar.push_back(w);
    combMM.push_back(mm);
    muComb.push_back(mu);
    muCombErr.push_back(muE);
    sigComb.push_back(sg);
    sigCombErr.push_back(sgE);
    NpComb.push_back(Np);
  }

  ccomb.Print((combPerPdf + "]").c_str());

  // combined summary pdf
  const std::string combSumPdf = std::string(resDir) + "/combined_summary.pdf";
  TCanvas ccs("ccs", "combined-summary", 900, 650);
  ccs.Print((combSumPdf + "[").c_str());

  CalibParams calComb;
  ResoParams resComb;

  ccs.Clear();
  make_calib_canvas(ccs, "Calibration (comb): E vs Mean(E_{comb})", E, muComb,
                    muCombErr, calComb);
  ccs.Print(combSumPdf.c_str());
  ccs.Clear();
  make_res_canvas(ccs, "Resolution (comb): #sigma(E)/E vs E", E, sigComb,
                  sigCombErr, calComb, resComb);
  ccs.Print(combSumPdf.c_str());

  // w*(E)
  {
    auto *grw = new TGraphErrors((int)E.size());
    for (int i = 0; i < (int)E.size(); ++i) {
      grw->SetPoint(i, E[i], wStar[i]);
      grw->SetPointError(i, 0.0, 0.0);
    }
    ccs.Clear();
    ccs.cd();
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.12);

    grw->SetTitle("Weight (comb): w*(E)");
    grw->SetMarkerStyle(20);
    grw->GetXaxis()->SetTitle("E [GeV]");
    grw->GetYaxis()->SetTitle("w*");
    grw->Draw("AP");

    ccs.Print(combSumPdf.c_str());
    delete grw;
  }

  ccs.Print((combSumPdf + "]").c_str());

  // summary txt
  const std::string parTxt = std::string(resDir) + "/fit_params.txt";
  write_all_results_txt(parTxt, calS, calC, calComb, resS, resC, resComb, E,
                        muS, sigS, muC, sigC, muComb, muCombErr, sigComb,
                        sigCombErr, wStar, combMM, NpComb, sigSerr, sigCerr);
}
