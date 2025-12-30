#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

static inline double quantile_sorted(const std::vector<double> &v_sorted,
                                     double q) {
  if (v_sorted.empty())
    return 0.0;
  if (q <= 0)
    return v_sorted.front();
  if (q >= 1)
    return v_sorted.back();
  const double pos = q * (v_sorted.size() - 1);
  const int i = (int)std::floor(pos);
  const double frac = pos - i;
  if (i + 1 >= (int)v_sorted.size())
    return v_sorted.back();
  return v_sorted[i] * (1 - frac) + v_sorted[i + 1] * frac;
}

// Freedman–Diaconis binning
static inline int fd_nbins(std::vector<double> v, double &xmin, double &xmax) {
  const int n = (int)v.size();
  if (n < 2) {
    xmin = (n ? v[0] : 0.0);
    xmax = xmin + 1.0;
    return 10;
  }
  std::sort(v.begin(), v.end());
  xmin = v.front();
  xmax = v.back();
  if (!(xmax > xmin)) {
    xmax = xmin + 1.0;
    return 10;
  }
  const double q1 = quantile_sorted(v, 0.25);
  const double q3 = quantile_sorted(v, 0.75);
  const double iqr = std::max(0.0, q3 - q1);

  double bw = 0.0;
  if (iqr > 0)
    bw = 2.0 * iqr * std::pow((double)n, -1.0 / 3.0);

  int nb = 0;
  if (bw > 0)
    nb = (int)std::ceil((xmax - xmin) / bw);
  if (nb < 5)
    nb = 5;
  return nb;
}
