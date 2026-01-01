#include "RodSD.hh"
#include "AnalysisManager.hh"
#include "CLHEP/Random/RandPoissonQ.h"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"

RodSD::RodSD(const G4String &name) : G4VSensitiveDetector(name) {}

void RodSD::Initialize(G4HCofThisEvent *) {
  // per-event init done in EventAction (below). Here do nothing.
}

static int FindTowerId(const G4VTouchable *touch) {
  for (int i = 0; i <= touch->GetHistoryDepth(); ++i) {
    auto pv = touch->GetVolume(i);
    if (pv && pv->GetName() == "TowerPV")
      return pv->GetCopyNo();
  }
  return -1;
}

static bool TouchHasVolumeName(const G4VTouchable *touch, const char *name) {
  for (int i = 0; i <= touch->GetHistoryDepth(); ++i) {
    auto pv = touch->GetVolume(i);
    if (pv && pv->GetName() == name)
      return true;
  }
  return false;
}

G4bool RodSD::ProcessHits(G4Step *step, G4TouchableHistory *) {
  auto pre = step->GetPreStepPoint();
  auto touch = pre->GetTouchable();

  const int towerId = FindTowerId(touch);
  if (towerId < 0 || towerId >= 16)
    return false;

  const bool inSFiber = TouchHasVolumeName(touch, "SFiberPV");
  const bool inCFiber = TouchHasVolumeName(touch, "CFiberPV");
  // if energy is in the blank space, then do not count
  if (!inSFiber && !inCFiber)
    return false;

  // Scintillation
  if (inSFiber) {
    const double yield = 10000.0 / MeV;
    const double eff = 0.9;
    double edep = step->GetTotalEnergyDeposit();
    double mean = edep * yield * eff;
    if (mean > 0) {
      long n = (long)CLHEP::RandPoissonQ::shoot(mean);
      AnalysisManager::Instance().AddS(towerId, n);
    }
  }

  // Cherenkov
  if (inCFiber) {
    auto track = step->GetTrack();
    double q = track->GetParticleDefinition()->GetPDGCharge();
    if (q == 0.0) return true;

    const double nref = 1.458;
    const double eff = 0.9;
    double beta = pre->GetBeta();
    if (beta > 0) {
      double denom = beta * beta * nref * nref;
      double factor = 1.0 - 1.0 / denom;
      if (factor > 0) {
        double dNdL = 369.0 * factor; // per cm
        double stepLcm = step->GetStepLength() / cm;
        double mean = dNdL * stepLcm * eff;
        if (mean > 0) {
          long n = (long)CLHEP::RandPoissonQ::shoot(mean);
          AnalysisManager::Instance().AddC(towerId, n);
        }
      }
    }
  }

  return true;
}

void RodSD::EndOfEvent(G4HCofThisEvent *) {}
