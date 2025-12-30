#pragma once
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class RodSD : public G4VSensitiveDetector {
public:
  explicit RodSD(const G4String &name);
  ~RodSD() override = default;

  void Initialize(G4HCofThisEvent *) override;
  G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;
  void EndOfEvent(G4HCofThisEvent *) override;
};
