#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction() = default;
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume *Construct() override;
  void ConstructSDandField() override;

  void SetRotationX(G4double v) { m_rotX = v; }
  void SetRotationY(G4double v) { m_rotY = v; }

private:
  G4double m_rotX = 0.0;
  G4double m_rotY = 0.0;
};

#endif
