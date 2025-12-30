#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"

#include "RodSD.hh"
#include <cmath>
#include <vector>

G4VPhysicalVolume *DetectorConstruction::Construct() {
  auto nist = G4NistManager::Instance();
  auto air = nist->FindOrBuildMaterial("G4_AIR");

  // ---------------- World ----------------
  auto worldS = new G4Box("World", 5 * m, 5 * m, 5 * m);
  auto worldLV = new G4LogicalVolume(worldS, air, "WorldLV");
  auto worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV", nullptr,
                                   false, 0, true);

  const G4double rodL = 2.0 * m;
  const G4double rodXY = 4.0 * mm;

  const G4double holeR = 1.25 * mm; // 2.5 mm diameter
  const G4double fibR = 0.40 * mm;  // 0.8 mm diameter

  const int nTowerX = 4, nTowerY = 4;
  const int nRodX = 16, nRodY = 16;

  const G4double rodPitch = rodXY;

  const G4double towerSizeX = nRodX * rodPitch;
  const G4double towerSizeY = nRodY * rodPitch;
  const G4double caloSizeX = nTowerX * towerSizeX;
  const G4double caloSizeY = nTowerY * towerSizeY;

  auto caloS = new G4Box("Calo", caloSizeX / 2, caloSizeY / 2, rodL / 2);
  auto caloLV = new G4LogicalVolume(caloS, air, "CaloLV");

  auto rot = new G4RotationMatrix();
  rot->rotateY(m_rotY * deg);
  rot->rotateX(m_rotX * deg);

  new G4PVPlacement(rot, {}, caloLV, "CaloPV", worldLV, false, 0, true);

  // Tower
  auto towerS = new G4Box("Tower", towerSizeX / 2, towerSizeY / 2, rodL / 2);
  auto towerLV = new G4LogicalVolume(towerS, air, "TowerLV");

  for (int ty = 0; ty < nTowerY; ++ty) {
    for (int tx = 0; tx < nTowerX; ++tx) {
      int towerId = ty * nTowerX + tx;
      G4double x = -caloSizeX / 2 + (tx + 0.5) * towerSizeX;
      G4double y = -caloSizeY / 2 + (ty + 0.5) * towerSizeY;

      new G4PVPlacement(nullptr, {x, y, 0}, towerLV, "TowerPV", caloLV, false,
                        towerId, true);
    }
  }

  // Rod
  auto Cu = nist->FindOrBuildMaterial("G4_Cu");
  auto rodS = new G4Box("Rod", rodXY / 2, rodXY / 2, rodL / 2);
  auto rodLV = new G4LogicalVolume(rodS, Cu, "RodLV");

  // Hole
  auto holeS = new G4Tubs("Hole", 0, holeR, rodL / 2, 0, 360 * deg);
  auto holeLV = new G4LogicalVolume(holeS, air, "HoleLV");
  new G4PVPlacement(nullptr, {}, holeLV, "HolePV", rodLV, false, 0, true);

  // Fibers
  auto polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  auto quartz = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  auto fibS = new G4Tubs("Fiber", 0, fibR, rodL / 2, 0, 360 * deg);
  auto sFibLV = new G4LogicalVolume(fibS, polystyrene, "SFiberLV");
  auto cFibLV = new G4LogicalVolume(fibS, quartz, "CFiberLV");

  const G4double a = 0.80 * mm;
  const G4double rt3 = std::sqrt(3.0);

  const std::vector<G4ThreeVector> pos7 = {
      {0, 0, 0},                 // 0 center
      {0, +a, 0},                // 1 top
      {+rt3 / 2 * a, +a / 2, 0}, // 2 upper-right
      {+rt3 / 2 * a, -a / 2, 0}, // 3 lower-right
      {0, -a, 0},                // 4 bottom
      {-rt3 / 2 * a, -a / 2, 0}, // 5 lower-left
      {-rt3 / 2 * a, +a / 2, 0}  // 6 upper-left
  };

  for (int i = 0; i < 7; ++i) {
    if (i == 1 || i == 3 || i == 5) {
      new G4PVPlacement(nullptr, pos7[i], sFibLV, "SFiberPV", holeLV, false, i,
                        true);
    } else {
      new G4PVPlacement(nullptr, pos7[i], cFibLV, "CFiberPV", holeLV, false, i,
                        true);
    }
  }

  // Place rods in tower
  for (int ry = 0; ry < nRodY; ++ry) {
    for (int rx = 0; rx < nRodX; ++rx) {
      int rodId = ry * nRodX + rx;
      G4double x = -towerSizeX / 2 + (rx + 0.5) * rodPitch;
      G4double y = -towerSizeY / 2 + (ry + 0.5) * rodPitch;

      new G4PVPlacement(nullptr, {x, y, 0}, rodLV, "RodPV", towerLV, false,
                        rodId, true);
    }
  }

  return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
  auto sdman = G4SDManager::GetSDMpointer();
  auto sd = new RodSD("RodSD");
  sdman->AddNewDetector(sd);

  SetSensitiveDetector("SFiberLV", sd);
  SetSensitiveDetector("CFiberLV", sd);
}
