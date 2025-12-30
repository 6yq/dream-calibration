#include "FTFP_BERT.hh"
#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"

#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"

#include "Randomize.hh"
#include <string>

// output
class AnaMessenger : public G4UImessenger {
public:
  AnaMessenger() {
    dir = new G4UIdirectory("/ana/");
    cmd = new G4UIcmdWithAString("/ana/output", this);
  }
  ~AnaMessenger() override {
    delete cmd;
    delete dir;
  }

  void SetNewValue(G4UIcommand *command, G4String newValue) override {
    if (command == cmd) {
      AnalysisManager::Instance().SetOutput(newValue);
    }
  }

private:
  G4UIdirectory *dir = nullptr;
  G4UIcmdWithAString *cmd = nullptr;
};

// primary gun
class PrimaryGun : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGun() {
    gun = new G4ParticleGun(1);
    SetParticle("pi-");
    SetEnergy(1.0 * GeV);

    gun->SetParticleMomentumDirection({0, 0, 1});
    gun->SetParticlePosition({0 * mm, 0 * mm, -2.0 * m});
  }

  ~PrimaryGun() override { delete gun; }

  void GeneratePrimaries(G4Event *evt) override {
    gun->GeneratePrimaryVertex(evt);
  }

  void SetParticle(const G4String &name) {
    auto *p = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if (p)
      gun->SetParticleDefinition(p);
  }

  void SetEnergy(G4double e) { gun->SetParticleEnergy(e); }

  G4double GetEnergy() const { return gun ? gun->GetParticleEnergy() : 0.0; }

private:
  G4ParticleGun *gun = nullptr;
};

// gun
class GunMessenger : public G4UImessenger {
public:
  explicit GunMessenger(PrimaryGun *pg) : gun(pg) {
    dir = new G4UIdirectory("/gun/");

    particle = new G4UIcmdWithAString("/gun/particle", this);
    particle->SetGuidance("Set particle name, e.g. /gun/particle pi-");

    energy = new G4UIcmdWithADoubleAndUnit("/gun/energy", this);
    energy->SetGuidance("Set kinetic energy, e.g. /gun/energy 20 GeV");
    energy->SetUnitCategory("Energy");
  }

  ~GunMessenger() override {
    delete energy;
    delete particle;
    delete dir;
  }

  void SetNewValue(G4UIcommand *command, G4String newValue) override {
    if (command == particle) {
      gun->SetParticle(newValue);
    } else if (command == energy) {
      gun->SetEnergy(energy->GetNewDoubleValue(newValue));
    }
  }

private:
  PrimaryGun *gun = nullptr;
  G4UIdirectory *dir = nullptr;
  G4UIcmdWithAString *particle = nullptr;
  G4UIcmdWithADoubleAndUnit *energy = nullptr;
};

class MyRunAction : public G4UserRunAction {
public:
  void BeginOfRunAction(const G4Run *) override {
    AnalysisManager::Instance().BeginRun();
  }
  void EndOfRunAction(const G4Run *) override {
    AnalysisManager::Instance().EndRun();
  }
};

class MyEventAction : public G4UserEventAction {
public:
  explicit MyEventAction(const PrimaryGun *pg) : gun(pg) {}

  void BeginOfEventAction(const G4Event *) override {
    const double eGeV = gun ? (gun->GetEnergy() / GeV) : 0.0;
    AnalysisManager::Instance().BeginEvent(eGeV);
  }

  void EndOfEventAction(const G4Event *) override {
    AnalysisManager::Instance().EndEvent();
  }

private:
  const PrimaryGun *gun = nullptr;
};

int main(int argc, char **argv) {
  std::string macro;
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "-m" && i + 1 < argc)
      macro = argv[++i];
  }

  // set random seed
  G4Random::setTheSeed(42);

  auto *runManager = new G4RunManager;
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new FTFP_BERT);

  auto *gun = new PrimaryGun();
  runManager->SetUserAction(static_cast<G4VUserPrimaryGeneratorAction *>(gun));
  runManager->SetUserAction(new MyRunAction());
  runManager->SetUserAction(new MyEventAction(gun));

  AnaMessenger anaMessenger;
  GunMessenger gunMessenger(gun);

  auto *visManager = new G4VisExecutive();
  visManager->Initialize();

  auto *ui = G4UImanager::GetUIpointer();
  if (!macro.empty()) {
    ui->ApplyCommand("/control/execute " + G4String(macro));
  } else {
    G4UIExecutive uiexec(argc, argv);
    uiexec.SessionStart();
  }

  delete visManager;
  delete runManager;
  return 0;
}
