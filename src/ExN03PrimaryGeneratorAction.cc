#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03PrimaryGeneratorMessenger.hh"
#include "ExN03DetectorConstruction.hh"
#include <math.h>
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4HEPEvtInterface.hh"
#include "G4UnitsTable.hh"
extern G4double zsource;

G4VPrimaryGenerator* HEPEvt;

ExN03PrimaryGeneratorAction::ExN03PrimaryGeneratorAction(ExN03DetectorConstruction* ExN03DC):ExN03Detector(ExN03DC)
{
    HEPEvt = new G4HEPEvtInterface("generator.data");
    G4cout << "Opened generator data file " << G4endl;
}

ExN03PrimaryGeneratorAction::~ExN03PrimaryGeneratorAction()
{
    delete HEPEvt;
}

void ExN03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    HEPEvt->SetParticlePosition(G4ThreeVector(0.*CLHEP::cm,0.*CLHEP::cm,zsource));
    HEPEvt->GeneratePrimaryVertex(anEvent);
}
