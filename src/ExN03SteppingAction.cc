#include "ExN03SteppingAction.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03EventAction.hh"
#include "ExN03RunAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4TrackVector.hh"

extern TFile *outhe3; 
extern TTree *the3;
extern G4int npulses;
extern std::vector<double> Ecounters,tcounters;
extern std::vector<int> idcounters;
//G4int npulses[4];
//G4double x0he3[3];
//G4double Ecounters[4][10];
//G4double tcounters[4][10];
extern G4double lasthittime[4],firsthittime[4];
extern G4int nhe3primary;
extern std::vector<double> E0he3;
extern std::vector<double> x0he3;
extern std::vector<double> y0he3;
extern std::vector<double> z0he3;
extern std::vector<double> deltathe3;
extern G4double timecut;  //If a capture in a detector occurs within timecut ns of the previous capture, it is treated as contributing to the same pulse


ExN03SteppingAction::ExN03SteppingAction(ExN03DetectorConstruction* det,
                                         ExN03EventAction* evt)
:detector(det), eventaction(evt)					 
{ }



ExN03SteppingAction::~ExN03SteppingAction()
{ }



void ExN03SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //We will use pre step point object for coordinates, time, and volume
  //We will use aStep object to get energy deposit
  G4StepPoint *start = aStep->GetPreStepPoint();
  G4StepPoint *end = aStep->GetPostStepPoint();
  G4VPhysicalVolume *volume = start->GetTouchableHandle()->GetVolume();
  G4Track *track = aStep->GetTrack();
  G4String pN = track->GetDefinition()->GetParticleName();
  G4double edep = aStep->GetTotalEnergyDeposit()/CLHEP::keV;
  G4double time=start->GetGlobalTime()/CLHEP::ns;
  //Fill quantities, assuming the only tracking is not going on simultaneously in two or more counters
  G4int dum1, dum2;
  G4double tmp1, tmp2;
  if(edep > 0){
    //    G4cout << "edep, volume: " << edep << " " << volume->GetName() << G4endl;
    if(volume == detector->GetCounter1()){
	if(time-lasthittime[0] > timecut){
	  npulses=npulses+1;
	  Ecounters.push_back(edep);
	  tcounters.push_back(time);
	  deltathe3.push_back(0.);
	  idcounters.push_back(1);
	  firsthittime[0]=time;
	  lasthittime[0]=time;
	}
	else{
	  Ecounters[npulses-1]=Ecounters[npulses-1]+edep;
	  tcounters[npulses-1]=time;
	  lasthittime[0]=time;
	  deltathe3[npulses-1]=time-firsthittime[0];
	}
      }
    if(volume == detector->GetCounter2()){
	if(time-lasthittime[1] > timecut){
	  npulses=npulses+1;
	  Ecounters.push_back(edep);
	  tcounters.push_back(time);
	  idcounters.push_back(2);
	  deltathe3.push_back(0.);
	  firsthittime[1]=time;
	  lasthittime[1]=time;
	}
	else{
	  Ecounters[npulses-1]=Ecounters[npulses-1]+edep;
	  tcounters[npulses-1]=time;
	  lasthittime[1]=time;
	  deltathe3[npulses-1]=time-firsthittime[1];
	}
      }
    if(volume == detector->GetCounter3()){
	if(time-lasthittime[2] > timecut){
	  npulses=npulses+1;
	  Ecounters.push_back(edep);
	  tcounters.push_back(time);
	  idcounters.push_back(3);
	  deltathe3.push_back(0.);
	  firsthittime[2]=time;
	  lasthittime[2]=time;
	}
	else{
	  Ecounters[npulses-1]=Ecounters[npulses-1]+edep;
	  tcounters[npulses-1]=time;
	  lasthittime[2]=time;
	  deltathe3[npulses-1]=time-firsthittime[2];
	}
      }
    if(volume == detector->GetCounter4()){
	if(time-lasthittime[3] > timecut){
	  npulses=npulses+1;
	  Ecounters.push_back(edep);
	  tcounters.push_back(time);
	  idcounters.push_back(4);
	  deltathe3.push_back(0.);
	  firsthittime[3]=time;
	  lasthittime[3]=time;
	}
	else{
	  Ecounters[npulses-1]=Ecounters[npulses-1]+edep;
	  tcounters[npulses-1]=time;
	  lasthittime[3]=time;
	  deltathe3[npulses-1]=time-firsthittime[3];
	}
      }
  }
  //register primary neutrons here
  if(track->GetParentID()==0) {//primary
    if(pN == "neutron"){
      if(track->GetCurrentStepNumber()==1) {
    G4double nke= (track->GetVertexKineticEnergy())/CLHEP::MeV;
    x0he3.push_back(track->GetVertexPosition().x()/CLHEP::cm);
    y0he3.push_back(track->GetVertexPosition().y()/CLHEP::cm);
    z0he3.push_back(track->GetVertexPosition().z()/CLHEP::cm);
    nhe3primary++;
    E0he3.push_back(nke);
    //    G4cout << "primary neutron generated in volume " << volume->GetName() << G4endl;
      }
    }
  }
}
  

