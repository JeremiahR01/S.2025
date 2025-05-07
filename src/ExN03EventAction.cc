#include "ExN03EventAction.hh"

#include "ExN03RunAction.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>


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


ExN03EventAction::ExN03EventAction(ExN03RunAction* run)
:runAct(run)
{

}


ExN03EventAction::~ExN03EventAction()
{
}


void ExN03EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtn = evt->GetEventID();
 if(!(evtn%500))
  G4cout<<"EventID: "<<evtn<<G4endl;

 // ECyl = 0;
 // Evol.clear();
 // for(int i=0;i<4;++i) {//for every hole 
 //  for(int j=0;j<3;++j)//for every module = {air_hole,stainless_wall,he-3_counter}
 //   Evol.push_back(0);
 // }//redo with enums for clarity

// procArr.clear();
// pnameArr.clear();
// EstepArr.clear();
// np=0;
// for(int i=0;i<3;++i) {x0he3[i]=0;};
 npulses=0;
 E0he3.clear();
 Ecounters.clear();
 tcounters.clear();
 idcounters.clear();
 E0he3.clear();
 x0he3.clear();
 y0he3.clear();
 z0he3.clear();
 deltathe3.clear();
 nhe3primary=0;
 for(int i=0;i<4;i++){
   lasthittime[i]=-9999.;
   firsthittime[i]=0;
 }
}

void ExN03EventAction::EndOfEventAction(const G4Event* evt)
{

  G4int evtn = evt->GetEventID();
//  G4int nsec = (G4int)pnameArr.size();

//   for(int i=0;i<nsec;++i) 
//    if( (!strcmp(pnameArr.at(i),"triton"))&&(!strcmp(procArr.at(i),"NeutronInelastic"))) {
//     np=1;
//     break;
//    };
 
//  runAct->fillPerEvent(ECyl, Evol, evtn, procArr,pnameArr,EstepArr,np,x0,E0,nsec);
the3->Fill();
 if (nhe3primary < 1){
 G4cout << "Event : " << evtn << " Number of primary neutrons: " << nhe3primary << G4endl;
 }
}  

