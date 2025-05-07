#include "ExN03RunAction.hh"
#include "ExN03EventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <vector>


#include <TFile.h>
#include <TTree.h>
#include <TString.h>

TFile *outhe3; 
TTree *the3;
G4int npulses;
std::vector<double> Ecounters,tcounters;
std::vector<int> idcounters;
//G4int npulses[4];
//G4double x0he3[3];
//G4double Ecounters[4][10];
//G4double tcounters[4][10];
G4double lasthittime[4],firsthittime[4];
G4int nhe3primary;
std::vector<double> E0he3;
std::vector<double> x0he3;
std::vector<double> y0he3;
std::vector<double> z0he3;
std::vector<double> deltathe3;


ExN03RunAction::ExN03RunAction()
{
 outhe3 = new TFile("outhe3.root","RECREATE");
 
 the3 = new TTree("the3","tree");
// t->Branch("ECyl",&ECyl,"ECyl/D");
// t->Branch("Eholes",Eholes,"Eholes[4]/D");
// t->Branch("Etubes",Etubes,"Eholes[4]/D");
 the3->Branch("Ecounters",&Ecounters);
 the3->Branch("npulses",&npulses,"npulses/I");
 the3->Branch("tcounters",&tcounters);
 the3->Branch("idcounters",&idcounters);
 // t->Branch("EventID",&EventID,"EventID/I");
 // t->Branch("procArr",&procArr,"procArr/I");
 // t->Branch("pnameArr",&pnameArr,"pnameArr/I");
 // t->Branch("EstepArr",&EstepArr,"EstepArr/I");
 // t->Branch("np",&np,"np/I");
 //  the3->Branch("x0he3",x0he3,"x0he3[3]/D");
 the3->Branch("E0he3",&E0he3);
 the3->Branch("nhe3primary",&nhe3primary,"nhe3primary/I");
 the3->Branch("x0he3",&x0he3);
 the3->Branch("y0he3",&y0he3);
 the3->Branch("z0he3",&z0he3);
 the3->Branch("deltathe3",&deltathe3);
}
  

ExN03RunAction::~ExN03RunAction()
{
 outhe3->cd();
 the3->Write();
 outhe3->Close();

}



void ExN03RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

}






void ExN03RunAction::EndOfRunAction(const G4Run* aRun)
{

 outhe3->cd();
 the3->Write();
 outhe3->Close();
}
