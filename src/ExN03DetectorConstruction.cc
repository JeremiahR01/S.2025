#include "ExN03DetectorConstruction.hh"

#include <fstream>
using namespace std;

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Tubs.hh"

#include "G4Transform3D.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include <vector>

G4double zsource;
G4double timecut;

ExN03DetectorConstruction::ExN03DetectorConstruction()
{
 world_phys=0;
 box_phys=0;
 crete_phys1a=0;
 crete_phys1b=0;
 crete_phys2a=0;
 crete_phys2b=0;
 crete_phys3=0;
 for(int i=0;i<4;++i) {
  holes_phys[i]=0;
  tubes_phys[i]=0;
  counters_phys[i]=0;
 };
//World should be at least bigger than the system
 worldsize=3*CLHEP::m;
}

ExN03DetectorConstruction::~ExN03DetectorConstruction() {}



G4VPhysicalVolume* ExN03DetectorConstruction::Construct()
{
  return ConstructCounter();
}



G4VPhysicalVolume* ExN03DetectorConstruction::ConstructCounter()
{

  // Clean old geometry, if any
  //
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

  // Import detector geometry and fill gas specifications
  ///////////////////////////////////////////////////////////
    ifstream geometry("geometry.dat");
    G4double boxh,countr,counth,countd,boxs,He3Pressure,CounterVolume,He3Moles,He3Density,R,pi;

    R = 0.08206;
    pi = 3.141596;

    G4double depthsource,detectorsetting,dumtimecut;
    geometry >> depthsource; // distance of source center from top of source tube
    geometry >> detectorsetting;  // distance of detector from edge of flange
    geometry >> dumtimecut; // Successive captures in the same counter within dumtimecut ns are treated as contributing to the same pulse

    G4cout << " Distance (in) of source center from top of source tube: " << depthsource << G4endl;
    G4cout << " Distance (in) of edge of detector holder from edge of center flange: " << detectorsetting << G4endl;
    G4cout << " Time cut (ns) for pile-up: " << dumtimecut << G4endl;
    timecut=dumtimecut;

    //Water parameters (cm)
    G4double WaterD, WaterH, WaterVolume;

    WaterD = (59 / 2) * 2.54 * CLHEP::cm;  // Diameter of water volume
    WaterH = (59 / 2) * 2.54 * CLHEP::cm;  // Height of water volume at nominal filling height

    G4double LeadD, LeadL;
    LeadD=(1.00/2)*2.54*CLHEP::cm;  //diameter  of lead weights at the bottom of each detector holder
    LeadL=(2.00/2)*2.54*CLHEP::cm;  //length of lead weights at the bottom of each detector holder

  //Detector parameters (cm)
    G4double SSTubeID, SSTubeOD, DetectorHolderL, DetectorWireL, DetectorTubeL, DetectorTubeD,DetectorHe3L, DetectorHe3D, DetectorHe3L_dead_lower, DetectorHe3L_dead_upper, DetectorwhiteElecL, DetectornonwhiteElecL;

    SSTubeID=0.5*(2.69875)*CLHEP::cm;  // ID of SS tubes used for detector holders and source holder
    SSTubeOD=1.5875 * CLHEP::cm; // OD of SS tubes used for detector holders and source holder
    DetectorHolderL = (31 / 2) * 2.54 * CLHEP::cm; //Length of Detector Holder tube
    DetectorWireL = (42.2656 / 2) * CLHEP::cm; //Length of Detector inner wire
    DetectorHe3L = (13.00 / 2) * 2.54 * CLHEP::cm; //Length of the cylinder of active he3 gas in each detector
    DetectorHe3D = ((1.00-0.032) / 2) * 2.54 * CLHEP::cm; //Diameter of the cylinder of active or dead he3 gas in each detector
    DetectorHe3L_dead_lower = (0.95 / 2) * 2.54 * CLHEP::cm; //Length of the cylinder of dead he3 gas at lower end of detector
    DetectorHe3L_dead_upper = (1.25 / 2) * 2.54 * CLHEP::cm; //Length of the cylinder of dead he3 gas at upper end of  detector
    DetectorTubeL = ((16.64-1.0) / 2) * 2.54 * CLHEP::cm; //Length of the detector tube minus overlap with electronics
    DetectornonwhiteElecL = (6.6675 / 2) * CLHEP::cm; //Length of the electronics for non-white colored detectors
    DetectorwhiteElecL = (9.525 / 2) * CLHEP::cm; // Length of the electronics for the white colored detector

    //Source holder tube parameters; OD and ID are same as for detector holder tubes
    G4double SourceHolderTubeL, SourceHolderTubeID, SourceHolderTubeOD, SourceHolderL, SourceHolderOD, SourceHolderAirL, SourceHolderAirD, SourceHolderBottomWallThickness;
    SourceHolderTubeL=0.5*36.0*2.54*CLHEP::cm; 
    SourceHolderTubeOD=SSTubeOD;
    SourceHolderTubeID=SSTubeID;
    SourceHolderL=(2.0/2)*2.54*CLHEP::cm; //outer length of PET source holder
    SourceHolderOD=(1.0/2)*2.54*CLHEP::cm; // outer diameer of PET source holder
    SourceHolderBottomWallThickness=0.25*2.54*CLHEP::cm; //thickness of bottom wall of source holder
    SourceHolderAirL=(1.0/2)*2.54*CLHEP::cm; // length of air cavity inside source holder
    SourceHolderAirD=(0.625/2)*2.54*CLHEP::cm; //diameter of air cavity inside source holder
    
    //Offset in z from top of water to tops of source and detector holders 
    G4double HolderZOffset;
    HolderZOffset=4.625*2.54*CLHEP::cm;
    
  // Calculate necessary parameters -- units are Liters and Atmospheres
  ////////////////////////////////////////////////////////////////////////
  //Detector Distances
    G4double DetSet,Detdist,DetHolderheight;
    DetSet = detectorsetting+4.25; // This is the distance setting in inches to the centers of the detectors (If offset is discovered, adjust to account for source flange)
    Detdist = (DetSet*2.54)*CLHEP::cm; // This is the real distance from the center of the tank in cm
    DetHolderheight = WaterH - DetectorHolderL + HolderZOffset; // This is the height difference from center of water to center of the detector holders.  It depends on the actual height of the water.  Here we use that tops of detector holder tubes and source holder tube are 4.625 inches above the top of the water.
  //Water Construction
    WaterVolume = pi*pow(WaterD,2)*2*WaterH/1000; // V in L

  //jb-20160519  CounterVolume = pi*pow(countr,2)*counth*1000;
  //make materials
  ///////////////////////////////////////////////////////////
    G4int nelements, natoms, ncomponents;
    G4double density,a,z;
    G4String name,symbol;


  G4NistManager* man = G4NistManager::Instance();
  //Air
  G4Material* air = man->FindOrBuildMaterial("G4_AIR");

  //Lead
  G4Material* lead = man->FindOrBuildMaterial("G4_Pb");
  

 
  //Aluminum
  G4Element* Al = new G4Element(name="Aluminum",symbol="Cu",z=13,a=26.982*CLHEP::g/CLHEP::mole);
  G4Material* Aluminum = new G4Material(name="Aluminum",density=2.7*CLHEP::g/CLHEP::cm3,nelements=1);
  Aluminum->AddElement(Al,100*CLHEP::perCent);

  //Water
  //G4Material* h2o =  man->FindOrBuildMaterial("G4_WATER");


  //Define water as needed for more precise treatment of interactions of thermal neutrons
  G4Element* elTSHW = new G4Element("TS_H_of_Water", "H_WATER", 1.0, 1.0079*CLHEP::g/CLHEP::mole);
  G4Element* elO = new G4Element(name="Oxygen",symbol="O",z=8,a=15.999*CLHEP::g/CLHEP::mole);
  G4State phase;
  G4double temperature, pressure;
  G4Material* h2o_TS = new G4Material("Water_TS", density=1.0*CLHEP::g/CLHEP::cm3, ncomponents=2,phase=kStateLiquid,temperature=293.15*CLHEP::kelvin, pressure=1.0*CLHEP::atmosphere);
  h2o_TS->AddElement(elTSHW,natoms=2);
  h2o_TS->AddElement(elO,natoms=1);


 //PET
  //  G4Material* pet = man->FindOrBuildMaterial("G4_POLYETHYLENE");


  
  //Define polyethylene as needed for more precise treatment of interactions of thermal neutrons
  G4Element *elTSHpet = new G4Element("TS_H_of_Polyethylene", "H_pet", 1., 1.0079*CLHEP::g/CLHEP::mole);
  G4Element *elC = new G4Element(name="Carbon",symbol="C",z=6,a=12.011*CLHEP::g/CLHEP::mole);
  G4Material *pet_TS = new G4Material("pet_TS",density=0.94*CLHEP::g/CLHEP::cm3, ncomponents=2,phase=kStateSolid,temperature=293.15*CLHEP::kelvin, pressure=1.0*CLHEP::atmosphere);
  pet_TS->AddElement(elTSHpet,natoms=2);
  pet_TS->AddElement(elC,natoms=1);

  //Brass
  G4Element* Cu = new G4Element(name="Copper",symbol="Cu",z=29,a=63.546*CLHEP::g/CLHEP::mole);
  G4Element* Zn = new G4Element(name="Zinc",symbol="Zn",z=30,a=65.38*CLHEP::g/CLHEP::mole);
  G4Material* Brass = new G4Material(name="Brass",density=8.5*CLHEP::g/CLHEP::cm3,nelements=2);
  Brass->AddElement(Cu,70*CLHEP::perCent);
  Brass->AddElement(Zn,30*CLHEP::perCent);

  //PTFE
  a=3.0160293*CLHEP::g/CLHEP::mole;
  G4Element* C = new G4Element(name="Carbon",symbol="C",z=6,a=12.011*CLHEP::g/CLHEP::mole);
  G4Element* F = new G4Element(name="Fluorine",symbol="F",z=9,a=18.998*CLHEP::g/CLHEP::mole);
  G4Material* PTFE = new G4Material(name="PTFE",density=2.20*CLHEP::g/CLHEP::cm3,nelements=2);
  PTFE->AddElement(C,24.0183*CLHEP::perCent);
  PTFE->AddElement(F,75.9817*CLHEP::perCent);

//Counter gas
  // He3
  a=3.0160293*CLHEP::g/CLHEP::mole;
  G4Isotope* He3 = new G4Isotope(name="He3",2,3,a);
  G4Element* myHe = new G4Element(name="myHe",symbol="myHe",nelements=1);
  myHe->AddIsotope(He3,100.*CLHEP::perCent);
  //Ar
  G4Element *Ar = new G4Element(name="Argon",symbol="Ar",z=18,a=35.95*CLHEP::g/CLHEP::mole);
  //Carbon dioxide
    G4Material *carbondioxide = man->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  //
    G4Material* CounterGas = new G4Material(name="CounterGas",density=0.004325*CLHEP::g/CLHEP::cm3,nelements=3,kStateGas,293.*CLHEP::kelvin,(170.5/14.7)*CLHEP::atmosphere);
  CounterGas->AddElement(myHe,28.3*CLHEP::perCent);
  CounterGas->AddElement(Ar,66.8*CLHEP::perCent);
  CounterGas->AddMaterial(carbondioxide,4.9*CLHEP::perCent);

  //SS for the tube walls
  G4Element* Fe = new G4Element(name="Iron",symbol="Fe",z=26,a=55.845*CLHEP::g/CLHEP::mole);
  G4Element* Ni = new G4Element(name="Nickel",symbol="Ni",28,58.693*CLHEP::g/CLHEP::mole);
  G4Element* Cr = new G4Element(name="Chromium",symbol="Cr",24,51.996*CLHEP::g/CLHEP::mole);
  G4Material* SS = new G4Material(name="Stainless",density=7.87*CLHEP::g/CLHEP::cm3,nelements=3);
  SS->AddElement(Fe,71.*CLHEP::perCent);
  SS->AddElement(Cr,19.*CLHEP::perCent);
  SS->AddElement(Ni,10.*CLHEP::perCent);

  /*Electronics are composed of:
  Hydrogen
  Carbon
  Copper

  The ratios of Total Guesses are provided as percentages
  */
  G4Element* H = new G4Element(name="Hydrogen",symbol="H",z=1,a=1.008*CLHEP::g/CLHEP::mole);
  G4Material* Elec = new G4Material(name="Electronics",density=3.00*CLHEP::g/CLHEP::cm3,nelements=3);
  Elec->AddElement(H,25.*CLHEP::perCent);
  Elec->AddElement(C,54.*CLHEP::perCent);
  Elec->AddElement(Cu,21.*CLHEP::perCent);
  /////////////////////////////////////////////////////////////


  //make objects
  //////////////////////////////////////////////////////////////

  //world

    G4Box* world_solid = new G4Box("world_solid",worldsize,worldsize,worldsize);
    G4LogicalVolume* world_log = new G4LogicalVolume(world_solid,air,"world_log");
    world_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),world_log,"world_phys",0,false,0);

  //Water
    G4VSolid* Water_solid = new G4Tubs("Water_solid", 0 * CLHEP::cm, WaterD, WaterH, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* Water_log = new G4LogicalVolume(Water_solid, h2o_TS, "Water_log");
    //    G4LogicalVolume* Water_log = new G4LogicalVolume(Water_solid, h2o, "Water_log");
    Water_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Water_log, "Water_phys", world_log, 0, false, 0);

          // Inner Radii set to 0 to allow daughter volumes to exist
          // All positions and dimensions are relative to object center

    //Detector holder 1
    G4VSolid* DetHolder_solid1 = new G4Tubs("DetHolder_solid1", 0 * CLHEP::cm, SSTubeOD, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_log1 = new G4LogicalVolume(DetHolder_solid1, SS, "DetHolder_log1");
    DetHolder_phys1 = new G4PVPlacement(0, G4ThreeVector(0, Detdist, DetHolderheight), DetHolder_log1, "DetHolder_phys1", Water_log, 0, false, 0);

    //Detector holder 2
    G4VSolid* DetHolder_solid2 = new G4Tubs("DetHolder_solid2", 0 * CLHEP::cm, SSTubeOD, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_log2 = new G4LogicalVolume(DetHolder_solid2, SS, "DetHolder_log2");
    DetHolder_phys2 = new G4PVPlacement(0, G4ThreeVector(Detdist, 0, DetHolderheight), DetHolder_log2, "DetHolder_phys2", Water_log, 0, false, 0);

    //Detector holder 3
    G4VSolid* DetHolder_solid3 = new G4Tubs("DetHolder_solid3", 0 * CLHEP::cm, SSTubeOD, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_log3 = new G4LogicalVolume(DetHolder_solid3, SS, "DetHolder_log3");
    DetHolder_phys3 = new G4PVPlacement(0, G4ThreeVector(0, -Detdist, DetHolderheight), DetHolder_log3, "DetHolder_phys3", Water_log, 0, false, 0);

    //Detector holder 4
    G4VSolid* DetHolder_solid4 = new G4Tubs("DetHolder_solid4", 0 * CLHEP::cm, SSTubeOD, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_log4 = new G4LogicalVolume(DetHolder_solid4, SS, "DetHolder_log4");
    DetHolder_phys4 = new G4PVPlacement(0, G4ThreeVector(-Detdist, 0, DetHolderheight), DetHolder_log4, "DetHolder_phys4", Water_log, 0, false, 0);

    // Construct air volume to make detector holders hollow tubes with ID 1.0625 inches (ignores steel bottom of tube and flanges at top of tube)
      
    //Air Volume-Detector holder 1
    G4VSolid* DetHolder_air_solid1 = new G4Tubs("DetHolder_air_solid1", 0 * CLHEP::cm, SSTubeID, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_air_log1 = new G4LogicalVolume(DetHolder_air_solid1, air, "DetHolder_air_log1");
    DetHolder_air_phys1 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), DetHolder_air_log1, "DetHolder_air_phys1", DetHolder_log1, 0, false, 0);

    //Air Volume-Detector holder 2
    G4VSolid* DetHolder_air_solid2 = new G4Tubs("DetHolder_air_solid2", 0 * CLHEP::cm, SSTubeID, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_air_log2 = new G4LogicalVolume(DetHolder_air_solid2, air, "DetHolder_air_log2");
    DetHolder_air_phys2 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), DetHolder_air_log2, "DetHolder_air_phys2", DetHolder_log2, 0, false, 0);

    //Air Volume-Detector holder 3
    G4VSolid* DetHolder_air_solid3 = new G4Tubs("DetHolder_air_solid3", 0 * CLHEP::cm, SSTubeID, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_air_log3 = new G4LogicalVolume(DetHolder_air_solid3, air, "DetHolder_air_log3");
    DetHolder_air_phys3 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), DetHolder_air_log3, "DetHolder_air_phys3", DetHolder_log3, 0, false, 0);


    //Air Volume-Detector holder 4
    G4VSolid* DetHolder_air_solid4 = new G4Tubs("DetHolder_air_solid4", 0 * CLHEP::cm, SSTubeID, DetectorHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHolder_air_log4 = new G4LogicalVolume(DetHolder_air_solid4, air, "DetHolder_air_log4");
    DetHolder_air_phys4 = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), DetHolder_air_log4, "DetHolder_air_phys4", DetHolder_log4, 0, false, 0);


    //Lead weight
    G4VSolid* DetLead_solid = new G4Tubs("DetLead_solid", 0 * CLHEP::cm, LeadD, LeadL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetLead_log = new G4LogicalVolume(DetLead_solid, lead, "DetLead_log");

    //Detector tubes
    G4VSolid* DetTube_solid1 = new G4Tubs("DetTube_solid1", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectorTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetTube_log1 = new G4LogicalVolume(DetTube_solid1, Aluminum, "DetTube_log1");

    G4VSolid* DetTube_solid2 = new G4Tubs("DetTube_solid2", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectorTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetTube_log2 = new G4LogicalVolume(DetTube_solid2, Aluminum, "DetTube_log2");
    
    G4VSolid* DetTube_solid3 = new G4Tubs("DetTube_solid3", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectorTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetTube_log3 = new G4LogicalVolume(DetTube_solid3, Aluminum, "DetTube_log3");

    G4VSolid* DetTube_solid4 = new G4Tubs("DetTube_solid4", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectorTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetTube_log4 = new G4LogicalVolume(DetTube_solid4, Aluminum, "DetTube_log4");

    //Nonwhite electronics (does not contain low voltage to high voltage converter)

    G4VSolid* DetElec_nonwhite_solid = new G4Tubs("DetElec_nonwhite_solid", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectornonwhiteElecL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetElec_nonwhite_log = new G4LogicalVolume(DetElec_nonwhite_solid, Elec, "DetElec_nonwhite_log");

    //White electronics (contains low voltage to high voltge converter)
    G4VSolid* White_DetElec_solid = new G4Tubs("White_DetElec_solid", 0 * CLHEP::cm, 1.27 * CLHEP::cm, DetectorwhiteElecL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* White_DetElec_log = new G4LogicalVolume(White_DetElec_solid, Elec, "White_DetElec_log");

    //He-3 regions
    G4VSolid* DetHe3_solid = new G4Tubs("DetHe3_solid", 0 * CLHEP::cm, DetectorHe3D, DetectorHe3L, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHe3_log = new G4LogicalVolume(DetHe3_solid, CounterGas, "DetHe3_log");

    G4VSolid* DetHe3_dead_lower_solid = new G4Tubs("DetHe3_dead_lower_solid", 0 * CLHEP::cm, DetectorHe3D, DetectorHe3L_dead_lower, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHe3_dead_lower_log = new G4LogicalVolume(DetHe3_dead_lower_solid, CounterGas, "DetHe3_dead_lower_log");

    G4VSolid* DetHe3_dead_upper_solid = new G4Tubs("DetHe3_dead_upper_solid", 0 * CLHEP::cm, DetectorHe3D, DetectorHe3L_dead_upper, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetHe3_dead_upper_log = new G4LogicalVolume(DetHe3_dead_upper_solid, CounterGas, "DetHe3_dead_upper_log");

    //Detector wire-not placed for now
    G4VSolid* DetWire_solid = new G4Tubs("DetWire_solid1", 0 * CLHEP::cm, 0.0025 * CLHEP::cm, DetectorWireL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* DetWire_log = new G4LogicalVolume(DetWire_solid, SS, "DetWire_log");

    //Stack lead, detector tube, and electronics in each detector holder and stack counter gas regions in detector tube; neglect wire for now
    //Non-white detector 1
    DetLead_phys1 = new G4PVPlacement(0,G4ThreeVector(0,0,-DetectorHolderL+LeadL),DetLead_log,"DetLead_phys1",DetHolder_air_log1,false,0,false);
    DetTube_phys1 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+DetectorTubeL), DetTube_log1, "DetTube_phys1", DetHolder_air_log1, false, 0, false);
    DetElec_phys1 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+2.0*DetectorTubeL+DetectornonwhiteElecL), DetElec_nonwhite_log, "DetElec_phys1", DetHolder_air_log1, false,0,false);
    DetHe3_dead_lower_phys1 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+DetectorHe3L_dead_lower),DetHe3_dead_lower_log,"DetHe_dead_lower_phys1",DetTube_log1,false,0,false);
    DetHe3_phys1 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+DetectorHe3L),DetHe3_log,"DetHe3_phys1",DetTube_log1,false,0,false);
    DetHe3_dead_upper_phys1 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+2.0*DetectorHe3L+DetectorHe3L_dead_upper),DetHe3_dead_upper_log,"DetHe3_dead_upper_phys1",DetTube_log1,false,0,false);
    //Non-white detector 2
    DetLead_phys2 = new G4PVPlacement(0,G4ThreeVector(0,0,-DetectorHolderL+LeadL),DetLead_log,"DetLead_phys1",DetHolder_air_log2,false,0,false);
    DetTube_phys2 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+DetectorTubeL), DetTube_log2, "DetTube_phys2", DetHolder_air_log2, false, 0, false);
    DetElec_phys2 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+2.0*DetectorTubeL+DetectornonwhiteElecL), DetElec_nonwhite_log, "DetElec_phys2", DetHolder_air_log2, false,0,false);
    DetHe3_dead_lower_phys2 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+DetectorHe3L_dead_lower),DetHe3_dead_lower_log,"DetHe_dead_lower_phys2",DetTube_log2,false,0,false);
    DetHe3_phys2 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+DetectorHe3L),DetHe3_log,"DetHe3_phys2",DetTube_log2,false,0,false);
    DetHe3_dead_upper_phys2 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+2.0*DetectorHe3L+DetectorHe3L_dead_upper),DetHe3_dead_upper_log,"DetHe3_dead_upper_phys2",DetTube_log2,false,0,false);

    //Non-white detector 3
    DetLead_phys3 = new G4PVPlacement(0,G4ThreeVector(0,0,-DetectorHolderL+LeadL),DetLead_log,"DetLead_phys3",DetHolder_air_log3,false,0,false);
    DetTube_phys3 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+DetectorTubeL), DetTube_log3, "DetTube_phys3", DetHolder_air_log3, false, 0, false);
    DetElec_phys3 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+2.0*DetectorTubeL+DetectornonwhiteElecL), DetElec_nonwhite_log, "DetElec_phys3", DetHolder_air_log3, false,0,false);
    DetHe3_dead_lower_phys3 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+DetectorHe3L_dead_lower),DetHe3_dead_lower_log,"DetHe_dead_lower_phys3",DetTube_log3,false,0,false);
    DetHe3_phys3 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+DetectorHe3L),DetHe3_log,"DetHe3_phys3",DetTube_log3,false,0,false);
    DetHe3_dead_upper_phys3 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+2.0*DetectorHe3L+DetectorHe3L_dead_upper),DetHe3_dead_upper_log,"DetHe3_dead_upper_phys3",DetTube_log3,false,0,false);

    //White detector
    DetLead_phys4 = new G4PVPlacement(0,G4ThreeVector(0,0,-DetectorHolderL+LeadL),DetLead_log,"DetLead_phys4",DetHolder_air_log4,false,0,false);
    DetTube_phys4 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+DetectorTubeL), DetTube_log4, "DetTube_phys4", DetHolder_air_log4, false, 0, false);
    DetElec_phys4 = new G4PVPlacement(0,G4ThreeVector(0, 0, -DetectorHolderL+2.0*LeadL+2.0*DetectorTubeL+DetectorwhiteElecL), White_DetElec_log, "DetElec_phys4", DetHolder_air_log4, false, 0, false);
    DetHe3_dead_lower_phys4 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+DetectorHe3L_dead_lower),DetHe3_dead_lower_log,"DetHe_dead_lower_phys4",DetTube_log4, false,0, false);
    DetHe3_phys4 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+DetectorHe3L),DetHe3_log,"DetHe3_phys4",DetTube_log4,false,0, false);
    DetHe3_dead_upper_phys4 = new G4PVPlacement(0, G4ThreeVector(0,0,-DetectorTubeL+2.0*DetectorHe3L_dead_lower+2.0*DetectorHe3L+DetectorHe3L_dead_upper),DetHe3_dead_upper_log,"DetHe3_dead_upper_phys4",DetTube_log4,false,0,false);


    //Source holder tube with two lead weights at the bottom and PET source holder 
    //SS steel cylinder
    G4VSolid* SourceHolderTube_solid = new G4Tubs("SourceHolderTube_solid",  0, SourceHolderTubeOD, SourceHolderTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* SourceHolderTube_log = new G4LogicalVolume(SourceHolderTube_solid, SS, "SourceHolderTube_log");
    G4double SourceHolderTubeheight=WaterH+HolderZOffset-SourceHolderTubeL;
    SourceHolderTube_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, SourceHolderTubeheight), SourceHolderTube_log, "SourceHolderTube_phys", Water_log, false,0,false);
    //Source holder tube air
    G4VSolid* SourceHolderTube_air_solid = new G4Tubs("SourceHolderTube_air_solid",  0, SourceHolderTubeID, SourceHolderTubeL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* SourceHolderTube_air_log = new G4LogicalVolume(SourceHolderTube_air_solid, air, "SourceHolderTube_air_log");
    SourceHolderTube_air_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), SourceHolderTube_air_log, "SourceHolderTube_air_phys", SourceHolderTube_log, false,0,false);
    //Source holder weights
    SourceHolderTube_lead1_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-SourceHolderTubeL+LeadL),DetLead_log,"SourceHolderTube_lead1_phys",SourceHolderTube_air_log,0,false,0);
    SourceHolderTube_lead2_phys = new G4PVPlacement(0,G4ThreeVector(0,0,-SourceHolderTubeL+2.0*LeadL+LeadL),DetLead_log,"SourceHolderTube_lead2_phys",SourceHolderTube_air_log,0,false,0);
    //Source holder PET
    G4VSolid* SourceHolder_pet_solid = new G4Tubs("SourceHolder_pet_solid",  0, SourceHolderOD, SourceHolderL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* SourceHolder_pet_log = new G4LogicalVolume(SourceHolder_pet_solid, pet_TS, "SourceHolder_pet_log");
    // G4LogicalVolume* SourceHolder_pet_log = new G4LogicalVolume(SourceHolder_pet_solid, pet, "SourceHolder_pet_log");

    //Source holder air
    G4VSolid* SourceHolder_air_solid = new G4Tubs("SourceHolder_air_solid",  0, SourceHolderAirD, SourceHolderAirL, 0. * CLHEP::deg, 360. * CLHEP::deg);
    G4LogicalVolume* SourceHolder_air_log = new G4LogicalVolume(SourceHolder_air_solid, air, "SourceHolder_air_log");

    //place source holder, with center of air cavity being source activity center

    SourceHolder_pet_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, SourceHolderTubeL-depthsource*2.54*CLHEP::cm-SourceHolderAirL-SourceHolderBottomWallThickness+SourceHolderL), SourceHolder_pet_log, "SourceHolder_pet_phys", SourceHolderTube_air_log, false,0,false);
    SourceHolder_air_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -SourceHolderL+SourceHolderBottomWallThickness+SourceHolderAirL), SourceHolder_air_log, "SourceHolder_air_phys", SourceHolder_pet_log, false,0,false);
    
    //Set z of source center in world coordinates
    //Input is distance of source in inches from top of source holder tube
    zsource=WaterH+HolderZOffset-depthsource*2.54*CLHEP::cm;
    G4cout << "Source z position in world coordinates (mm): " << zsource << G4endl;
    G4double zhe3active=DetHolderheight-DetectorHolderL+2.0*LeadL+DetectorTubeL-DetectorTubeL+2.0*DetectorHe3L_dead_lower+DetectorHe3L;
    G4cout << "z of center of active He-3 region in world coordinates (mm): " << zhe3active << G4endl;
    G4double zairsource = SourceHolderTubeheight+SourceHolderTubeL-depthsource*2.54*CLHEP::cm-SourceHolderAirL-SourceHolderBottomWallThickness+SourceHolderL-SourceHolderL+SourceHolderBottomWallThickness+SourceHolderAirL;
    G4cout << "z of center of source air volume in world coordinates (mm): " << zairsource << G4endl;



// jkb-20200704-commented followiing code for source tube pole and capsule, including visualization setup
//    Source tube pole
//    G4VSolid* SourcePole_solid = new G4Tubs("SourcePole_solid", 0 * CLHEP::cm, 0.555625 * CLHEP::cm, 81.28 / 2 * CLHEP::cm, 0. * CLHEP::deg, 360. * CLHEP::deg);
//    G4LogicalVolume* SourcePole_log = new G4LogicalVolume(SourcePole_solid, SS, "SourcePole_log");
//    SourcePole_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 43.14), SourcePole_log, "SourcePole_phys", SourceHolderTube_log, 0, false, 0);



  //Source Capsule (Placed at center but insignificant. Change z value of G4ThreeVector if necessary)
//    G4VSolid* SourceCaps_solid = new G4Tubs("SourceCaps_solid", 0 * CLHEP::cm, 1.35 * CLHEP::cm, 2.5 * CLHEP::cm, 0. * CLHEP::deg, 360. * CLHEP::deg);
//    G4LogicalVolume* SourceCaps_log = new G4LogicalVolume(SourceCaps_solid, SS, "SourceCaps_log");
//    SourceCaps_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), SourceCaps_log, "SourceCaps_phys", SourceHolderTube_log, 0, false, 0);

    //stopping point Aug 2

  // Visualization attributes
   // world_log->SetVisAttributes (G4VisAttributes::Invisible);
    G4VisAttributes* waterVisAtt = new G4VisAttributes(G4Colour(0, 0, 1));
    waterVisAtt->SetForceAuxEdgeVisible(true);

    G4VisAttributes* cylinderVisAtt= new G4VisAttributes(G4Colour(0,1,0));
    cylinderVisAtt->SetForceAuxEdgeVisible(true);

    G4VisAttributes* counterVisAtt= new G4VisAttributes(G4Colour(1,0,0));
    counterVisAtt->SetForceSolid(true);

    G4VisAttributes* tubeVisAtt= new G4VisAttributes(G4Colour(1,1,1));
    tubeVisAtt->SetForceAuxEdgeVisible(true);

    G4VisAttributes* ElecVisAtt= new G4VisAttributes(G4Colour(1,1,0));
    ElecVisAtt->SetForceSolid(true);

    G4VisAttributes* SourceHolderTubeAirVisAtt = new G4VisAttributes(G4Colour::White());
    SourceHolderTubeAirVisAtt->SetForceAuxEdgeVisible(true);
    SourceHolderTube_air_log->SetVisAttributes(SourceHolderTubeAirVisAtt);

    G4VisAttributes* SourceHolderPETVisAtt = new G4VisAttributes(G4Colour::Brown());
    SourceHolderPETVisAtt->SetForceAuxEdgeVisible(true);
    SourceHolder_pet_log->SetVisAttributes(SourceHolderPETVisAtt);

    G4VisAttributes* SourceHolderAirVisAtt = new G4VisAttributes(G4Colour::Red());
    SourceHolderAirVisAtt->SetForceSolid(true);
    SourceHolder_air_log->SetVisAttributes(SourceHolderAirVisAtt);
    
    Water_log->SetVisAttributes(waterVisAtt);

    //    DetHolder_air_log1->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetHolder_air_log2->SetVisAttributes(G4VisAttributes::Invisible);
    //   DetHolder_air_log3->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetHolder_air_log4->SetVisAttributes(G4VisAttributes::Invisible);

    DetHolder_air_log1->SetVisAttributes(cylinderVisAtt);
    DetHolder_air_log2->SetVisAttributes(cylinderVisAtt);
    DetHolder_air_log3->SetVisAttributes(cylinderVisAtt);
    DetHolder_air_log4->SetVisAttributes(cylinderVisAtt);

    DetHolder_log1->SetVisAttributes(G4VisAttributes::Invisible);
    DetHolder_log2->SetVisAttributes(G4VisAttributes::Invisible);
    DetHolder_log3->SetVisAttributes(G4VisAttributes::Invisible);
    DetHolder_log4->SetVisAttributes(G4VisAttributes::Invisible);

    //    DetTube_log1->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetTube_log2->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetTube_log3->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetTube_log4->SetVisAttributes(G4VisAttributes::Invisible);

    DetTube_log1->SetVisAttributes(tubeVisAtt);
    DetTube_log2->SetVisAttributes(tubeVisAtt);
    DetTube_log3->SetVisAttributes(tubeVisAtt);
    DetTube_log4->SetVisAttributes(tubeVisAtt);

    DetElec_nonwhite_log->SetVisAttributes(ElecVisAtt);
    White_DetElec_log->SetVisAttributes(ElecVisAtt);

    DetHe3_log->SetVisAttributes(counterVisAtt);
    
    //    DetHe3_log_det1->SetVisAttributes(counterVisAtt);
    //    DetHe3_log_det2->SetVisAttributes(counterVisAtt);
    //    DetHe3_log_det3->SetVisAttributes(counterVisAtt);
    //    DetHe3_log_det4->SetVisAttributes(counterVisAtt);

    //    DetWire_log1->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetWire_log2->SetVisAttributes(G4VisAttributes::Invisible);
    //    DetWire_log3->SetVisAttributes(G4VisAttributes::Invisible);
    //    White_DetWire_log->SetVisAttributes(G4VisAttributes::Invisible);

    //    SourceHolderTube_log->SetVisAttributes(counterVisAtt);
//    SourcePole_log->SetVisAttributes(counterVisAtt);
//    SourceCaps_log->SetVisAttributes(counterVisAtt);
    //stopping point July 31


//List the materials defined
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world_phys;
}
