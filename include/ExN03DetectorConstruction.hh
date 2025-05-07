
#ifndef ExN03DetectorConstruction_h
#define ExN03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;

//keele
class G4Box;

class ExN03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    ExN03DetectorConstruction();
   ~ExN03DetectorConstruction();

  public:
     
  
     
     G4VPhysicalVolume* Construct();
     
    //keele
     const G4VPhysicalVolume* GetBox() {return box_phys;};
     const G4VPhysicalVolume* GetWall1a() {return crete_phys1a;};    
     const G4VPhysicalVolume* GetWall1b() {return crete_phys1b;};    
     const G4VPhysicalVolume* GetWall2a() {return crete_phys2a;};    
     const G4VPhysicalVolume* GetWall2b() {return crete_phys2b;};    
     const G4VPhysicalVolume* GetWall3() {return crete_phys3;};    
     const G4VPhysicalVolume* GetphysiWorld() {return world_phys;};           
     const G4VPhysicalVolume* GetCylinder()   {return cylinder_phys;};
     const G4VPhysicalVolume* GetHole(int n)           {return holes_phys[n];};
     const G4VPhysicalVolume* GetTube(int n)           {return tubes_phys[n];};
     const G4VPhysicalVolume* GetCounter(int n)        {return counters_phys[n];};
     const G4VPhysicalVolume* GetCounter1()        {return DetHe3_phys1;};
     const G4VPhysicalVolume* GetCounter2()        {return DetHe3_phys2;};
     const G4VPhysicalVolume* GetCounter3()        {return DetHe3_phys3;};
     const G4VPhysicalVolume* GetCounter4()        {return DetHe3_phys4;};
  //     const G4VPhysicalVolume* GetCounter4()        {return White_DetHe3_phys;};
     
     G4double GetWorldSize() {return worldsize;};
    
      
  private:
          G4VPhysicalVolume *holes_phys[4], *tubes_phys[4], *counters_phys[4], *world_phys, *cylinder_phys, *box_phys, *crete_phys1a, *crete_phys1b, *crete_phys2a, *crete_phys2b, *crete_phys3;
          G4double worldsize;
          G4VPhysicalVolume *Water_phys, *DetHolder_phys1, *DetHolder_phys2, *DetHolder_phys3, *DetHolder_phys4;
          G4VPhysicalVolume *DetHolder_air_phys1, *DetHolder_air_phys2, *DetHolder_air_phys3, *DetHolder_air_phys4;
          G4VPhysicalVolume *DetLead_phys1, *DetTube_phys1, *DetElec_phys1, *DetHe3_dead_lower_phys1, *DetHe3_phys1, *DetHe3_dead_upper_phys1; 
          G4VPhysicalVolume *DetLead_phys2, *DetTube_phys2, *DetElec_phys2, *DetHe3_dead_lower_phys2, *DetHe3_phys2, *DetHe3_dead_upper_phys2; 
          G4VPhysicalVolume *DetLead_phys3, *DetTube_phys3, *DetElec_phys3, *DetHe3_dead_lower_phys3, *DetHe3_phys3, *DetHe3_dead_upper_phys3; 
          G4VPhysicalVolume *DetLead_phys4, *DetTube_phys4, *DetElec_phys4, *DetHe3_dead_lower_phys4, *DetHe3_phys4, *DetHe3_dead_upper_phys4;
  G4VPhysicalVolume *SourceHolderTube_phys, *SourceHolderTube_air_phys, *SourceHolderTube_lead1_phys, *SourceHolderTube_lead2_phys, *SourceHolder_pet_phys, *SourceHolder_air_phys;


          G4VPhysicalVolume* ConstructCounter();     
};


#endif

