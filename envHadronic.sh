#
# for neutronHP
#
unset G4NEUTRONHP_SKIP_MISSING_ISOTOPES 
unset G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE
unset G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION
unset G4NEUTRONHP_NELECT_DOPPLER
unset G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS

G4NEUTRONHP_SKIP_MISSING_ISOTOPES=1
export G4NEUTRONHP_SKIP_MISSING_ISOTOPES

G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE=1
export G4NEUTRONHP_DO_NOT_ADJUST_FINAL_STATE
   
#### G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1
#### export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION
  
#### G4NEUTRONHP_NELECT_DOPPLER=1
#### export G4NEUTRONHP_NELECT_DOPPLER
  
#### G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS=1
#### export G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS
#
# for Bertini cascade
#
unset G4CASCADE_USE_PRECOMPOUND  
#### G4CASCADE_USE_PRECOMPOUND=1
#### export G4CASCADE_USE_PRECOMPOUND             

env |grep G4NEUTRONHP
env |grep G4CASCADE
