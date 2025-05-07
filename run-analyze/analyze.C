#define analyze_cxx
// The class definition in analyze.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("analyze.C")
// root> T->Process("analyze.C","some options")
// root> T->Process("analyze.C+")
//


#include "analyze.h"
#include <TH2.h>
#include <TStyle.h>

// Define the histograms we want to fill here

TH1F *h1 = new TH1F("h1","primary neutron multiplicity",11,-0.5,10.5);
TH1F *h2 = new TH1F("h2","primary neutron energy (MeV)",2000,0.,20.);
TH1F *h3 = new TH1F("h3","vertex-x (cm)",100,-10.,10.);
TH1F *h4 = new TH1F("h4","vertex-y (cm)",100,-10.,10.);
TH1F *h5 = new TH1F("h5","vertex-z (cm)",200,-100.,100.);
TH1F *h11 = new TH1F("h11","Number of pulses per event",11,-0.5,10.5);
TH1F *h12 = new TH1F("h12","Energy deposited (keV)",200,0.,2000.);
TH1F *h13 = new TH1F("h13","Counter ID",4,0.5,4.5);
TH1F *h14 = new TH1F("h14","Capture time (microseconds)",250,0.,1000.);
TH1F *h15 = new TH1F("h15","Pulse duration (ns)",100,0.,50.);

//Initialize pulse counts to zero
Int_t total_pulses=0;
Int_t total_pulses_above_threshold=0;
Int_t number_generated_neutrons=0;

// Set threshold in keV
Double_t threshold=735.;

void analyze::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void analyze::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t analyze::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
      //   std::cout << "Number of primary neutrons is " << *nhe3primary << std::endl;
   h1->Fill(*nhe3primary);
   number_generated_neutrons=number_generated_neutrons+*nhe3primary;
   for(int i=0;i<*nhe3primary;i++){
        h2->Fill(E0he3[i]);
        h3->Fill(x0he3[i]);
        h4->Fill(y0he3[i]);
        h5->Fill(z0he3[i]);
   }
   total_pulses=total_pulses+*npulses;
   h11->Fill(*npulses);
   if(*npulses > 0){
   for(int j=0;j<*npulses;j++){
        h12->Fill(Ecounters[j]);
        h13->Fill(idcounters[j]);
        h14->Fill(tcounters[j]/1000.);
	h15->Fill(deltathe3[j]);
        if(Ecounters[j] > threshold){
           total_pulses_above_threshold++;
        }
   }
   }


   return kTRUE;
}

void analyze::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void analyze::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   TFile *out = new TFile("outhe3.hist","RECREATE");
   h1->Write();
   h2->Write();
   h3->Write();
   h4->Write();
   h5->Write();
   h11->Write();
   h12->Write();
   h13->Write();
   h14->Write();
   std::cout << "***********" << std::endl;
   std::cout << "Threshold is " << threshold << " keV" << std::endl;
   std::cout << "Total number of generated neutrons: " << number_generated_neutrons << std::endl;
   std::cout << "Total number of pulses: " << total_pulses << std::endl;
   std::cout << "Total number of pulses above threshold: " << total_pulses_above_threshold << std::endl;
Double_t eff,unceff;
   eff=1./((Double_t)number_generated_neutrons/(Double_t)total_pulses_above_threshold);
   unceff=sqrt(eff*(1-eff)/(Double_t)number_generated_neutrons);
   std::cout << "Estimated efficiency: " << eff << " +/- " << unceff << std::endl;
   std::cout << "***********" << std::endl;   
   std::cout << "abcd" << eff << " " << unceff << std::endl;
   out->Close();

}
