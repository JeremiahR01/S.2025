//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 25 00:10:05 2020 by ROOT version 6.14/04
// from TTree the3/tree
// found on file: outhe3.root
//////////////////////////////////////////////////////////

#ifndef analyze_h
#define analyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class analyze : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<double> Ecounters = {fReader, "Ecounters"};
   TTreeReaderValue<Int_t> npulses = {fReader, "npulses"};
   TTreeReaderArray<double> tcounters = {fReader, "tcounters"};
   TTreeReaderArray<int> idcounters = {fReader, "idcounters"};
   TTreeReaderArray<double> E0he3 = {fReader, "E0he3"};
   TTreeReaderValue<Int_t> nhe3primary = {fReader, "nhe3primary"};
   TTreeReaderArray<double> x0he3 = {fReader, "x0he3"};
   TTreeReaderArray<double> y0he3 = {fReader, "y0he3"};
   TTreeReaderArray<double> z0he3 = {fReader, "z0he3"};
   TTreeReaderArray<double> deltathe3 = {fReader, "deltathe3"};


   analyze(TTree * /*tree*/ =0) { }
   virtual ~analyze() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(analyze,0);

};

#endif

#ifdef analyze_cxx
void analyze::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t analyze::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef analyze_cxx
