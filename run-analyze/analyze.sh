#!/bin/sh
echo $1
cd $1
root -b <<!
TFile test1("outhe3.root");
the3->Process("../analyze.C");
.q!
!
cd ..
