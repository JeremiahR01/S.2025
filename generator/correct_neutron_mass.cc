//Reads in generated event data in HEPEvt format and outputs the same data except that neutron mass is set to 0.939566 GeV 

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

int ngen,pdgid,i,idum0,idum1,idum2;
double px,py,pz,mass;
double neutronmass=0.939566;

using namespace std;

int main()
{
 ifstream infile;
 infile.open("generator.data");
 while(!infile.eof()) {
   infile >> ngen;
   if(infile.eof()) break;
   cout << ngen << endl;
   for(i=0;i<ngen;i++){
     infile >> idum0 >>pdgid >> idum1 >> idum2 >> px >> py >> pz >> mass;
     if (pdgid == 2112) {
	 mass=0.939566;
       }
     cout << idum0 << " " << pdgid << " " << idum1 << " " << idum2 << " " << px << " " << py << " " << pz << " " << mass << " " << endl;
   }
 }
 return 0;
}
