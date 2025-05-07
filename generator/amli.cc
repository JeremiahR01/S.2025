//J. Busenitz, based on modification of Am-Be generator by Ryan McKeown

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>

double energy[999],cumulative_prob[999];
int ndat=0;

void amli(double* rnd,int &ngen,double &mass1,int &idpdg1,
     double* p1,double &mass2,int &idpdg2,double* p2)
{

/*
 *** this routine generates the final state neutron for an Am-Li source

 *** references:
       1) NIM, volume 198, page 599 (1982).
 *** known limitations: 
      1) Gammas accompanying Am decay are not generated.
      2) Spectrum is approximated by generating a mixture of mono-energetic neutrons ranging in energy
         from 50 keV to 1.5 MeV.

 *** rnd(3) (real*4): array of 3 random numbers U(0,1)
 *** ngen (integer): number of particles generated
 *** mass1,mass2 (real*4): masses in GeV of particles generated
 *** idpdg1,idpdg2 (integer*4): PDG codes of particles generated
 *** p1(3), p2(3) (real*4): momentum in GeV/c of particles generated
 */

  int i;
  float probdum, edum;
  
    static bool amlifirst=true;
    //    double neumass=0.939566;
    int j1,j2;
    double esel,nke,gke,costh,sinth,phi,ptot,sum;
    double twopi=6.283185;

    std::ifstream spectrum;
    //  float energy[999],cumulative_prob[999];


  if(amlifirst){
      amlifirst=false;
  // Read in energies and respective weights of monoenergetic neutrons
      spectrum.open("amli_spectrum.dat");
      for(i=0;i<999;i++){
	energy[i]=0.;
	cumulative_prob[i]=0.;
      }
      while(!spectrum.eof()){
	ndat++;
	spectrum >> probdum >> edum;
	energy[ndat-1]=edum*1.e-6;
	if(ndat==1){
	  cumulative_prob[ndat-1]=probdum;
	}
	else{
	  cumulative_prob[ndat-1]=cumulative_prob[ndat-2]+probdum;
	}
      }
      // Normalize weights
      for(i=0;i<ndat;i++){
	cumulative_prob[i]=cumulative_prob[i]/cumulative_prob[ndat-1];
	//       	std::cout << energy[i] << " " << cumulative_prob[i] << std::endl;
      }
  }
  ngen=1;
  //  mass1=neumass;
  //  mass2=0;
  //  idpdg1=2112;
  //  idpdg2=0;

  for(int i=0;i<3;i++)
    {
      p1[i]=0;
      p2[i]=0;
     }
  //Throw neutron energy
  //  std::cout << "Ready to start loop to throw energy" << std::endl;
  for(int j=0;j<ndat;j++){
  //      std::cout << "At top of energy check.  j, cumulative probability and rnd[0]: " << j << " " << cumulative_prob[j] << "  " << rnd[0] << std::endl;
    if(cumulative_prob[j] >= rnd[0]){
  //          std::cout << " j, cumulative probablity, and rnd[0] are " << j << " " << cumulative_prob[j] << " " << rnd[0] << std::endl;
	nke=energy[j];
        break;
      }
  }
  //Throw polar and azimuthal angles for isotropic distribution
  costh=2*rnd[1]-1;
  phi=twopi*rnd[2];
  //3-momentum
  ptot=sqrt(2.*mass1*nke);
  // std::cout << "Energy and total momentum are " << nke << " " << ptot << std::endl;
  p1[0]=ptot*sqrt(1.-pow(costh,2.0))*cos(phi);
  p1[1]=ptot*sqrt(1.-pow(costh,2.0))*sin(phi);
  p1[2]=ptot*costh;
  //Done
}
  


int main(int argc, char**argv)
{

  int nevent=100000;
  unsigned int user_sd=0;  

  if(argc > 1)
    nevent=atoi(argv[1]);

  if(argc > 2)
    user_sd = (unsigned int) atoi(argv[2]);
  else 
    user_sd = (unsigned int) time(0);
  srand(user_sd);

  double rnd[3],p1[3],p2[3];
  
  int ngen,idpdg1,idpdg2;
  double mass1,mass2;

  for (int j=0;j<nevent;j++){
//    seedRand(static_cast<unsigned int>(j) );
    for (int i=0;i<3;i++)//need 3 random numbers to pass to ambe
      rnd[i]=rand()/(RAND_MAX+1.0);
    //    std::cout << rnd[0] << " " << rnd[1] << " " << rnd[2] << std::endl;
    ngen=1;
    mass1=0.939566;
    idpdg1=2112;
    mass2=0.;
    idpdg2=22;
    amli(rnd,ngen,mass1,idpdg1,p1,mass2,idpdg2, p2);

    std::cout << ngen <<std::endl ; //number of particles gen

    std::cout << "1 " <<idpdg1<<" 0 0 " << p1[0] << " "<< p1[1]
	   << " " << p1[2]<<" " << mass1 << std::endl ;
    if(ngen>1){
      std::cout << "1 " << idpdg2<< " 0 0 " << p2[0] << " " << p2[1]
	     <<" " << p2[2] << " " <<mass2<<std::endl;
    }
  }
  return 0;
}
