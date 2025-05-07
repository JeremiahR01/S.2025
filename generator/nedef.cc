//Questions email rwm33@drexel.edu Ryan McKeown

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>

void emon(double* rnd,int &ngen,double &mass1,int &idpdg1, double *p1, double &ekemev)
{

/*
 *** this routine generates a monoenergetic particle of mass1, ID id

 *** references:
       1) A.D. Vijaya and Arun Kumar, Nucl. Inst. and Meth. 111 (1973) 435.
       2) K.W. Geiger and L. van der Zwan, Nucl Inst. and Meth. 131 (1975) 315.
 *** known limitations: 
      1) energy spectrum of neutrons accompanied by a 4.4 MeV gamma simplified
         near endpoints, namely 1.8 MeV and 6.2 MeV.

 *** rnd(5) (real*4): array of 5 random numbers U(0,1)
 *** ngen (integer): number of particles generated
 *** mass1,mass2 (real*4): masses in GeV of particles generated
 *** idpdg1,idpdg2 (integer*4): PDG codes of particles generated
 *** p1(3), p2(3) (real*4): momentum in GeV/c of particles generated
 */

    double neumass=0.939566;
    double deneu=0.00010989;
    int j1,j2;
    double esel,nke,gke,costh,sinth,phi,ptot,sum;
    double twopi=6.283185;
    
    ngen=1;
    mass1=neumass;
    idpdg1=2112;

    for(int i=0;i<3;i++)
      {
	p1[i]=0;
      }
    costh=(1-2.0*rnd[0]);
    sinth=sqrt(1.-pow(costh,2.0));
    phi=twopi*rnd[1];
    ptot=sqrt(2.0*neumass*ekemev/1000.);
    p1[0]=ptot*sinth*cos(phi);
    p1[1]=ptot*sinth*sin(phi);
    p1[2]=ptot*costh;
}


int main(int argc, char**argv)
{

  int nevent=100000;
  unsigned int user_sd=0;  
  double ekemev;

  if(argc > 1)
    nevent=atoi(argv[1]);

  if(argc > 2)
    ekemev=atof(argv[2]);
//  std::cout << " kinetic energy in MeV is " << ekemev << std::endl;

  if(argc > 3)
    user_sd = (unsigned int) atoi(argv[3]);
  else 
    user_sd = (unsigned int) time(0);
  srand(user_sd);

  double rnd[5],p1[3],p2[3];
  
  int ngen,idpdg1,idpdg2;
  double mass1,mass2;

  for (int j=0;j<nevent;j++){
//    seedRand(static_cast<unsigned int>(j) );
    for (int i=0;i<5;i++)//need 5 random numbers to pass to ambe
      rnd[i]=rand()/(RAND_MAX+1.0);
    
    emon(rnd,ngen,mass1,idpdg1,p1,ekemev);
    
    std::cout << ngen <<std::endl ; //number of particles gen

    std::cout << "1 " <<idpdg1<<" 0 0 " << p1[0] << " "<< p1[1]
	   << " " << p1[2]<<" " << mass1 << std::endl ;
    }
}
