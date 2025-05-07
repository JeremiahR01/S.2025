//Questions email rwm33@drexel.edu Ryan McKeown

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <ctime>

void ambe(double* rnd,int &ngen,double &mass1,int &idpdg1,
     double* p1,double &mass2,int &idpdg2,double* p2)
{

/*
 *** this routine generates the final state neutron and gamma, if any,
 *** for an Am-Be source

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

  const int numbin=102;
  
  double contents1[numbin-2];
  static double cusum[numbin-2];

  double contents[numbin]={
    0.,13.,14.,13.5,12.,10.,8.8,8.,7.8,7.,
    6.5,6.3,6.,5.8,5.2,4.5,4.3,5.5,6.3,7.,
    7.5,7.6,7.3,7.1,7.8,10.,11.,14.,16.,18.5,
    18.,16.5,15.5,15.,14.3,14.,13.8,13.6,13.3,13.2,
    13.2,13.7,14.,14.5,15.,15.5,15.3,14.8,14.,13.,
    12.5,12.3,11.8,11.,10.,8.5,7.8,7.7,8.,8.2,
    7.9,7.7,6.6,6.3,6.1,6.2,6.5,7.0,7.5,8.0,
    8.3,8.35,8.2,7.8,7.3,6.8,5.8,5.,4.2,3.5,
    2.8,2.4,1.8,1.7,1.8,2.2,2.7,3.1,3.5,3.7,
    3.5,3.3,3.1,2.8,2.7,2.3,1.8,1.6,1.1,0.8,
    0.4,0.
  };
      
    static bool ambefirst=true;
    double neumass=0.939566;
    double deneu=0.00010989;
    int j1,j2;
    double esel,nke,gke,costh,sinth,phi,ptot,sum;
    double twopi=6.283185;
    
    ngen=0;
    mass1=0;
    mass2=0;
    idpdg1=0;
    idpdg2=0;

    for(int i=0;i<3;i++)
      {
	p1[i]=0;
	p2[i]=0;
      }

      if(ambefirst)
	{
	  sum=0;
	  for(int i=0;i<numbin;i++)
	    sum+=contents[i];
	  for(int i=1;i<numbin-1;i++)
            contents1[i-1]=contents[i]/sum;
	  
	  sum=0;
	  
	  for(int i=0;i<numbin-2;i++){
            sum+=contents1[i];
            cusum[i]=sum;
	  }
	  
	  ambefirst=false;
	}
      // here i decremented all array values to count for the fact
      // c arrays start at 0
      j1=1;
      for(int i=0;i<numbin-2;i++){
	
	if(rnd[0] > cusum[i]) 
	  j1=i;
      }
      
      if(j1 == 1) 
	nke=deneu*rnd[0]/cusum[0];
      else
	{
	  j2=j1+1;
	  nke=deneu*(float(j1)+(rnd[0]-cusum[j1])/(cusum[j2]-cusum[j1]));
	}
      
      if(nke > 0.0062) 
	gke=0;
      else if(nke > 0.0018) 
	gke=0.00443;
      else if(nke > 0.0005) 
	gke=0;
      else
	gke=0;
      
      if(nke > 0)
	{
	  ngen=1;
	  mass1=neumass;
	  idpdg1=2112;
	  costh=1.-2.*rnd[1];
	  sinth=sqrt(1. - pow(costh,2) );
	  phi=twopi*rnd[2];
	  ptot=sqrt(2.*neumass*nke);
	  p1[0]=ptot*sinth*cos(phi);
	  p1[1]=ptot*sinth*sin(phi);
	  p1[2]=ptot*costh;
	}

      if(gke > 0.)
	{
	  ngen=2;
	  mass2=0;
	  idpdg2=22;
	  costh=1.-2.*rnd[3];
	  sinth=sqrt(1.- pow(costh,2) );
	  phi=twopi*rnd[4];
	  ptot=gke;
	  p2[0]=ptot*sinth*cos(phi);
	  p2[1]=ptot*sinth*sin(phi);
	  p2[2]=ptot*costh;
	}
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

  double rnd[5],p1[3],p2[3];
  
  int ngen,idpdg1,idpdg2;
  double mass1,mass2;

  for (int j=0;j<nevent;j++){
//    seedRand(static_cast<unsigned int>(j) );
    for (int i=0;i<5;i++)//need 5 random numbers to pass to ambe
      rnd[i]=rand()/(RAND_MAX+1.0);
    
    ambe(rnd,ngen,mass1,idpdg1,p1,mass2,idpdg2, p2);
    
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
