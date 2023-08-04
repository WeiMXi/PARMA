#include <iostream>
using namespace std;

double getHPcpp(int,int,int);
double getrcpp(double,double);
double getdcpp(double,double);
double getSpecCpp(int,double,double,double,double,double);
double getSpecAngFinalCpp(int,double,double,double,double,double,double);

int main()
{
const int npart = 33;
int IangPart[npart+1] = {1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6};

// Set condition 
int ip = 0; // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
double e = 100; // Energy (MeV/n)
int iyear = 2019; // Year
int imonth = 2;   // Month
int iday = 1;    // Day
double glat = 30.5; // Latitude (deg), -90 =< glat =< 90
double glong = -76.2; // Longitude (deg), -180 =< glong =< 180
double alti = 0.0; // Altitude (km)
double g = 0.15; // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
double ang = -0.5; // cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)

// calculate parameters
double s = getHPcpp(iyear,imonth,iday); // W-index (solar activity)
double r = getrcpp(glat,glong);  // Vertical cut-off rigidity (GV)
double d = getdcpp(alti,glat);   // Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

double Flux,DifFlux;

Flux=getSpecCpp(ip,s,r,d,e,g);
cout << "Angular Integrated Flux(/cm2/s/(MeV/n))= " << Flux << "\n";

if(IangPart[ip] > 0) {
 DifFlux=Flux*getSpecAngFinalCpp(IangPart[ip],s,r,d,e,g,ang);
 cout << "Angular Differential Flux(/cm2/s/(MeV/n)/sr)= " << DifFlux << "\n";
}

}
