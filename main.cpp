#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

double getHPcpp(int,int,int);
double getrcpp(double,double);
double getdcpp(double,double);
double getSpecCpp(int,double,double,double,double,double);
double getSpecAngFinalCpp(int,double,double,double,double,double,double);
double get511fluxCpp(double,double,double);

int main()
{
const int nebin=140; // number of energy bin
const int npart = 33;
const int nang=21; // from -1 to 1, 0.2 step
int IangPart[npart+1] = {1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6};

double emid[nebin+2],ewid[nebin+2];
double dcc[npart+1][nebin+1] = {}; // dose conversion coefficient (pSv*cm^2)
double Flux[npart+1][nebin+2] = {}; // energy differential flux, ie=nebin+1 for energy integrated
double dose[npart+2] = {};    // dose from each particle, npart+1 for total
double ang[nang+1]; // angle, 0 for angular integrated
double FluxAngEdif[nebin+2][nang+1] = {}; // Energy differntial & Angular differential flux

int ie511 = 78; // energy bin for 511 keV photon

string dname,str,uname;

string pname[npart+2] = {"neutro","proton","he---4","li---7","be---9","b---11","c---12","n---14","o---16","f---19" 
                        ,"ne--20","na--23","mg--24","al--27","si--28","p---31","s---32","cl--35","ar--40","k---39" 
                        ,"ca--40","sc--45","ti--48","v---51","cr--52","mn--55","fe--56","co--59","ni--59", 
                         "muplus","mumins","electr","positr","photon","total-"};

string condition = "Ang-EnergyDep"; // default input condition
string dccname = "ICRP116"; // Effective dose for ISO irradation (uSv/h)
// string dccname = "NM06tub"; // Count rate of 6 tube neutron monitor (100 count per hour)
// string dccname = "Air-pGy"; // Dose rate in air (uGy/h)
// string dccname = "h10ICRP"; // H*(10) (uSv/h)

int ip,iyear,imonth,iday,isout,istyle,ia,ie,i;
double e,glat,glong,alti,s,r,d,unitconv,tmp,dori;
double g; // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin

dname="Condition/"+condition+".inp";
ifstream conf(dname,ios::in);
getline(conf, str); 
istringstream conf1(str);
conf1 >> isout >> istyle;
// isout=0: No output for flux data for each condition, 
//      =1: Output Angular integrated flux in "SpecOut" folder
//      =2: Output Angular differential flux in "Angout" folder as a function of energy
//      =3: Output Angular differential flux in "Angout" folder as a function of cos(theta)
//      =4: Output Angular differential flux in "Angout" folder as a function of cos(theta) (absolute value)
// istyle = 0:force field poteintal (MV), cutoff rigidity (GV), atmospheric detph(g/cm^2), g(direct input)
//        = 1:year, month, day, latitude, longitude, atmospheric detph(g/cm^2), g(direct input)
//        = 2:year, month, day, latitude, longitude, altitude (m), g(direct input) 
//        = 3:year, month, day, latitude, longitude, altitude (ft), g(direct input)
//        = 4:year, month, day, cutoff rigidity, altitude (m), g(direct input)

if(dccname=="NM06tub") { // for neutron monitor count rate
 uname="100cph";
 unitconv=3600.0/100.0;
} else { //for dose
 uname="uSv/h ";
 unitconv=1.0e-6*3600;
}

// open output files
ofstream sf("SpecOut/"+condition+".out",ios::out);
sf << " Summary of calculated flux in (/(MeV/n)/cm2/sr/s) \n";
for (i=1;i<=npart+2;i++) {sf << setw(14) << i;} sf<<"\n";

ofstream df("DoseOut/"+condition+"."+dccname,ios::out);
df << " Summary of calculated dose in " << uname << "\n";
for (i=1;i<=npart+9;i++) {df << setw(15) << i;} df<<"\n";
df << "        W-index         Rc(Gv)   Depth(g/cm2)              g";
for (ip=0;ip<=npart+1;ip++) {df << "         " << pname[ip];} 
if(istyle>=1) df << "           Year          Month            Day";
df << "\n";

ofstream af("AngOut/"+condition+".out",ios::out);
if(isout==2) {for (i=1;i<=nang+1;i++) {af << setw(14) << i;} af<<"\n";}


//     Read Dose Conversion Coefficient Data
ifstream dccf("dcc/"+dccname+".inp",ios::in);
getline(dccf, str);
getline(dccf, str);
for(ie=1;ie<=nebin;ie++) {
 getline(dccf, str);
 istringstream dccf1(str);
 dccf1 >> emid[ie] >> ewid[ie];
 for(ip=0;ip<=npart;ip++) {dccf1 >> dcc[ip][ie];}
}
emid[nebin+1]=0.0;
ewid[nebin+1]=1.0;

// Set angle
ang[1]=-1.0;
for(ia=2;ia<=nang;ia++) {
 ang[ia]=ang[ia-1]+2.0/(nang-1);
 if(abs(ang[ia]) < 1.0e-10) ang[ia]=0.0;
} 

// Main calculation start here
while(true) {
 getline(conf, str); 
 if(conf.eof()) {break;}
 istringstream conf2(str);
 if(istyle==0) {
  conf2 >> s >> r >> d >> g;
 } else if(istyle>=1 && istyle <=3) {
  conf2 >> iyear >> imonth >> iday >> glat >> glong >> dori >> g;
  s = getHPcpp(iyear,imonth,iday);
  r = getrcpp(glat,glong);
  if(istyle==1) {
   d=dori; // direct input (g/cm2)
  } else if(istyle==2){
   alti=dori*0.001; // dori in m, alti in km
   d = getdcpp(alti,glat);
  } else {
   alti=dori*0.3048*0.001; // dori in ft, alti in km
   d = getdcpp(alti,glat);
  }
 } else if(istyle==4) {
  conf2 >> iyear >> imonth >> iday >> r >> d >> g;
  s = getHPcpp(iyear,imonth,iday);
 }
 for(ie=1;ie<=nebin;ie++) {  
  e=emid[ie];
  for(ip=0;ip<=npart;ip++) {
   Flux[ip][ie]=getSpecCpp(ip,s,r,d,e,g);
//   cout << tmp << " " << ip << " " << e << "\n";
   if(ip==npart && ie==ie511) {Flux[ip][ie]=Flux[ip][ie]+get511fluxCpp(s,r,d)/ewid[ie];} // add 511 keV flux
   Flux[ip][nebin+1]=Flux[ip][nebin+1]+Flux[ip][ie]; // energy integrated flux
   dose[ip]=dose[ip]+Flux[ip][ie]*dcc[ip][ie]*ewid[ie]; // Each Particle Dose (pSv/s)
   dose[npart+1]=dose[npart+1]+Flux[ip][ie]*dcc[ip][ie]*ewid[ie];  // Total Dose (pSv/s)
  }
 }
 if(isout==1) {
  sf << " *****************  New condition ***********************************************************************************************\n";
  if(istyle==0) {sf << "W-index= " << s << " Rc(GV)= " << r << " depth(g/cm2)= " << d << " g= " << g << "\n";}
  else if(istyle==1) {sf << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " depth(g/cm2)= " << dori << " g= " << g << "\n";}
  else if(istyle==2) {sf << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " Alti(m)= " << dori << " g= " << g << "\n";}
  else if(istyle==3) {sf << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " Alti(ft)= " << dori << " g= " << g << "\n";}
  else if(istyle==4) {sf << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Rc(GV)= " << r << " Depth(g/cm2)= " << d << " g= " << g << "\n";}
  sf << " Total dose in " << uname << " " << dose[npart+1]*unitconv << "\n";
  sf << " Dose rate    "; for(ip=0;ip<=npart;ip++) {sf << setw(14) << scientific << dose[ip]*unitconv;} sf << "\n";
  sf << " Dose contri  "; for(ip=0;ip<=npart;ip++) {sf << setw(14) << scientific << dose[ip]/dose[npart+1];} sf << "\n";
  sf << " Energy(MeV/n)"; for(ip=0;ip<=npart;ip++) {sf << "        " << pname[ip];} sf << " (/(MeV/n)/cm2/s) \n";
  for(ie=1;ie<=nebin;ie++) {
   sf << setw(14) << scientific << emid[ie];
   for(ip=0;ip<=npart;ip++) {sf << setw(14) << Flux[ip][ie];}
   sf << "\n";
  }
 }
 df << setw(15) << scientific << s;
 df << setw(15) << r;
 df << setw(15) << d;
 df << setw(15) << g;
 for(ip=0;ip<=npart+1;ip++) {df << setw(15) << dose[ip]*unitconv;}
 if(istyle>=1) {
  tmp=iyear+(imonth-1)/12.0+(iday-1)/365.0; // rough year
  df << setw(15) << tmp << setw(15) << imonth << setw(15) << iday;  
 }
 df << "\n";
 if(isout>=2) { // angular differential calculation
  af << " *****************  New condition ***********************************************************************************************\n";
  if(istyle==0) {af << defaultfloat << "W-index= " << s << " Rc(GV)= " << r << " depth(g/cm2)= " << d << " g= " << g << "\n";}
  else if(istyle==1) {af << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " depth(g/cm2)= " << dori << " g= " << g << "\n";}
  else if(istyle==2) {af << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " Alti(m)= " << dori << " g= " << g << "\n";}
  else if(istyle==3) {af << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Lat(deg)= " << glat << " Lon(deg)= " << glong << " Alti(ft)= " << dori << " g= " << g << "\n";}
  else if(istyle==4) {af << " Year= " << iyear << " Month= " << imonth << " Day= " << iday << " Rc(GV)= " << r << " Depth(g/cm2)= " << d << " g= " << g << "\n";}
  for(ip=0;ip<=npart;ip++) {
   if(IangPart[ip] > 0) {
    af << " Particle name = " << pname[ip] << "\n";
    int iemin=61;
    if(ip==0) {iemin=1;} // only for neutron, down to 1e-8 MeV
    if(isout>=3) {for (i=1;i<=nebin-iemin+2;i++) {af << defaultfloat << setw(14) << i;} af << "\n";} // output line number for each particle
    for(ie=iemin;ie<=nebin+1;ie++) {
     e=emid[ie];
     for(ia=1;ia<=nang;ia++) {
      FluxAngEdif[ie][ia]=Flux[ip][ie]*getSpecAngFinalCpp(IangPart[ip],s,r,d,e,g,ang[ia]);
     }
    }
    if(isout==2) { // output energy dependence
     af << " E(MeV/n) /cos";
     for(ia=1;ia<=nang;ia++) {af << defaultfloat << setw(14) << ang[ia];}
     af << scientific << " (/(MeV/n)/cm2/s/sr)\n";
     for(ie=iemin;ie<=nebin;ie++) {
      af << setw(14) << emid[ie];
      for(ia=1;ia<=nang;ia++) {af << setw(14) << FluxAngEdif[ie][ia];}
      af << "\n";
     }
    }
    else { // output angular dependence
     af << scientific << " cos /E(MeV/n)";
     for(ie=iemin;ie<=nebin;ie++) {af << setw(14) << emid[ie];}
     if(isout==3) {af << " (/sr)\n";}
     else {af << " (/(MeV/n)/cm2/s/sr)\n";}
     for(ia=1;ia<=nang;ia++) {
      af << defaultfloat << setw(14) << ang[ia];
      if(isout==3) {for(ie=iemin;ie<=nebin;ie++) {af << scientific << setw(14) << FluxAngEdif[ie][ia]/Flux[ip][ie];}} // Relative value
      else {for(ie=iemin;ie<=nebin;ie++) {af << scientific << setw(14) << FluxAngEdif[ie][ia];}} // Absolute value 
      af << "\n";
     }
    }
   }
  }
 }
}

} // end of program