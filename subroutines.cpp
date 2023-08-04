#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
using namespace std;

// **********************************************************
double getRfromECpp(int iz, int ia, double Ek, double Em) // get Rigidity in MV from Kinetic Energy (MeV/n)
// **********************************************************
{
     double getRfromE;
     getRfromE = sqrt(pow(ia * Ek, 2) + 2 * Ek * ia * Em) / iz;
     return getRfromE;
}

// *******************************************************
double getFFPfromWCpp(double s)
// *******************************************************
{
     double getFFPfromW;
     if (s >= 0)
     {
          getFFPfromW = 370.0 + 3.0e-1 * pow(s, 1.45); // FFP in MV
     }
     else
     {
          getFFPfromW = 370.0 - 3.0e-1 * pow(abs(s), 1.45); // FFP in MV
     }
     return getFFPfromW;
}

// *******************************************************
double getPowCpp(int ip, double d, double r)
// *******************************************************
{
     double getPow;
     const int npart = 13; // number of particle type, neutron, proton, alpha, photon, electron positron, mu+, mu-, Li-O
     const int nBdata = 2; // B(1) - B(4) : Pow = b1 + b2*d
     const int nAdata = 5; // ! A(1) - A(5) : B = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
     string chatmp;
     static string pname[npart + 1] = {"neutro", "proton", "alphaa", "elemag", "elemag", "elemag", "muon--", "muon--", "ions  ", "ions  ", "ions  ", "ions  ", "ions  ", "ions  "};
     static double A[npart + 1][nBdata + 1][nAdata + 1];
     static double B[nBdata + 1];
     static int ifirst = 0;

     int i, ia, ib;

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (i = 0; i <= 7; i++)
          {
               if (i <= 2)
               { // neutron, proton alpha
                    dname = "input/" + pname[i] + "/solar-dep.inp";
               }
               else if (i == 3)
               {
                    dname = "input/" + pname[i] + "/solar-dep-EL.inp";
               }
               else if (i == 4)
               {
                    dname = "input/" + pname[i] + "/solar-dep-PO.inp";
               }
               else if (i == 5)
               {
                    dname = "input/" + pname[i] + "/solar-dep-PH.inp";
               }
               else if (i == 6)
               {
                    dname = "input/" + pname[i] + "/solar-dep.plus";
               }
               else if (i == 7)
               {
                    dname = "input/" + pname[i] + "/solar-dep.mins";
               }
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (ib = 1; ib <= nBdata; ib++)
               {
                    getline(ifs, str); // read data
                    istringstream s(str);
                    for (ia = 1; ia <= nAdata; ia++)
                    {
                         s >> A[i][ib][ia];
                    }
               }
          }
          dname = "input/ions/solar-dep.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          for (i = npart - 5; i <= npart; i++)
          {
               for (ib = 1; ib <= nBdata; ib++)
               {
                    getline(ifs1, str); // read data
                    istringstream s1(str);
                    for (ia = 1; ia <= nAdata; ia++)
                    {
                         s1 >> A[i][ib][ia];
                    }
               }
          }
     }

     for (ib = 1; ib <= nBdata; ib++)
     {
          B[ib] = A[ip][ib][1] + A[ip][ib][2] * r + A[ip][ib][3] / (1.0 + exp((r - A[ip][ib][4]) / A[ip][ib][5]));
     }

     getPow = B[1] + B[2] * d;

     return getPow;
}

// *******************************************************
double CorrNeutCpp(double s, double r, double d, double e)
// *******************************************************
{
     double CorrNeut;
     const int nhensu = 9; // number of hensu
     const int mpara = 5;  // number of parameters to COR dependence
     const int ndep = 26;  // number of depth
     const int nsol = 2;   // solar minimum and maximum
     static double ainp[mpara + 1][nhensu + 1][ndep + 1][nsol + 1];
     static double dep[ndep + 1];
     static double spot[nsol + 1] = {0.0, 0.0, 150.0};
     double b[nhensu + 1], c[nsol + 1]; // temporary dimension

     static int ifirst = 0;

     string dname;
     string str;

     string chatmp;
     int ih, id, is, ip, itmp;
     double ratio, d1, d2, rc, powidx, A1, A2;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/neutro/correction-depth-rigid.inp";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          for (ih = 1; ih <= nhensu; ih++)
          {
               for (id = 1; id <= ndep; id++)
               {
                    getline(ifs, str); // read data
                    istringstream s(str);
                    s >> itmp >> dep[id];
                    for (is = 1; is <= nsol; is++)
                    {
                         for (ip = 1; ip <= mpara; ip++)
                         {
                              s >> ainp[ip][ih][id][is];
                         }
                    }
               }
          }
     }

     // find depth ID
     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
               break;
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     rc = max(1.0, r);
     for (is = 1; is <= nsol; is++)
     {
          for (ih = 1; ih <= nhensu; ih++)
          { // determine 9 hensu used in the correction equation
               d1 = ainp[1][ih][id - 1][is] + ainp[2][ih][id - 1][is] * rc + ainp[3][ih][id - 1][is] / (1 + exp((rc - ainp[4][ih][id - 1][is]) / ainp[5][ih][id - 1][is]));
               d2 = ainp[1][ih][id - 0][is] + ainp[2][ih][id - 0][is] * rc + ainp[3][ih][id - 0][is] / (1 + exp((rc - ainp[4][ih][id - 0][is]) / ainp[5][ih][id - 0][is]));
               b[ih] = d1 + (d2 - d1) * ratio;
          }
          c[is] = pow(10, (b[1] + (b[2] * log10(e) + b[3]) * (1 - tanh(b[4] * log10(e / b[5]))) + (b[6] * log10(e) + b[7]) * (1 + tanh(b[8] * log10(e / b[9])))));
     }

     ip = 0;                       // always neutron
     powidx = getPowCpp(ip, d, r); // get Power index
     A2 = (c[1] - c[2]) / (pow(getFFPfromWCpp(spot[1]), powidx) - pow(getFFPfromWCpp(spot[2]), powidx));
     A1 = c[1] - A2 * pow(getFFPfromWCpp(spot[1]), powidx);
     CorrNeut = A1 + A2 * pow(getFFPfromWCpp(s), powidx);

     return CorrNeut;
}

// *******************************************************
void getGparaCpp(double g, double *geo)
// *******************************************************
{
     static double p[25] = {};
     string dname;
     string str;
     string chatmp;

     if (p[1] == 0)
     { // Initialization
          dname = "input/neutro/Geo-Dep.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Read p(1)-p(3)
          istringstream s1(str);
          s1 >> p[1] >> p[2] >> p[3];
          dname = "input/neutro/Water-Dep.inp";
          ifstream ifs2(dname, ios::in);
          getline(ifs2, str); // Kara-yomi
          getline(ifs2, str); // Read data
          istringstream s20(str);
          s20 >> chatmp >> p[4] >> p[5] >> p[6];
          getline(ifs2, str); // Read data
          istringstream s21(str);
          s21 >> chatmp >> p[7] >> p[8] >> p[9];
          getline(ifs2, str); // Read data
          istringstream s22(str);
          s22 >> chatmp >> p[10] >> p[11] >> p[12] >> p[13] >> p[14];
          dname = "input/neutro/Aircraft-Dep.inp";
          ifstream ifs3(dname, ios::in);
          getline(ifs3, str); // Kara-yomi
          getline(ifs3, str); // Read data
          istringstream s30(str);
          s30 >> p[15] >> p[16] >> p[17] >> p[18] >> p[19];
          getline(ifs3, str); // Read data
          istringstream s31(str);
          s31 >> p[20] >> p[21] >> p[22] >> p[23] >> p[24];
     }

     if (g >= 10.0)
     {                  // in semi-infite atmosphere
          geo[3] = 1.0; // if geo(3)=0, the value should be in NaN
          geo[5] = 1.0; // if geo(5)=0, the value should be in NaN
     }
     else if (g >= 0.0)
     { // for normal ground case
          geo[1] = p[1];
          geo[2] = p[2];
          geo[3] = pow(10.0, (p[4] + p[5] / (p[6] + g)));
          geo[4] = p[3];
          geo[5] = p[7] + p[8] * g + p[9] * pow(g, 2);
          geo[6] = (p[10] + p[11] * exp(-p[12] * g)) / (1 + p[13] * exp(-p[14] * g));
     }
     else
     { // pilot or cabin location
          int is = 14;
          if (g <= -10.0)
          {
               is = is + 5;
          } // for passenger & small aircraft configuration, skip 5 more data
          for (int i = 1; i <= 5; i++)
          {
               geo[i] = p[is + i];
          }
          geo[6] = 0.0; // no thermal component
     }

     return;
}

// *******************************************************
double getA4Cpp(double r, double d)
// *******************************************************
{
     double getA4 = 0.0;
     const int nBdata = 4; // B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
     const int nAdata = 6; // A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
     static double A[nAdata + 1] = {};
     static double B[nBdata + 1] = {};

     string dname;
     string str;

     double tmp1, tmp2;
     string chatmp;
     int ia;

     if (B[2] == 0.0)
     { // first time
          dname = "input/neutro/Depth-Dep-mid.out";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          getline(ifs, str); // Kara-yomi
          getline(ifs, str);
          istringstream s(str);
          s >> tmp1 >> tmp2 >> B[2] >> B[3] >> B[4];
          dname = "input/neutro/Rigid-Dep.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // read B(5) data
          istringstream s1(str);
          s1 >> chatmp;
          for (ia = 1; ia <= nAdata; ia++)
          {
               s1 >> A[ia];
          }
     }

     B[1] = A[1] + A[2] * r + A[3] / (1 + exp((r - A[4]) / A[5]));
     getA4 = B[1] + B[2] * d / (1 + B[3] * exp(B[4] * d));

     return getA4;
}

// *******************************************************
double getA12Cpp(double r, double d)
// *******************************************************
{
     double getA12 = 0.0;
     const int nBdata = 4; // B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
     const int nAdata = 6; // A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
     static double A[nBdata + 1][nAdata + 1] = {};
     static double B[nBdata + 1] = {};

     string dname;
     string str;

     double tmp0, tmp1, tmp2, tmp3;
     string chatmp;
     int ia, ib;

     if (B[4] == 0.0)
     { // first time
          dname = "input/neutro/Depth-Dep-hig.out";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          getline(ifs, str); // Kara-yomi
          getline(ifs, str);
          istringstream s(str);
          s >> tmp0 >> tmp1 >> tmp2 >> tmp3 >> B[4];
          dname = "input/neutro/Rigid-Dep.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi B(1)
          getline(ifs1, str); // Kara-yomi B(2)
          getline(ifs1, str); // Kara-yomi B(3)
          getline(ifs1, str); // Kara-yomi B(4)
          getline(ifs1, str); // Kara-yomi B(5)
          for (ib = 1; ib <= 3; ib++)
          {
               getline(ifs1, str); // read B data
               istringstream s1(str);
               s1 >> chatmp;
               for (ia = 1; ia <= nAdata; ia++)
               {
                    s1 >> A[ib][ia];
               }
          }
     }

     for (ib = 1; ib <= 3; ib++)
     {
          B[ib] = A[ib][1] + A[ib][2] * r + A[ib][3] / (1 + exp((r - A[ib][4]) / A[ib][5]));
     }

     getA12 = B[1] * (exp(-B[2] * d) + B[3] * exp(-B[4] * d));

     return getA12;
}

// *******************************************************
double getBestRCpp(int is, double r, double d, int ip)
// *******************************************************
{
     double getBestR;
     const int nBdata = 6;
     const int ndep = 26;
     const int npart = 3; // high-altitude correction is necessary only for electron, positron, photon
     static double A[npart + 1][nBdata + 1][ndep + 1] = {};
     static double B[nBdata + 1];
     static double dep[ndep + 1];
     string chatmp;
     static string pname[npart + 1] = {"  ", "EL", "PO", "PH"};
     static int ifirst = 0;

     int i, ib, id;
     double ratio;

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (i = 0; i <= npart; i++)
          {
               if (i == 0)
               {
                    dname = "input/neutro/bestR.inp";
               }
               else
               {
                    dname = "input/elemag/bestR-" + pname[i] + ".inp";
               }
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (id = 1; id <= ndep; id++)
               {
                    getline(ifs, str); // read data
                    istringstream s(str);
                    s >> dep[id];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s >> A[i][ib][id];
                    }
               }
          }
     }

     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
          {
               break;
          }
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     for (ib = 1; ib <= nBdata; ib++)
     {
          B[ib] = A[ip][ib][id - 1] + (A[ip][ib][id] - A[ip][ib][id - 1]) * ratio;
     }

     if (is == 1)
     { // solar minimum
          getBestR = pow(10, B[1] + B[2] * r + B[3] / r);
     }
     else
     {
          getBestR = pow(10, B[4] + B[5] * r + B[6] / r);
     }

     return getBestR;
}

// *******************************************************
double getFlCpp(int ip, double s, double r, double d)
// *******************************************************
{
     double getFl = 0.0;
     const int npart = 11;  // number of particle type, 0:neutron, 1:proton, 2:alpha, 3:electron, 4:positron, 5:photon, 6-11: Li-O
     const int nBdata = 4;  // B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
     const int nAdata = 10; // A(1) - A(10) : Bmin = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5))), Bmin = A(6)+A(7)*r+A(8)/(1+exp((r-A(9))/A(10)))
     const int nsor = 2;    //  solar minimum & maximum
     string chatmp;
     static string pname[npart + 1] = {"neutro", "proton", "alphaa", "elemag", "elemag", "elemag", "ions  ", "ions  ", "ions  ", "ions  ", "ions  ", "ions  "};
     static double A[npart + 1][nBdata + 1][nAdata + 1], B[nBdata + 1];
     static double spot[nsor + 1] = {0.0, 0.0, 150.0};
     static double Fl[nsor + 1];
     static int ifirst = 0;
     int i, ia, ib, is, ipidx;
     double r1, powidx, A1, A2;

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (i = 0; i <= 5; i++)
          {
               if (i <= 2)
               {
                    dname = "input/" + pname[i] + "/Rigid-Dep.inp"; // neutron, proton, alpha
               }
               else if (i == 3)
               {
                    dname = "input/" + pname[i] + "/Rigid-Dep-EL.inp"; // electron
               }
               else if (i == 4)
               {
                    dname = "input/" + pname[i] + "/Rigid-Dep-PO.inp"; // positron
               }
               else if (i == 5)
               {
                    dname = "input/" + pname[i] + "/Rigid-Dep-PH.inp"; // photon
               }
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (ib = 1; ib <= nBdata; ib++)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    s >> chatmp;
                    for (ia = 1; ia <= nAdata; ia++)
                    {
                         s >> A[i][ib][ia];
                    } // A(4)&A(12) is s,r,d-dependence, so will be changed
               }
          }
          dname = "input/ions/Rigid-Dep.inp";
          ifstream ifs1(dname, ios::in);
          for (i = npart - 5; i <= npart; i++)
          {
               getline(ifs1, str);
               for (ib = 1; ib <= nBdata; ib++)
               {
                    getline(ifs1, str);
                    istringstream s1(str);
                    s1 >> chatmp;
                    for (ia = 1; ia <= nAdata; ia++)
                    {
                         s1 >> A[i][ib][ia];
                    } // A(4)&A(12) is s,r,d-dependence, so will be changed
               }
          }
     }

     for (is = 1; is <= nsor; is++)
     { // solar minimum and maximum
          if (ip == 1 || ip == 2 || ip >= npart - 5)
          { // ! need not correction
               r1 = r;
          }
          else
          { // for neutron, electron, positron and photon, Rc for high-altitude should be corrected
               if (ip == 0)
               {
                    ipidx = ip; // for neutron
               }
               else
               {
                    ipidx = ip - 2; // 1:electron, 2:positron, 3:photon
               }
               r1 = r * getBestRCpp(is, r, d, ipidx);
          }
          for (ib = 1; ib <= nBdata; ib++)
          {
               if (is == 1)
               { // solar minimum
                    B[ib] = A[ip][ib][1] + A[ip][ib][2] * r1 + A[ip][ib][3] / (1 + exp((r1 - A[ip][ib][4]) / A[ip][ib][5]));
               }
               else
               {
                    B[ib] = A[ip][ib][6] + A[ip][ib][7] * r1 + A[ip][ib][8] / (1 + exp((r1 - A[ip][ib][9]) / A[ip][ib][10]));
               }
          }
          Fl[is] = B[1] * (exp(-B[2] * d) - B[3] * exp(-B[4] * d));
     }

     if (ip <= 5)
     {
          powidx = getPowCpp(ip, d, r); // get Power index
     }
     else
     {
          powidx = getPowCpp(ip + 2, d, r); // in getPow, ip should be +2 for ions
     }

     A2 = (Fl[1] - Fl[2]) / (pow(getFFPfromWCpp(spot[1]), powidx) - pow(getFFPfromWCpp(spot[2]), powidx));
     A1 = Fl[1] - A2 * pow(getFFPfromWCpp(spot[1]), powidx);
     getFl = A1 + A2 * pow(getFFPfromWCpp(s), powidx);

     return getFl;
}

// *******************************************************
double getNeutSpecCpp(double s, double r1, double d1, double g, double e)
//     s:Force Field Potential (MV)
//     r:cut off rigidity (GV)
//     d:air depth (g/cm^2)
//     g:local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
//     e:neutron energy (MeV)
// *******************************************************
{
     double getNeutSpec = 0.0;
     const int nA = 12;       // number of basic spectrum parameter
     const int nG = 6;        // number of geometry parameter
     static double a[nA + 1]; // basic spectrum parameter
     double geo[nG + 1] = {}; // geometry parameter, static in Fortran version, but probably not necessary to be remembered

     static int ifirst = 0;
     static double airbus = 2.45; // weight of airbus340 (100t)
     static double Eth = 2.5e-8;  // themal energy

     double basic, fG, gtmp, geofactor, ther;

     string dname;
     string str;

     int ia;

     if (ifirst == 0)
     {    // first time, get universal parameter (i.e: independent of all parameters)
          //      Read A parameter
          dname = "input/neutro/fitting-lowspec.inp";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          getline(ifs, str);
          istringstream s(str);
          for (ia = 1; ia <= nA; ia++)
          {
               s >> a[ia];
          } // A(4)&A(12) is s,r,d-dependence, so will be changed
          ifirst = 1;
     }

     double r = max(1.0, r1);  // secondary particle fluxes are the same for Rc<1GV
     double d = max(0.15, d1); // secondary particle fluxes are the same for d < 0.15 g/cm2

     //     get condition dependent parameters
     double Fl = getFlCpp(0, s, r, d);
     a[12] = getA12Cpp(r, d);
     a[4] = getA4Cpp(r, d);
     getGparaCpp(g, geo); // obtain G parameters

     //     calculate flux (/cm^2/s/lethargy)
     double x = e; // x is energy

     double evap = a[1] * pow((x / a[2]), a[3]) * exp(-x / a[2]);
     double gaus = a[4] * exp(-pow((log10(x) - log10(a[5])), 2) / (2 * pow(log10(a[6]), 2)));
     double conti = a[7] * log10(x / a[8]) * (1 + tanh(a[9] * log10(x / a[10]))) * (1 - tanh(a[11] * log10(x / a[12])));

     basic = conti + evap + gaus;

     basic = basic * CorrNeutCpp(s, r, d, e); // correction for high altitdue

     fG = geo[1] + geo[2] * log10(x / geo[3]) * (1 - tanh(geo[4] * log10(x / geo[5])));
     if (g < 0.0)
     {
          if (g < -10.0)
          { // cabin
               gtmp = g + 10.0;
          }
          else
          { // pilot
               gtmp = g;
          }
          fG = fG * (abs(gtmp) - 10 * int(abs(gtmp) / 10.0)) / airbus; // consider size of aircraft
     }

     geofactor = pow(10.0, fG);
     ther = geo[6] * pow(x / Eth, 2) * exp(-(x / Eth));

     getNeutSpec = Fl * (basic * geofactor + ther) / e;

     return getNeutSpec;
}

// *******************************************************
double getEfromRCpp(int iz, double Rm, double COR) // get Kinetic Energy (MeV) from Rigidity (MV)
// *******************************************************
{
     double getEfromR;
     getEfromR = sqrt(pow(iz * COR, 2) + pow(Rm, 2)) - Rm;
     return getEfromR;
}

// *******************************************************

// **********************************************************
double getTOAspecCpp(int iz, int ia, double Ek, double Spot) // get TOA spectrum in (/(MeV/n)/s/m^2/sr)
//    Ek: Kinetic Energy in MeV/n
//    Spot: Force Field Potential
// **********************************************************
{
     double getTOAspec;
     const int npart = 28;
     static double Dpara[npart + 1] = {0.0, 1.85e4, 3.69e3, 19.50, 17.70, 49.20, 103.00, 36.70, 87.40, 3.19, 16.40, 4.43, 19.30, 4.17, 13.40, 1.15, 3.06, 1.30, 2.33, 1.87, 2.17, 0.74, 2.63, 1.23, 2.12, 1.14, 9.32, 0.10, 0.48};
     static double alpha[npart + 1] = {0.0, 2.85, 3.12, 3.41, 4.30, 3.93, 3.18, 3.77, 3.11, 4.05, 3.11, 3.14, 3.65, 3.46, 3.00, 4.04, 3.30, 4.40, 4.33, 4.49, 2.93, 3.78, 3.79, 3.50, 3.28, 3.29, 3.01, 4.25, 3.52};
     static double gamma[npart + 1] = {0.0, 2.74, 2.77, 2.82, 3.05, 2.96, 2.76, 2.89, 2.70, 2.82, 2.76, 2.84, 2.70, 2.77, 2.66, 2.89, 2.71, 3.00, 2.93, 3.05, 2.77, 2.97, 2.99, 2.94, 2.89, 2.74, 2.63, 2.63, 2.63};
     static double Emp = 938.27; // mass of proton, nucleus mass is simply assumed to be A*Emp

     double R, beta, dR2dE, SpecLIS, R0, delta;

     // ***** Determine LIS spectra based on DLR Model *****************
     R = getRfromECpp(iz, ia, Ek, Emp * ia) * 0.001; // rigidity in GV
     beta = sqrt(1 - pow(Emp * ia / (Emp * ia + Ek * ia), 2));
     dR2dE = 0.001 / iz / beta * ia; // convert GV to MV, MeV to MeV/n
     SpecLIS = Dpara[iz] * pow(beta, alpha[iz]) / pow(R, gamma[iz]) * dR2dE;
     // ****************************************************************

     // consider solar modulation
     R0 = getFFPfromWCpp(Spot) * 0.001; // FFP in GV from W, taken from Matthia-ASR2013
     delta = 0.02 * Spot + 4.7;         // taken from Matthia-ASR2013
     getTOAspec = SpecLIS * pow(R / (R + R0), delta);

     return getTOAspec;
}

// *************************************************************
double getPrimaryCpp(int iz, int ia, double s, double d, double e) // get Primary Flux
// *************************************************************
{
     double getPrimary;
     const int nAdata = 3;
     const int npart = 28;
     const int nepoint = 6;
     const int ngroup = 6;
     const int ndep = 26;
     static double a[nAdata + 1][ngroup + 1][ndep + 1];
     static double dEdxTable[npart + 1][nepoint + 1], epoint[nepoint + 1]; // dE/dx for each particle
     static double dep[ndep + 1];                                          // depth (g/cm2)
     static double down = 0.0;
     static int igidx[npart + 1] = {0, 1, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6}; // group index
     double b[nAdata + 1];                                                                                                  // temporary used dimension
     static string gname[ngroup + 1] = {"  ", "H-", "He", "Be", "N-", "Si", "Fe"};

     string dname;
     string str;

     int i, ie, ig, id, ip, itmp;
     double ratio, tmp, dEdx, Eini;

     if (down == 0.0)
     { // Initialization
          down = 1.720313e0;
          for (ig = 1; ig <= ngroup; ig++)
          {
               dname = "input/ions/primary-" + gname[ig] + ".inp";
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (id = 1; id <= ndep; id++)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    s >> dep[id];
                    for (i = 1; i <= nAdata; i++)
                    {
                         s >> a[i][ig][id];
                    }
               }
          }
          dname = "input/ions/dEdx-table.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str); // Kara-yomi
          getline(ifs1, str);
          istringstream s11(str);
          for (ie = 1; ie <= nepoint; ie++)
          {
               s11 >> epoint[ie];
          }
          for (ip = 1; ip <= npart; ip++)
          {
               getline(ifs1, str);
               istringstream s12(str);
               s12 >> itmp >> itmp;
               for (ie = 1; ie <= nepoint; ie++)
               {
                    s12 >> dEdxTable[ip][ie];
               }
          }
     }

     ig = igidx[iz];

     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
          {
               break;
          }
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     for (i = 1; i <= nAdata; i++)
     {
          b[i] = a[i][ig][id - 1] + (a[i][ig][id] - a[i][ig][id - 1]) * ratio;
     }

     // find dEdx
     for (ie = 2; ie <= nepoint - 1; ie++)
     {
          if (e <= epoint[ie])
          {
               break;
          }
     }
     ratio = (log(e) - log(epoint[ie - 1])) / (log(epoint[ie]) - log(epoint[ie - 1])); // log-logg interpolation
     ratio = max(0.0, min(1.0, ratio));
     tmp = log(dEdxTable[iz][ie - 1]) + (log(dEdxTable[iz][ie]) - log(dEdxTable[iz][ie - 1])) * ratio; // log-log interpolation
     dEdx = exp(tmp);

     Eini = e + dEdx * d; // Energy at the TOA
     getPrimary = getTOAspecCpp(iz, ia, Eini, s) * (b[1] * exp(-b[2] * d) + (1.0 - b[1]) * exp(-b[3] * d));
     getPrimary = getPrimary * 4.0 * acos(-1.0) * 1.0e-4 / down; // convert (/(MeV/n)/s/m^2/sr) to (/(MeV/n)/s/cm^2)

     return getPrimary;
}

// *******************************************************
double getSecondaryCpp(int ip, double s, double r1, double d1, double e)
// *******************************************************
{
     double getSecondary = 0.0;
     const int nBdata = 8;
     const int ndep = 26;
     const int npart = 11; // proton, alpha, electron, positron, photon, Li, Be, B, C, N, O
     static double a[nBdata + 1][npart + 1][ndep + 1], b[nBdata + 1];
     static double dep[ndep + 1];
     static int ifirst = 0;

     string dname;
     string str;

     int i, id, ib;
     double r, d, ratio;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (i = 1; i <= 5; i++)
          {
               if (i == 1)
               {
                    dname = "input/proton/fitting-lowspec.inp";
               }
               else if (i == 2)
               {
                    dname = "input/alphaa/fitting-lowspec.inp";
               }
               else if (i == 3)
               {
                    dname = "input/elemag/fitting-lowspec-EL.inp";
               }
               else if (i == 4)
               {
                    dname = "input/elemag/fitting-lowspec-PO.inp";
               }
               else if (i == 5)
               {
                    dname = "input/elemag/fitting-lowspec-PH.inp";
               }
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (id = 1; id <= ndep; id++)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    s >> dep[id];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s >> a[ib][i][id];
                    }
               }
          }
          dname = "input/ions/fitting-lowspec.inp";
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str); // Kara-yomi
          for (i = npart - 5; i <= npart; i++)
          {
               getline(ifs1, str);
               istringstream s1(str);
               for (ib = 1; ib <= nBdata; ib++)
               {
                    s1 >> a[ib][i][1];
               }
               for (ib = 1; ib <= nBdata; ib++)
               {
                    for (id = 1; id <= ndep; id++)
                    {
                         a[ib][i][id] = a[ib][i][1]; // for Li to O, depth independent
                    }
               }
          }
     }

     r = max(1.0, r1);  // to get Fl, minimum Rc=1GV
     d = max(0.15, d1); // secondary particle fluxes are the same for d < 0.15 g/cm2

     if (e < 1.0e-2 || ip > npart)
     { // for lower energy of heavy ions, no output
          getSecondary = 0.0;
          return getSecondary;
     }

     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
          {
               break;
          }
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     for (ib = 1; ib <= nBdata; ib++)
     {
          b[ib] = a[ib][ip][id - 1] + (a[ib][ip][id] - a[ib][ip][id - 1]) * ratio;
     }

     if (ip == 1 || ip == 2)
     { // proton or alpha
          getSecondary = getFlCpp(ip, s, r, d) * (b[1] * pow(e, b[2])) / (1 + b[3] * pow(e, b[4])) / (1 + b[5] * pow(e, b[6])) * (1 + exp(-b[7] * (log(e) + b[8])));
     }
     else if (ip == 3 || ip == 4)
     { // electron or positron
          getSecondary = getFlCpp(ip, s, r, d) * (b[1] * pow(e, b[2])) / (1 + b[3] * pow(e, b[4])) / (1 + b[5] * pow(e, b[6]));
     }
     else if (ip == 5)
     { // photon
          getSecondary = getFlCpp(ip, s, r, d) * (b[1] * pow(e, b[2])) * (1 + b[3] * pow(e, b[4])) / (1 + b[5] * pow(e, b[6])) / (1 + exp(-b[7] * (log(e) + b[8])));
     }
     else
     { // Li,Be,B,C,N,O
          getSecondary = getFlCpp(ip, s, r, d) * (b[1] * pow(e, b[2])) / (1 + b[3] * pow(e, b[4])) / (1 + b[5] * pow(e, b[6])) / (1 + exp(-b[7] * (log(e) + b[8])));
     }

     return getSecondary;
}

// *******************************************************
double getIonSpecCpp(int iz, double s, double r, double d, double e)
//     s:Force Field Potential (MV)
//     r:cut off rigidity (GV)
//     d:air depth (g/cm^2)
//     e:ion energy (MeV/n)
// *******************************************************
{
     double getIonSpec = 0.0;
     const int nAdata = 6;
     const int npart = 28;
     const int ngroup = 6;

     string chatmp;

     static double a[nAdata + 1][ngroup + 1]; // parameter used in combine.for
     static int iAnum[npart + 1] = {0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51, 52, 55, 56, 59, 59};
     static int ifirst = 0;

     static int igidx[npart + 1] = {0, 1, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6}; // group index

     static double restmass = 938.27e0;

     string dname;
     string str;

     int i, ia, ig, ip;
     double x, Ecut, EcPri, EcSec, tmp;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/ions/Combine.inp";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          for (i = 1; i <= ngroup; i++)
          {
               getline(ifs, str);
               istringstream s(str);
               for (ia = 1; ia <= nAdata; ia++)
               {
                    s >> a[ia][i];
               } // A(4)&A(12) is s,r,d-dependence, so will be changed
          }
     }

     if (e < 1.0e-2)
     { // for lower energy, no output
          getIonSpec = 0.0;
          return getIonSpec;
     }

     x = e; // x is energy
     ig = igidx[iz];

     tmp = restmass * iAnum[iz];

     Ecut = getEfromRCpp(iz, tmp, r * 1000.0) / iAnum[iz] - a[6][ig] * d; // MeV/n
     EcPri = max(a[1][ig], Ecut * a[3][ig]);
     EcSec = max(a[2][ig], Ecut * a[3][ig]);

     if (iz <= 2)
     {
          ip = iz;
     }
     else
     {
          ip = iz + 3; // in getsecondary, ip=iz+3 (electron, positron, photon) for ions
     }

     getIonSpec = getPrimaryCpp(iz, iAnum[iz], s, d, x) * 0.5 * (tanh(a[4][ig] * (x / EcPri - 1)) + 1.0) + getSecondaryCpp(ip, s, r, d, x) * 0.5 * (tanh(a[5][ig] * (1 - x / EcSec)) + 1.0);

     return getIonSpec;
}

// *******************************************************
void getAmuon(int ip, int is, double d, double r, double *acurr)
// *******************************************************
{
     const int npart = 2;   // Muon+ or Muon-
     const int nAdata = 7;  // ! number of A parameter
     const int nBdata = 10; // A_min = B1+B2*r+B3/(1+exp((r-B4)/B5), A_max = B6+B7*r+B8/(1+exp((r-B9)/B10)
     const int ndep = 26;   // number of depth

     static double Bdata[nAdata + 1][nBdata + 1][npart + 1][ndep + 1];
     static double dep[ndep + 1];

     static string charge[npart + 1] = {"    ", "plus", "mins"};
     static int ifirst = 0;

     string dname;
     string str;

     int ip2, id2, ib, ia, id;
     double rc, B1, B2, ratio;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (ip2 = 1; ip2 <= npart; ip2++)
          {
               dname = "input/muon--/final135." + charge[ip2];
               ifstream ifs1(dname, ios::in);
               getline(ifs1, str); // Kara-yomi
               for (id2 = 1; id2 <= ndep; id2++)
               {
                    getline(ifs1, str);
                    istringstream s1(str);
                    s1 >> dep[id2] >> Bdata[1][1][ip2][id2] >> Bdata[3][1][ip2][id2] >> Bdata[5][1][ip2][id2];
               }
               dname = "input/muon--/final2467." + charge[ip2];
               ifstream ifs2(dname, ios::in);
               getline(ifs2, str); // Kara-yomi
               for (id2 = 1; id2 <= ndep; id2++)
               {
                    getline(ifs2, str);
                    istringstream s21(str);
                    s21 >> dep[id2];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s21 >> Bdata[2][ib][ip2][id2];
                    }
               }
               getline(ifs2, str); // Kara-yomi
               for (id2 = 1; id2 <= ndep; id2++)
               {
                    getline(ifs2, str);
                    istringstream s22(str);
                    s22 >> dep[id2];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s22 >> Bdata[4][ib][ip2][id2];
                    }
               }
               getline(ifs2, str); // Kara-yomi
               for (id2 = 1; id2 <= ndep; id2++)
               {
                    getline(ifs2, str);
                    istringstream s23(str);
                    s23 >> dep[id2];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s23 >> Bdata[6][ib][ip2][id2];
                    }
               }
               getline(ifs2, str); // Kara-yomi
               for (id2 = 1; id2 <= ndep; id2++)
               {
                    getline(ifs2, str);
                    istringstream s24(str);
                    s24 >> dep[id2];
                    for (ib = 1; ib <= nBdata; ib++)
                    {
                         s24 >> Bdata[7][ib][ip2][id2];
                    }
               }
          }
     }

     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
          {
               break;
          }
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     rc = max(1.0, r); // minimum Rc = 1.0GV

     for (ia = 1; ia <= nAdata; ia++)
     {
          if (ia == 1)
          {
               acurr[ia] = log(Bdata[ia][1][ip][id - 1]) + (log(Bdata[ia][1][ip][id]) - log(Bdata[ia][1][ip][id - 1])) * ratio; // log-interpolation
               acurr[ia] = exp(acurr[ia]);
          }
          else if (ia == 3 || ia == 5)
          {
               acurr[ia] = Bdata[ia][1][ip][id - 1] + (Bdata[ia][1][ip][id] - Bdata[ia][1][ip][id - 1]) * ratio;
          }
          else
          {
               if (is == 1)
               {
                    B1 = Bdata[ia][1][ip][id - 1] + Bdata[ia][2][ip][id - 1] * rc + Bdata[ia][3][ip][id - 1] / (1 + exp((rc - Bdata[ia][4][ip][id - 1]) / Bdata[ia][5][ip][id - 1]));
                    B2 = Bdata[ia][1][ip][id] + Bdata[ia][2][ip][id] * rc + Bdata[ia][3][ip][id] / (1 + exp((rc - Bdata[ia][4][ip][id]) / Bdata[ia][5][ip][id]));
               }
               else
               {
                    B1 = Bdata[ia][6][ip][id - 1] + Bdata[ia][7][ip][id - 1] * rc + Bdata[ia][8][ip][id - 1] / (1 + exp((rc - Bdata[ia][9][ip][id - 1]) / Bdata[ia][10][ip][id - 1]));
                    B2 = Bdata[ia][6][ip][id] + Bdata[ia][7][ip][id] * rc + Bdata[ia][8][ip][id] / (1 + exp((rc - Bdata[ia][9][ip][id]) / Bdata[ia][10][ip][id]));
               }
               acurr[ia] = B1 + (B2 - B1) * ratio;
          }
     }

     return;
}
// *******************************************************
double getMuonSpecCpp(int ip, double s, double r1, double d1, double e)
// *******************************************************
{
     double getMuonSpec = 0.0;
     const int nsol = 2;   // solar minimum and maximum
     const int nAdata = 7; // ! number of A parameter

     double acurr[nAdata + 1];
     double Fl[nsol + 1];
     static double spot[nsol + 1] = {0.0, 0.0, 150.0};
     static double restmass = 105.6; // rest mass of muon
     static double ethre = 3.0e5;    // threshold energy for high-energy muon correction

     double r, d, beta, tmp, powidx, A1, A2;
     int iptmp, is;

     if (e < 1.0e-2)
     {
          getMuonSpec = 0.0;
          return getMuonSpec;
     }

     r = max(1.0, r1);  // muon fluxes are the same for Rc < 1 GV
     d = max(0.15, d1); // secondary particle fluxes are the same for d < 0.15 g/cm2

     beta = sqrt(1.0 - pow(restmass / (restmass + e), 2));
     tmp = max(2.0, log10(e)); // below 100 MeV, this value should be constant

     for (is = 1; is <= nsol; is++)
     {
          iptmp = ip; // 1 for mu+, 2 for mu+
          getAmuon(iptmp, is, d, r, acurr);
          if (e > ethre)
          {                               // high energy correction
               acurr[5] = acurr[5] + 0.4; // high energy correction, see fit/muon/PowerCorrection
               acurr[1] = acurr[1] * pow(ethre, 0.4);
          }
          Fl[is] = acurr[1] * pow(e + (acurr[2] + acurr[4] * tmp) / pow(beta, acurr[3]), -acurr[5]) * (1 + exp(-acurr[6] * (log(e) + acurr[7])));
     }

     iptmp = ip + 5;                  // 6 for mu+, 7 for mu-
     powidx = getPowCpp(iptmp, d, r); // get Power index

     A2 = (Fl[1] - Fl[2]) / (pow(getFFPfromWCpp(spot[1]), powidx) - pow(getFFPfromWCpp(spot[2]), powidx));
     A1 = Fl[1] - A2 * pow(getFFPfromWCpp(spot[1]), powidx);
     getMuonSpec = A1 + A2 * pow(getFFPfromWCpp(s), powidx);

     return getMuonSpec;
}

// *******************************************************
double getSpecCpp(int ip, double s, double r, double d, double e, double g)
// ip: particle ID
// s: W index
// r: cut-off rigidity in GV
// d: atmospheric depth in g/cm2
// e: energy in MeV/n
// g: local geometry effect
// *******************************************************
{
     double getSpec = 0.0;
     int iptmp;

     if (ip == 0)
     { // neutron
          getSpec = getNeutSpecCpp(s, r, d, g, e);
     }
     else if (ip >= 1 && ip <= 28)
     { // proton to Ni
          getSpec = getIonSpecCpp(ip, s, r, d, e);
     }
     else if (ip == 29 || ip == 30)
     { // Muon
          iptmp = ip - 28;
          getSpec = getMuonSpecCpp(iptmp, s, r, d, e);
     }
     else
     {
          iptmp = ip - 28; // 3:electron, 4:positron, 5:photon
          getSpec = getSecondaryCpp(iptmp, s, r, d, e);
     }

     return getSpec;
}

// *******************************************************
double funcAngCpp(double x, double *a)
// *******************************************************
{
     double funcAng;

     static double atlow = -0.05; // lower threshold angle for 90 degree interpolation
     static double athig = 0.05;  // higher threshold angle for 90 degree interpolation

     double tmp1, tmp2;

     if (x <= atlow)
     { // backward
          funcAng = a[1] + a[2] * pow(abs(x), a[3]);
     }
     else if (x < athig)
     { // intermediate
          tmp1 = a[1] + a[2] * pow(abs(atlow), a[3]);
          tmp2 = a[4] + a[5] * pow(athig, a[6]);
          funcAng = tmp1 + (tmp2 - tmp1) * (x - atlow) / (athig - atlow); // simply interpolate
     }
     else if (x <= a[7])
     { // forward
          funcAng = a[4] + a[5] * pow(x, a[6]);
     }
     else
     { // after peak
          tmp1 = a[4] + a[5] * pow(a[7], a[6]);
          tmp2 = a[8];
          funcAng = tmp1 + (tmp2 - tmp1) * (x - a[7]) / (1.0 - a[7]);
     }

     return funcAng;
}

// *******************************************************
double getintCpp(double x0, double x1, double a1, double a2, double a3)
// *******************************************************
{
     double getint;
     getint = a1 * (x1 - x0) + a2 * (pow(x1, (a3 + 1)) - pow(x0, (a3 + 1))) / (a3 + 1);
     return getint;
}

// *******************************************************
void adjustParaAdepCpp(double *a)
// *******************************************************
{
     static double one = 1.0;
     static double atlow = -0.05; // lower threshold angle for 90 degree interpolation
     static double athig = 0.05;  // higher threshold angle for 90 degree interpolation

     double sum1, sum2, sum3, sum4, sum;

     sum1 = max(0.0, getintCpp(abs(atlow), one, a[1], a[2], a[3]));                          // -1.0 to atlow
     sum2 = max(0.0, (funcAngCpp(atlow, a) + funcAngCpp(athig, a)) * 0.5 * (athig - atlow)); // atlow to athig
     sum3 = max(0.0, getintCpp(athig, a[7], a[4], a[5], a[6]));                              // athig to a(7)
     if (a[7] < one)
     {
          sum4 = max(0.0, (funcAngCpp(a[7], a) + funcAngCpp(one, a)) * 0.5 * (one - a[7])); // a(7) to 1.0
     }
     else
     {
          sum4 = 0.0;
     }
     sum = max(1.0e-20, (sum1 + sum2 + sum3 + sum4) * 2.0 * acos(-1.0));

     a[1] = a[1] / sum;
     a[2] = a[2] / sum;
     a[4] = a[4] / sum;
     a[5] = a[5] / sum;
     a[8] = a[8] / sum;

     return;
}

// *******************************************************
double getParaAdepCpp(int i, double x, double *a)
// i: parameter index
// x: energy
// A: ParaEdep
// *******************************************************
{
     double getParaAdep;
     double a4, a7, a10;

     a4 = max(0.01, a[4]);
     a7 = max(0.01, a[7]);
     a10 = max(0.01, a[10]);

     getParaAdep = a[1] + a[2] / (1 + exp((a[3] - log10(x)) / a4)) + a[5] / (1 + exp((a[6] - log10(x)) / a7)) + a[8] / (1 + exp((a[9] - log10(x)) / a10));

     if (i == 3 || i == 6)
     {
          getParaAdep = pow(10.0, getParaAdep);
     } // a(3) & a(6) are determined by log10 fit
     if (i == 7)
     {
          getParaAdep = max(0.1, min(1.0, getParaAdep));
     } // a(7) should be 0.1 < a7 1.0
     if (i == 1 || i == 4 || i == 8)
     {
          getParaAdep = max(0.0, getParaAdep);
     } // a(1), a(4) & a(8) should be positive

     return getParaAdep;
}

// *******************************************************
double getGneutCpp(double emid, int id)
// *******************************************************
{
     double getGneut;
     const int maxfit = 8;  // number of fitting parameter for angular dependence
     const int maxEfit = 3; // number of fitting parameter for energy dependence
     static double Gneut[maxEfit + 1][maxfit + 1];
     static int ifirst = 0;

     string dname;
     string str;

     int i, ie;
     double elog;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/angle/NeutronGround.out";
          ifstream ifs(dname, ios::in);
          for (i = 1; i <= maxfit; i++)
          {
               getline(ifs, str);
               istringstream s(str);
               for (ie = 1; ie <= maxEfit; ie++)
               {
                    s >> Gneut[ie][i];
               }
          }
     }

     elog = max(-8.0, log10(emid)); // parameters are effective only above 0.01 eV
     getGneut = Gneut[1][id] / (1 + exp((elog - Gneut[2][id]) / Gneut[3][id]));

     return getGneut;
}

// *******************************************************
double BHfactorCpp(int ip, double e, double ang)
// *******************************************************
{
     double BHfactor;
     const int npart = 6;   // number of particle type (1:neutron, 2:proton, 3:heavy ion, 4:muon, 5:electron&positron, 6:photon
     const int nBHpara = 3; // number of parameter
     const int nBHeach = 7; // number of parameter to represent each parameter

     static double BHpara[nBHeach + 1][nBHpara + 1][npart + 1];
     double dimtmp[nBHpara + 1]; // temporary used dimension

     static string pname[npart + 1] = {"      ", "neutro", "proton", "he---4", "muon--", "elepos", "photon"};

     static int ifirst = 0;

     string dname;
     string str;

     int i, ii, ip2;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (ip2 = 1; ip2 <= npart; ip2++)
          {
               dname = "input/angle/BkH-" + pname[ip2] + ".inp";
               ifstream ifs(dname, ios::in);
               getline(ifs, str); // Kara-yomi
               for (i = 1; i <= nBHpara; i++)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    for (ii = 1; ii <= nBHeach; ii++)
                    {
                         s >> BHpara[ii][i][ip2];
                    }
               }
          }
     }

     if (ang < 0.0)
     {
          BHfactor = 0.0; // backward is always 0
     }
     else
     {
          for (i = 1; i <= nBHpara; i++)
          {
               dimtmp[i] = BHpara[1][i][ip] + BHpara[2][i][ip] / (1 + exp((BHpara[3][i][ip] - log10(e)) / BHpara[4][i][ip])) + BHpara[5][i][ip] / (1 + exp((BHpara[6][i][ip] - log10(e)) / BHpara[7][i][ip]));
          }
          BHfactor = dimtmp[1] + (dimtmp[2] - dimtmp[1]) * pow(ang, dimtmp[3]);
     }

     return BHfactor;
}

// *******************************************************
double getGmuonCpp(double emid, double ang)
// *******************************************************
{
     double getGmuon;
     double elog, a2, a3, a4, x, Scal;

     elog = max(4.062, log10(emid)); // parameters are effective only between 10GeV to 1 TeV
     a2 = 9.9873006E-01 + 2.9141114E+00 / (1.0 + exp((5.8030900E+00 - elog) / 2.4585039E-01));
     a3 = 2.4226042e1 - 1.5142933e1 * elog + 3.2012346 * pow(elog, 2) - 2.2325286e-1 * pow(elog, 3);
     elog = min(6.0, log10(emid));
     a4 = 1.4970401e1 - 5.3110524 * elog + 4.7458357e-1 * pow(elog, 2);

     if (ang < 0.0)
     {
          getGmuon = 0.0;
     }
     else
     {
          x = 1.0 / max(0.001, ang);
          Scal = a2 * sqrt(1 - exp(-(pow(1.0 / a2, 2))));
          getGmuon = (a2 * sqrt(1 - exp(-(pow(x / a2, 2)))) - (1 - exp(-a4 * pow(x - 1, a3)))) / Scal;
     }

     return getGmuon;
}

// *******************************************************
double getSpecAngCpp(int ip, double s, double r, double d, double e, double g, double ang)
// *******************************************************
{
     double getSpecAng;
     const int npart = 6;   // number of particle type (1:neutron, 2:proton, 3:heavy ion, 4:muon, 5:electron&positron, 6:photon
     const int nsur = 18;   // number of surface, upto 52 km
     const int ncor = 7;    // number of cut-off ridigity upto 20 GV
     const int maxfit = 8;  // number of parameters to express angular distribution
     const int mpeach = 10; // number of parameters to express energy dependence of each parameter
     const int ifour = 4;   // 4

     static double ParaEdep[mpeach + 1][maxfit + 1][ncor + 1][nsur + 1][npart + 1]; // parameter for expressing energy-differential energy dependence
     static double ParaEint[maxfit + 1][ncor + 1][nsur + 1][npart + 1];             // parameter for expressing energy-integrated energy dependence
     static double depth[nsur + 1];                                                 // depth (g/cm2) for parameters
     static double cor[ncor + 1];                                                   // Rc (GV) for parameters
     static double ParaAdep[maxfit + 1][ifour + 1];                                 // parameter for expressing angular distribution, 1-4 is for each depth & Rc condition
     static double ratio1, ratio2;                                                  // ratio must be saved
     static double emin[npart + 1] = {0.0, 1.0e-7, 1.0e0, 1.0e0, 1.0e1, 1.0e-1, 1.0e-2};
     static double emax[npart + 1] = {0.0, 1.0e4, 1.0e4, 1.0e4, 1.0e5, 1.0e4, 1.0e4};
     static double phimin = 1.0e-3;

     static string pname[npart + 1] = {"      ", "neutro", "proton", "he---4", "muon--", "elepos", "photon"};

     static int ifirst = 0;
     static int ipold = 0;
     static double sold = 0;
     static double rold = 0;
     static double dold = 0;
     static double eold = 0;
     static double gold = 0;

     double dimtmp[100]; // temporary used dimension

     string dname;
     string str;

     int i, ii, ip1, is, ic, is1, ic1, idx, itmp;
     double B1, B2, B3, B4, C1, C2, tmp, ene;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          for (ip1 = 1; ip1 <= npart; ip1++)
          {
               dname = "input/angle/" + pname[ip1] + ".out";
               ifstream ifs1(dname, ios::in);
               getline(ifs1, str); // Kara-yomi
               for (i = 1; i <= maxfit; i++)
               {
                    for (is = 1; is <= nsur; is++)
                    {
                         for (ic = 1; ic <= ncor; ic++)
                         {
                              getline(ifs1, str);
                              istringstream s1(str);
                              s1 >> itmp >> depth[is] >> cor[ic];
                              for (ii = 1; ii <= mpeach; ii++)
                              {
                                   s1 >> ParaEdep[ii][i][ic][is][ip1];
                              }
                         }
                    }
               }
               dname = "input/angle/" + pname[ip1] + "-Eint.out";
               ifstream ifs2(dname, ios::in);
               getline(ifs2, str); // Kara-yomi
               for (is = 1; is <= nsur; is++)
               {
                    for (ic = 1; ic <= ncor; ic++)
                    {
                         getline(ifs2, str);
                         istringstream s2(str);
                         s2 >> tmp >> tmp;
                         for (i = 1; i <= maxfit; i++)
                         {
                              s2 >> ParaEint[i][ic][is][ip1];
                         }
                    }
               }
          }
     }

     // check previous condition
     if (ip == ipold && s == sold && r == rold && d == dold && e == eold && g == gold)
     {
     } // same condition
     else
     {
          ipold = ip;
          sold = s;
          rold = r;
          dold = d;
          eold = e;
          gold = g;
          // find depth
          for (is = 1; is <= nsur; is++)
          {
               if (d < depth[is])
                    break;
          }
          if (is == 1)
          {
               ratio1 = 0.0;
               is = 2;
          }
          else if (is == nsur + 1)
          {
               ratio1 = 1.0;
               is = nsur;
          }
          else
          {
               ratio1 = (d - depth[is - 1]) / (depth[is] - depth[is - 1]);
          }
          // find Rc
          for (ic = 1; ic <= ncor; ic++)
          {
               if (r < cor[ic])
                    break;
          }
          if (ic == 1)
          {
               ratio2 = 0.0;
               ic = 2;
          }
          else if (ic == ncor + 1)
          {
               ratio2 = 1.0;
               ic = ncor;
          }
          else
          {
               ratio2 = (r - cor[ic - 1]) / (cor[ic] - cor[ic - 1]);
          }
          idx = 0;
          for (is1 = is - 1; is1 <= is; is1++)
          {
               for (ic1 = ic - 1; ic1 <= ic; ic1++)
               {
                    idx = idx + 1;
                    for (i = 1; i <= maxfit; i++)
                    {
                         if (e == 0.0)
                         { // energy integrated
                              ParaAdep[i][idx] = ParaEint[i][ic1][is1][ip];
                         }
                         else
                         {
                              ene = max(min(e, emax[ip]), emin[ip]);
                              //     ParaAdep[i][idx]=getParaAdepCpp(i,ene,ParaEdep[1][i][ic1][is1][ip]);
                              for (ii = 1; ii <= mpeach; ii++)
                              {
                                   dimtmp[ii] = ParaEdep[ii][i][ic1][is1][ip];
                              }
                              ParaAdep[i][idx] = getParaAdepCpp(i, ene, dimtmp);
                              if (ip == 1 && g >= 0.0 && g <= 1.0)
                                   ParaAdep[i][idx] = ParaAdep[i][idx] + getGneutCpp(e, i); // Ground level neutron correction, e instead of ene is used because data are down to 1e-8
                              if (i == 2 && (ParaAdep[1][idx] + ParaAdep[2][idx]) < phimin)
                                   ParaAdep[2][idx] = -ParaAdep[1][idx] + min(ParaAdep[1][idx], phimin); // avoid negative value at cos(theta)=-1.0
                              if (i == 5 && (ParaAdep[4][idx] + ParaAdep[5][idx]) < phimin)
                                   ParaAdep[5][idx] = -ParaAdep[4][idx] + min(ParaAdep[4][idx], phimin); // avoid negative value at cos(theta)=1.0
                         }
                    }
               }
          }
          for (idx = 1; idx <= ifour; idx++)
          {
               for (i = 1; i <= maxfit; i++)
               {
                    dimtmp[i] = ParaAdep[i][idx];
               }
               adjustParaAdepCpp(dimtmp); // integration value is adjusted to 1
               for (i = 1; i <= maxfit; i++)
               {
                    ParaAdep[i][idx] = dimtmp[i];
               }
          }
     }

     for (i = 1; i <= maxfit; i++)
     {
          dimtmp[i] = ParaAdep[i][1];
     }
     B1 = funcAngCpp(ang, dimtmp);
     for (i = 1; i <= maxfit; i++)
     {
          dimtmp[i] = ParaAdep[i][2];
     }
     B2 = funcAngCpp(ang, dimtmp);
     for (i = 1; i <= maxfit; i++)
     {
          dimtmp[i] = ParaAdep[i][3];
     }
     B3 = funcAngCpp(ang, dimtmp);
     for (i = 1; i <= maxfit; i++)
     {
          dimtmp[i] = ParaAdep[i][4];
     }
     B4 = funcAngCpp(ang, dimtmp);

     C1 = B1 + (B2 - B1) * ratio2; // Rc interpolation
     C2 = B3 + (B4 - B3) * ratio2; // Rc interpolation

     getSpecAng = C1 + (C2 - C1) * ratio1; // Depth interpolation

     return getSpecAng;
}

// *******************************************************
double getSpecAngFinalCpp(int ip, double s, double r, double d, double e, double g, double ang)
// ip: particle ID '1:neutro','2:proton','3:he---4','4:muon','5:elepos','6:photon'
// s: force field potential
// r: cut-off rigidity in GV
// d: atmospheric depth in g/cm2
// e: energy in MeV/n
// g: local geometry effect
// ang: cos(theta)
// *******************************************************
{
     double getSpecAngFinal = 0.0;
     double emin, ratio;

     if (ip == 4 && g >= 0 && g <= 1.0)
     { // ground level muon
          emin = 1.1535e4;
          if (e >= emin)
          {
               getSpecAngFinal = getSpecAngCpp(ip, s, r, d, e, g, 1.0) * getGmuonCpp(e, ang);
          }
          else
          {
               ratio = (getSpecAngCpp(ip, s, r, d, emin, g, 1.0) * getGmuonCpp(emin, ang)) / getSpecAngCpp(ip, s, r, d, emin, g, ang);
               getSpecAngFinal = getSpecAngCpp(ip, s, r, d, e, g, ang) * ratio;
          }
     }
     else
     {
          getSpecAngFinal = getSpecAngCpp(ip, s, r, d, e, g, ang);
     }

     if (g >= 100.0)
          getSpecAngFinal = getSpecAngFinal * BHfactorCpp(ip, e, ang); // Consider black hole

     return getSpecAngFinal;
}

// ************************************************
double getHPcpp(int iy0, int im0, int id0)
// ************************************************
{
     const int nmonth = 12;
     const int nday = 31;
     const int iymax = 2100;
     const int iymin = 1614;

     static double FFP[iymax + 1][nmonth + 1][nday + 1] = {}; // Daily force field potential
     static double FFPuso[iymax + 1] = {};                    // Annual force field potential

     static int ifirst = 0;
     static int iystart, iyend, iysUs, iyeUs;
     double getHP;
     int itmp;
     int ic; // originally used in main routine, but not used in c++ version
     //  ic=1: obtain FFP from daily FFP
     //  ic=2: obtain FFP from annual FFP
     //  ic=3: suspected ground level event
     //  ic=4: Too long time ago or future
     //  ic=5: no such date

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/FFPtable.day"; // Read daily force potential
          ifstream ifs(dname, ios::in);
          getline(ifs, str);
          istringstream s(str);
          s >> iystart >> iyend;
          for (int im = 1; im <= nmonth; im++)
          {
               for (int id = 1; id <= nday; id++)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    s >> itmp >> itmp;
                    for (int iy = iystart; iy <= iyend; iy++)
                    {
                         s >> FFP[iy][im][id];
                    }
                    // for (int iy=iystart; iy <= iyend; iy++) {cout << FFP[iy][im][id] << " ";}
               }
          }
          dname = "input/FFPtable.uso"; // Read annual force field potential
          ifstream ifs1(dname, ios::in);
          getline(ifs1, str);
          istringstream s1(str);
          s1 >> iysUs >> iyeUs;
          for (int iy = iysUs; iy <= iyeUs; iy++)
          {
               getline(ifs1, str);
               istringstream s1(str);
               s1 >> itmp >> FFPuso[iy];
               // cout << iy << " " << FFPuso[iy] << "\n";
          }
     }

     // **** Check Year, month, day *****************
     if ((iy0 < iymin) || (iy0 > iymax))
     {
          ic = 4;
          getHP = 0.0;
          cout << "Out of period"
               << "\n";
          return getHP;
     }
     if ((im0 < 1) || (im0 > nmonth) || (id0 < 1) || (id0 > nday))
     {
          ic = 5;
          getHP = 0.0;
          cout << "No such day"
               << "\n";
          return getHP;
     }

     // **********Determine FFP from Neutron Monitor *************
     if (FFP[iy0][im0][id0] > -99.0)
     { // data exist
          if ((iy0 >= iystart) && (FFP[iy0][im0][id0] == -1000.0))
          { // no such date
               ic = 5;
               getHP = 0.0;
          }
          else
          { // neutron monitor data exist
               getHP = FFP[iy0][im0][id0];
               if (getHP > 1000.0)
               {
                    ic = 3; // GLE occurred
                    getHP = getHP - 10000.0;
               }
               else
               {
                    ic = 1; // normal data
               }
               return getHP;
          }
     }

     //  ***** Determine FFP from Usoskin's data  *****
     if (FFPuso[iy0] != 0.0)
     {
          getHP = FFPuso[iy0];
          ic = 2; // determine FFP from Usoskin's data
          return getHP;
     }

     //  **** No data **************************
     ic = 4;
     getHP = 0.0;
     if (ic == 4)
     {
          cout << "No FFP data is available"
               << "\n";
     }
     return getHP;
}

// ******************************************************
double getrcpp(double cido, double ckei) // cido and ckei are center of ido&keido of each grid
// ******************************************************
{
     double getr;
     static double cordata[182][362];      // maximum 1 deg step, +1 mean
     static double dpido[182], dpkei[362]; // data point ido & keido
     static int mkei, mido;
     static double skei, sido;
     static int ifirst = 0;

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/CORdata.inp"; // Read vertical cut-off rigidity database
          ifstream ifs(dname, ios::in);
          getline(ifs, str);
          istringstream s(str);
          s >> mkei >> mido;         // read step size
          getline(ifs, str);         // skip 1 line
          skei = 360.0 / (mkei - 1); // keido step
          sido = 180.0 / (mido - 1); // ido step
          for (int id = mido; id >= 1; id--)
          {
               for (int ik = mkei; ik >= 1; ik--)
               {
                    getline(ifs, str);
                    istringstream s(str);
                    s >> dpkei[ik] >> dpido[id] >> cordata[id][ik]; // read vertical cut-off rigidity data
               }
          }
     }

     // **** Determine Cut-off Rigidity ***********
     if ((ckei > 180) && (ckei < 360))
     {
          ckei = ckei - 360.0;
     }                                                         // should be given in western longitude
     int id = min(mido - 1, (int)((cido + 90.0) / sido) + 1);  // lower ido bin
     int ik = min(mkei - 1, (int)((ckei + 180.0) / skei) + 1); // lower keido bin
     double cor1 = cordata[id][ik] * (dpido[id + 1] - cido) / (dpido[id + 1] - dpido[id]) + cordata[id + 1][ik] * (cido - dpido[id]) / (dpido[id + 1] - dpido[id]);
     double cor2 = cordata[id][ik + 1] * (dpido[id + 1] - cido) / (dpido[id + 1] - dpido[id]) + cordata[id + 1][ik + 1] * (cido - dpido[id]) / (dpido[id + 1] - dpido[id]);
     getr = cor1 * (dpkei[ik + 1] - ckei) / (dpkei[ik + 1] - dpkei[ik]) + cor2 * (ckei - dpkei[ik]) / (dpkei[ik + 1] - dpkei[ik]);

     return getr;
}

// ******************************************************
double getdcpp(double alti, double cido) // getd in g/cm^2, alti in km
// ******************************************************
{
     double getd;
     int iMSIS = 0;           // 0:US standard atmosphere, 1:NRLMSISE database
     const int maxUS = 75;    // number of altitude bin for US-Standard Air
     const int maxMSIS = 129; // number of altitude bin for NRLMSISE-00
     const int maxlat = 36;   // number of latitude bin for NRLMSISE-00

     static double altUS[maxUS + 1];                 // altitude data for US-Standard Air 1976
     static double altMSIS[maxMSIS + 1];             // altitude data for NRLMSISE-00
     static double depUS[maxUS + 1];                 // atmospheric depth data for US-Standard Air 1976
     static double depMSIS[maxMSIS + 1][maxlat + 1]; // atmospheric depth data for each altitude & latitude for NRLMSISE-00
     static double glat[maxlat + 1];                 // latitude data
     static int ifirst = 0;

     double ratio, ratio1, dep1, dep2;
     int ia, ido;

     string dname;
     string str;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/AtomDepth.inp"; // Read vertical cut-off rigidity database
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // kara-yomi
          for (ia = 1; ia <= maxUS; ia++)
          {
               getline(ifs, str);
               istringstream s(str);
               s >> altUS[ia] >> depUS[ia];
          }
          getline(ifs, str); // kara-yomi
          getline(ifs, str); // Ido data
          istringstream s(str);
          for (ido = 1; ido <= maxlat; ido++)
          {
               s >> glat[ido];
          }
          for (ia = 1; ia <= maxMSIS; ia++)
          { // read NRLMSISE-00 data
               getline(ifs, str);
               istringstream s(str);
               s >> altMSIS[ia];
               for (ido = 1; ido <= maxlat; ido++)
               {
                    s >> depMSIS[ia][ido];
               }
          }
     }

     if ((cido < -90.01 || cido > 90.01) || iMSIS == 0)
     { //! US standard atmosphere 1976 mode
          for (ia = 1; ia <= maxUS; ia++)
          {
               if (altUS[ia] > alti)
               {
                    break;
               }
          }
          if (ia == 1)
          { // out of range
               cout << "Error in function getd"
                    << "\n";
               cout << "Altitude = " << alti << " (km) should be higher than " << altUS[1] << " (km)"
                    << "\n";
               exit(1);
          }
          if (ia == maxUS + 1)
          { // ! out of range
               cout << "Warning in function getd"
                    << "\n";
               cout << "Altitude = " << alti << " (km) is too high. It is assumed to be " << altUS[maxUS] << " (km)"
                    << "\n";
               getd = depUS[maxUS];
               return getd;
          }
          ratio = (alti - altUS[ia - 1]) / (altUS[ia] - altUS[ia - 1]);
          getd = depUS[ia - 1] + ratio * (depUS[ia] - depUS[ia - 1]);
     }
     else
     { //  latitude is specified, so use NRLMSISE-00 data
          for (ia = 1; ia <= maxMSIS; ia++)
          {
               if (altMSIS[ia] > alti)
               {
                    break;
               }
          }
          if (ia == 1)
          { // out of range
               cout << "Error in function getd"
                    << "\n";
               cout << "Altitude = " << alti << " (km) should be higher than " << altMSIS[1] << " (km)"
                    << "\n";
               exit(1);
          }
          if (ia == maxMSIS + 1)
          { // ! out of range
               cout << "Warning in function getd"
                    << "\n";
               cout << "Altitude = " << alti << " (km) is too high. It is assumed to be " << altMSIS[maxMSIS] << " (km)"
                    << "\n";
               alti = altMSIS[maxMSIS];
               ia = maxMSIS;
          }
          ratio = (alti - altMSIS[ia - 1]) / (altMSIS[ia] - altMSIS[ia - 1]);
          for (ido = 2; ido <= maxlat - 1; ido++)
          {
               if (glat[ido] > cido)
               {
                    break;
               }
          }
          ratio1 = min(1.0, max(0.0, (cido - glat[ido - 1]) / (glat[ido] - glat[ido - 1])));
          dep1 = depMSIS[ia - 1][ido - 1] + ratio * (depMSIS[ia][ido - 1] - depMSIS[ia - 1][ido - 1]);
          dep2 = depMSIS[ia - 1][ido] + ratio * (depMSIS[ia][ido] - depMSIS[ia - 1][ido]);
          getd = dep1 + ratio1 * (dep2 - dep1);
     }
     return getd;
}

// ******************************************************
double get511fluxCpp(double s, double r, double d) // get 511 keV photon flux in (/cm2/s)
// ******************************************************
{
     double get511flux;
     const int ndep = 26;
     static double F511[ndep + 1], dep[ndep + 1];
     static int ifirst = 0;

     string dname;
     string str;

     int id, iptmp;
     double ratio, ene, Fratio;

     if (ifirst == 0)
     { // Initialization
          ifirst = 1;
          dname = "input/elemag/flux511keV.inp";
          ifstream ifs(dname, ios::in);
          getline(ifs, str); // Kara-yomi
          for (id = 1; id <= ndep; id++)
          {
               getline(ifs, str); // read data
               istringstream s(str);
               s >> dep[id] >> F511[id];
          }
     }

     // find depth ID
     for (id = 1; id <= ndep; id++)
     {
          if (d < dep[id])
               break;
     }
     if (id == 1)
     {
          ratio = 0.0;
          id = 2;
     }
     else if (id == ndep + 1)
     {
          ratio = 1.0;
          id = ndep;
     }
     else
     {
          ratio = (d - dep[id - 1]) / (dep[id] - dep[id - 1]);
     }

     Fratio = F511[id - 1] + (F511[id] - F511[id - 1]) * ratio;

     iptmp = 5;   // photon index
     ene = 0.511; // 511 keV
     get511flux = getSecondaryCpp(iptmp, s, r, d, ene) * Fratio;

     return get511flux;
}