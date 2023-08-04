#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

double getHPcpp(int, int, int);
double getrcpp(double, double);
double getdcpp(double, double);
double getSpecCpp(int, double, double, double, double, double);
double getSpecAngFinalCpp(int, double, double, double, double, double, double);
double getGenerationCpp(double, double, int *, int, double *, double *);

int main()
{
    const int npart = 33;
    static int IangPart[npart + 1] = {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 6};
    const int nebin = 100;                           // number of energy mesh (divided by log)
    const int nabin = 100;                           // number of angle mesh (divided by linear)
    static double ehigh[nebin + 1], emid[nebin + 1]; // higher and middle point of energy bin
    static double ahigh[nabin + 1], amid[nabin + 1]; // higher and middle point of angular bin
    static double etable[nebin + 1] = {};            // probability table (0.0 for 0, 1.0 for nebin)
    static double atable[nabin + 1][nebin + 1] = {}; // probability table (0.0 for 0, 1.0 for nabin)
    double atable2[nabin + 1];                       // temporary used dimension for anguluar probability
    static double Flux[nabin + 1][nebin + 1] = {};   // Monte Carlo generated flux

    int ia, ie, i, ia2;
    double e, phi, u, v, w, sx, cx, xd, yd, zd, x, y, z;

    // Set condition
    int nevent = 1000;    // number of particles to be generated
    int ip = 1;           // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
    int iyear = 2019;     // Year
    int imonth = 2;       // Month
    int iday = 1;         // Day
    double glat = 30.5;   // Latitude (deg), -90 =< glat =< 90
    double glong = -76.2; // Longitude (deg), -180 =< glong =< 180
    double alti = 0.0;    // Altitude (km)
    double g = 0.15;      // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
    double radi = 100.0;  // radius of the target area in cm (put your target inside this area)

    // Set energy and angle ranges for generation
    double emin = 1.0e0; // Minimum energy of particle
    double emax = 1.0e5; // Maximum energy of particle
    double amin = -1.0;  // Minimum cosine of particle
    double amax = 1.0;   // Maximum cosine of particle

    // calculate parameters
    double s = getHPcpp(iyear, imonth, iday); // solar modulation potential
    double r = getrcpp(glat, glong);          // Vertical cut-off rigidity (GV)
    double d = getdcpp(alti, glat);           // Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

    if (IangPart[ip] == 0)
    {
        cout << "Angular distribution is not available for the particle";
        exit(1);
    }

    if (ip == 0 && emin < 1.0e-8)
    {
        emin = 1.0e-8;
    } // Minimum energy for neutron is 10 meV
    if (ip != 0 && emin < 1.0e-2)
    {
        emin = 1.0e-2;
    } // Minimum energy for other particle is 10 keV

    // Make energy and angle mesh
    double elog = log10(emin);
    double estep = (log10(emax) - log10(emin)) / nebin;
    for (ie = 0; ie <= nebin; ie++)
    {
        ehigh[ie] = pow(10, elog);
        if (ie != 0)
            emid[ie] = sqrt(ehigh[ie] * ehigh[ie - 1]);
        elog = elog + estep;
    }

    double astep = (amax - amin) / nabin;
    for (ia = 0; ia <= nabin; ia++)
    {
        ahigh[ia] = amin + astep * ia;
        if (ia != 0)
            amid[ia] = (ahigh[ia] + ahigh[ia - 1]) * 0.5;
    }

    // Make probability table (absolute value)
    for (ie = 1; ie <= nebin; ie++)
    {
        for (ia = 1; ia <= nabin; ia++)
        {
            atable[ia][ie] = atable[ia - 1][ie] + getSpecCpp(ip, s, r, d, emid[ie], g) * getSpecAngFinalCpp(IangPart[ip], s, r, d, emid[ie], g, amid[ia]) * (2.0 * acos(-1.0)) * (ahigh[ia] - ahigh[ia - 1]); // angular integrated value
        }
    }
    for (ie = 1; ie <= nebin; ie++)
    {
        etable[ie] = etable[ie - 1] + atable[nabin][ie] * (ehigh[ie] - ehigh[ie - 1]); // energy integrated value
    }
    double TotalFlux = etable[nebin]; // Total Flux (/cm2/s), used for normalization

    // Make probability table (normalized to 1)
    for (ie = 1; ie <= nebin; ie++)
    {
        etable[ie] = etable[ie] / etable[nebin];
        for (ia = 1; ia <= nabin; ia++)
        {
            atable[ia][ie] = atable[ia][ie] / atable[nabin][ie];
        }
    }

    // Particle Generation
    mt19937 engine;                                     // Mersenne Twister
    uniform_real_distribution<double> rand01(0.0, 1.0); // random number between 0 to 1

    ofstream sf("GeneOut/generation.out", ios::out);
    sf << "ip= " << ip << " ,W-index= " << s << " ,Rc(GV)= " << r << " ,depth(g/cm2)= " << d << " ,g= " << g << "\n";
    sf << "Total Flux (/cm2/s)= " << TotalFlux << "\n";
    sf << "  Energy(MeV/n)              u              v              w              x              y              z\n";
    for (i = 1; i <= nevent; i++)
    {
        e = getGenerationCpp(rand01(engine), rand01(engine), &ie, nebin, ehigh, etable); // energy
        phi = 2.0 * acos(-1.0) * (rand01(engine) - 0.5);                                 // azimuth angle (rad)
        for (ia2 = 0; ia2 <= nabin; ia2++)
        {
            atable2[ia2] = atable[ia2][ie];
        }
        cx = getGenerationCpp(rand01(engine), rand01(engine), &ia, nabin, ahigh, atable2); // z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
        do
        {
            xd = (rand01(engine) - 0.5) * 2.0 * radi;
            yd = (rand01(engine) - 0.5) * 2.0 * radi;
        } while (sqrt(xd * xd + yd * yd) > radi);
        zd = radi;

        sx = sqrt(1 - cx * cx); // sin(theta)

        x = xd * cx * cos(phi) - yd * sin(phi) + zd * sx * cos(phi);
        y = xd * cx * sin(phi) + yd * cos(phi) + zd * sx * sin(phi);
        z = -xd * sx + zd * cx;
        u = -sx * cos(phi);
        v = -sx * sin(phi);
        w = -cx;
        sf << scientific << setw(15) << e << setw(15) << u << setw(15) << v << setw(15) << w << setw(15) << x << setw(15) << y << setw(15) << z << "\n";
        Flux[ia][ie] = Flux[ia][ie] + TotalFlux / (ehigh[ie] - ehigh[ie - 1]) / ((ahigh[ia] - ahigh[ia - 1]) * 2.0 * acos(-1.0)); // /cm2/s/sr/MeV
        Flux[0][ie] = Flux[0][ie] + TotalFlux / (ehigh[ie] - ehigh[ie - 1]);                                                      // /cm2/s/MeV
    }

    ofstream of("GeneOut/flux.out", ios::out);
    of << "ip= " << ip << " ,W-index= " << s << " ,Rc(GV)= " << r << " ,depth(g/cm2)= " << d << " ,g= " << g << "\n";
    of << "Total Flux (/cm2/s)= " << TotalFlux << "\n";
    of << "Angular and energy differential fluxes in /cm2/s/(MeV/n)/sr\n";
    of << "   E_low(MeV/n) /cm2/s/(MeV/n)";
    for (ia = 1; ia <= nabin; ia++)
    {
        of << setw(15) << ahigh[ia - 1];
    }
    of << "\n";
    for (ie = 1; ie <= nebin; ie++)
    {
        of << scientific << setw(15) << ehigh[ie - 1];
        for (ia = 0; ia <= nabin; ia++)
        {
            of << setw(15) << Flux[ia][ie] / nevent;
        }
        of << "\n";
    }

} // End of main

// ***********************************************************
double getGenerationCpp(double rand1, double rand2, int *ibin, int nbin, double *high, double *table)
// ***********************************************************
{
    double getGeneration = 0.0;

    int i;

    for (i = 1; i <= nbin - 1; i++)
    {
        if (rand1 <= table[i])
        {
            break;
        }
    }
    *ibin = i; // bin ID

    getGeneration = high[*ibin - 1] * rand2 + high[*ibin] * (1.0 - rand2);

    return getGeneration;
}
