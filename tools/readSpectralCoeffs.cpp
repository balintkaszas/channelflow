/**
 * This file is a part of channelflow version 2.0 https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"
#include <fstream>
#include <string>
#include <sstream>
#include <string>
#include <iomanip>
using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
    string purpose("Read spectral coefficients from file and produce a flowfield");
   
    ArgList args(argc, argv, purpose);

    const string input = args.getstr(1, "<coeffs>", "initial condition");
//    const string outdir = args.getpath("-o", "--outdir", "spectral/", "output directory");


    const int Nx = args.getint("-Nx", "--Nx", "# x gridpoints");
    const int Ny = args.getint("-Ny", "--Ny", "# y gridpoints");
    const int Nz = args.getint("-Nz", "--Nz", "# z gridpoints");
    const Real alpha = args.getreal("-a", "--alpha", 0, "Lx = 2 pi/alpha");
    const Real gamma = args.getreal("-g", "--gamma", 0, "Lz = 2 pi/gamma");
        const Real lx = (alpha == 0.0) ? args.getreal("-lx", "--lx", 0.0, "Lx = 2 pi lx") : 1 / alpha;
        const Real lz = (gamma == 0.0) ? args.getreal("-lz", "--lz", 0.0, "Lz = 2 pi lz") : 1 / gamma;


    const Real Lx = (lx == 0.0) ? args.getreal("-Lx", "--Lx", "streamwise (x) box length") : 2 * pi * lx;
    const Real Lz = (lz == 0.0) ? args.getreal("-Lz", "--Lz", "spanwise   (z) box length") : 2 * pi * lz;
    const Real ymin = args.getreal("-ymin", "--ymin", -1, "lower wall height (y is wallnormal) ");
    const Real ymax = args.getreal("-ymax", "--ymax", +1, "upper wall height (y is wallnormal) ");

    const string symmstr = args.getstr("-symms", "--symmetries", "", "file of symmetries to satisfy");
    //const int numberofentries = 3*Ny*Nz*Nx;
    const int numberofmodes = Nx*Ny*Nz;
    std::string base_filename = input.substr(input.find_last_of("/\\") + 1);
    std::cout << base_filename << '\n';

   ifstream infile;
   std::string output = "FlowField_" + base_filename + ".nc";
   infile.open(input);
   //FlowField u(input);
    FlowField u(Nx, Ny, Nz, 3, Lx, Lz, ymin, ymax); //generate empty flowfield
    u.setPadded(true);

    SymmetryList s;
    if (symmstr.length() > 0) {
        s = SymmetryList(symmstr);
        cout << "Restricting random field to invariant subspace generated by symmetries" << endl;
        cout << s << endl;
    }
   //u.makeSpectral();
   std::string line;
    int readLines = 0;
    for(int i=0; i< numberofmodes; ++i){
        readLines++;
        int kx, my, kz;
        Real RealPart, ImagPart;
        infile >> kx >> my >> kz >> std::setprecision(16) >> RealPart >> ImagPart;
        if(i == 996){
            std::cout << i << '\t' <<  line << '\n';
            std::cout << u.kx(u.mx(kx)) <<  '\t' << my << '\t' << '\t' << u.kz(u.mz(kz))<< ' ' << std::setprecision(16) << RealPart << '\t' << ImagPart << '\n';
        } 
        if(kz>=0){
            u.cmplx(u.mx(kx), my, u.mz(kz), 0) = RealPart + ImagPart*I;
        }
    }
    for(int i=0; i< numberofmodes; ++i){
        readLines++;
        int kx, my, kz;
        Real RealPart, ImagPart;
        infile >> kx >> my >>  kz >>std::setprecision(16) >>RealPart >> ImagPart;
        if(kz>=0){
            u.cmplx(u.mx(kx), my, u.mz(kz), 1) = RealPart + ImagPart*I;
        }
    }

    for(int i=0; i< numberofmodes; ++i){
        readLines++;
        int kx, my, kz;
        Real RealPart, ImagPart;
        infile >> kx >> my >> kz >> std::setprecision(16) >> RealPart >> ImagPart;
        if(readLines == 3018){
            std::cout << i << '\t' <<  line << '\n';
            std::cout << u.kx(u.mx(kx)) <<  '\t' << my << '\t' << '\t' << u.kz(u.mz(kz))<< ' ' << std::setprecision(16) << RealPart << '\t' << ImagPart << '\n';
        } 


        if(kz>=0){
            u.cmplx(u.mx(kx), my, u.mz(kz), 2) = RealPart + ImagPart*I;
        }
    }
	// fill coefficients 
	// only the chebyshev position counts
	std::cout << "Read " << readLines << " lines \n";
    std::cout << "Read " << numberofmodes << " lines \n";
    std::cout << u.cmplx(u.mx(-1), 2, u.mz(0), 2) << '\n';
    //u.setPadded(true);

    u.makePhysical();
    u.makeSpectral();
    std::cout << u.cmplx(u.mx(-1), 2, u.mz(0), 2) << '\n';
    infile.close();
    u.save(output);
 }
    cfMPI_Finalize();
}
