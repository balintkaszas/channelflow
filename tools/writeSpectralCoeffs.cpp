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
#include <iomanip>

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
    string purpose("Output spectral coefficients of the given flowfield");
   
    ArgList args(argc, argv, purpose);

    const string input = args.getstr(1, "<flowfield>", "initial condition");
//    const string outdir = args.getpath("-o", "--outdir", "spectral/", "output directory");

    std::string base_filename = input.substr(input.find_last_of("/\\") + 1);
    std::cout << base_filename << '\n';

   ofstream outfile;
   std::string output = "Spectral_u_" + base_filename + ".dat";
   outfile.open(output);
   FlowField u(input);
   u.makeSpectral();

   for (int i=0; i<u.Nd(); ++i)
     for (int mx=0; mx<u.Mx(); ++mx) {
        int kx = u.kx(mx);
        for (int my=0; my<u.My(); ++my)
           for (int mz=0; mz<u.Mz(); ++mz) {
              int kz = u.kz(mz);
                outfile << kx << ' ' << my << ' '<< kz << ' '<< std::setprecision(16) <<  Re(u.cmplx(mx,my,mz,i)) << ' ' << Im(u.cmplx(mx, my, mz, i)) << '\n';
           }
     }
	// fill coefficients 
	// only the chebyshev position counts
	
    outfile.close();

 }
    cfMPI_Finalize();
}
