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

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {



   ofstream outfileLin;
   outfileLin.open("SpectralCoeffs_Lin_v_p_Re300.dat");
   ofstream outfileNLin;

   outfileNL.open("SpectralCoeffs_NLin_v_p_Re300.dat");


   FlowField u("ubest.nc");
   FlowField p("p0.nc");
   u.makeSpectral();
   p.makeSpectral();
   

  //evaluate the linear term first.
  // need to add the base flow
  //
  //i
  //
  //
  //
  string purpose("Calculate spectral coefficients of the NS equations");

  ArgList args(argc, argv, purpose);
  DNSFlags flags;
  string Uname, Wname;
  DNSFlags baseflags = setBaseFlowFlags(args, Uname, Wname); 

  FlowField uBase;
  uBase.makeSpectral();

  vector<ChebyCoeff> base_Flow = baseFlow(u.Ny(), u.a(), u.b(), flags, Uname, Wname);
  uBase += base_Flow;

  FlowField LinearTerm(u);
  lapl(u, LinearTerm);
  LinearTerm *= nu;

  LinearTerm.makeSpectral();

  FlowField nonlin1;
  FlowField nonlin2;
  FlowField tmp(u);
  dotgrad(u, uBase,nonlin1,tmp);

  FlowField tmp2(u);
  dotgrad(uBase, u, nonlin2, tmp2);
  FlowField gradP(u);
  grad(p, gradP);
  linearTerm += -1*nonlin1;
  linearTerm += -1*nonlin2;
  linearTerm += -1*gradP;  

  FlowField Nonlin(u);
  FlowField tmp3(u);
  Nonlin.makeSpectral();
  dotgrad(u, u, Nonlin, tmp3);

  
   for (int i=0; i<Nonlin.Nd(); ++i)
     for (int mx=0; mx<Nonlin.Mx(); ++mx) {
        int kx = Nonlin.kx(mx);
        for (int my=0; my<Nonlin.My(); ++my)
           for (int mz=0; mz<Nonlin.Mz(); ++mz) {
              int kz = Nonlin.kz(mz);
outfileNL << kx << ' ' << my << ' '<< kz << ' '<< Re(Nonlin.cmplx(mx,my,mz,i)) << ' ' << Im(Nonlin.cmplx(mx, my, mz, i)) << '\n';
           }
     }
	// fill coefficients 
	// only the chebyshev position counts
	
    outfileNL.close();
    }
     for (int i=0; i<LinearTerm.Nd(); ++i)
     for (int mx=0; mx<LinearTerm.Mx(); ++mx) {
        int kx = Nonlin.kx(mx);
        for (int my=0; my<LinearTerm.My(); ++my)
           for (int mz=0; mz<LinearTerm.Mz(); ++mz) {
              int kz = LinearTerm.kz(mz);
outfileLin << kx << ' ' << my << ' '<< kz << ' '<< Re(LinearTerm.cmplx(mx,my,mz,i)) << ' ' << Im(LinearTerm.cmplx(mx, my, mz, i)) << '\n';
           }
     }
	// fill coefficients 
	// only the chebyshev position counts
	
    outfileLin.close();
    }
   cfMPI_Finalize();
}
