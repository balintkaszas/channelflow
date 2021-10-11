/**
 * This file is a part of channelflow version 2.0 https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"
#include "channelflow/diffops.h"
#include "channelflow/poissonsolver.h"

#include <fstream>
#include <string>

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
    string purpose("Check the residual of a given solution, when substituted into the NS equations");

    ArgList args(argc, argv, purpose);
    const Real Re = args.getreal("-R", "--Reynolds", 400, "Reynolds number");

    const string input = args.getstr(1, "<flowfield>", "initial condition");

 // Add the possibility to change Reynolds number and input file with command line arguments


   FlowField u(input); // this is taken from the channelflow database
   //FlowField p("p0.nc");
   u.makeSpectral();
   //p.makeSpectral();
   

  //evaluate the linear term first.
  // need to add the base flow
  //
  //i
  //
  //
  //

  DNSFlags flags;
  string Uname, Wname;
  DNSFlags baseflags = setBaseFlowFlags(args, Uname, Wname); 
   //BaseFlow b = LinearBase;
   //flags.b
  FlowField uBase(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b()); //Initialize to 0

  uBase.makeSpectral();
  Real Reyn = Re;
  Real nu = 1./Reyn;
  vector<ChebyCoeff> base_Flow = baseFlow(u.Ny(), u.a(), u.b(), baseflags, Uname, Wname); //Create the base flow from cmd line argument
  uBase += base_Flow;
  std::cout << "Generated base flow \n";
  FlowField Residual(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());
    Residual.makeSpectral();
  FlowField uSum = u + uBase; // u = u_pert + u_base
  lapl(uSum, Residual);
  Residual *= nu; //1/Re * Lapl(u)

 //Compute the pressure explicitly here: copied from pressure.cpp
      NonlinearMethod nl = Convection; //poissonsolver says that this nonlinearity works best for the pressure. 


      // compute pressure
      PressureSolver poisson(uSum, nu, baseflags.Vsuck,nl);

      FlowField q = poisson.solve(uSum);

      q.setPadded(true);
     // q.save(pname);
  FlowField nonlinear(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());
  nonlinear.makeSpectral();


  FlowField tmp(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());
  tmp.makeSpectral();

  dotgrad(uSum, uSum,nonlinear,tmp); // (u \cdot \nabla) u

  FlowField gradP(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b());
  gradP.makeSpectral();

  grad(q, gradP); //nabla P

  Residual += -1*nonlinear;
  Residual += -1*gradP;  

   std::cout << "Generated nonlin terms  \n";
    std::cout << "|f(u)| = " << L2Norm(Residual, true) << '\n';
    }
   cfMPI_Finalize();
}


