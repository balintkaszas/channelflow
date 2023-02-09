/**
 * This file is a part of channelflow version 2.0 https://channelflow.ch.
 * License is GNU GPL version 2 or later: ./LICENCE
 */

#include <channelflow/laurettedsi.h>
#include "channelflow/cfdsi.h"
#include "nsolver/dsi.h"

#include "channelflow/flowfield.h"
#include "nsolver/nsolver.h"

using namespace std;
using namespace Eigen;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        ArgList args(argc, argv, "find an invariant solution using Newton-Krylov-hookstep algorithm");

        /** Choose the Newton algorithm to be used. Currently, two options are available: simple Newton without any
         * trust region optimization, and Newton with Hookstep (default). For the simple Newton, you can choose either a
         * full-space algorithm to solve the Newton equations (-solver "eigen") or between the two iterative algorithms
         * GMRES and BiCGStab. Newton-Hookstep requires GMRES. Note that the available parameters depend on your choice
         * of the algorithm.
         */

        unique_ptr<Newton> N;
        NewtonSearchFlags searchflags(args);
        searchflags.save(searchflags.outdir);
        N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));
        WriteProcessInfo(argc, argv);
        //dnsflags.save();

        //CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);


        //create a lambda: function to evaluate 
        //VectorXd lambd = [] (const VectorXd & x) {return 0.*x;};
        /** Construct the dynamical-systems interface object depending on the given parameters. Current options are
         * either standard (f(u) via forward time integration) or Laurette (f(u) via Laurettes method)
         */
        unique_ptr<DSI> dsi;
        //dsi = unique_ptr<DSI>(new DSI());

      VectorXd x0(10);
      setToZero(x0);
      Real residual = 0;

      N->solve(*dsi, x0, residual);
    }

    cfMPI_Finalize();

    return 0;
}
