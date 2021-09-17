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
#include <string>

#include <algorithm>

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
    std::vector<FlowField> UB = readHeteroclinicOrbit("../exactHeteroclinic/LBtoUB/u", 0, 3000, 10);
    std::vector<FlowField> UC = readHeteroclinicOrbit("../exactHeteroclinic/LBtoLam/u", 0, 2600, 10);
   std::ofstream file("data/distFromHet.dat");
   for(int i =0; i<4000; i++) {
      FlowField u("data/u" + std::to_string(i));
      Real distUB = distFromHeteroclinic(UB, u);
      Real distUC = distFromHeteroclinic(UC, u);
      file << i << '\t' << std::min(distUB, distUC) << std::endl;
      std::cout << i << std::endl;
   }
   file.close();
    }
    cfMPI_Finalize();
}
