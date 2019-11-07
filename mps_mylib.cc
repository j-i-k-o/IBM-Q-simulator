#include "itensor/all.h"
#include "itensor/mps/siteset.h"
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "ibmq.h"

using namespace itensor;
using namespace ibmq;

int main(int argc, char const* argv[]){

    IBMQPeps circuit;

    circuit.PrintMat();
    circuit.PrintCursor();

    std::cout << "---------------" << std::endl;

    circuit.shift_to(6, {"Cutoff", 1E-8});
    circuit.PrintMat();
    circuit.PrintCursor();

    std::cout << "---------------" << std::endl;

    circuit.decomposePsi({"Cutoff", 1E-8});
    circuit.PrintMat();
    circuit.PrintCursor();

	return 0;
}
