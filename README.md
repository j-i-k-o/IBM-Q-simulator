#IBM Q simulator 

An IBM Q simulator of 53-qubit quantum processor (rochester) by using tensor contraction (tensor network)

##Installation

Install [ITensor](https://itensor.org/) and change `LIBRARY_DIR` in Makefile to your ITensor directory.

##Example

```c++
#include "itensor/all.h"
#include "itensor/mps/siteset.h"
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "ibmq.h"

using namespace itensor;
using namespace ibmq;

int main(int argc, char const* argv[]){

    std::array<std::pair<double, double>, IBMQPeps::NUM_BITS> init_param;

    //start with |0> (53qubit)
    for(auto& elem : init_param){
        elem = std::make_pair(1.0, 0.0);
    }

    IBMQPeps circuit(init_param);


    //Below is the demonstration of genearating bell state

    //the default cursor is located on qubit number 1 and 2
    //shift the cursor to qubit number 5
    circuit.shift_to(3);
    circuit.shift_to(4);
    circuit.shift_to(5);

    //apply Hadamard to gate 4
    circuit.apply(H(circuit.site(4))*Id(circuit.site(5)));
    //apply CNOT to gate 4 and 5
    circuit.apply(CNOT(circuit.site(4), circuit.site(5)));
    //the result should be bell state (1/sqrt(2))(|00> + |11>)
    


    //to show that there is bell state, calc the overlap between |0...00....0> and |0...11....0> where 00 and 11 are located on the qubit number 4 and 5 respectively.
    
    //|0...00....0>
    IBMQPeps circuit00(init_param, circuit.site());

    //|0...11....0>
    IBMQPeps circuit11(init_param, circuit.site());
    //shift the cursor to qubit number 5
    circuit11.shift_to(3);
    circuit11.shift_to(4);
    circuit11.shift_to(5);
    //flip the qubit number 4 and 5
    circuit11.apply(X(circuit.site(4))*X(circuit.site(5)));
    
    myvector<ITensor> op(IBMQPeps::NUM_BITS);
    for(auto i : range1(IBMQPeps::NUM_BITS)){
        op[i] = Id(circuit.site(i));
    }

    Print(overlap(circuit, op, circuit00)); //result should be 1/sqrt(2)
    Print(overlap(circuit, op, circuit11)); //result should be 1/sqrt(2)

	return 0;
}
```


