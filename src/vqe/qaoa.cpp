#include "xacc.hpp"

#include <iostream>

void run_dummy_qaoa(){
    xacc::Initialize();
    xacc::set_verbose(true);
    xacc::setOption("quest-verbose", "true");
    auto qpu = xacc::getAccelerator("quest", {
        std::make_pair("quest-visitor", "exatn-mps"),
        std::make_pair("shots", 10),
    });
    //xacc::getAccelerator("quest");
    // Allocate a register of 40 qubits
    auto qubitReg = xacc::qalloc(3);

    // Create a Program
    auto xasmCompiler = xacc::getCompiler("xasm");
    auto ir = xasmCompiler->compile(R"(__qpu__ void test(qbit q) {
        H(q[0]);
		X(q[1]);
        Measure(q[0]);
        Measure(q[1]);
    })", qpu);

    // Request the quantum kernel representing
    // the above source code
    auto program = ir->getComposite("test");
    // Execute!
    qpu->execute(qubitReg, program);

    for(auto &m : qubitReg->getMeasurementCounts())//print();
    	std::cout << "oo: " << std::get<0>(m) <<" " <<std::get<1>(m)<< std::endl;

    // Finalize the XACC Framework
    xacc::Finalize();

}
/*#include <xacc.hpp>
#include "xacc_observable.hpp"
#include "xacc_service.hpp"
#include <random>

#include "../logger.h"

// Use XACC built-in QAOA to solve a QUBO problem
// QUBO function:
// y = -5x1 - 3x2 - 8x3 - 6x4 + 4x1x2 + 8x1x3 + 2x2x3 + 10x3x4

void run_dummy_qaoa(){
	//int main(int argc, char **argv) {
	//xacc::Initialize(argc, argv);
	xacc::Initialize();
	std::cout<<"init\n";

 // Use the Qpp simulator as the accelerator
   auto acc = xacc::getAccelerator("quest");

   auto buffer = xacc::qalloc(4);
   // The corresponding QUBO Hamiltonian is:
   auto observable = xacc::quantum::getObservable(
         "pauli",
         std::string("-5.0 - 0.5 Z0 - 1.0 Z2 + 0.5 Z3 + 1.0 Z0 Z1 + 2.0 Z0 Z2 + 0.5 Z1 Z2 + 2.5 Z2 Z3"));

   const int nbSteps = 1;//12;
   const int nbParams = nbSteps*11;
   std::vector<double> initialParams;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<> dis(-2.0, 2.0);

   // Init random parameters
   for (int i = 0; i < nbParams; ++i)
   {
      initialParams.emplace_back(dis(gen));
   }

   auto optimizer = xacc::getOptimizer("nlopt",
      xacc::HeterogeneousMap {
         std::make_pair("initial-parameters", initialParams),
         std::make_pair("nlopt-maxeval", nbParams*100) });

   auto qaoa = xacc::getService<xacc::Algorithm>("QAOA");

   const bool initOk = qaoa->initialize({
                           std::make_pair("accelerator", acc),
                           std::make_pair("optimizer", optimizer),
                           std::make_pair("observable", observable),
                           // number of time steps (p) param
                           std::make_pair("steps", nbSteps)
                        });
   qaoa->execute(buffer);
   std::cout << "Min QUBO: " << (*buffer)["opt-val"].as<double>() << "\n";
}*/
