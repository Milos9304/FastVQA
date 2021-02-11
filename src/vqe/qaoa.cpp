#include <xacc.hpp>
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
   auto acc = xacc::getAccelerator("qpp");

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
}
