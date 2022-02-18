/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#include "optimizer.h"
#include "nlopt.hpp"
#include "logger.h"

#include <iostream>
using namespace std::placeholders;

double c_wrapper(const std::vector<double> &x, std::vector<double> &grad,
                 void *extra) {
  auto e = reinterpret_cast<fastVQA::ExtraNLOptData *>(extra);
  return e->f(x, grad);
}

const std::string fastVQA::NLOptimizer::get_algorithm() const {
  std::string optimizerAlgo = "cobyla";
  /*if (options.stringExists("algorithm")) {
    optimizerAlgo = options.getString("algorithm");
  }
  if (options.stringExists("nlopt-optimizer")) {
    optimizerAlgo = options.getString("nlopt-optimizer");
  }*/
  return optimizerAlgo;
}

const bool fastVQA::NLOptimizer::isGradientBased() const {

  std::string optimizerAlgo = "cobyla";
  /*if (options.stringExists("algorithm")) {
    optimizerAlgo = options.getString("algorithm");
  }
  if (options.stringExists("nlopt-optimizer")) {
    optimizerAlgo = options.getString("nlopt-optimizer");
  }*/

  if (optimizerAlgo == "l-bfgs") {
    return true;
  } else {
    return false;
  }
}

fastVQA::OptResult fastVQA::NLOptimizer::optimize(OptFunction &function, std::vector<double> x, double tol, int maxeval, std::vector<double> lowerBounds, std::vector<double> upperBounds) {

  auto dim = function.dimensions();
  nlopt::algorithm algo = nlopt::algorithm::LN_COBYLA;

  bool maximize = false;

  std::string optimizerAlgo;
  /*if (options.stringExists("algorithm")) {
    optimizerAlgo = options.getString("algorithm");
  }
  if (options.stringExists("nlopt-optimizer")) {
    optimizerAlgo = options.getString("nlopt-optimizer");
  }

  if (options.stringExists("nlopt-optimizer") || options.stringExists("algorithm")) {
    if (optimizerAlgo == "cobyla") {
      algo = nlopt::algorithm::LN_COBYLA;
    } else if (optimizerAlgo == "nelder-mead") {
      algo = nlopt::algorithm::LN_NELDERMEAD;
    } else if (optimizerAlgo == "l-bfgs") {
      algo = nlopt::algorithm::LD_LBFGS;
    } else {
      loge("Invalid optimizer at this time: " + optimizerAlgo);
    }
  }*/

    /*if (options.keyExists<std::vector<double>>("initial-parameters")) {
    x = options.get_with_throw<std::vector<double>>("initial-parameters");
  } else if (options.keyExists<std::vector<int>>("initial-parameters")) {
    auto tmpx = options.get<std::vector<int>>("initial-parameters");
    x = std::vector<double>(tmpx.begin(), tmpx.end());
  }*/

  ExtraNLOptData data;
  data.f = function;
  auto d = reinterpret_cast<void *>(&data);

  nlopt::opt _opt(algo, dim);
  if (maximize) {
    _opt.set_max_objective(c_wrapper, d);
  } else {
    _opt.set_min_objective(c_wrapper, d);
  }
  // Default lower bounds
  //std::vector<double> lowerBounds(dim, -3.1415926);
  /*if (options.keyExists<std::vector<double>>("lower-bounds")) {
    lowerBounds = options.get<std::vector<double>>("lower-bounds");
  }
  if (options.keyExists<std::vector<double>>("nlopt-lower-bounds")) {
    lowerBounds = options.get<std::vector<double>>("nlopt-lower-bounds");
  }*/

  // Default upper bounds
  //std::vector<double> upperBounds(dim, 3.1415926);
  /*if (options.keyExists<std::vector<double>>("upper-bounds")) {
    upperBounds = options.get<std::vector<double>>("upper-bounds");
  }
  if (options.keyExists<std::vector<double>>("nlopt-upper-bounds")) {
    upperBounds = options.get<std::vector<double>>("nlopt-upper-bounds");
  }*/

  _opt.set_lower_bounds(lowerBounds);
  _opt.set_upper_bounds(upperBounds);
  _opt.set_maxeval(maxeval);
  _opt.set_ftol_rel(tol);

  /*if (options.keyExists<double>("stopval")) {
    const double stopVal = options.get<double>("stopval");
    logi("[NLOpt] function stopval set to " + std::to_string(stopVal));
    _opt.set_stopval(stopVal);
  }
  if (options.keyExists<double>("nlopt-stopval")) {
    const double stopVal = options.get<double>("nlopt-stopval");
    logi("[NLOpt] function stopval set to " + std::to_string(stopVal));
    _opt.set_stopval(stopVal);
  }*/

  if (dim != x.size()) {
    loge("Invalid optimization configuration: function dim == " +
                std::to_string(dim) +
                ", param_size == " + std::to_string(x.size()));
  }
  double optF;
  nlopt::result r;
  try {
    r = _opt.optimize(x, optF);
  } catch (std::exception &e) {
    loge("NLOpt failed with error code = " + std::to_string(r) + ", " +
                std::string(e.what()));
  }

  return OptResult{optF, x};
}
