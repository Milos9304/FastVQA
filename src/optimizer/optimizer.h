/*
 * optimizer.h
 *
 *  Created on: Apr 15, 2021
 *      Author: Milos Prokop
 */

#ifndef SRC_OPTIMIZER_H_
#define SRC_OPTIMIZER_H_

#include <functional>
#include <string>
#include "../io/logger.h"

using OptimizerFunctorNoGrad =
    std::function<double(const std::vector<double> &)>;
using OptimizerFunctor =
    std::function<double(const std::vector<double> &, std::vector<double> &)>;
using OptResult = std::pair<double, std::vector<double>>;

using OptFunctionPtr = double (*)(const std::vector<double> &,
                                  std::vector<double> &, void *);
class OptFunction {
protected:
  OptimizerFunctor _function;
  int _dim = 0;

public:
  OptFunction() = default;

  // Standard constructor, takes function that takes params as
  // first arg and gradient as second arg
  OptFunction(OptimizerFunctor f, const int d) : _function(f), _dim(d) {}
  OptFunction(OptimizerFunctorNoGrad f, const int d)
      : _function([&](const std::vector<double> &x, std::vector<double> &) {
          return f(x);
        }),
        _dim(d) {}
  // OptFunction(OptimizerFunctorNoGradValue f, const int d)
  //     : _function([&](const std::vector<double> &x, std::vector<double> &) {
  //         return f(x);
  //       }),
  //       _dim(d) {}
  virtual const int dimensions() const { return _dim; }
  virtual double operator()(const std::vector<double> &x,
                            std::vector<double> &dx) {
    return _function(x, dx);
  }
  virtual double operator()(const std::vector<double> &&x) {
    std::vector<double> dx;
    return _function(x, dx);
  }
};


class Optimizer{
public:
	virtual OptResult optimize(OptFunction &function) = 0;
	virtual const std::string get_algorithm() const { return ""; }
	virtual const bool isGradientBased() const { return false; }
};

struct ExtraNLOptData {
    std::function<double(const std::vector<double>&, std::vector<double>&)> f;
};

class NLOptimizer : public Optimizer {
public:
  OptResult optimize(OptFunction &function) override;
  const bool isGradientBased() const override;
  virtual const std::string get_algorithm() const;
};

#endif /* SRC_OPTIMIZER_H_ */
