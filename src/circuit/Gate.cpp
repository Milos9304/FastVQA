#include "Gate.h"

const std::string Parameter::blank_param_name="";
const std::shared_ptr<Parameter> Gate::blankParam=std::shared_ptr<Parameter>(new Parameter(Parameter::blank_param_name, 0.));

