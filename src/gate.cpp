#include "gate.h"
#include "logger.h"

namespace FastVQA{

const std::string Parameter::blank_param_name="";
const std::shared_ptr<Parameter> Gate::blankParam=std::shared_ptr<Parameter>(new Parameter(Parameter::blank_param_name, 0.));

Gate::Gate(GateCode code, int qubit1, std::shared_ptr<Parameter> param){
		if(code & two_qubit){
		  loge(std::to_string(code) + "is a two qubit gate");
		  throw;
		}

		this->code=code;
		this->qubit1=qubit1;
		this->qubit2=-1;

		this->param=param;
}

Gate::Gate(GateCode code, int qubit1, int qubit2, std::shared_ptr<Parameter> param){
		if((!(code & two_qubit)) && qubit2 != -1){
		  loge(std::to_string(code) + "is a one qubit gate");
		  throw;
		}

		this->code=code;
		this->qubit1=qubit1;
		this->qubit2=qubit2;
		this->param=param;
}
}
