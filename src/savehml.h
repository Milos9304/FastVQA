
void saveHml(Expression* qubo, std::string filename){
	ofstream hmlfile(filename);
	logw(filename);
	hmlfile<<qubo->variables.size()-1<<"\n";
	for(auto &var:qubo->variables){
		if(var->name != "id")
			hmlfile<<var->name<<"\n";
	}

	for(auto &term:qubo->polynomial){

		if(term.second == 0)
			continue;

		std::string id1 = qubo->idMap[term.first.first]->name;

		if(term.first.second == -1){ //id
			hmlfile<<term.second<<"\n";
			continue;
		}
		std::string id2 = qubo->idMap[term.first.second]->name;
		if(id1 == "id")
			hmlfile<<term.second<<" "<<id2<<"\n";
		else
			hmlfile<<term.second<<" "<<id1<<" "<<id2<<"\n";

	}


	hmlfile.close();

}
