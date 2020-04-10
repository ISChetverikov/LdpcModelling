#include "../include/SimulationType.h"

simulationType simulationTypeFromString(std::string str){
	std::unordered_map<std::string, simulationType> simulationTypeResolver = {
		{"MC", MC},
		{"FFH", FFH},
		{"LFH", LFH}
	};

	if (simulationTypeResolver.count(str) > 0)
		return simulationTypeResolver[str];
	
	return UnknownSimulation;
}