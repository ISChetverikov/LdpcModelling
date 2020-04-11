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

std::string simulationTypeToString(simulationType type) {
	std::unordered_map<simulationType, std::string> simulationTypeStringResolver = {
		{MC, "MC"},
		{FFH, "FFH"},
		{LFH, "LFH"},
		{UnknownSimulation, "UnknownSimulation"}
	};

	return simulationTypeStringResolver[type];
}

