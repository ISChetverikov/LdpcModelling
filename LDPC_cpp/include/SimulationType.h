#pragma once

#include <unordered_map>

enum simulationType { MC, FFH, LFH, UnknownSimulation };

simulationType simulationTypeFromString(std::string str);
std::string simulationTypeToString(simulationType);