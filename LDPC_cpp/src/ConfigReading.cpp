#include <filesystem>
#include <unordered_map>
#include "../include/SimulationParameters.h"
#include "../include/ConfigReading.h"
#include "../lib/pugixml-1.10/src/pugixml.hpp"
#include "../lib/pugixml-1.10/src/pugiconfig.hpp"
#include "../include/Exceptions.h"
#include "../include/SimulationType.h"
#include "../include/DecoderType.h"

// if not exist then exception
pugi::xml_node GetChildNode(pugi::xml_node node, std::string childName) {
	auto child = node.child(childName.c_str());
	if (child == NULL)
		throw ConfigParseException("Parse config error: can not find \"" + childName + "\" section");

	return child;
}

// if not exists then returns NULL
pugi::xml_node TryGetChildNode(pugi::xml_node node, std::string childName) {
	auto child = node.child(childName.c_str());
	
	return child;
}

pugi::xml_attribute TryGetAttribute(pugi::xml_node node, std::string name) {
	auto attr = node.attribute(name.c_str());
	if (attr == NULL)
		throw ConfigParseException("Parse config error: can not find \"" + name + "\" atribute");

	return attr;
}

SimulationParams ReadSimulationSection(pugi::xml_node simulation_node) {
	SimulationParams params;

	params.H_MatrixFilename = GetChildNode(simulation_node, "MatrixFilename").child_value();
	params.resultsFilename = GetChildNode(simulation_node, "ResultsFilename").child_value();

	pugi::xml_node simulator_node = GetChildNode(simulation_node, "Simulator");

	std::string simulationTypeStr = TryGetAttribute(simulator_node, "type").value();
	simulationType simulationType_ = simulationTypeFromString(simulationTypeStr);
	if (simulationType_ == simulationType::UnknownSimulation)
		throw ConfigParseException("Parse config error: unsupported value of simulation type: " + simulationTypeStr);

	params.type = simulationType_;

	std::unordered_map <std::string, std::string > simulationParams;
	for (pugi::xml_node simulationParam_node : simulator_node.children("SimulationParam"))
	{
		simulationParams[TryGetAttribute(simulationParam_node, "name").value()] = simulationParam_node.child_value();
	}
	params.simulationTypeParams = simulationParams;

	pugi::xml_node decoder_node = GetChildNode(simulation_node, "Decoder");

	std::string decoderTypeStr = TryGetAttribute(decoder_node, "type").value();
	decoderType decoderType_ = decoderTypeFromString(decoderTypeStr);
	if (decoderType_ == decoderType::UnknownDecoder)
		throw ConfigParseException("Parse config error: unsupported value of decoder type: " + decoderTypeStr);

	params.decoder = decoderTypeFromString(decoderTypeStr);

	std::unordered_map <std::string, std::string > decoderParams;
	for (pugi::xml_node decoderParam_node : decoder_node.children("DecoderParam"))
	{
		decoderParams[TryGetAttribute(decoderParam_node, "name").value()] = decoderParam_node.child_value();
	}
	params.decoderParams = decoderParams;

	pugi::xml_node snrRange_node = TryGetChildNode(simulation_node, "SnrRange");
	bool isSnrPresent = false;
	if (snrRange_node != NULL) {
		for (double i = std::stod(TryGetAttribute(snrRange_node, "start").value());
			i < std::stod(TryGetAttribute(snrRange_node, "stop").value());
			i += std::stod(TryGetAttribute(snrRange_node, "step").value()))
			params.snrArray.push_back(i);
		isSnrPresent = true;
	}

	pugi::xml_node snrArray_node = TryGetChildNode(simulation_node, "SnrArray");
	if (snrArray_node != NULL) {
		for (pugi::xml_node snr_node : snrArray_node.children("Snr"))
		{
			params.snrArray.push_back(std::stod(snr_node.child_value()));
		}
		isSnrPresent = true;
	}

	if (!isSnrPresent)
		throw ConfigParseException("Parse config error: either \"SnrArray\" and \"SnrRange\" sections are absent.");

	return params;
}

std::vector<SimulationParams> ReadConfig(std::string configFilename) {
	std::vector<SimulationParams> paramsArray;

	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(configFilename.c_str());

	if (!result)
	{
		throw ConfigParseException("Error during open config xml file:" + std::string(result.description()));
	}

	pugi::xml_node config_node = GetChildNode(doc, "Config");
	for (pugi::xml_node simulation_node : config_node.children("Simulation"))
	{
		paramsArray.push_back(ReadSimulationSection(simulation_node));
	}

	if (paramsArray.empty())
		throw ConfigParseException("Parse config error: any \"Simulation\" section is absent.");

	return paramsArray;
}