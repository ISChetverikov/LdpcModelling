#include <filesystem>
#include <unordered_map>
#include "../include/SimulationParameters.h"
#include "../include/ConfigReading.h"
#include "../lib/pugixml-1.10/src/pugixml.hpp"
#include "../lib/pugixml-1.10/src/pugiconfig.hpp"
#include "../include/Exceptions.h"
#include "../include/SimulationType.h"
#include "../include/DecoderType.h"


pugi::xml_node GetChildNode(pugi::xml_node node, std::string childName) {
	auto child = node.child(childName.c_str());
	if (child == NULL)
		throw new ConfigParseException("Parse config error: can not find \"" + childName + "\" section");

	return child;
}

pugi::xml_attribute GetAttribute(pugi::xml_node node, std::string name) {
	auto attr = node.attribute(name.c_str());
	if (attr == NULL)
		throw new ConfigParseException("Parse config error: can not find \"" + name + "\" atribute");

	return attr;
}


SimulationParams ReadConfig(std::string configFilename) {
	SimulationParams params;

	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(configFilename.c_str());

	if (!result)
	{
		throw new ConfigParseException("Error during open config xml file:" + std::string(result.description()));
	}

	pugi::xml_node config_node = GetChildNode(doc, "Config");
	params.H_MatrixFilename = GetChildNode(config_node, "MatrixFilename").child_value();
	params.resultsFilename = GetChildNode(config_node, "ResultsFilename").child_value();
	
	pugi::xml_node simulation_node = GetChildNode(config_node, "Simulation");

	std::string simulationTypeStr = GetAttribute(simulation_node, "type").value();
	simulationType simulationType_ = simulationTypeFromString(simulationTypeStr);
	if (simulationType_ == simulationType::UnknownSimulation)
		throw new ConfigParseException("Parse config error: unsupported value of simulation type: " + simulationTypeStr);

	params.type = simulationType_;

	std::unordered_map <std::string, std::string > simulationParams;
	for (pugi::xml_node simulationParam_node : simulation_node.children("SimulationParam"))
	{
		simulationParams[GetAttribute(simulationParam_node, "name").value()] = simulationParam_node.child_value();
	}
	params.simulationTypeParams = simulationParams;

	pugi::xml_node decoder_node = GetChildNode(config_node, "Decoder");

	std::string decoderTypeStr = GetAttribute(decoder_node, "type").value();
	decoderType decoderType_ = decoderTypeFromString(decoderTypeStr);
	if (decoderType_ == decoderType::UnknownDecoder)
		throw new ConfigParseException("Parse config error: unsupported value of decoder type: " + decoderTypeStr);

	params.decoder = decoderTypeFromString(decoderTypeStr);

	std::unordered_map <std::string, std::string > decoderParams;
	for (pugi::xml_node decoderParam_node : decoder_node.children("DecoderParam"))
	{
		decoderParams[GetAttribute(decoderParam_node, "name").value()] = decoderParam_node.child_value();
	}
	params.decoderParams = decoderParams;

	pugi::xml_node snrRange_node = GetChildNode(config_node, "SnrRange");
	for (double i = std::stod(GetAttribute(snrRange_node, "start").value());
		i < std::stod(GetAttribute(snrRange_node, "stop").value());
		i += std::stod(GetAttribute(snrRange_node, "step").value()))
		params.snrArray.push_back(i);

	pugi::xml_node snrArray_node = GetChildNode(config_node, "SnrArray");
	for (pugi::xml_node snr_node : snrArray_node.children("Snr"))
	{
		params.snrArray.push_back(std::stod(snr_node.child_value()));
	}
	
	return params;
}