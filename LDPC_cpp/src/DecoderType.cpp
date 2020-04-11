#include "../include//DecoderType.h"

decoderType decoderTypeFromString(std::string str) {
	std::unordered_map<std::string, decoderType> decoderTypeResolver = {
		{"ONMS", decoderType::ONMS},
		{"MS", decoderType::MS},
		{"BF", decoderType::BF},
		{"SP", decoderType::SP}
	};

	if (decoderTypeResolver.count(str) > 0)
		return decoderTypeResolver[str];

	return UnknownDecoder;
}

std::string decoderTypeToString(decoderType type) {

	std::unordered_map<decoderType, std::string> decoderTypeStringResolver = {
		{decoderType::ONMS, "ONMS"},
		{decoderType::MS, "MS"},
		{decoderType::BF, "BF"},
		{decoderType::SP, "SP"}
	};

	return decoderTypeStringResolver[type];
}