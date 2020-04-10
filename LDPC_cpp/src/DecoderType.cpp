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