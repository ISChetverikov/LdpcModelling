#pragma once

#include <string>
#include <unordered_map>

enum decoderType { ONMS, MS, BF, SP, UnknownDecoder };

decoderType decoderTypeFromString(std::string str);
std::string decoderTypeToString(decoderType type);