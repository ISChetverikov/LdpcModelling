#pragma once

#include <exception>
#include <string>

class IncorrectMatrixDimensionsException : public std::exception
{
public:
	IncorrectMatrixDimensionsException(const std::string err) : std::exception(err.c_str()) {} ;
};

class IncorrectCodewordException : public std::exception
{
public:
	IncorrectCodewordException(const std::string err) : std::exception(err.c_str()) {};
};

class IncorrectDimensionsException : public std::exception
{
public:
	IncorrectDimensionsException(const std::string err) : std::exception(err.c_str()) {};
};