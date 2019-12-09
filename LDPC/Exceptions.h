#include <exception>
#include <string>

class IncorrectMatrixDimensionsException : public std::exception
{
public:
	IncorrectMatrixDimensionsException(const string err) : std::exception(err.c_str()) {} ;
};

class IncorrectCodewordException : public std::exception
{
public:
	IncorrectCodewordException(const string err) : std::exception(err.c_str()) {};
};