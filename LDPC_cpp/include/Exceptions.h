#pragma once

#include <exception>
#include <string>

class IncorrectMatrixDimensionsException : public std::exception
{
private:
	std::string m_error;
public:
	IncorrectMatrixDimensionsException(const std::string err) : m_error(err) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class IncorrectCodewordException : public std::exception
{
private:
	std::string m_error;
public:
	IncorrectCodewordException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class IncorrectDimensionsException : public std::exception
{
private:
	std::string m_error;
public:
	IncorrectDimensionsException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class NotBinaryMatrixException : public std::exception
{
private:
	std::string m_error;
public:
	NotBinaryMatrixException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class MatrixRowSkippedException : public std::exception
{
private:
	std::string m_error;
public:
	MatrixRowSkippedException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class ExtensionException : public std::exception
{
private:
	std::string m_error;
public:
	ExtensionException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};

class ConfigParseException : public std::exception
{
private:
	std::string m_error;
public:
	ConfigParseException(const std::string err) : m_error(err.c_str()) {};
	const char* what() const noexcept { return m_error.c_str(); }
};
