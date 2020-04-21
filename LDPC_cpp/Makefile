all: bin main

bin:
	mkdir -p bin 

main: bin/main.o bin/BaseSimulator.o bin/simulate.o bin/MonteCarloSimulator.o bin/FastFlatHistSimulator.o bin/LocalFlatHistSimulator.o bin/Base_decoder.o bin/ONMS_decoder.o bin/BF_decoder.o bin/SP_decoder.o bin/MathOperations.o bin/MatrixReading.o bin/ConfigReading.o bin/SimulationType.o bin/DecoderType.o
	g++ bin/main.o bin/BaseSimulator.o bin/simulate.o bin/MonteCarloSimulator.o bin/FastFlatHistSimulator.o bin/LocalFlatHistSimulator.o bin/Base_decoder.o bin/ONMS_decoder.o bin/SP_decoder.o bin/BF_decoder.o bin/MathOperations.o bin/MatrixReading.o bin/ConfigReading.o bin/SimulationType.o bin/DecoderType.o -pthread -std=c++17 -ltbb -Wall -W -o main

bin/main.o: src/main.cpp include/simulate.h include/SimulationParameters.h include/MatrixReading.h
	g++ -std=c++17 -c src/main.cpp -Wall -W -o bin/main.o

bin/BaseSimulator.o: src/BaseSimulator.cpp include/BaseSimulator.h include/Base_decoder.h
	g++ -std=c++17 -c src/BaseSimulator.cpp -Wall -W -o bin/BaseSimulator.o

bin/simulate.o: src/simulate.cpp include/simulate.h include/SimulationParameters.h include/Base_decoder.h include/ONMS_decoder.h include/BF_decoder.h include/MonteCarloSimulator.h include/FastFlatHistSimulator.h include/LocalFlatHistSimulator.h include/BaseSimulator.h include/MatrixReading.h include/ConfigReading.h
	g++ -std=c++17 -c src/simulate.cpp -Wall -W -o bin/simulate.o

bin/MonteCarloSimulator.o: src/MonteCarloSimulator.cpp include/MonteCarloSimulator.h include/Base_decoder.h include/BaseSimulator.h
	g++ -std=c++17 -c src/MonteCarloSimulator.cpp -Wall -W -o bin/MonteCarloSimulator.o

bin/FastFlatHistSimulator.o: src/FastFlatHistSimulator.cpp include/FastFlatHistSimulator.h include/Base_decoder.h include/BaseSimulator.h
	g++ -std=c++17 -c src/FastFlatHistSimulator.cpp -Wall -W -o bin/FastFlatHistSimulator.o

bin/LocalFlatHistSimulator.o: src/LocalFlatHistSimulator.cpp include/LocalFlatHistSimulator.h include/Base_decoder.h include/BaseSimulator.h
	g++ -std=c++17 -c src/LocalFlatHistSimulator.cpp -Wall -W -o bin/LocalFlatHistSimulator.o

bin/Base_decoder.o: src/Base_decoder.cpp 
	g++ -std=c++17 -c src/Base_decoder.cpp -Wall -W -o bin/Base_decoder.o

bin/ONMS_decoder.o: src/ONMS_decoder.cpp include/ONMS_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -std=c++17 -c src/ONMS_decoder.cpp -Wall -W -o bin/ONMS_decoder.o

bin/BF_decoder.o: src/BF_decoder.cpp include/BF_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -std=c++17 -c src/BF_decoder.cpp -Wall -W -o bin/BF_decoder.o

bin/SP_decoder.o: src/SP_decoder.cpp include/SP_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -std=c++17 -c src/SP_decoder.cpp -Wall -W -o bin/SP_decoder.o

bin/MathOperations.o: src/MathOperations.cpp include/MathOperations.h include/Exceptions.h
	g++ -std=c++17 -c src/MathOperations.cpp -Wall -W -o bin/MathOperations.o

bin/MatrixReading.o: src/MatrixReading.cpp include/MatrixReading.h
	g++ -std=c++17 -c src/MatrixReading.cpp -Wall -W -o bin/MatrixReading.o

bin/SimulationType.o: src/SimulationType.cpp include/SimulationType.h
	g++ -std=c++17 -c src/SimulationType.cpp -Wall -W -o bin/SimulationType.o

bin/DecoderType.o: src/DecoderType.cpp include/DecoderType.h
	g++ -std=c++17 -c src/DecoderType.cpp -Wall -W -o bin/DecoderType.o

bin/ConfigReading.o: src/ConfigReading.cpp include/ConfigReading.h include/SimulationParameters.h include/Exceptions.h include/SimulationType.h include/DecoderType.h lib/pugixml-1.10/src/pugiconfig.hpp lib/pugixml-1.10/src/pugixml.hpp
	g++ -std=c++17 -c src/ConfigReading.cpp -Wall -W -o bin/ConfigReading.o

clean:
	rm -rf bin main