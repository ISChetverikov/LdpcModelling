all: bin main

bin:
	mkdir -p bin 

main: bin/main.o bin/BaseSimulator.o bin/simulate.o bin/FFH_simulate.o bin/MonteCarloSimulator.o bin/Base_decoder.o bin/ONMS_decoder.o bin/BF_decoder.o bin/SP_decoder.o bin/MathOperations.o bin/MatrixReading.o
	g++ bin/main.o bin/BaseSimulator.o bin/simulate.o bin/FFH_simulate.o bin/MonteCarloSimulator.o bin/Base_decoder.o bin/ONMS_decoder.o bin/SP_decoder.o bin/BF_decoder.o bin/MathOperations.o bin/MatrixReading.o -pthread -std=c++17 -ltbb -Wall -W -o main

bin/main.o: src/main.cpp include/simulate.h include/FFH_simulate.h include/SimulationParameters.h
	g++ -std=c++17 -c src/main.cpp -Wall -W -o bin/main.o

bin/BaseSimulator.o: src/BaseSimulator.cpp include/BaseSimulator.h include/Base_decoder.h
	g++ -std=c++17 -c src/BaseSimulator.cpp -Wall -W -o bin/BaseSimulator.o

bin/simulate.o: src/simulate.cpp include/simulate.h include/SimulationParameters.h include/Base_decoder.h include/ONMS_decoder.h include/BF_decoder.h include/MonteCarloSimulator.h include/BaseSimulator.h include/MatrixReading.h
	g++ -std=c++17 -c src/simulate.cpp -Wall -W -o bin/simulate.o

bin/FFH_simulate.o: src/FFH_simulate.cpp include/FFH_simulate.h include/ONMS_decoder.h include/MatrixReading.h
	g++ -std=c++17 -c src/FFH_simulate.cpp -Wall -W -o bin/FFH_simulate.o

bin/MonteCarloSimulator.o: src/MonteCarloSimulator.cpp include/MonteCarloSimulator.h include/Base_decoder.h include/BaseSimulator.h
	g++ -std=c++17 -c src/MonteCarloSimulator.cpp -Wall -W -o bin/MonteCarloSimulator.o

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

clean:
	rm -rf bin main