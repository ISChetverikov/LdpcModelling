all: bin simulate

bin:
	mkdir -p bin 

simulate: bin/simulate.o bin/FFH_simulate.o bin/Base_decoder.o bin/ONMS_decoder.o bin/BF_decoder.o bin/SP_decoder.o bin/MathOperations.o bin/MatrixReading.o
	g++ bin/simulate.o bin/FFH_simulate.o bin/Base_decoder.o bin/ONMS_decoder.o bin/SP_decoder.o bin/BF_decoder.o bin/MathOperations.o bin/MatrixReading.o -pthread -o simulate

bin/simulate.o: src/simulate.cpp include/ONMS_decoder.h include/BF_decoder.h
	g++ -c src/simulate.cpp -o bin/simulate.o

bin/FFH_simulate.o: src/FFH_simulate.cpp include/FFH_simulate.h include/ONMS_decoder.h include/MatrixReading.h
	g++ -c src/FFH_simulate.cpp -o bin/FFH_simulate.o

bin/Base_decoder.o: src/Base_decoder.cpp 
	g++ -c src/Base_decoder.cpp -o bin/Base_decoder.o

bin/ONMS_decoder.o: src/ONMS_decoder.cpp include/ONMS_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -c src/ONMS_decoder.cpp -o bin/ONMS_decoder.o

bin/BF_decoder.o: src/BF_decoder.cpp include/BF_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -c src/BF_decoder.cpp -o bin/BF_decoder.o

bin/SP_decoder.o: src/SP_decoder.cpp include/SP_decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -c src/SP_decoder.cpp -o bin/SP_decoder.o

bin/MathOperations.o: src/MathOperations.cpp include/MathOperations.h include/Exceptions.h
	g++ -c src/MathOperations.cpp -o bin/MathOperations.o

bin/MatrixReading.o: src/MatrixReading.cpp include/MatrixReading.h
	g++ -c src/MatrixReading.cpp -o bin/MatrixReading.o

clean:
	rm -rf bin simulate
