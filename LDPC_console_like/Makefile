all: bin simulate

bin:
	mkdir -p bin 

simulate: bin/simulate.o bin/Decoder.o bin/MathOperations.o
	g++ bin/simulate.o bin/Decoder.o bin/MathOperations.o -pthread -o simulate

bin/simulate.o: src/simulate.cpp include/Decoder.h
	g++ -c src/simulate.cpp -o bin/simulate.o

bin/Decoder.o: src/Decoder.cpp include/Decoder.h include/MathOperations.h include/Exceptions.h 
	g++ -c src/Decoder.cpp -o bin/Decoder.o

bin/MathOperations.o: src/MathOperations.cpp include/MathOperations.h include/Exceptions.h
	g++ -c src/MathOperations.cpp -o bin/MathOperations.o

clean:
	rm -rf bin simulate
