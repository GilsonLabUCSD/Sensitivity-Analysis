objects =  sensitivity.o 

all:$(objects)
	nvcc -arch=sm_30 $(objects) -o sensitivity

%.o: %.cpp
	nvcc -x cu -arch=sm_30 -I. -dc $< -o $@

clean:
	rm -f *.o sensitivity
