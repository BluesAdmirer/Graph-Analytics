main : main.cu
	nvcc main.cu -o main -arch=sm_61
