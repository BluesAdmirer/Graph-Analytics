# finalalgorithm : finalalgorithm.cpp
# 			g++ finalalgorithm.cpp -fopenmp -o finalalgorithm
# 			./finalalgorithm

pagerank : pagerank.cu
		nvcc pagerank.cu -o pagerank -arch=sm_60
		./pagerank

#test: test.cu
#	nvcc test.cu -o test -arch=sm_60
#	./test
# 	
