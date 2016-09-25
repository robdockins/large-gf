rungf : rungf.c ../gf2_16.c ../gf2_64.c
	cc rungf.c ../gf2_16.c ../gf2_64.c -o rungf -O3 \
           -lOpenCL -L/opt/cuda/lib64 \
           -I${HOME}/code/gf-complete/include \
           ${HOME}/code/gf-complete/src/.libs/libgf_complete.a

#hello : hello.c
#	cc hello.c -o hello -lOpenCL -L/opt/cuda/lib64

#hello: hello.cpp
#	g++ -std=c++0x -o hello hello.cpp -lOpenCL -L/opt/cuda/lib64
