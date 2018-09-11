all : main main64 main128 main_gfc main_gfc128

main : main.c gf2_32.c gf2_32.h gf2_16.h gf2_16.c
	gcc -O3 -o main main.c gf2_32.c gf2_16.c

main64 : main64.c gf2_64.c gf2_64.h gf2_16.h gf2_16.c
	gcc -O3 -o main64 main64.c gf2_64.c gf2_16.c

main128 : main128.c gf2_128.c gf2_128.h gf2_16.h gf2_16.c 
	gcc -O3 -o main128 main128.c gf2_128.c gf2_16.c

main_gfc : main_gfc.c
	gcc -O3 -o main_gfc main_gfc.c -I${HOME}/code/gf-complete/include \
           ${HOME}/code/gf-complete/src/.libs/libgf_complete.a

main_gfc128 : main_gfc128.c
	gcc -O3 -o main_gfc128 main_gfc128.c -I${HOME}/code/gf-complete/include \
           ${HOME}/code/gf-complete/src/.libs/libgf_complete.a

.phony : all
