all : main main64

main : main.c gf2_32.c gf2_32.h gf2_16.h gf2_16.c
	gcc -O3 -o main main.c gf2_32.c gf2_16.c

main64 : main64.c gf2_64.c gf2_64.h gf2_16.h gf2_16.c
	gcc -O3 -o main64 main64.c gf2_64.c gf2_16.c

.phony : all
