makevis.so: makevis.c ./Makefile
	gcc -O3 -Wall -fPIC -shared -I/usr/include/python2.4/ -lpython2.4 \
	-o makevis.so makevis.c
