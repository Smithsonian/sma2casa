makevis.so: makevis.c ./Makefile
	gcc -fPIC -shared -I/usr/include/python2.4/ -lpython2.4 \
	-o makevis.so makevis.c
