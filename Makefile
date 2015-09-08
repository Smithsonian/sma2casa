makevis.so: makevis.c ./Makefile
	gcc -O3 -Wall -fPIC -shared -I/sma/python/anaconda/pkgs/python-2.7.6-1/include/python2.7/ /sma/python/anaconda/pkgs/python-2.7.6-1/lib/libpython2.7.so \
	-o makevis.so makevis.c
