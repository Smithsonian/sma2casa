#include "/sma/SMAusers/taco/anaconda/pkgs/python-2.7.6-1/include/python2.7/Python.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#define TRUE (1)
#define FALSE (0)

int fD;
const unsigned char *mapping;

static PyObject *makevis_open(PyObject *self, PyObject *args)
{
  int i;
  long long sts;
  struct stat s;
  size_t size;
  const char *fileName;

  if (!PyArg_ParseTuple(args, "s", &fileName))
    return NULL;
  fD = open(fileName, O_RDONLY);
  if (fD < 0) {
    fprintf(stderr, "makevis.open: error returned by open()\n");
    perror(fileName);
    exit(-1);
  }
  i = fstat(fD, &s);
  if (i != 0) {
    fprintf(stderr, "makevis.open: error returned by fstat()\n");
    perror(fileName);
    exit(-1);
  }
  size = s.st_size;
  mapping = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fD, 0);
  if (mapping == MAP_FAILED) {
    fprintf(stderr, "makevis.open: error returned by mmap()\n");
    perror(fileName);
    exit(-1);
  }
  sts = size;
  return Py_BuildValue("l", sts);
}

static PyObject *makevis_scanno(PyObject *self, PyObject *args)
{
  long int offset;
  int newFormat, scanNo;

  if (!PyArg_ParseTuple(args, "li", &offset, &newFormat))
    return NULL;
  if (newFormat)
    scanNo = *((unsigned int *)&mapping[offset]);
  else
    scanNo = (mapping[offset]<<24) + (mapping[offset+1]<<16) + (mapping[offset+2]<<8) + mapping[offset+3];
  return Py_BuildValue("i", scanNo);
}

static PyObject *makevis_recsize(PyObject *self, PyObject *args)
{
  long int offset;
  int newFormat, recSize;

  if (!PyArg_ParseTuple(args, "li", &offset, &newFormat))
    return NULL;
  if (newFormat)
    recSize = *((unsigned int *)&mapping[offset+4]);
  else
    recSize = (mapping[offset+8]<<24) + (mapping[offset+9]<<16) + (mapping[offset+10]<<8) + mapping[offset+11];
  return Py_BuildValue("i", recSize);
}

static PyObject *makevis_scaleexp(PyObject *self, PyObject *args)
{
  int scaleExp;
  long int offset;
  int newFormat;

  if (!PyArg_ParseTuple(args, "li", &offset, &newFormat))
    return NULL;
  if (newFormat)
    scaleExp = *((unsigned short *)&mapping[offset]);
  else
    scaleExp = (mapping[offset+8]<<8) + mapping[offset+9];
  if (scaleExp > 32767)
    scaleExp -= 65536;
  return Py_BuildValue("i", scaleExp);
}

static PyObject *makevis_convert(PyObject *self, PyObject *args)
{
  static int firstCall = TRUE;
  long int offset, offsetInc;
  int nPoints, newFormat, i, trim, first, last, reverse;
  int real, imag;
  double scale, weight, fReal, fImag;
  static PyObject *list;
  PyObject *num, *pyWeight;

  if (firstCall)
    firstCall = FALSE;
  else
    Py_DECREF(list);

  if (!PyArg_ParseTuple(args, "idlidiiii",
			&nPoints, &scale, &offset, &newFormat, &weight, &trim, &first, &last, &reverse))
    return NULL;
  if (newFormat)
    offset += 2;
  else
    offset += 10;
  list = PyList_New(3*nPoints);
  if (!list)
    return NULL;
  pyWeight = PyFloat_FromDouble(weight);
  if (!pyWeight) {
    Py_DECREF(list);
    return NULL;
  }
  if (reverse) {
    offset += (nPoints-1)*4;
    offsetInc = -4;
  } else {
    offsetInc = 4;
  }
  for (i = 0; i < nPoints; i++) {
    if ((weight <= 0.0) || (trim && ((i < first) || (i > last)))) {
      fReal = fImag = 0.0;
    } else {
      if (newFormat) {
	real = *((unsigned short *)&mapping[offset]);
	imag = *((unsigned short *)&mapping[offset+2]);
      } else {
	real = (mapping[offset]<<8) + mapping[offset+1];
	imag = (mapping[offset+2]<<8) + mapping[offset+3];
      }
      if (real > 32767)
	real -= 65536;
      if (imag > 32767)
	imag -= 65536;
      fReal = ((double)real)*scale;
      fImag = ((double)-imag)*scale;
    }
    num = PyFloat_FromDouble(fReal);
    if (!num) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, 3*i, num);
    num = PyFloat_FromDouble(fImag);
    if (!num) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, 3*i + 1, num);
    PyList_SET_ITEM(list, 3*i + 2, pyWeight);
    offset += offsetInc;
  }
  return Py_BuildValue("O", list);
}

static PyMethodDef MakevisMethods[] = {
  {"scaleexp",  makevis_scaleexp, METH_VARARGS,
   "Return the scaling exponent"},
  {"recsize",  makevis_recsize, METH_VARARGS,
   "Return the record size"},
  {"scanno",  makevis_scanno, METH_VARARGS,
   "Return the scan number"},
  {"open",  makevis_open, METH_VARARGS,
   "Open and mmap the sch_read file"},
  {"convert",  makevis_convert, METH_VARARGS,
   "Convert raw data to signed integers."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initmakevis(void)
{
  PyObject *m;
  
  m = Py_InitModule("makevis", MakevisMethods);
  if (m == NULL)
    return;
}
