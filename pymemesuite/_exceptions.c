#include <stdarg.h>
#include <stdbool.h>
#include <Python.h>

#define ERROR_BUFFER_SIZE 1024

void die
  (char *format,
   ...)
{
  va_list argp;
  char    buffer[ERROR_BUFFER_SIZE];

  va_start(argp, format);
  vsnprintf(&buffer[0], ERROR_BUFFER_SIZE, format, argp);
  va_end(argp);

  PyErr_SetString(PyExc_RuntimeError, &buffer[0]);
  abort();
}

void myassert
  (bool die_on_error,
   bool test,
   char * const    format,
   ...)
{
  va_list  argp;
  char    buffer[ERROR_BUFFER_SIZE];

  if (!test) {

    va_start(argp, format);
    vsnprintf(&buffer[0], ERROR_BUFFER_SIZE, format, argp);
    va_end(argp);

    if (die_on_error) {
      PyErr_SetString(PyExc_AssertionError, &buffer[0]);
      abort();
    } else {
      PyErr_WarnEx(PyExc_RuntimeWarning, &buffer[0], 1);
    }
  }
}
