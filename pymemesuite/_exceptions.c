#include <stdarg.h>
#include <Python.h>

void die
  (char *format,
   ...)
{
  va_list argp;
  char    buffer[1024];

  va_start(argp, format);
  vsnprintf(&buffer[0], 1024, format, argp);
  va_end(argp);

  PyErr_SetString(PyExc_RuntimeError, &buffer[0]);
  abort();
}
