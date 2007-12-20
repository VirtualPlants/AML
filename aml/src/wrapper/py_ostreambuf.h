#ifndef PYTHON_PY_OSTREAMBUF_HPP
#define PYTHON_PY_OSTREAMBUF_HPP

#include <Python.h>
#include <string>
#include <streambuf>

namespace python {


class py_ostreambuf : public std::streambuf
{
 private:
    PyObject* file;
    bool have_flush;

 public:
    py_ostreambuf(PyObject* file_object)
      : file( file_object)
    {
        have_flush = PyObject_HasAttrString( file, "flush");
    }

 protected:
    virtual int_type overflow( int_type c)
    {
	PyObject* write_method = PyString_FromString( "write" );
	char cc = static_cast<char>(c);
	const char * str = std::string(&cc,1).c_str();
	PyObject* pystr = PyString_FromString(str); 
        PyObject_CallMethodObjArgs(file, write_method, pystr, NULL);
	return c;
    }

    virtual int xsputn(const char_type* s, std::streamsize n)
    {
        const char * str = std::string( s, n).c_str();
        PyObject* pystr = PyString_FromString(str); 
 	PyObject* write_method = PyString_FromString( "write" ); 
	PyObject_CallMethodObjArgs(file, write_method, pystr, NULL); 
        return n;
    }
    
    virtual int sync()
    {
        if (have_flush) {
	  
	  PyObject* write_method = PyString_FromString( "flush" );
	  PyObject_CallMethodObjArgs(file, write_method, NULL);
        }
        return 0;
    }
};

} // !namespace python

#endif // !defined
