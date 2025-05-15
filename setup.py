#!/usr/bin/env python
# coding: utf-8

import collections
import configparser
import glob
import io
import os
import platform
import re
import socket
import subprocess
import sys
from unittest import mock

import setuptools
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.build_clib import build_clib as _build_clib
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Extension, Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Utils ------------------------------------------------------------------

def _split_multiline(value):
    value = value.strip()
    sep = max('\n,;', key=value.count)
    return list(filter(None, map(lambda x: x.strip(), value.split(sep))))

def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

# --- `setup.py` commands ----------------------------------------------------

class build_ext(_build_ext):
    """A `build_ext` that disables optimizations if compiled in debug mode.
    """

    def finalize_options(self):
        _build_ext.finalize_options(self)
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.force = self.force
        self._clib_cmd.debug = self.debug

    def run(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError("Cython is required to run `build_ext` command") from cythonize

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {},
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
            }
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["gdb_debug"] = True
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
            cython_args["compiler_directives"]["profile"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False
            cython_args["compiler_directives"]["cdivision"] = True

        # cythonize and patch the extensions
        self.extensions = cythonize(self.extensions, **cython_args)
        for ext in self.extensions:
            ext._needs_stub = False

        # check the libraries have been built already
        if not self.distribution.have_run["build_clib"]:
            self._clib_cmd.run()

        # build the extensions as normal
        _build_ext.run(self)

    def build_extension(self, ext):
        # update compile flags if compiling in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-Og")
                ext.extra_compile_args.append("--coverage")
                ext.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Od")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # use GNU89 extensions if
        if platform.system() == "Linux" and self.compiler.compiler_type == "unix":
            ext.extra_compile_args.append("-fgnu89-inline")

        # update link and include directories
        ext.include_dirs.append(self._clib_cmd.build_clib)
        ext.library_dirs.append(self._clib_cmd.build_clib)

        # use static linking by directly using the static library path
        # instead of a `-l` flag, and remove the library so that the
        # compiler doesn't try to link against the system libraries
        for name in ext.libraries.copy():
            lib = self._clib_cmd.get_library(name)
            if lib is not None:
                ext.libraries.remove(name)
                ext.include_dirs.extend(lib.include_dirs)
                ext.extra_objects.append(self.compiler.library_filename(
                    lib.name, output_dir=self._clib_cmd.build_clib
                ))

        # build the rest of the extension as normal
        _build_ext.build_extension(self, ext)


class configure(_build_clib):
    """A ``./configure``-style command to generate C configuration header.

    Derives from the `build_clib` command from setuptools to be able to use
    the same configuration values (``self.build_temp``, ``self.build_clib``,
    etc).

    """

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Autotools-like helpers ---

    def _has_header(self, headername):
        _eprint('checking whether <{}> can be included'.format(headername), end="... ")

        slug = re.sub("[./-]", "_", headername)
        testfile = os.path.join(self.build_temp, "have_{}.c".format(slug))
        objects = []

        with open(testfile, "w") as f:
            f.write('#include "{}"\n'.format(headername))
        try:
            with mock.patch.object(sys, "stdout", new=io.StringIO()):
                objects = self.compiler.compile([testfile], debug=self.debug)
        except CompileError as err:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    _function_headers = {
        "floor": ["math.h"],
        "fork": ["unistd.h"],
        "getcwd": ["unistd.h"],
        "gethostbyname": ["netdb.h"],
        "gmtime": ["time.h"],
        "isascii": ["ctype.h"],
        "isinf": ["math.h"],
        "isnan": ["math.h"],
        "localtime": ["time.h"],
        "malloc": ["stdlib.h"],
        "memset": ["string.h"],
        "pow": ["math.h"],
        "realloc": ["stdlib.h"],
        "rint": ["math.h"],
        "socket": ["sys/socket.h"],
        "sqrt": ["math.h"],
        "stat": ["sys/stat.h"],
        "strchr": ["string.h"],
        "strcspn": ["string.h"],
        "strdup": ["string.h"],
        "strftime": ["time.h"],
        "strlcpy": ["string.h"],
        "strspn": ["string.h"],
        "strstr": ["string.h"],
        "vfork": ["unistd.h"],
    }

    _function_args = {
        "floor": ["3.14"],
        "getcwd": ["0", "0"],
        "gethostbyname": ['"localhost"'],
        "gmtime": ["0"],
        "isascii": ["0"],
        "isinf": ["0.0"],
        "isnan": ["0.0"],
        "localtime": ["0"],
        "malloc": ["0"],
        "memset": ["0", "0", "0"],
        "pow": ["1.0", "0.0"],
        "realloc": ["0", "0"],
        "rint": ["0.0"],
        "socket": ["0", "0", "0"],
        "sqrt": ["0"],
        "stat": ["0", "0"],
        "strchr": ["0", "0"],
        "strcspn": ["0", "0"],
        "strdup": ["0"],
        "strftime": ["0", "0", '"%a"', "0"],
        "strlcpy": ["0", "0", "0"],
        "strspn": ["0", "0"],
        "strstr": ["0", "0"],
    }

    def _has_function(self, funcname):
        _eprint('checking whether function', repr(funcname), 'is available', end="... ")

        testfile = os.path.join(self.build_temp, "have_{}.c".format(funcname))
        binfile = os.path.join(self.build_temp, "have_{}.bin".format(funcname))
        objects = []

        with open(testfile, "w") as f:
            for header in self._function_headers.get(funcname, []):
                f.write('#include <{}>\n'.format(header))
            args = ", ".join(self._function_args.get(funcname, []))
            f.write('int main() {{ {}({}); return 0;}}\n'.format(funcname, args))
        try:
            with mock.patch.object(sys, "stdout", new=io.StringIO()):
                objects = self.compiler.compile([testfile], debug=self.debug, extra_preargs=["-Wno-unused-value", "-Wno-nonnull"])
                self.compiler.link_executable(objects, binfile)
        except:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    # --- Required interface for `setuptools.Command` ---

    def build_libraries(self, libraries):
        # read `setup.cfg`, which should be next to this file, so we can
        # use the `__file__` magic variable to locate it
        _cfg = configparser.ConfigParser()
        _cfg.read([__file__.replace(".py", ".cfg")])

        # ensure the output directory exists, otherwise create it
        self.mkpath(self.build_clib)

        # run the `configure_library` method sequentially on each library,
        # unless the header already exists and the configuration has not
        # been edited
        for library in self.distribution.libraries:
            section = "configure.{}".format(library.name)

            try:
                config = dict(_cfg.items(section))
            except configparser.NoSectionError:
                continue

            headers = _split_multiline(config.get("headers", ""))
            functions = _split_multiline(config.get("functions", ""))
            constants = dict(
                map(str.strip, line.split("="))
                for line in _split_multiline(config.get("constants", ""))
            )

            if "filename" in config:
                self.mkpath(os.path.join(self.build_clib, library.name))
                self.make_file(
                    [__file__.replace(".py", ".cfg")],
                    os.path.join(self.build_clib, library.name, config["filename"]),
                    self.configure_library,
                    (
                        library,
                        config["filename"],
                        constants,
                        headers,
                        functions,
                    ),
                    exec_msg='generating "{}" for {} library'.format(
                        config["filename"], library.name
                    )
                )
                library.include_dirs.insert(0, os.path.join(self.build_clib, library.name))

    def configure_library(self, library, filename, constants=None, headers=None, functions=None):
        # create a mapping to store the defines, and make sure it is ordered
        # (this is to keep compatibility with Python 3.5, a dict would do fine)
        defines = collections.OrderedDict()

        # fill the defines with the constants
        constants = constants or {}
        for k, v in constants.items():
            defines[k] = v

        # check endianness
        if sys.byteorder == "big":
            defines["WORDS_BIGENDIAN"] = 1

        # fill the defines if headers are found
        headers = headers or []
        for header in headers:
            if self._has_header(header):
                slug = re.sub("[./-]", "_", header).upper()
                defines["HAVE_{}".format(slug)] = 1

        # fill the defines if functions are found
        functions = functions or []
        for func in functions:
            if self._has_function(func):
                defines["HAVE_{}".format(func.upper())] = 1

        # write the header file
        slug = re.sub("[./-]", "_", filename).upper()
        with open(os.path.join(self.build_clib, library.name, filename), "w") as f:
            f.write("#ifndef {}_INCLUDED\n".format(slug))
            f.write("#define {}_INCLUDED\n".format(slug))
            for k, v in defines.items():
                f.write("#define {} {}\n".format(k, '' if v is None else v))
            f.write("#endif\n".format(slug))


class build_clib(_build_clib):
    """A custom `build_clib` that compiles out of source.
    """

    # --- Distutils command interface ---

    def finalize_options(self):
        _build_clib.finalize_options(self)
        self._configure_cmd = self.get_finalized_command("configure")
        self._configure_cmd.force = self.force

    # --- Compatibility with base `build_clib` command ---

    def check_library_list(self, libraries):
        pass

    def get_library_names(self):
        return [ lib.name for lib in self.libraries ]

    def get_source_files(self):
        return [ source for lib in self.libraries for source in lib.sources ]

    # --- Helpers ---

    def get_library(self, name):
        return next((lib for lib in self.libraries if lib.name == name), None)

    # --- Build code ---

    def run(self):
        # make sure the C headers were generated already
        if not self.distribution.have_run["configure"]:
            self._configure_cmd.run()
        # build the libraries normally
        _build_clib.run(self)

    def build_libraries(self, libraries):
        # check platform
        if platform.system() == "Linux":
            self.compiler.define_macro("Linux")
        # get hostname
        self.compiler.define_macro("HOSTNAME", '"{}"'.format(socket.gethostname()))
        # compile each library
        self.mkpath(self.build_clib)
        for library in libraries:
            if platform.system() == "Linux" and self.compiler.compiler_type == "unix":
                library.extra_compile_args.append("-fgnu89-inline")
            self.make_file(
                library.sources,
                self.compiler.library_filename(library.name, output_dir=self.build_clib),
                self.build_library,
                (library,)
            )

    def build_library(self, library):
        # update compile flags if compiling in debug or release mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                library.extra_compile_args.append("-Og")
                library.extra_compile_args.append("--coverage")
                library.extra_link_args.append("--coverage")
            elif self.compiler.compiler_type == "msvc":
                library.extra_compile_args.append("/Od")

        # record extra include dirs: the compile needs to use the `config.h`
        # file that we generated in the `self.build_clib` folder, but some
        # C sources try to include it as `../config.h` so we need to add a
        # dummy include folder as well
        include_dirs = library.include_dirs.copy()
        include_dirs.append(os.path.join(self.build_clib, library.name))
        include_dirs.append(os.path.join(self.build_clib, library.name, "vendor"))

        # build objects and create a static library
        objects = self.compiler.compile(
            library.sources,
            output_dir=os.path.join(self.build_temp, library.name),
            include_dirs=include_dirs,
            macros=library.define_macros,
            debug=self.debug,
            depends=library.depends,
            extra_preargs=library.extra_compile_args,
        )
        self.compiler.create_static_lib(
            objects,
            library.name,
            output_dir=self.build_clib,
            debug=self.debug,
        )


class clean(_clean):

    def remove_file(self, filename):
        if os.path.exists(filename):
            _eprint("removing", repr(filename))
            os.remove(filename)
        else:
            _eprint(repr(filename), "does not exist -- can't clean it")

    def run(self):
        _clean.run(self)

        _build_cmd = self.get_finalized_command("build_ext")
        _build_cmd.inplace = True

        for ext in self.distribution.ext_modules:
            filename = _build_cmd.get_ext_filename(ext.name)
            if self.all:
                self.remove_file(filename)
            basename = _build_cmd.get_ext_fullname(ext.name).replace(".", os.path.sep)
            for ext in ["c", "html"]:
                filename = os.path.extsep.join([basename, ext])
                self.remove_file(filename)


# --- C static libraries -----------------------------------------------------

# fmt: off
libraries = [
    Library(
        "xml2",
        include_dirs=[
            os.path.join("vendor", "meme", "src"),
            os.path.join("vendor", "meme", "src", "libxml2", "include"),
        ],
        sources=[
            os.path.join("vendor", "meme", "src", "libxml2", "{}.c".format(src))
            for src in (
                "chvalid", "debugXML", "dict", "encoding", "entities",
                "error", "globals", "hash", "HTMLparser", "HTMLtree", "list",
                "parserInternals", "parser", "pattern", "relaxng", "tree",
                "SAX2", "threads", "uri", "valid", "xmlIO", "xmlmodule",
                "xmlreader", "xmlregexp", "xmlsave", "xmlstring", "xmlschemas",
                "xmlschemastypes", "xmlunicode", "xmlwriter", "xmlmemory",
                "xpath",
            )
        ],
    ),
    Library(
        "xslt",
        libraries=["xml2"],
        include_dirs=[
            os.path.join("vendor", "meme", "src"),
            os.path.join("vendor", "meme", "src", "libxml2", "include"),
        ],
        sources=[
            os.path.join("vendor", "meme", "src", "libxslt", "{}.c".format(src))
            for src in (
            	"attributes", "attrvt", "documents", "extensions",
                "extra",
            	"functions", "imports", "keys", "namespaces", "numbers",
                "pattern", "preproc", "security", "templates", "transform",
                "variables", "xslt", "xsltlocale", "xsltutils",
            )
        ],
    ),
    Library(
        "meme",
        libraries=["xml2", "xslt"],
        define_macros=[
            # ("MT_GENERATE_CODE_IN_HEADER", "0"),
        ],
        include_dirs=[
            "include",
            os.path.join("vendor", "meme", "src"),
            os.path.join("vendor", "meme", "src", "libxml2", "include"),
        ],
        sources=[
            os.path.join("vendor", "meme", "src", basename)
            for basename in [
                "alphabet.c",
                "alph-in.c",
                "array.c",
                "array-list.c",
                "binary-search.c",
                "cisml.c",
                "cisml-sax.c",
                "dreme-sax.c",
                "hash_table.c",
                "heap.c",
                "html-data.c",
                "io.c",
                "fasta-get-markov.c",
                "json-checker.c",
                "json-reader.c",
                "json-writer.c",
                "linked-list.c",
                "matrix.c",
                "meme-sax.c",
                "motif.c",
                "motif-db.c",
                "motif-in.c",
                "motif-in-common.c",
                "motif-in-dreme-xml.c",
                "motif-in-meme-html.c",
                "motif-in-meme-json.c",
                "motif-in-meme-text.c",
                "motif-in-meme-xml.c",
                "motif-in-streme-xml.c",
                "mtwist.c",
                "parser-message.c",
                "pssm.c",
                "prior-dist.c",
                "qvalue.c",
                "red-black-tree.c",
                "regex-utils.c",
                "reservoir.c",
                "sax-parser-utils.c",
                "scanned-sequence.c",
                "seq.c",
                "streme-sax.c",
                "string-builder.c",
                "string-list.c",
                "string-match.c",
                "ushuffle.c",
                "utils.c",
                "xlate-in.c",
                "xml-out.c",
                "xml-util.c",
            ]
        ],
    )
]


# --- Cython extensions ------------------------------------------------------

# fmt: off
extensions = [
    Extension(
        "pymemesuite.errors",
        sources=[
            os.path.join("pymemesuite", "errors.pyx"),
            os.path.join("pymemesuite", "_exceptions.c"),
        ],
        include_dirs=[os.path.join("vendor", "meme", "src")],
    ),
    Extension(
        "pymemesuite.common",
        sources=[
            os.path.join("pymemesuite", "common.pyx"),
            os.path.join("pymemesuite", "_globals.c"),
            os.path.join("pymemesuite", "_exceptions.c"),
        ],
        libraries=["m", "xml2", "meme"],
        include_dirs=[
            os.path.join("vendor", "meme", "src"),
        ],
    ),
    Extension(
        "pymemesuite.cisml",
        sources=[
            os.path.join("pymemesuite", "cisml.pyx"),
            os.path.join("pymemesuite", "_globals.c"),
            os.path.join("pymemesuite", "_exceptions.c"),
        ],
        libraries=["m", "xml2", "xslt", "meme"],
        include_dirs=[
            os.path.join("vendor", "meme", "src"),
            os.path.join("vendor", "meme", "src", "libxml2", "include")
        ],
    ),
    Extension(
        "pymemesuite.fimo",
        sources=[
            os.path.join("pymemesuite", "fimo.pyx"),
            os.path.join("pymemesuite", "_globals.c"),
            os.path.join("pymemesuite", "_exceptions.c"),
        ],
        libraries=["m", "xml2", "xslt", "meme"],
        include_dirs=[
            os.path.join("vendor", "meme", "src"),
            os.path.join("vendor", "meme", "src", "libxml2", "include")
        ],
    ),
]


# --- Setup ------------------------------------------------------------------

setuptools.setup(
    ext_modules=extensions,
    libraries=libraries,
    cmdclass=dict(
        build_ext=build_ext,
        build_clib=build_clib,
        clean=clean,
        configure=configure,
    ),
)
