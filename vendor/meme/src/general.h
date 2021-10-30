#ifndef GENERAL_H
#define GENERAL_H

/* 
	The macro DO_STANDARD_COMMAND_LINE(n, stmts) processes argv[1..argc-1]
	according to the (semicolon-separated) statements in stmts.
	In the following, X has the following meaning:
		X=0, ignore statement
		X=1, execute statement
		X=2, execute statement but omit from usage message
		     unless EXP_USAGE_MESSAGE has been used as an ACTION
	The possible statements are as follows.
	  SIMPLE_FLAG_OPTN(X,s,desc,var), which sets var to true on
	    encountering the flag -s (s is an unquoted string of
	    letters and digits).
	  FLAG_OPTN(X,s,desc,stmt), which executes stmt on encountering -s.
	    Desc is used as descriptive text in a "USAGE" message.
	  SIMPLE_CFLAG_OPTN(X,s,var) and CFLAG_OPTN(X,s,stmt), which are like
	    SIMPLE_FLAG_OPTN or FLAG_OPTN but require s to be a single
	    character and allow multiple flags to be concatenated 
	    (as in ls -AF).
	  DATA_OPTN(X,s,desc1,desc2,stmt), which assumes a switch of the form
	    -s data, and executes stmt with _OPTION_ set to the data string.  
	    desc1 and desc2 are used as descriptive text in a "USAGE" message.
	  NON_SWITCH(X,desc,stmt), which executes stmt for any argument
	    that does not begin with "-".  
	    _OPTION_ is set to the argument.
	    Here, desc is used as descriptive text in a "USAGE" message.  
	  ANY_OPTION(X,desc,stmt), which executes stmt for any argument.
	    _OPTION_ is set to the argument.
	    Here, desc is used as descriptive text in a "USAGE" message.  
	  USAGE(desc), which prints desc in case of a command line error,
	    and has on effect otherwise.

	If there are fewer than <n> arguments on the command line, an
	error occurs.

	If argv[0][0] == \0, the "USAGE: " header is not printed.

	In case an error is detected, an error message is constructed from
	the flag names and descriptive text, and the program exits
	abnormally after printing the message.  Within one of the statements,
	the programmer can cause such an error message and termination with
	the macro COMMAND_LINE_ERROR.

	A usage message can be printed (eg, for -h) with USAGE_MESSAGE.
*/

#define SIMPLE_FLAG_OPTN(X,string,desc,var) \
  FLAG_OPTN(X, string, desc, (var) = true)

#define FLAG_OPTN(X,string,desc,stmt) \
  if (X) { \
    if (__ACTION__ == 0) { \
      if (strcmp(_OPTION_, "-" #string) == 0) { \
	{stmt;} \
	continue; \
      } \
    } \
    else if ((__EXP__ || X != 2) && (__ACTION__ == 1 || __ACTION__ == 4)) { \
      fprintf(stderr, "\t[-" #string "]" #desc "\n"); \
    } \
  }

#define CFLAG_OPTN(X,char,stmt) \
  if (X) { \
    if (__ACTION__ == 0 && _OPTION_[1] == (#char)[0]) { \
      __ACTION__ = 2; \
      _OPTION_ += 1; \
    }  \
    if (__ACTION__ == 2 && _OPTION_[0] == (#char)[0]) { \
      {stmt;} \
      continue; \
    } else if (X != 2 && __ACTION__ == 1) { \
      fprintf(stderr,"\t-" #char); \
      __ACTION__ = 3; \
    } else if (X != 2 && __ACTION__ == 3) { \
      fprintf(stderr, #char "\n"); \
    } else if (__ACTION__ == 5) { \
      __ACTION__ = 4; \
    } \
  }
  
#define SIMPLE_CFLAG_OPTN(X,char,var) \
  CFLAG_OPTN(X, char, (var) = 1)

#define DATA_OPTN(X,string,desc1,desc2,stmt) \
  if (X) { \
    if (__ACTION__ == 0 \
	&& strcmp(_OPTION_, "-" #string) == 0) { \
      __i__ += 1; \
      if (__i__ >= argc) __ACTION__ = 1; \
      _OPTION_ = argv[__i__]; \
      if (__ACTION__ != 1) {stmt;} \
      continue; \
    } else if (X != 2 && (__ACTION__ == 1 || __ACTION__ == 4)) { \
      fprintf(stderr, "\t[-" #string " " #desc1 "]" #desc2 "\n"); \
    } \
  }

#define NON_SWITCH(X,desc,stmt) \
  if (X) { \
    if (__ACTION__ == 0 && _OPTION_[0] != '-') { \
      stmt; \
      continue; \
    } \
    else if (X != 2 && (__ACTION__ == 1 || __ACTION__ == 4)) \
      fprintf(stderr, "\t" #desc); \
  }

#define ANY_OPTION(X,desc,stmt) \
  if (X) { \
    if (__ACTION__ == 0) { \
      stmt; \
      continue; \
    } \
    else if (X != 2 && (__ACTION__ == 1 || __ACTION__ == 4)) \
      fprintf(stderr, " " #desc "\n"); \
  }

#define USAGE(desc) \
  if (__ACTION__ == 1 || __ACTION__ == 4) \
    fprintf(stderr, " " #desc "\n"); \

#define COMMAND_LINE_ERROR { __ACTION__ = 1; continue; }

#define USAGE_MESSAGE { __ACTION__ = 6; continue; }

#define EXP_USAGE_MESSAGE { __EXP__ = 1; __ACTION__ = 6; continue; }

#define DO_STANDARD_COMMAND_LINE(n, stmts) \
  { \
    int __i__; \
    int __ACTION__ = 0; \
    int __EXP__ = 0; \
    char *_OPTION_ = NULL; \
    __i__ = 0; \
    while (1) { \
      if (__ACTION__ == 0 && argc < n+1) COMMAND_LINE_ERROR; \
      if (__ACTION__ == 2) { \
        _OPTION_ += 1; \
        if (_OPTION_[0] == '\0') __ACTION__ = 0; \
      } \
      if (__ACTION__ == 0) { \
        __i__ += 1; \
        if (__i__ >= argc) break; \
        _OPTION_ = argv[__i__]; \
      } \
      if (__ACTION__ == 1) { \
        if (__i__) fprintf(stderr, "error at: %s\n", argv[__i__]); \
      } \
      if (__ACTION__ == 6) { __ACTION__ = 1; } \
      if (__ACTION__ == 1) { \
        if (argv[0][0] != '\0') fprintf(stderr, "\nUSAGE:\n\t%s", argv[0]); \
      } \
      stmts; \
      if (__ACTION__ == 1 || __ACTION__ == 4) { \
        (void) putc('\n', stderr); \
        exit(1); \
      } \
      if (__ACTION__ == 3) { \
        (void) putc(']', stderr); \
        __ACTION__ = 5; \
      } else { \
         __ACTION__ = 1; \
      } \
    } \
  }

#endif
