#include <stdbool.h>

#include "macros.h"
#include "logs.h"
#include "utils.h"

VERBOSE_T verbosity = INVALID_VERBOSE;

double log_table[2*((const int) log_precision)+2];
double exp_table[((const int) BITS) * ((const int) exp_precision) + 2];

bool TRACE = false;
bool VERBOSE = false;
bool PRINTALL = false;
bool PRINT_FASTA = false;
bool PRINT_LL = false;
bool PRINT_W = false;
bool PRINT_Z = false;
bool PRINT_STARTS = false;
// bool NO_STATUS = false;
bool NO_STATUS = true;
bool TEXT_ONLY = false;
bool DOC = false;
char* OFFSET_FILE = NULL;
int TIMER = 0;

// bool TRACE = true;
// bool VERBOSE = true;
// bool PRINTALL = true;
// bool PRINT_FASTA = true;
// bool PRINT_LL = true;
// bool PRINT_W = true;
// bool PRINT_Z = true;
// bool PRINT_STARTS = true;
// bool NO_STATUS = true;
// bool TEXT_ONLY = true;
// bool DOC = true;
// char* OFFSET_FILE = NULL;
// int TIMER = 0;

char* __crash_x__ = NULL;
