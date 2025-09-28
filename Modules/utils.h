/**
 * @file utils.h
 * @author David Ramírez Betancourth
 * @brief Utility functions and macros
 */

#ifndef UTILS_H
#define UTILS_H

#include <complex.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "datatypes.h"

#define TO_MHZ(x) ((x) * 1000000.0)

signal_iq_t* load_cs8(const char* filename);

Paths_t get_paths(void);

// Nueva firma: ahora recibe el tiempo de captura en segundos
int instantaneous_capture(BackendParams_t* params, Paths_t* paths);


int hackrf_is_connected(void);

#endif // UTILS_H
