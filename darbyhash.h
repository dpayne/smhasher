#pragma once
#include <x86intrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

uint64_t darbyhash(const void *data, size_t len, uint64_t seed);

#ifdef __cplusplus
}
#endif
