/* A global RNG for use throughout GCMC */

#ifndef __RNG_H__
#define __RNG_H__

#include <random>

// Get random bits from system RNG
size_t sysrandom(void*, size_t);

namespace gcmc {
// RNG type used throughout code
typedef std::mt19937 rng_type;

// Actual instantiated rng
extern rng_type rng;

// Uniform distribution [0, 1)
double rand_uniform();

// Initialize the global rng
void init_rng();
}

#endif  // __RNG_H__
