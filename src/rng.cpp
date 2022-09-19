/* A global RNG for use throughout GCMC */

#include <stdexcept>
#include <array>
#include <random>
#include <unistd.h>

#include "rng.h"


// Check if getrandom syscall is defined
#if defined(__linux__) || defined(linux) || defined(__linux)
#include <linux/version.h>
// 'getrandom' defined in version Kernel 3.17 and later
#if LINUX_VERSION_CODE >= KERNEL_VERSION(3,17,0)
#define HAVE_GETRANDOM
#endif // version >3.17
#endif // linux


// System specific random numbers
#ifdef HAVE_GETRANDOM
#include <sys/syscall.h>

// Get random bits with the getrandom syscall
size_t sysrandom(void *dest, size_t dest_len) {
    int bytes = syscall(SYS_getrandom, dest, dest_len, 0);
    if (bytes != dest_len) {
        throw std::runtime_error("Unable to read bits from system RNG");
    }
    return dest_len;
}

#else
#include <fcntl.h>

// Get random bits from /dev/urandom
size_t sysrandom(void *dest, size_t dest_len) {
    // Read unbuffered
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd == -1) {
        throw std::runtime_error("Could not open /dev/urandom");
    }
    if (read(fd, dest, dest_len) != dest_len) {
        close(fd);
        throw std::runtime_error("Could not read from system RNG");
    }
    close(fd);
    return dest_len;
}

#endif


namespace gcmc {
rng_type rng;


void init_rng() {
    std::array<int, rng_type::state_size> seed_data;
    sysrandom(seed_data.data(), seed_data.size() * sizeof(decltype(seed_data)::value_type));

    // Seed the actual rng
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq);

    // Seed the old C rng
    srand(rng());
}


double rand_uniform() {
    std::uniform_real_distribution<double> dist(0, 1);
    return dist(rng);
}

}
