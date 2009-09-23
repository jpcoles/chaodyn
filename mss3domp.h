#ifndef MSS3DOMP_H
#define MSS3DOMP_H

#include <assert.h>
#include <stdio.h>
#include <time.h>

#include <inttypes.h>

//==============================================================================
#define WITH_SOFTENING 1
#define WITH_INTEGERS 1
#define WITH_LONGS 0
#define WITH_FLOATS 0
#if (WITH_INTEGERS+WITH_LONGS+WITH_FLOATS) > 1
#error "Only one data type allowed."
#endif 
//==============================================================================

//==============================================================================
// Debugging, verbosity, and logging macros
//==============================================================================
#define DBG_LEVEL 0
#define DBG(__lvl) if (DBG_LEVEL >= __lvl)
#define VL(__lvl) if (env->opt.verbosity >= __lvl)

#define LOG(...) do {  \
    if (env->opt.logfp)  { fprintf(env->opt.logfp,  __VA_ARGS__); fflush(env->opt.logfp); } \
    if (env->opt.stdout) { fprintf(env->opt.stdout, __VA_ARGS__); fflush(env->opt.stdout); } \
} while (0)

#define eprintf(...) do { fprintf (stderr, __VA_ARGS__); } while (0)


//==============================================================================
#define CL "\033[2K\r"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


//==============================================================================
//                             Physical Constants
//==============================================================================
#define C_G (1.0)                      /* Gravitationl Constant */
//#define G (4.3e-3)                  /* pc/Msun(km/s)^2  */

//==============================================================================
//                                   Types
//==============================================================================

typedef __uint128_t uint128_t; // __attribute__((mode(TI)));
//typedef unsigned int uint128_t __attribute__((mode(TI)));
typedef __int128_t int128_t; // __attribute__((mode(TI)));
//typedef long long int128_t;

#define L_TIMET "%"PRId128
#define L_DISTT "%"PRId128
#define L_POST  "%"PRId128
#define L_VELT  "%"PRId128
#define L_ACCT  "%"PRId128
#define L_RHOT  "%e"
#define L_MASST "%e"
#define L_ENGYT "%e"
#define L_FORCET "%e"
//#define L_FORCET "%"PRId64
#define L_SOFTT  "%e"

#define I_TIMET "%"PRId64
#define I_DISTT "%"PRId64
#define I_POST  "%"PRId64
#define I_VELT  "%"PRId64
#define I_ACCT  "%"PRId64
#define I_RHOT  "%e"
#define I_MASST "%e"
#define I_ENGYT "%.20e"
#define I_FORCET "%e"
//#define I_FORCET "%"PRId64
#define I_SOFTT  "%e"

#define F_TIMET "%e"
#define F_DISTT "%e"
#define F_POST  "%e"
#define F_VELT  "%e"
#define F_ACCT  "%e"
#define F_RHOT  "%e"
#define F_MASST "%e"
#define F_ENGYT "%e"
#define F_FORCET "%e"
#define F_SOFTT  "%e"

#if WITH_LONGS
#define PRId128 "x"

typedef int128_t tyme_t;
typedef int128_t pos_t;
typedef int128_t dist_t;
typedef int128_t vel_t;
typedef int128_t acc_t;
typedef double  force_t;
#define TIMET  L_TIMET 
#define POST   L_POST  
#define DISTT  L_DISTT  
#define VELT   L_VELT  
#define ACCT   L_ACCT  
#define RHOT   L_RHOT  
#define MASST  L_MASST 
#define ENGYT  L_ENGYT 
#define FORCET L_FORCET
#define SOFTT  L_SOFTT
#endif

#if WITH_INTEGERS
typedef int64_t tyme_t;
typedef int64_t pos_t;
typedef int64_t dist_t;
typedef int64_t vel_t;
typedef int64_t acc_t;
typedef double  force_t;
#define TIMET  I_TIMET 
#define POST   I_POST  
#define DISTT  I_DISTT  
#define VELT   I_VELT  
#define ACCT   I_ACCT  
#define RHOT   I_RHOT  
#define MASST  I_MASST 
#define ENGYT  I_ENGYT 
#define FORCET I_FORCET
#define SOFTT  I_SOFTT
#endif

#if WITH_FLOATS
typedef double pos_t;
typedef double dist_t;
typedef double vel_t;
typedef double acc_t;
typedef double tyme_t;
typedef double force_t;
#define TIMET  F_TIMET 
#define DISTT  F_DISTT 
#define POST   F_POST  
#define VELT   F_VELT  
#define ACCT   F_ACCT  
#define RHOT   F_RHOT  
#define MASST  F_MASST 
#define ENGYT  F_ENGYT 
#define FORCET F_FORCET
#define SOFTT  F_SOFTT
#endif

typedef double mass_t;             /* Msun             */
typedef double energy_t;           /* */
typedef double soft_t;
typedef unsigned char class_t;

typedef struct
{
    pos_t x[3];
    vel_t v[3];
    class_t class;
} particle_t;

typedef struct
{
    uint32_t nr;
    uint32_t nc;
    unsigned char *image;
} image_t;

typedef struct
{
    double  G;
    double  M;
    double  L;
    double T;

    double  Msun;
    double  kpc;
    double Myr;
} units_t;

typedef struct 
{
    char     save;
    char     save_image;
    char     save_path_image;
    char     dump_sim;

    size_t   save_every;
    size_t   reverse_at;
    size_t   save_image_every;
    size_t   save_path_every;
    size_t   start_step;
    uint32_t modify_mode;

    char    *tag;
    char    *inputfile;
    size_t   verbosity; 

    size_t   Nsteps;
    size_t   Nclasses;
    size_t   output_every;

    int32_t img_rows;
    int32_t img_cols;

    FILE *logfp;
    FILE *stdout;

    char *ic_gen_name;

    struct 
    {
        char dup_ic;
    } X;

    char restart;

} options_t;

typedef struct
{
    size_t   N;
    particle_t * restrict p0;
    particle_t * restrict p;
    force_t *    restrict F[3];

    force_t *** restrict Ft;

    soft_t eps;

    mass_t   M;
    tyme_t  dt;
    tyme_t  t;
    size_t step;
    size_t end_step;
    time_t seed;
    pos_t    radius;

    char **class_color;

    options_t opt;

    image_t snapshot;
    image_t path;

    units_t units;

} env_t;

typedef struct
{
    char version[16];
    uint32_t with_integers;
    size_t   N;
    mass_t   M;
    soft_t   eps;
    tyme_t  dt;
    size_t   Nclasses;
    size_t   step;
    size_t   dummy0;
    pos_t    radius;
    size_t   dummy1;
    time_t   seed;
    uint32_t img_rows;
    uint32_t img_cols;
    units_t  units;

} header_t;

#endif

