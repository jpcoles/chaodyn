#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ic.h"

struct iclist_s ics[] = 
{ {"random",  ic_random},
  {"rshell",  ic_uniform_random_shell},
  {"2shells", ic_two_shells},
  {"shell",   ic_one_shell},
  {"8",       ic_figure8},
  {"2",       ic_two_points},
  {NULL, NULL}
};

//==============================================================================
//                                  get_ics
//==============================================================================
struct iclist_s *iclist()
{
    return ics;
}

//==============================================================================
//                                 ic_random
//==============================================================================
void ic_random(env_t *env)
{
    size_t i;
    size_t N = env->N;
    for (i=0; i < N; i++)
    {
        env->p[i].x[0] = (pos_t)(env->radius * (2*drand48()-1));
        env->p[i].x[1] = 0; //(pos_t)(env->radius * (2*drand48()-1));
        env->p[i].x[2] = 0; //(pos_t)(env->radius * (2*drand48()-1));
        env->p[i].v[0] = 0;
        env->p[i].v[1] = 0;
        env->p[i].v[2] = 0;
    }
}

//==============================================================================
//                                   sphere
//==============================================================================
void sphere(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1, class_t class)
{
    size_t i,j;
    double x,y,z, t;
    int good;

    for (i=N0; i < N1; i++)
    {
        do
        {
            z = 2.0 * drand48() - 1.0;
            t = 2.0 * M_PI * drand48();
            x = sqrt(1-(z*z)) * cos(t);
            y = sqrt(1-(z*z)) * sin(t);

            x = x0 + x*R;
            y = y0 + y*R;
            z = z0 + z*R;

            good = 1;
            for (j=0; j < i; j++)
            {
                double r2 = pow(env->p[i].x[0] - x, 2)
                          + pow(env->p[i].x[1] - y, 2)
                          + pow(env->p[i].x[2] - z, 2);
                if (r2 < 10000) { good = 0; break; }
            }
        } while (!good);

        env->p[i].x[0] = x;
        env->p[i].x[1] = y;
        env->p[i].x[2] = z;
        env->p[i].v[0] = 0;
        env->p[i].v[1] = 0;
        env->p[i].v[2] = 0;
        env->p[i].class = class;

        DBG(2) eprintf("%ld] %e %e %e\n", i, x, y, z);
    }
}

//==============================================================================
//                                   shell
//==============================================================================
void shell(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1, class_t class)
{
    size_t i,j;
    double x,y,z, t,w;

    assert(N1 % 2 == 0);
    for (i=N0; i < N1; i+=2)
    {
redo:   {
            z = 2.0 * drand48() - 1.0;
            t = 2.0 * M_PI * drand48();
            w = sqrt(1 - z*z);
            x = w * cos(t);
            y = w * sin(t);

            x *= R;
            y *= R;
            z *= R;

            x += x0;
            y += y0;
            z += z0;

            for (j=0; j < i; j++)
            {
                double r = sqrt(pow(env->p[i].x[0] - x, 2)
                              + pow(env->p[i].x[1] - y, 2)
                              + pow(env->p[i].x[2] - z, 2));
                if (r < 4*env->eps) goto redo;
            }
        } 

        env->p[i].x[0] = x;
        env->p[i].x[1] = y;
        env->p[i].x[2] = z;
        env->p[i].v[0] = 0;
        env->p[i].v[1] = 0;
        env->p[i].v[2] = 0;
        env->p[i].class = class;

        env->p[i+1].x[0] = -x;
        env->p[i+1].x[1] = -y;
        env->p[i+1].x[2] = -z;
        env->p[i+1].v[0] = 0;
        env->p[i+1].v[1] = 0;
        env->p[i+1].v[2] = 0;
        env->p[i+1].class = class;

        DBG(2) eprintf("%ld] %e %e %e\n", i, x, y, z);
    }
}

//==============================================================================
//                          ic_uniform_random_shell
//==============================================================================
void ic_uniform_random_shell(env_t *env)
{
    env->N      = MAX(env->N, 100);
    env->M      = 1e22 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));
    env->opt.Nclasses = 1;

    shell(env, 0,0,0, env->radius, 0, env->N, 0);
}

//==============================================================================
//                               ic_two_shells
//==============================================================================
void ic_two_shells(env_t *env)
{
    env->N      = MAX(env->N, 100);
    env->M      = 1e22 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));
    env->opt.Nclasses = 2;

    shell(env, -env->radius*(1.-1./5),0,0, env->radius/5, 0,env->N/2, 0);
    shell(env,  env->radius*(1.-1./5),0,0, env->radius/5, env->N/2,env->N, 1);

    fflush(stdout);
}

//==============================================================================
//                                ic_one_shell
//==============================================================================
void ic_one_shell(env_t *env)
{
    assert(env->N != 0);

    double Myr  = env->units.Myr  = 1e6 * 365 * 24 * 60 * 60;   // [s]
    double kpc  = env->units.kpc  = 3.08568025e19;              // [m]
    double Msun = env->units.Msun = 1.98892e30;                 // [kg]

    double G = 6.67300e-11;                                     // [m^3 kg^-1 s^-2]
    double L = env->units.L = 1e-15 * kpc;
    double T = env->units.T = 1     * Myr;
    double M = env->units.M = 1e-8  * Msun;
               env->units.G = G * (pow(L,-3) * M * pow(T,2));

    env->N      = env->N;
    env->M      = 1e12*Msun / M / env->N;
    env->radius = 80*kpc / L;
    env->eps    = 10*kpc / L;

#if 1
    LOG("unit G = %e g\n", env->units.G);
    LOG("env->M = %e\n", env->M);
    LOG("1 M = %e\n", env->M / Msun);
    LOG("1 T = %e\n", T / Myr);
    LOG("1 L = %e\n", L / kpc);
    LOG("1 V = %e\n", (L/T) / (kpc/Myr));
#else
    LOG("1 M = %e Msun\n", M / Msun);
    LOG("1 T = %e Myr\n", T / Myr);
    LOG("1 L = %e kpc\n", L / kpc);
    LOG("1 V = %e kpc/Myr\n", (L/T) / (kpc/Myr));
    LOG("1 kpc = %e L\n", kpc / L);
    LOG("env->M      = "MASST" [%e Msun]\n", env->M, env->M/Msun);
    LOG("env->radius = "POST"\n", env->radius);
    LOG("env->eps    = "SOFTT"\n", env->eps);
#endif

    env->p = malloc(env->N * sizeof(*(env->p)));
    env->opt.Nclasses = 1;

    shell(env, 0,0,0, 40*kpc / L, 0,env->N, 0);

    fflush(stdout);
}

//==============================================================================
//                             ic_colliding_halos
//==============================================================================
void ic_colliding_halos(env_t *env)
{
    /* This assumes a halo has already been loaded */

    size_t i;
    env->p = realloc(env->p, 2 * env->N * sizeof(*(env->p)));
    memcpy(env->p + env->N, env->p, env->N * sizeof(*(env->p)));

    for (i=0; i < env->N; i++)
        env->p[i].x[0] -= 50*env->units.kpc / env->units.L;

    for (i=env->N; i < 2*env->N; i++)
        env->p[i].x[0] += 50*env->units.kpc / env->units.L;

    //env->radius *= 2;

    env->N *= 2;
}

//==============================================================================
//                                 ic_sphere
//==============================================================================
void ic_sphere(env_t *env)
{
    env->N      = env->N;
    env->M      = 1e20 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));
    env->opt.Nclasses = 3;

    sphere(env, 0,0,0, env->radius, 0, env->N, 0);
}

//==============================================================================
//                                 ic_figure8
//==============================================================================
void ic_two_points(env_t *env)
{
    double Myr  = env->units.Myr  = 1e6 * 365 * 24 * 60 * 60;   // [s]
    double kpc  = env->units.kpc  = 3.08568025e19;              // [m]
    double Msun = env->units.Msun = 1.98892e30;                 // [kg]

    double G = 6.67300e-11;                                     // [m^3 kg^-1 s^-2]
    double L = env->units.L = 1e-12 * kpc;
    double T = env->units.T = 1    * Myr;
    double M = env->units.M = 1e-8  * Msun;
               env->units.G = G * (pow(L,-3) * M * pow(T,2));

    env->N      = 3;
    env->M      = 1e12*Msun / M / env->N;
    env->radius = 400*kpc / L;
    env->eps    = 10*kpc / L;
    env->opt.Nclasses = 1;
    env->p            = malloc(env->N * sizeof(*(env->p)));

    eprintf("unit G = %e g\n", env->units.G);
    eprintf("env->M = %e\n", env->M);
    eprintf("1 M = %e\n", env->M / Msun);
    eprintf("1 T = %e\n", T / Myr);
    eprintf("1 L = %e\n", L / kpc);
    eprintf("1 V = %e\n", (L/T) / (kpc/Myr));

    env->p[0].x[0] =  10 * kpc / L;
    env->p[0].x[1] =  0.0;
    env->p[0].x[2] =  0.0;
    env->p[0].v[0] =  0.0;
    env->p[0].v[1] =  0.0;
    env->p[0].v[2] =  0.0;
    env->p[0].class = 0;

    env->p[1].x[0] = -10 * kpc / L;
    env->p[1].x[1] =  0.0;
    env->p[1].x[2] =  0.0;
    env->p[1].v[0] =  0.0;
    env->p[1].v[1] =  0.0;
    env->p[1].v[2] =  0.0;
    env->p[1].class = 0;

#if 1
    env->p[1].x[0] =  0.0 * kpc / L;
    //env->p[1].x[1] =  0.0;
    env->p[1].x[1] = -17.32 * kpc / L;
    env->p[1].x[2] =  0.0;
    env->p[1].v[0] =  0.0;
    env->p[1].v[1] =  0.0;
    env->p[1].v[2] =  0.0;
    env->p[1].class = 0;
#endif
}

//==============================================================================
//                                 ic_figure8
//==============================================================================
void ic_figure8(env_t *env)
{
    double Myr  = env->units.Myr  = 1e6 * 365 * 24 * 60 * 60;   // [s]
    double kpc  = env->units.kpc  = 3.08568025e16;              // [km]
    double Msun = env->units.Msun = 1.98892e33;                 // [g]

//  double G = env->units.G = 1;
//  double L = env->units.L = 1e-2 * kpc;
//  double T = env->units.T = 1e-6 * Myr;
//  double M = env->units.M = pow(L,3) / G / pow(T,2);

    double G = 1;
    double L = env->units.L = 1e-3 * kpc;
    double T = env->units.T = 1e-12 * Myr;
    double M = env->units.M = pow(L,3) / G / pow(T,2);

    env->N            = 3;
    env->M            = 1* M;
    env->radius       = 2* L;
    env->eps          = 1e-8* L;
    env->opt.Nclasses = 3;
    env->p            = malloc(env->N * sizeof(*(env->p)));

    double vfac = L/T;
    double rfac = L;

    VL(1) LOG("M %e\n", M);
    VL(1) LOG("L %e\n", L);
    VL(1) LOG("T %e\n", T);
    VL(1) LOG("V %e\n", L/T);
    VL(1) LOG("vfac %e\n", vfac);
    VL(1) LOG("rvac %e\n", rfac);

    env->p[0].x[0] =  1.07614373351092    * rfac;
    env->p[0].x[1] =  0.0                 * rfac;
    env->p[0].x[2] =  0.0                 * rfac;
    env->p[0].v[0] =  0.0                 * vfac;
    env->p[0].v[1] = -0.468266218409064   * vfac;
    env->p[0].v[2] =  0.0                 * vfac;
    env->p[0].class = 0;

    env->p[1].x[0] = -0.538071866755461   * rfac;
    env->p[1].x[1] =  0.343706827758244   * rfac;
    env->p[1].x[2] =  0.0                 * rfac;
    env->p[1].v[0] = -1.09960375207519    * vfac;
    env->p[1].v[1] =  0.23413310920453    * vfac;
    env->p[1].v[2] =  0.0                 * vfac;
    env->p[1].class = 1;

    env->p[2].x[0] = -0.538071866755461   * rfac;
    env->p[2].x[1] = -0.343706827758244   * rfac;
    env->p[2].x[2] =  0.0                 * rfac;
    env->p[2].v[0] =  1.09960375207519    * vfac;
    env->p[2].v[1] =  0.23413310920453    * vfac;
    env->p[2].v[2] =  0.0                 * vfac;
    env->p[2].class = 2;
}

#if 0
void ic_figure8(env_t *env)
{
    double yr   = env->units.Myr  = 365;                        // [s]
    double kpc  = env->units.kpc  = 3.08568025e16;              // [km]
    double Msun = env->units.Msun = 1.98892e33;                 // [g]

    double G = env->units.G = 1;
    double L = env->units.L = 1 * kpc;
    double T = env->units.T = 1 * yr;
    double M = env->units.M = pow(L,3) / G / pow(T,2);

    env->N            = 3;
    env->M            = M;
    env->radius       = 2 * L;
    env->eps          = 1e-3 * L;
    env->opt.Nclasses = 3;
    env->p            = malloc(env->N * sizeof(*(env->p)));

    double vfac = L/T;
    double rfac = L;

    VL(1) LOG("M %e\n", M);
    VL(1) LOG("L %e\n", L);
    VL(1) LOG("T %e\n", T);
    VL(1) LOG("V %e\n", L/T);
    VL(1) LOG("vfac %e\n", vfac);
    VL(1) LOG("rvac %e\n", rfac);

    env->p[0].x[0] =  0.97000436    * rfac;
    env->p[0].x[1] = -0.24308753    * rfac;
    env->p[0].x[2] =  0.0           * rfac;
    env->p[0].v[0] =  0.93240737/2  * vfac;
    env->p[0].v[1] =  0.86473146/2  * vfac;
    env->p[0].v[2] =  0.0           * vfac;
    env->p[0].class = 0;

    env->p[1].x[0] =  0.0           * rfac;
    env->p[1].x[1] =  0.0           * rfac;
    env->p[1].x[2] =  0.0           * rfac;
    env->p[1].v[0] = -0.93240737    * vfac;
    env->p[1].v[1] = -0.86473146    * vfac;
    env->p[1].v[2] =  0.0           * vfac;
    env->p[1].class = 1;

    env->p[2].x[0] = -0.97000436    * rfac;
    env->p[2].x[1] =  0.24308753    * rfac;
    env->p[2].x[2] =  0.0           * rfac;
    env->p[2].v[0] =  0.93240737/2  * vfac;
    env->p[2].v[1] =  0.86473146/2  * vfac;
    env->p[2].v[2] =  0.0           * vfac;
    env->p[2].class = 2;
}
#endif
