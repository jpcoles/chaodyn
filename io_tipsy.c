#include "chaodyn.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <arpa/inet.h>

struct tipsy_header
{
    double   time;
    uint32_t nBodies;
    uint32_t nDims;
    uint32_t nSph;
    uint32_t nDark;
    uint32_t nStar;
    uint32_t pad;
};

struct tipsy_dark_particle
{
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
};

float htonf(float f)
{
    union { float f; uint32_t i; } ret;
    ret.f = f;
    ret.i = htonl(ret.i);
    return ret.f;
}

void save_ic_tipsy(env_t *env)
{
    size_t i;
    char *fname = (char *)malloc(strlen(env->opt.tag)+1+2+1+5+1);
    sprintf(fname, "%s.ic.tipsy", env->opt.tag);

    //--------------------------------------------------------------------------
    // Setup our unit conversions. We can choose the mass and length unit but
    // PKDGRAV will assume that G=1, so our time unit is constrained.  Below,
    // the mass unit is chosen such that the total mass will sum to 1.
    //--------------------------------------------------------------------------


    double Mt = 1.0 / (env->M * env->N);
    double Lt = 1.0 / (env->units.L / env->radius);
    double Tt = sqrt(env->units.G * pow(Lt,-3) * Mt);

    eprintf("%ld %e %e\n", env->N, env->M, env->M * Mt);
    LOG("%e [M] = 1 [Mt] = %e Msun = %e\n", Mt, 1*env->units.Msun / env->units.M / Mt, Mt);
    LOG("%e [L] = 1 [Lt] = %e kpc  = %e\n", Lt, 1*env->units.kpc  / env->units.L / Lt, Lt);
    LOG("%e [T] = 1 [Tt] = %e Myr  = %e\n", Tt, 1*env->units.Myr  / env->units.T / Tt, Tt);
    

    //--------------------------------------------------------------------------
    // Write the header and all the particles. The particles are assumed to be
    // dark matter particles. Unit conversions are performed on the fly.
    //--------------------------------------------------------------------------


    LOG("Writing %s\n", fname);

    struct tipsy_header h;
    struct tipsy_dark_particle d;

    FILE *fp = fopen(fname, "wb");

    h.time     = htonf(0.0);
    h.nBodies  = htonl((uint32_t)env->N);
    h.nDims    = htonl(3);
    h.nDark    = htonl((uint32_t)env->N);
    h.nStar    = htonl(0);
    h.nSph     = htonl(0);

    fwrite(&h, sizeof(h), 1, fp);

    for (i=0; i < env->N; i++)
    {
        d.mass   = htonf(env->M   / Mt);
        d.eps    = htonf(env->eps / Lt);
        d.phi    = htonf(0.0);

        d.pos[0] = htonf(env->p[i].x[0] / Lt);
        d.pos[1] = htonf(env->p[i].x[1] / Lt);
        d.pos[2] = htonf(env->p[i].x[2] / Lt);

        d.vel[0] = htonf(env->p[i].v[0] / (Lt / Tt));
        d.vel[1] = htonf(env->p[i].v[1] / (Lt / Tt));
        d.vel[2] = htonf(env->p[i].v[2] / (Lt / Tt));
        
        fwrite(&d, sizeof(d), 1, fp);
    }

    fclose(fp);

    free(fname);
}

