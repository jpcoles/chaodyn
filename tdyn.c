#include <math.h>
#include <stdlib.h>
#include "chaodyn.h"
#include "io.h"

#define PI (3.1415926535)

typedef struct
{
    mass_t sumr[3];
    mass_t N;
    pos_t x[3];
    pos_t R;
} cm_t;

typedef struct
{
    pos_t x[3];
    pos_t R;
    vel_t V;
    size_t id;
} tdyn_particle_t;

void usage()
{
    eprintf("usage: tdyn <sim-file>\n");
    exit(2);
}

void cm(size_t N, tdyn_particle_t *p, pos_t x[3], pos_t R, cm_t *CM)
{
    size_t i;

    CM->sumr[0] = 0;
    CM->sumr[1] = 0;
    CM->sumr[2] = 0;
    CM->N       = 0;
    pos_t maxR  = 0;
    for (i=0; i < N; i++)
    {
        double r = sqrt(pow(p[i].x[0] - x[0],2)
                      + pow(p[i].x[1] - x[1],2)
                      + pow(p[i].x[2] - x[2],2));

        if (r < R)
        {
            CM->sumr[0] += p[i].x[0];
            CM->sumr[1] += p[i].x[1];
            CM->sumr[2] += p[i].x[2];
            CM->N++;
            if (r > maxR) maxR = r;
        }
    }
    CM->x[0] = CM->sumr[0] / CM->N;
    CM->x[1] = CM->sumr[1] / CM->N;
    CM->x[2] = CM->sumr[2] / CM->N;
    CM->R    = maxR;
}

int main(int argc, char **argv)
{
    size_t i;
    size_t N;

    cm_t CM;

    env_t *env = calloc(1, sizeof(env_t)); assert(env != NULL);

    if (argc != 2) usage();

    load(env, argv[1]);
    N = env->N;
    tdyn_particle_t *p = malloc(env->N * sizeof(tdyn_particle_t)); assert(p != NULL);

#define MAG2(x) pow((x)[0],2) + pow((x)[1],2) + pow((x)[2],2)
#define MAG(x) sqrt(MAG2((x)))

    //--------------------------------------------------------------------------
    // Largest simulation radius.
    //--------------------------------------------------------------------------
    double R2=0;
    double v;
    vel_t V=0;
    for (i=0; i < env->N; i++)
    {
        p[i].id   = i;
        p[i].x[0] = env->p[i].x[0];
        p[i].x[1] = env->p[i].x[1];
        p[i].x[2] = env->p[i].x[2];
        double r2 = MAG2(p[i].x);
        p[i].R    = sqrt(r2);

        p[i].V = v = MAG(env->p[i].v);

        if (r2 > R2) R2 = r2;
        if (v  > V)  V = v;
    }

    tyme_t t = 100;
    eprintf("Maximum speed is %ld\n", V); 
    eprintf("For t="TIMET" M should be %e\n", t, t*pow(V,3)/(2*PI));
    dist_t R = sqrt(R2);
    eprintf("R is "DISTT"\n", R);
    eprintf("R is %e kpc\n", R * env->units.L / env->units.kpc);
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Initial CM guess.
    //--------------------------------------------------------------------------
    pos_t center[3] = {0,0,0};
    cm(N,p, center, R, &CM);
    //--------------------------------------------------------------------------

    char converged = 0;
    double dR = CM.R * 0.001;
    double tol0 = CM.R * 0.000001;
    double tol1 = CM.R * 0.000005;
    while (!converged)
    {
        cm_t CMt = CM;

        CMt.R -= dR;

        cm(N,p, CM.x, CM.R, &CMt);

        double
        shift = sqrt( pow(CM.x[0]-CMt.x[0],2)
                    + pow(CM.x[1]-CMt.x[1],2)
                    + pow(CM.x[2]-CMt.x[2],2));

        if (shift < tol0) 
            converged = 1;
        else if (shift > tol1)
            dR = CMt.R * 0.01;

        CM = CMt;
    }

    eprintf("Center of mass: ");
    eprintf(POST" "POST" "POST"", CM.x[0], CM.x[1], CM.x[2]);
    eprintf(" [%e %e %e]\n", ((double)CM.x[0])/R, ((double)CM.x[1])/R, ((double)CM.x[2])/R);

    //--------------------------------------------------------------------------
    // Recenter.
    //--------------------------------------------------------------------------
    for (i=0; i < env->N; i++)
    {
        p[i].x[0] -= CM.x[0];
        p[i].x[1] -= CM.x[1];
        p[i].x[2] -= CM.x[2];
        p[i].R = MAG(p[i].x);
    }
    //--------------------------------------------------------------------------

    eprintf("Particle 0 is at "POST" " POST" "POST"\n", p[0].x[0], p[0].x[1], p[0].x[2]);
    eprintf("Particle 0 is at R=%e %e %e %e\n", 
        ((double)p[0].R)/R,
        ((double)p[0].x[0])/R, ((double)p[0].x[1])/R, ((double)p[0].x[2])/R);

    //--------------------------------------------------------------------------
    // Sort.
    //--------------------------------------------------------------------------
    int compar(const void *a0, const void *b0)
    {
        const tdyn_particle_t *a = (tdyn_particle_t *)a0;
        const tdyn_particle_t *b = (tdyn_particle_t *)b0;
        
        const double aR = a->R;
        const double bR = b->R;

        return (aR == bR) + (2*(aR > bR) - 1);
    }
    qsort(p, env->N, sizeof(*p), compar);
    //--------------------------------------------------------------------------

    FILE *fp = fopen("tdyn.MvR", "w"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        fprintf(fp, "%e %e\n", R, ((double)(i+1)) / env->N);
    }
    fclose(fp);

#define RHO(N,R)      ((N)*env->M / (4./3.*PI*pow((R),3)))
#define RHOa(N,R0,R1) ((N)*env->M / (4./3.*PI*(pow((R1),3) - pow((R0),3))))

    fp = fopen("tdyn.dencum", "w"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        fprintf(fp, "%e %e\n", R, RHO(i+1,R));
    }
    fclose(fp);

    fp = fopen("tdyn.den", "w"); assert(fp != NULL);
    size_t j=0;
    double R0 = 0;
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        if (R0 != 0)
            fprintf(fp, "%e %e\n", R, RHOa(i-j+1,R0,R));
        R0 = R;
        j = i;
    }
    fclose(fp);

    double frac = 0.5;
    size_t Mfrac_index = env->N * frac - 1;
    if (!(0 <= Mfrac_index && Mfrac_index <= env->N)) Mfrac_index = 0;
    double R_Mfrac  = p[Mfrac_index].R;

    eprintf("\n%.1f%% of the mass, %e Msun, is within r=%e kpc\n", 
        100.*frac, env->M * (Mfrac_index+1) * env->units.M/env->units.Msun,
        R_Mfrac*env->units.L/env->units.kpc);

    double t_dynamical = 1/sqrt(env->units.G * RHO(Mfrac_index+1,R_Mfrac)
                       * env->units.T / env->units.Myr);
    //double t_crossing  = sqrt(env->units.G * RHO(Mfrac_index+1,R_Mfrac))
    double t_crossing  = R_Mfrac / p[Mfrac_index].V
                       * env->units.T / env->units.Myr;
    double t_relax     = t_crossing * env->N / log(env->N) / 8
                       * env->units.T / env->units.Myr;

    eprintf("The dynamical time 1/sqrt(G * rho(< r)) is ......... %e Myr\n", t_dynamical);
    eprintf("The crossing time r/v is ........................... %e Myr\n", t_crossing);
    eprintf("The relaxation time 0.1*t_crossing * N / ln N is ... %e Myr\n", t_relax);
        
#if 0
    eprintf("fractional mass is %e\n", env->M * (Mfrac_index+1));
    eprintf("fractional mass is %e Msun\n", env->M * (Mfrac_index+1) * env->units.M/env->units.Msun);
    eprintf("fractional mass (%.2f) radius is "POST" [%e]\n", frac, (pos_t)R_Mfrac, R_Mfrac / R);
    eprintf("fractional mass (%.2f) radius is %e kpc [%e]\n", 
            frac, R_Mfrac*env->units.L/env->units.kpc, R_Mfrac / R);
#endif

#if 0
    eprintf("tdyn ~ %e\n", 1/sqrt(RHO(Mfrac_index+1,R_Mfrac)));
    eprintf("tdyn ~ %e Myr\n", 1/sqrt(RHO(Mfrac_index+1,R_Mfrac)) * env->units.T / env->units.Myr);
#endif

    fp = fopen("tdyn.tdyn", "w"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        t_dynamical = 1/sqrt(env->units.G * RHO(i+1,R)
                    * env->units.T / env->units.Myr);
        fprintf(fp, "%e %e\n", R, t_dynamical);
        //fprintf(fp, "%ld %e\n", i, 1/sqrt(RHO(i+1,R)));
    }
    fclose(fp);

    for (i=0; i < env->N; i++)
    {
        double R = p[i].R;
        //if (0.99*R_Mfrac <= R && R <= 1.01*R_Mfrac)
        if (0.9999*R_Mfrac <= R && R <= 1.0001*R_Mfrac)
            printf("%ld\n", p[i].id);
    }


    free(env);

    return EXIT_SUCCESS;
}

