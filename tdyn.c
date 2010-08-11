#define _GNU_SOURCE
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "chaodyn.h"
#include "io.h"

#define PI (3.1415926535)

#define RHO(N,R)      ((N)*env->M / (4./3.*PI*pow((R),3)))
#define RHOa(N,R0,R1) ((N)*env->M / (4./3.*PI*(pow((R1),3) - pow((R0),3))))

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
    char *tag = "tdyn";
    char fname[1024];

    cm_t CM;

    env_t *env = calloc(1, sizeof(env_t)); assert(env != NULL);

    if (argc < 2) usage();

    while(1)
    {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"tag", 1, 0, 0},
            {0, 0, 0, 0}
        };

        int c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) 
        {
            #define OPTSTR(s) (!strcmp(s, long_options[option_index].name))
            case 0:
                if OPTSTR("tag")
                    tag = optarg;
                break;
        }
    }

    if (optind == argc) usage();

    load(env, argv[optind]);
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
#if 0
    for (i=0; i < env->N; i++)
    {
        p[i].x[0] -= CM.x[0];
        p[i].x[1] -= CM.x[1];
        p[i].x[2] -= CM.x[2];
        p[i].R = MAG(p[i].x);
    }
#endif
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


    //--------------------------------------------------------------------------
    // Output important times
    //--------------------------------------------------------------------------


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
        

    //--------------------------------------------------------------------------
    // Mass vs. Radius
    //--------------------------------------------------------------------------


    snprintf(fname, sizeof(fname)-1, "%s.MvR", tag);
    FILE *fp = fopen(fname, "wt"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R * env->units.L / env->units.kpc;
        fprintf(fp, "%e %e\n", R, ((double)(i+1)) / env->N);
    }
    fclose(fp);

    snprintf(fname, sizeof(fname)-1, "%s.dencum", tag);
    fp = fopen(fname, "wt"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        fprintf(fp, "%e %e\n", R, RHO(i+1,R));
    }
    fclose(fp);


    //--------------------------------------------------------------------------
    // Density vs. Radius
    //--------------------------------------------------------------------------


    snprintf(fname, sizeof(fname)-1, "%s.den", tag);
    fp = fopen(fname, "wt"); assert(fp != NULL);
    size_t j=0;
    double R0 = 0;
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R * env->units.L / env->units.kpc;
        if (R0 != 0)
            fprintf(fp, "%e %e\n", R, RHOa(i-j+1,R0,R) * env->units.M / env->units.Msun);
        R0 = R;
        j = i;
    }
    fclose(fp);


    //--------------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // (dln rho / dln R) vs. Radius
    //--------------------------------------------------------------------------
    snprintf(fname, sizeof(fname)-1, "%s.dln", tag);
    fp = fopen(fname, "wt"); assert(fp != NULL);
    j=0;
    R0 = 0;
    double rho0 = 0;
    double rho  = 0;

#if 1
    typedef struct
    {
        double R;
        int64_t Nenc;
    } bin_t;

    const int nbins = 100;
    bin_t bin[nbins+1];
    R0 = p[env->N-1].R;
    for (i=0; i <= nbins; i++)
    {
        bin[i].R = R0 * i / nbins;
        bin[i].Nenc = 0;
    }

    for (i=0; i < env->N; i++)
        for (j=0; j <= nbins; j++)
            if (p[i].R <= bin[j].R) bin[j].Nenc++;

    double rho_frac = 1; //RHO(Mfrac_index, p[Mfrac_index].R * env->units.L / env->units.kpc);

    R0 = 0;
    for (j=1; j <= nbins; j++)
    {
        double R = bin[j].R * env->units.L / env->units.kpc;
        //eprintf("%ld %e\n", bin[j].Nenc, R);

        rho = RHO(bin[j].Nenc,R) * env->units.M / env->units.Msun;
        rho /= rho_frac;

        eprintf("%e %e %e %e\n", rho, rho0, R, R0);
        if (rho == 0 || rho0 == 0)
            fprintf(fp, "%e 0\n", R);
        else
            fprintf(fp, "%e %e\n", R, (log(rho)-log(rho0)) / (log(R)-log(R0)));
        R0 = R;
        rho0 = rho;
    }

#else

    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R * env->units.L / env->units.kpc;
        if (R0 != 0)
        {
#if 0
            double drho = RHOa(i-j+1,R0,R) * env->units.M / env->units.Msun;
            double dR = R - R0;
#else
            rho = RHO(i+1,R) * env->units.M / env->units.Msun;
#endif
            fprintf(fp, "%e %e\n", R, (log(rho)-log(rho0)) / (log(R)-log(R0)));
        }
        R0 = R;
        rho0 = rho;
        j = i;
    }
#endif
    fclose(fp);

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


    //--------------------------------------------------------------------------
    // Dynamical time vs. Radius
    //--------------------------------------------------------------------------


    snprintf(fname, sizeof(fname)-1, "%s.tdyn", tag);
    fp = fopen(fname, "wt"); assert(fp != NULL);
    for (i=0; i < env->N; i+=env->N*0.01)
    {
        double R = p[i].R;
        t_dynamical = 1/sqrt(env->units.G * RHO(i+1,R)
                    * env->units.T / env->units.Myr);
        fprintf(fp, "%e %e\n", R, t_dynamical);
        //fprintf(fp, "%ld %e\n", i, 1/sqrt(RHO(i+1,R)));
    }
    fclose(fp);


    //--------------------------------------------------------------------------



//  for (i=0; i < env->N; i++)
//  {
//      double R = p[i].R;
//      //if (0.99*R_Mfrac <= R && R <= 1.01*R_Mfrac)
//      if (0.9999*R_Mfrac <= R && R <= 1.0001*R_Mfrac)
//          printf("%ld\n", p[i].id);
//  }


    free(env);

    return EXIT_SUCCESS;
}

