/*
 * A 64-bit unsigned integer has 20 decimal digits of which 19 can range
 * from 0 to 9. The scales for a given simulation must be chosen so that
 * there is enough precision to the "right" of the decimal point and that
 * the "left" side can also range over the size of the simulation. At
 * least 10 decimal places is a good start, but then more then better
 * to conserve energy.
 */

#define _GNU_SOURCE
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <signal.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "chaodyn.h"
#include "io.h"
#include "ic.h"
#include "io_tipsy.h"


//==============================================================================
//                                 integrate
//
// Integrate using a time-reversible leap-frog DKD algorithm.
//==============================================================================
// Plummer softening
#define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
#define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
void integrate(env_t *env)
{
    size_t i,j; 
    size_t   N   = env->N;
    mass_t   M   = env->M;
    const tyme_t   dt  = env->dt;
    force_t       *Fx  = env->F[0];
    force_t       *Fy  = env->F[1];
    force_t       *Fz  = env->F[2];
    particle_t    *p   = env->p;

    force_t Fx0, Fy0, Fz0;

    force_t ***Ft = env->Ft;

    double x,y,z;
    int128_t dx,dy,dz;


    double r2,rinv,re;
    const double eps = env->eps;
    const double e   = 1. / eps; 
    const double G   = env->units.G;

//  eprintf("%ld] "POST" "POST" "POST"\n", 0, p[0].v[0], p[0].v[1], p[0].v[2]);
//  eprintf("%ld] "POST" "POST" "POST"\n", 1, p[1].v[0], p[1].v[1], p[1].v[2]);
//  eprintf("%ld] "POST" "POST" "POST"\n", 2, p[2].v[0], p[2].v[1], p[2].v[2]);

    //--------------------------------------------------------------------------
    // Drift
    //--------------------------------------------------------------------------
    #pragma omp parallel for firstprivate(N) shared(p)
    for (i=0; i < N; i++) 
    {
        p[i].x[0] += (pos_t)(p[i].v[0]);
        p[i].x[1] += (pos_t)(p[i].v[1]);
        p[i].x[2] += (pos_t)(p[i].v[2]);
    }

    //--------------------------------------------------------------------------
    // Acceleration / Potential
    //--------------------------------------------------------------------------
    #pragma omp parallel for \
                schedule(dynamic, 2) \
                firstprivate(M,N) \
                private(i,j) \
                private(x,y,z, dx,dy,dz,r2) \
                private(Fx0,Fy0,Fz0) \
                private(rinv,re) \
                shared(Fx,Fy,Fz) \
                shared(Ft) \
                shared(p) 
    for (i=0; i < N-1; i++)
    {
        x = p[i].x[0];
        y = p[i].x[1];
        z = p[i].x[2];

        for (j=i+1; j < N; j++)
        {
            dx = p[j].x[0] - x;
            dy = p[j].x[1] - y;
            dz = p[j].x[2] - z;
            r2 = dx*dx + dy*dy + dz*dz;

            rinv = 1.0 / sqrt(r2);

            re = sqrt(r2) * e;

            DBG(1) eprintf("r2=%e %ld,%ld rinv=%e re=%e gradphi=%e\n", r2,i,j, rinv, re, GRADPHI(re));
            DBG(1) eprintf("\n**\n%e\n**\n", ((M*e)*e)*(dx*rinv)*GRADPHI(re));
            
            // round'ing is important. Conserves momentum perfectly.
            double gradphi = GRADPHI(re);
            #define Fhat(x) round((G*((M*e)*e)*(((x)*rinv)*gradphi)))

            Fx0 = Fhat(dx);
            Fy0 = Fhat(dy);
            Fz0 = Fhat(dz);

            Ft[i][j][0] =  Fx0; Ft[i][j][1] =  Fy0; Ft[i][j][2] =  Fz0;
            Ft[j][i][0] = -Fx0; Ft[j][i][1] = -Fy0; Ft[j][i][2] = -Fz0;
        }
    }

    force_t Fxi, Fyi, Fzi;
    #pragma omp parallel for \
            private(j, Fxi, Fyi, Fzi) \
            shared(Fx,Fy,Fz) \
            firstprivate(N)
    for (i=0; i < N; i++) 
    {
        Fxi = Fyi = Fzi = 0;
        for (j=0; j < N; j++) 
        {
            Fxi += Ft[i][j][0];
            Fyi += Ft[i][j][1];
            Fzi += Ft[i][j][2];
        }

        Fx[i] = Fxi;
        Fy[i] = Fyi;
        Fz[i] = Fzi;

        DBG(1) eprintf("F[%ld] "FORCET" "FORCET" "FORCET"\n", i, Fx[i], Fy[i], Fz[i]);
        DBG(1) eprintf("%e\n", pow(Fxi,2)+pow(Fyi,2)+pow(Fzi,2));
    }

    //--------------------------------------------------------------------------
    // Kick
    //--------------------------------------------------------------------------
    #pragma omp parallel for firstprivate(N) shared(Fx,Fy,Fz,p)
    for (i=0; i < N; i++) 
    {
        p[i].v[0] += (vel_t)(Fx[i] * dt);
        p[i].v[1] += (vel_t)(Fy[i] * dt);
        p[i].v[2] += (vel_t)(Fz[i] * dt);
        DBG(1) eprintf("a[%ld] "ACCT" "ACCT" "ACCT"\n", i, (acc_t)Fx[i], (acc_t)Fy[i], (acc_t)Fz[i]);
        DBG(1) eprintf("v[%ld] "VELT" "VELT" "VELT"\n", i, p[i].v[0], p[i].v[1], p[i].v[2]);
    }

    //--------------------------------------------------------------------------
    // Drift
    //--------------------------------------------------------------------------
    #pragma omp parallel for firstprivate(N) private(x,y,z) shared(p) 
    for (i=0; i < N; i++) 
    {
        p[i].x[0] += (pos_t)(p[i].v[0]);
        p[i].x[1] += (pos_t)(p[i].v[1]);
        p[i].x[2] += (pos_t)(p[i].v[2]);

        DBG(1) eprintf("dx[%ld] "VELT"*"TIMET"=%ld\n", i, p[i].v[0],dt,p[i].v[0]*dt);
        //DBG(1) eprintf("dx[%ld] "POST" "POST" "POST"\n", i, dx, dy, dz);
        DBG(1) eprintf("x[%ld] "POST" "POST" "POST"\n", i, p[i].x[0], p[i].x[1], p[i].x[2]);
    }
}

//==============================================================================
//                                  capture
//==============================================================================
void capture(env_t *env, particle_t *p, image_t *img, int clear, int with_gradient)
{
    double scale = 1;

    if (clear) memset(img->image,
                      0, 
                      3*img->nc*img->nr*sizeof(*(img->image)));

    if (with_gradient)
        scale = pow(((double)env->step) / env->opt.Nsteps, 3);

    size_t i;
    for (i=0; i < env->N; i++)
    {
        int32_t c = ( p[i].x[0] + env->radius) / (2.0*env->radius) * env->opt.img_cols;
        int32_t r = (-p[i].x[1] + env->radius) / (2.0*env->radius) * env->opt.img_rows;
        if (!(0 <= r && r < env->opt.img_rows)) continue;
        if (!(0 <= c && c < env->opt.img_cols)) continue;
        img->image[3*(r*img->nc + c) + 0] = round(env->class_color[p[i].class][0] * scale);
        img->image[3*(r*img->nc + c) + 1] = round(env->class_color[p[i].class][1] * scale);
        img->image[3*(r*img->nc + c) + 2] = round(env->class_color[p[i].class][2] * scale);
    }
}

//==============================================================================
//                             reverse_velocities
//==============================================================================
void reverse_velocities(env_t *env, size_t skip)
{
    size_t i;
    for (i=0; i < env->N; i++) 
    {
        if (skip != -1 && skip == i) continue;
        env->p[i].v[0] *= -1;
        env->p[i].v[1] *= -1;
        env->p[i].v[2] *= -1;
    }
}

//==============================================================================
//                           add_phase_space_noise
//==============================================================================
void add_phase_space_noise(env_t *env, double eps)
{
    size_t i;
    for (i=0; i < env->N; i++) 
    {
#if 0
        env->p[i].x[0] *= 1.0 + (2*drand48() - 1) * eps;
        env->p[i].x[1] *= 1.0 + (2*drand48() - 1) * eps;
        env->p[i].x[2] *= 1.0 + (2*drand48() - 1) * eps;
        env->p[i].v[0] *= 1.0 + (2*drand48() - 1) * eps;
        env->p[i].v[1] *= 1.0 + (2*drand48() - 1) * eps;
        env->p[i].v[2] *= 1.0 + (2*drand48() - 1) * eps;
#else 
        env->p[i].x[0] += copysign(env->p[i].x[0]*eps, drand48()-.5);
        env->p[i].x[1] += copysign(env->p[i].x[1]*eps, drand48()-.5);
        env->p[i].x[2] += copysign(env->p[i].x[2]*eps, drand48()-.5);
        env->p[i].v[0] += copysign(env->p[i].v[0]*eps, drand48()-.5);
        env->p[i].v[1] += copysign(env->p[i].v[1]*eps, drand48()-.5);
        env->p[i].v[2] += copysign(env->p[i].v[2]*eps, drand48()-.5);
#endif
    }
}

#if 1
//==============================================================================
//                              periodic_report
//==============================================================================
env_t *g_env = NULL;
void periodic_report(int sig)
{
    static char buf[100];
    static size_t last_step=0;
    static size_t seconds=0;
    static float  rate = 0;
    static size_t eta = 0; 

    assert(sig == SIGALRM);

    if (g_env == NULL) return;
    if (last_step == 0) last_step = g_env->step;
    
    if (++seconds == 20)
    {
        rate = (g_env->step - last_step) / (float)seconds;
        last_step = g_env->step;
        seconds = 0;
    }

    if (rate)
    {
        eta = ceil((g_env->end_step - g_env->step) / rate);
        sprintf((char *)&buf, CL "Step %lu/%lu   [%.2f/s ETA %lus]", 
                g_env->step, g_env->end_step, rate, eta);
    }
    else
    {
        sprintf((char *)&buf, CL "Step %lu/%lu", g_env->step, g_env->end_step);
    }

    write(STDOUT_FILENO, buf, strlen(buf));
    fsync(STDOUT_FILENO);

    alarm(1);
}
#endif

//==============================================================================
//                                   stats
//==============================================================================
void stats(env_t *env, FILE *fpE)
{
    size_t i,j; 
    size_t   N   = env->N;
    mass_t   M   = env->M;

    force_t    *Fx  = env->F[0];
    force_t    *Fy  = env->F[1];
    force_t    *Fz  = env->F[2];
    particle_t *p   = env->p;

    energy_t T = 0;                 /* Kinetic energy   */
    energy_t U = 0;                 /* Potential energy */
    energy_t Tavg = 0;              /* Virial Theorem   */
    vel_t    Px=0,Py=0,Pz=0;        /* Momentum         */
    pos_t    CMx=0, CMy=0, CMz=0;
    double   Lx =0,Ly =0,Lz =0;        /* Angular momentum */
    double   Lxb=0,Lyb=0,Lzb=0;     /* Angular momentum */

    double L,P,CM;

    double r2,re;
    double x,y,z;
    double vx,vy,vz;
    int128_t dx,dy,dz;
    const double eps = env->eps;
    const double e   = 1. / eps; 
    const double G   = env->units.G;

    #pragma omp parallel for \
                private(x,y,z,dx,dy,dz) \
                shared(p) \
                reduction(+:T)   reduction(+:Tavg) \
                reduction(+:Lx)  reduction(+:Ly)  reduction(+:Lz)  \
                reduction(+:Lxb) reduction(+:Lyb) reduction(+:Lzb) \
                reduction(+:Px)  reduction(+:Py)  reduction(+:Pz)  \
                reduction(+:CMx) reduction(+:CMy) reduction(+:CMz)  
    for (i=0; i < N; i++) 
    {
        x  = p[i].x[0]; y  = p[i].x[1]; z  = p[i].x[2];
        vx = p[i].v[0]; vy = p[i].v[1]; vz = p[i].v[2];

        T    += pow(vx,2) + pow(vy,2) + pow(vz,2);
        Tavg += x*Fx[i] + y*Fy[i] + z*Fz[i];

        Lx += y*vz - z*vy;
        Ly += z*vx - x*vz;
        Lz += x*vy - y*vx;

        Lxb += fabs(y*vz) + fabs(z*vy);
        Lyb += fabs(z*vx) + fabs(x*vz);
        Lzb += fabs(x*vy) + fabs(y*vx);

        Px += p[i].v[0];
        Py += p[i].v[1];
        Pz += p[i].v[2];

        CMx += p[i].x[0];
        CMy += p[i].x[1];
        CMz += p[i].x[2];
    }

    T /= 2;
    T *= M;
    Tavg /= 2;

    #pragma omp parallel for \
                schedule(dynamic, 2) \
                firstprivate(M,N) \
                private(i,j) \
                private(x,y,z, dx,dy,dz,r2) \
                private(re) \
                shared(p) \
                reduction(-:U)
    for (i=0; i < N-1; i++)
    {
        x = p[i].x[0];
        y = p[i].x[1];
        z = p[i].x[2];

        for (j=i+1; j < N; j++)
        {
            dx = (int128_t)p[j].x[0] - x;
            dy = (int128_t)p[j].x[1] - y;
            dz = (int128_t)p[j].x[2] - z;
            r2 = dx*dx + dy*dy + dz*dz;

            assert(r2 != 0);
            re = sqrt(r2) * e;

            U -= (G*M*M) / sqrt(pow(eps,2) + r2);
        }
    }

    Lx /= (Lxb == 0 ? 1 : Lxb);
    Ly /= (Lyb == 0 ? 1 : Lyb);
    Lz /= (Lzb == 0 ? 1 : Lzb);
    L  = sqrtl(Lx*Lx + Ly*Ly + Lz*Lz);
    //P  = round(sqrtl(((int128_t)Px)*Px + Py*Py + Pz*Pz));
    P  = sqrtl(pow(Px,2)  + pow(Py,2)  + pow(Pz,2));
    CM = sqrtl(pow(CMx,2) + pow(CMy,2) + pow(CMz,2)) / N / (env->units.kpc);

    double vir = (Tavg == 0) ? 0 : T/Tavg/M;

    fprintf(fpE, "%5ld "ENGYT" "ENGYT" "ENGYT" % .6e % .6e % .6e % .6e "VELT" "VELT" "VELT" % .6e "ENGYT" "POST" "POST" "POST" % .6e\n", 
            env->step, T,U, (T+U), Lx,Ly,Lz,L, Px,Py,Pz,P, vir, CMx,CMy,CMz,CM);

    fflush(fpE);
}

//==============================================================================
//                                  new_env
//==============================================================================
env_t *new_env()
{
    env_t *env;
    env = malloc(sizeof(*env)); 
    assert(env != NULL);

    //--------------------------------------------------------------------------
    // Some default values.
    //--------------------------------------------------------------------------
    env->N                    = 0;
    env->M                    = 0;
    env->dt                   = 2;
    env->radius               = 0;
    env->p0                   = NULL;
    env->p                    = NULL;
    env->F[0]                 = 0;
    env->F[1]                 = 0;
    env->F[2]                 = 0;
    env->opt.Nclasses         = 0;
    env->class_color          = NULL;
    env->step                 = 0;
    env->opt.Nsteps           = 0;
    env->opt.output_every     = 1;
    env->seed                 = time(NULL);
    env->opt.reverse_at       = 0;
    env->opt.save_image       = 0;
    env->opt.save_path_image  = 0;
    env->opt.save_image_every = 0;
    env->opt.save_path_every  = 1;
    env->opt.modify_mode      = 1;
    env->opt.tag              = strcpy(malloc(9), "chaodyn");
    env->opt.inputfile        = NULL;
    env->opt.verbosity        = 1;
    env->opt.dump_sim         = 0;
    env->opt.start_step       = 0;
    env->opt.logfp            = NULL;
    env->opt.stdout           = stdout;
    env->opt.img_cols         = 1024;
    env->opt.img_rows         = 1024;
    env->opt.save             = 1;
    env->opt.save_every       = 1;
    env->opt.ic_gen_name      = NULL;
    env->opt.restart          = 0;

    env->units.Myr            = 1;
    env->units.kpc            = 1;
    env->units.Msun           = 1;

    env->units.G              = 1;
    env->units.M              = 1;
    env->units.T              = 1;
    env->units.L              = 1;

    env->opt.X.dup_ic         = 0;

    env->opt.tipsy            = 0;

    return env;
}

//==============================================================================
//                                   usage
//==============================================================================
void usage()
{
    fprintf(stderr, 
    "Usage: chaodyn [OPTION]... [-i inputfile | --gen-ic=<ic-name>]\n"
    "\n"
    "where OPTION can be\n"
    "  -N N                 Use N particles.\n"
    "  -v                   Increase verbose output. This can be used repeatedly. (Default 1)\n"
    "  -q                   Decrease verbose output. This can be used repeatedly.\n"
    "  --seed=N             Initialize random number generator (srand48) with N.\n"
    "  --steps=N            Run simulation for N steps.\n"
    "  --reverse-at=WHEN    When N=WHEN reverse velocities.\n"
    "  --save-image=[EVERY] Save an image every EVERY frames.\n"
    "  --save-path-image    Create a single image that shows the particles at each\n"
    "                       time step. Use --save-image with a value for EVERY to\n"
    "                       control which time steps are plotted.\n"
    "  --start-step=N       Set the first simulation step to N.\n"
    "  --modify-mode=MODE   MODE can be one of [1,2,3]. When reversing,\n"
    "                       (1) Don't reverse the first particle.\n"
    "                       (2) Flip the least significant bit in the velocity\n"
    "                           for the first particle.\n"
    "                       (3) Remove the last particle\n"
    "  --tag=NAME           Prefix output files with NAME. Defaults to 'chaodyn'.\n"
    "  --dump-sim           Write information about inputfile to the screen before\n"
    "                       continuing.\n"
    "\n"
    "Send bug reports, comments, and questions to Jonathan Coles <jonathan@physik.uzh.ch>\n"
    "\n");
    exit(2);
}

int main(int argc, char **argv)
{
    size_t i,j;

    env_t *env = new_env();

    struct iclist_s *ics = iclist();

    while (1) 
    {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"seed",            1, 0, 0},
            {"steps",           1, 0, 0},
            {"reverse-at",      1, 0, 0},
            {"save-image",      2, 0, 0},
            {"save-path-image", 2, 0, 0},
            {"dump-sim",        2, 0, 0},
            {"start-step",      1, 0, 0},
            {"run-backward",    0, 0, 0},
            {"modify-mode",     1, 0, 0},
            {"tag",             1, 0, 0},
            {"verbosity",       1, 0, 0},
            {"dt",              1, 0, 0},
            {"no-save",         0, 0, 0},
            {"save",            1, 0, 0},
            {"gen-ic",          1, 0, 0},
            {"restart",         0, 0, 0},
            {"list-ic",         0, 0, 0},
            {"tipsy",           0, 0, 0},
            {0, 0, 0, 0}
        };

        int c = getopt_long(argc, argv, "N:vqi:X:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) 
        {
            #define OPTSTR(s) (!strcmp(s, long_options[option_index].name))
            case 0:
                     if OPTSTR("seed")            env->seed                = atol(optarg);
                else if OPTSTR("steps")           env->opt.Nsteps          = atol(optarg);
                else if OPTSTR("reverse-at")      env->opt.reverse_at      = atol(optarg);
                else if OPTSTR("save")
                {
                    env->opt.save = 1;
                    if (optarg)
                        env->opt.save_every = atol(optarg);
                    else
                        env->opt.save_every = 1;
                }
                else if OPTSTR("save-image")
                {
                    env->opt.save_image = 1;
                    if (optarg)
                        env->opt.save_image_every = atol(optarg);
                }
                else if OPTSTR("save-path-image") 
                {
                    env->opt.save_path_image = 1;
                    if (optarg)
                        env->opt.save_path_every = atol(optarg);
                }
                else if OPTSTR("start-step")
                {
                    env->opt.start_step = atol(optarg);
                    env->step           = env->opt.start_step;
                }
                else if OPTSTR("modify-mode")     env->opt.modify_mode     = atol(optarg);
                else if OPTSTR("verbosity")       env->opt.verbosity       = atol(optarg);
                else if OPTSTR("tag") 
                {
                    env->opt.tag = malloc(strlen(optarg)+1);
                    strcpy(env->opt.tag, optarg);
                }
                else if OPTSTR("dt")              env->dt = atof(optarg);
                else if OPTSTR("dump-sim")
                {
                    env->opt.dump_sim = 1;
                    if (optarg)
                    {
                        if (!strcmp("header", optarg)) env->opt.dump_sim = 1;
                        if (!strcmp("all", optarg))    env->opt.dump_sim = 2;
                    }
                }
                else if OPTSTR("no-save") env->opt.save = env->opt.save_every = 0;
                else if OPTSTR("gen-ic")
                {
                    if (env->opt.ic_gen_name != NULL) usage();
                    env->opt.ic_gen_name = malloc(strlen(optarg)+1);
                    strcpy(env->opt.ic_gen_name, optarg);
                }
                else if OPTSTR("restart")
                {
                    env->opt.restart = 1;
                }
                else if OPTSTR("list-ic")
                {
                    for (i=0; ics[i].name != NULL; i++)
                        eprintf("%15s - %s\n", ics[i].name, ics[i].desc);
                    exit(EXIT_SUCCESS);
                }
                else if OPTSTR("tipsy")
                {
                    env->opt.tipsy = 1;
                }
                break;

            case 'i':
                env->opt.inputfile = malloc(strlen(optarg)+1);
                strcpy(env->opt.inputfile, optarg);
                break;
            case 'v':
                env->opt.verbosity++;
                break;
            case 'q':
                if (env->opt.verbosity > 0) env->opt.verbosity--;
                break;
            case 'N':
                if (!optarg) usage();
                env->N = atol(optarg);
                break;
            case 'X':
                if (!optarg) usage();
                if (!strcmp("dup-ic", optarg)) env->opt.X.dup_ic = 1;
                else usage();

                break;
            default:
                usage();
                break;
        }
    }

    if ((env->step > env->opt.Nsteps)
    ||  (env->opt.modify_mode > 7)
    ||  (env->opt.inputfile != NULL && env->opt.ic_gen_name != NULL)
    ||  (env->opt.inputfile == NULL && env->opt.ic_gen_name == NULL)
    )
    {
         eprintf("SADFASDFS\n");
         usage();
    }

    //--------------------------------------------------------------------------
    // Open the log file.
    //--------------------------------------------------------------------------
    char *logfile = malloc(strlen(env->opt.tag)+1+4+1);
    time_t now = time(NULL);
    sprintf(logfile, "%s.log", env->opt.tag);
    env->opt.logfp = fopen(logfile, "a+");
    fprintf(env->opt.logfp, "New log created  %s\n", ctime(&now));
    assert(env->opt.logfp != NULL);
    free(logfile);

    //--------------------------------------------------------------------------
    // Generate initial conditions and save them for later comparison.
    //--------------------------------------------------------------------------
    if (env->opt.inputfile)
    {
        load(env, env->opt.inputfile);

        if (!env->opt.restart)
        {
            env->step = 0;
        }

        if (env->opt.X.dup_ic)
        {
            ic_colliding_halos(env);
            save_ic(env);
            exit(EXIT_SUCCESS);
        }
    }
    else if (env->opt.ic_gen_name)
    {
        for (i=0; ics[i].name != NULL; i++)
        {
            if (!strcmp(ics[i].name, env->opt.ic_gen_name))
            {
                ics[i].f(env);
                if (env->opt.tipsy)
                {
                    save_ic_tipsy(env);
                }
                else
                {
                    save_ic(env);
                }
                exit(EXIT_SUCCESS);
            }
        }
        eprintf("Unknown initial condition specified.\n");
        usage();
    }

    //--------------------------------------------------------------------------
    // Open the Energy file.
    //--------------------------------------------------------------------------
    char *efile = malloc(strlen(env->opt.tag)+1+strlen("ENERGY")+1);
    sprintf(efile, "%s.ENERGY", env->opt.tag);
    FILE *fpE = fopen(efile, "w");
    assert(fpE != NULL);
    free(efile);

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    VL(1) LOG("Random seed %ld\n", env->seed);
    srand48(env->seed);

    VL(1) LOG("sizeof(pos_t)=%ld\n", sizeof(pos_t));
    VL(1) LOG("sizeof(int128_t)=%ld\n", sizeof(int128_t));

    //--------------------------------------------------------------------------
    // Memory to store the initial conditions. We compare the final state
    // to this in the case we are checking for reversibility.
    //--------------------------------------------------------------------------
    assert(env->N > 0);
    assert(env->p != NULL);
    env->p0 = malloc(env->N * sizeof(*env->p0)); assert(env->p0 != NULL);
    memcpy(env->p0, env->p, env->N * sizeof(*env->p0));

    //--------------------------------------------------------------------------
    // Now allocate space for the force calculations.
    //--------------------------------------------------------------------------
    env->F[0] = malloc(env->N * sizeof(*env->F[0])); assert(env->F[0] != NULL);
    env->F[1] = malloc(env->N * sizeof(*env->F[1])); assert(env->F[1] != NULL);
    env->F[2] = malloc(env->N * sizeof(*env->F[2])); assert(env->F[2] != NULL);
    env->Ft   = malloc(env->N * sizeof(*env->Ft));   assert(env->Ft   != NULL);
    for (i=0; i < env->N; i++)
    {
        env->Ft[i]    = malloc(env->N * sizeof(**env->Ft));
        assert(env->Ft[i] != NULL);
        for (j=0; j < env->N; j++)
            env->Ft[i][j] = malloc(3 * sizeof(***env->Ft));
        env->Ft[i][i][0] = 
        env->Ft[i][i][1] = 
        env->Ft[i][i][2] = 0;
    }

    //--------------------------------------------------------------------------
    // Setup variables for saving images.
    //--------------------------------------------------------------------------
    env->snapshot.nr = env->opt.img_rows;
    env->snapshot.nc = env->opt.img_cols;
    env->path.nr     = env->opt.img_rows;
    env->path.nc     = env->opt.img_cols;

    env->snapshot.image = calloc(env->snapshot.nr, env->snapshot.nc * 3);
    env->path.image     = calloc(env->path.nr,     env->path.nc     * 3);
    env->class_color    = malloc(env->opt.Nclasses * sizeof(*env->class_color));
    for (i=0; i < env->opt.Nclasses; i++)
    {
        env->class_color[i]    = malloc(3 * sizeof(**env->class_color));
        env->class_color[i][0] = 255;
        env->class_color[i][1] = 255;
        env->class_color[i][2] = 255;
    }

#if 0
    env->p[0].class = 0;
    env->p[1].class = 1;
    env->p[2].class = 2;
#endif

    env->class_color[0][0] = 255;
    env->class_color[0][1] = 0;
    env->class_color[0][2] = 0;

#if 0
    env->class_color[1][0] = 0;
    env->class_color[1][1] = 255;
    env->class_color[1][2] = 0;

    env->class_color[2][0] = 0;
    env->class_color[2][1] = 0;
    env->class_color[2][2] = 255;
#endif

    //--------------------------------------------------------------------------
    // Print a status update every second. Using a signal provides a regular
    // output regardless of how long the integration takes.
    //--------------------------------------------------------------------------
    g_env = env;
    VL(1)
    {
        signal(SIGALRM, periodic_report);
        periodic_report(SIGALRM);
    }

    //--------------------------------------------------------------------------
    // Run the simulation. Possibly reverse the velocities at a given time.
    //--------------------------------------------------------------------------

    env->end_step = env->step + env->opt.Nsteps;

    int forward = 1;

    while (env->step != env->end_step)
    {
        env->step++;

        if (env->step == env->opt.reverse_at)
        {
            forward = 0;

            switch (env->opt.modify_mode)
            {
                case 1:
                    LOG("Reversing all velocities.\n");
                    reverse_velocities(env, -1);
                    break;
                case 2:
                    LOG("Not reversing velocity for particle %ld.\n", 0L);
                    reverse_velocities(env, 0);
                    break;
                case 3:
#if WITH_INTEGERS
                    LOG("Flipping 1 velocity bit for particle %ld.\n", 0L);
                    env->p[i].v[0] ^= 1<<0;
                    env->p[i].v[1] ^= 1<<0;
                    env->p[i].v[2] ^= 1<<0;
#else
                    env->p[i].v[0] *= 1+1e-10;
                    env->p[i].v[1] *= 1+1e-10;
                    env->p[i].v[2] *= 1+1e-10;
#endif
                    reverse_velocities(env, -1);
                    break;
                case 4:
                    LOG("Removing particle %ld.\n", env->N-1);
                    env->N--;
                    reverse_velocities(env, -1);
                    break;
                case 5:
                    LOG("Changing dt from "TIMET" to "TIMET"\n", env->dt, env->dt*2);
                    env->dt *= 2;
                    break;
                case 6:
                    LOG("Changing dt from "TIMET" to "TIMET"\n", env->dt, env->dt/2);
                    env->dt /= 2;
                    break;
                case 7:
                    LOG(CL "Reversing all velocities.\n");
                    reverse_velocities(env, -1);
                    LOG(CL "Added 0.1%% errors to all phase-space coordinates\n");
                    add_phase_space_noise(env, 0.001);
            }
        }

        if (forward) stats(env, fpE);

        //----------------------------------------------------------------------
        integrate(env);
        //----------------------------------------------------------------------

        if (!forward) stats(env, fpE);

        //----------------------------------------------------------------------
        //----------------------------------------------------------------------
        if (env->opt.save_path_every && (env->step % env->opt.save_path_every) == 0)
            capture(env, env->p, &env->path, 0, 1);

        if (env->opt.save_image_every && (env->step % env->opt.save_image_every) == 0)
        {
            capture(env, env->p, &env->snapshot, 1, 0);
            save_snapshot(env);
        }

        if (env->opt.save_every && (env->step % env->opt.save_every) == 0)
            save(env);
    }

    if (env->opt.save_path_image)
        save_path_image(env);

    if (env->opt.save_image)
    {
#if 0
        pos_t R = env->radius;
        for (i=0; i < env->N; i++)
        {
            pos_t r  = sqrt(pow(env->p[i].x[0],2) 
                          + pow(env->p[i].x[1],2) 
                          + pow(env->p[i].x[2],2));

            pos_t r0 = sqrt(pow(env->p0[i].x[0],2)
                          + pow(env->p0[i].x[1],2)
                          + pow(env->p0[i].x[2],2));

            R = MAX(R, MAX(r, r0));
        }
#endif
        //env->radius *= 3;

        env->class_color[0][0] = 255; 
        env->class_color[0][1] = 0; 
        env->class_color[0][2] = 0;
        // env->class_color[1][0] = 0;   
        // env->class_color[1][1] = 0; 
        // env->class_color[1][2] = 255;
        capture(env, env->p0, &env->snapshot, 1,0);

        env->class_color[0][0] = 0;
        env->class_color[0][1] = 255;
        env->class_color[0][2] = 0;
        // env->class_color[1][0] = 255;
        // env->class_color[1][1] = 255;
        // env->class_color[1][2] = 255;
        capture(env, env->p, &env->snapshot, 0,0);

        save_comparison_image(env);
    }

    if (env->opt.Nsteps != 0)
    {
        //save(env);

        //--------------------------------------------------------------------------
        // Check that our final conditions are the same as the initial ones.
        //--------------------------------------------------------------------------
        if (env->opt.reverse_at)
        {
            size_t n_bad = 0;
            for (i=0; i < env->N; i++)
            {
                if ((env->p[i].x[0] != env->p0[i].x[0])
                ||  (env->p[i].x[1] != env->p0[i].x[1])
                ||  (env->p[i].x[2] != env->p0[i].x[2])
                ||  (env->p[i].v[0] != -env->p0[i].v[0])
                ||  (env->p[i].v[1] != -env->p0[i].v[1])
                ||  (env->p[i].v[2] != -env->p0[i].v[2]))
                {
                    n_bad++;
                    VL(2)
                    {
                        LOG("IC[%ld] x("POST" "POST" "POST") v("VELT" "VELT" "VELT")\n"
                            "  [%ld] x("POST" "POST" "POST") v("VELT" "VELT" "VELT")\n",
                            i, env->p0[i].x[0], env->p0[i].x[1], env->p0[i].x[2],
                               env->p0[i].v[0], env->p0[i].v[1], env->p0[i].v[2],
                            i, env->p[i].x[0],  env->p[i].x[1],  env->p[i].x[2],
                               env->p[i].v[0],  env->p[i].v[1],  env->p[i].v[2]);
                    }
                }
            }
            
            VL(1)
            {
                if (n_bad == 0)
                    LOG("PERFECT REVERSAL!\n");
                else
                {
                    LOG("%ld/%ld particles are not in the right place. "
                        "Did you reverse at the right place? Should be at middle step + 1.\n", n_bad, env->N);
                }
            }
        }
    }

    fclose(fpE);
    fclose(env->opt.logfp);

#define FREE(ptr) do { if (ptr) free(ptr); } while(0)
    FREE(env->p0);
    FREE(env->p);
    FREE(env->F[0]);
    FREE(env->F[1]);
    FREE(env->F[2]);
    FREE(env->snapshot.image);
    FREE(env->path.image);
    FREE(env->opt.tag);
    FREE(env->opt.inputfile);
    FREE(env);

    return 0;
}


#if 0
#if WITH_SOFTENING
            const double e  = 1. / env->eps; 
            const double re = sqrt(r2) * e;
            //printf("re=%e\n", re);
            if (re <= 1)
            {
#define K1 ((175./16. - ((147./8. - (135./16.*re)*re)*re)*re)*re)
                Fx0 = round(G*((M*e)*e)*((dx*rinv)*K1));
                Fy0 = round(G*((M*e)*e)*((dy*rinv)*K1));
                Fz0 = round(G*((M*e)*e)*((dz*rinv)*K1));
            }
            else
#endif
            {
                Fx0 = round(G*M*(((dx * rinv) * rinv) * rinv));
                Fy0 = round(G*M*(((dy * rinv) * rinv) * rinv));
                Fz0 = round(G*M*(((dz * rinv) * rinv) * rinv));
            }
#endif

