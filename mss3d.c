#define _GNU_SOURCE
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <png.h>
#include <time.h>

#define WITH_SOFTENING 1

#define DBG_LEVEL 0
#define DBG(__lvl) if (DBG_LEVEL >= __lvl)

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#if 1
typedef int64_t thyme_t;            /* Myr              */
typedef int64_t pos_t;              /* pc               */
typedef int64_t vel_t;              /* km/s             */
typedef int64_t acc_t;              /* km/s^2           */
typedef int64_t force_t;            /* Msun km/s^2      */
typedef double mass_t;             /* Msun             */
typedef double energy_t;           /* */
typedef double soft_t;

#define TIMET "%ld"
#define POST  "%ld"
#define VELT  "%ld"
#define ACCT  "%ld"
#define RHOT  "%e"
#define MASST "%e"
#define ENGYT "%e"
#define FORCET "%ld"

#else
typedef double pos_t;              /* pc               */
typedef double vel_t;              /* km/s             */
typedef double   acc_t;              /* km/s^2           */
typedef double   mass_t;             /* Msun             */
typedef double thyme_t;            /* Myr              */
typedef double energy_t;           /* */
typedef double soft_t;
#define TIMET "%e"
#define POST  "%e"
#define VELT  "%e"
#define ACCT  "%e"
#define RHOT  "%e"
#define MASST "%e"
#define ENGYT "%e"
#endif

double G = 1;

//#define G (4.3e-3)                  /* pc/Msun(km/s)^2  */

typedef struct
{
    pos_t x[3];
    vel_t v[3];
    char color[3];
} particle_t;

typedef struct
{
    size_t   N;
    particle_t * restrict p0;
    particle_t * restrict p;
    force_t *    restrict F[3];

    soft_t eps;

    mass_t   M;
    thyme_t  dt;
    thyme_t  t;
    size_t  Nsteps;

    pos_t    radius;

    size_t step;
    size_t output_every;

    unsigned char *image;
    uint32_t img_rows;
    uint32_t img_cols;
} env_t;

FILE *fpE;  /* Energy statistics */

//------------------------------------------------------------------------------
// Integrate using a time-reversible leap-frog DKD algorithm.
//------------------------------------------------------------------------------
void integrate(env_t *env, int last_step)
{
    size_t i,j; 
    const size_t   N   = env->N;
    const mass_t   M   = env->M;
    const thyme_t  dt  = env->dt;
    force_t       *Fx  = env->F[0];
    force_t       *Fy  = env->F[1];
    force_t       *Fz  = env->F[2];
    particle_t    *p   = env->p;

    energy_t   U    = 0;         /* Potential energy */
    energy_t   T    = 0;         /* Kinetic energy   */
    int64_t    L[3] = {0,0,0};   /* Angular momentum */
    int64_t    P[3] = {0,0,0};   /* Momentum         */

    //--------------------------------------------------------------------------
    // Drift
    //--------------------------------------------------------------------------
    for (i=0; i < N; i++) 
    {
        const pos_t  x = p[i].x[0];
        const pos_t  y = p[i].x[1];
        const pos_t  z = p[i].x[2];
        const pos_t dx = (pos_t)((p[i].v[0] * dt)/2);
        const pos_t dy = (pos_t)((p[i].v[1] * dt)/2);
        const pos_t dz = (pos_t)((p[i].v[2] * dt)/2);
        p[i].x[0] = x + dx;
        p[i].x[1] = y + dy;
        p[i].x[2] = z + dz;
    }

    if (last_step) return;
    
    memset(Fx, 0, N*sizeof(*Fx)); 
    memset(Fy, 0, N*sizeof(*Fy)); 
    memset(Fz, 0, N*sizeof(*Fz)); 

    //--------------------------------------------------------------------------
    // Acceleration / Potential
    //--------------------------------------------------------------------------
    const double e  = 1. / env->eps; 
    for (i=0; i < N-1; i++)
    {
        const pos_t x = p[i].x[0];
        const pos_t y = p[i].x[1];
        const pos_t z = p[i].x[2];

        force_t Fxi=0, Fyi=0, Fzi=0;

        for (j=i+1; j < N; j++)
        {
            const pos_t dx = x - p[j].x[0];
            const pos_t dy = y - p[j].x[1];
            const pos_t dz = z - p[j].x[2];
            const double r2 = ((double)dx)*((double)dx) 
                            + ((double)dy)*((double)dy) 
                            + ((double)dz)*((double)dz);

            DBG(2) printf("r2=%e %ld,%ld\n", r2,i,j);
            assert(r2 != 0);
            //DBG(2) printf("r=%e\n", r);
            const double rinv = 1 / sqrt(r2);

            const double re = sqrt(r2) * e;

            // Plummer softening
#define PHI(r)     (1./sqrt(1.+(r)*(r)))
#define GRADPHI(r) (pow(sqrt(1. + (r)*(r)), -3) * (r))
            
            const force_t Fx0 = round(G*((M*e)*e)*((dx*rinv)*GRADPHI(re)));
            const force_t Fy0 = round(G*((M*e)*e)*((dy*rinv)*GRADPHI(re)));
            const force_t Fz0 = round(G*((M*e)*e)*((dz*rinv)*GRADPHI(re)));

            U -= G*M*M*e * PHI(re);

            Fxi   += Fx0; Fyi   += Fy0; Fzi   += Fz0;
            Fx[j] += Fx0; Fy[j] += Fy0; Fz[j] += Fz0;
        }

        Fx[i] -= Fxi; Fy[i] -= Fyi; Fz[i] -= Fzi;
    }

    //--------------------------------------------------------------------------
    // Kick
    //--------------------------------------------------------------------------
    for (i=0; i < N; i++) 
    {
        p[i].v[0] += (vel_t)(Fx[i] * dt); // should divide F by M here
        p[i].v[1] += (vel_t)(Fy[i] * dt);
        p[i].v[2] += (vel_t)(Fz[i] * dt);
        DBG(1) printf("a[%ld] "FORCET" "FORCET" "FORCET"\n", i, Fx[i], Fy[i], Fz[i]);
        DBG(1) printf("v[%ld] "VELT" "VELT" "VELT"\n", i, p[i].v[0], p[i].v[1], p[i].v[2]);
    }

    //--------------------------------------------------------------------------
    // Drift
    //--------------------------------------------------------------------------
    for (i=0; i < N; i++) 
    {
        const pos_t  x = p[i].x[0];
        const pos_t  y = p[i].x[1];
        const pos_t  z = p[i].x[2];
        const pos_t dx = (pos_t)((p[i].v[0] * dt)/2);
        const pos_t dy = (pos_t)((p[i].v[1] * dt)/2);
        const pos_t dz = (pos_t)((p[i].v[2] * dt)/2);
        p[i].x[0] = x + dx;
        p[i].x[1] = y + dy;
        p[i].x[2] = z + dz;

        DBG(1) printf("dx[%ld] "POST" "POST" "POST"\n", i, dx, dy, dz);
        DBG(1) printf("x[%ld] "POST" "POST" "POST"\n", i, p[i].x[0], p[i].x[1], p[i].x[2]);

        T += 0.5 * M * (pow(p[i].v[0],2) + pow(p[i].v[1],2) + pow(p[i].v[2],2));

        L[0] += y*p[i].v[2] - z*p[i].v[1];
        L[1] += z*p[i].v[0] - x*p[i].v[2];
        L[2] += x*p[i].v[1] - y*p[i].v[0];

        P[0] += p[i].v[0];
        P[1] += p[i].v[1];
        P[2] += p[i].v[2];
    }

    double Lmag = sqrt(pow(L[0],2) + pow(L[1],2) + pow(L[2],2));
    double Pmag = sqrt(pow(P[0],2) + pow(P[1],2) + pow(P[2],2));

    //fprintf(fpE, "%5ld %23ld %23ld %23e\n", env->step, T, U, L);
    fprintf(fpE, "%5ld "ENGYT" "ENGYT" %23e "ENGYT" %23e ""\n", env->step, T, U, Lmag, T+U, Pmag);
}

void sphere(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1)
{
    size_t i,j;
    double x,y,z, t,w;
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

        DBG(2) printf("%ld] %e %e %e\n", i, x, y, z);
    }
}

void ic_sphere(env_t *env)
{
    env->N      = env->N;
    env->M      = 1e22 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));

    sphere(env, 0,0,0, env->radius, 0, env->N);
}

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


void shell(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1)
{
    size_t i,j;
    double x,y,z, t,w;

    for (i=N0; i < N1; i++)
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
                double r2 = pow(env->p[i].x[0] - x, 2)
                          + pow(env->p[i].x[1] - y, 2)
                          + pow(env->p[i].x[2] - z, 2);
                if (r2 < 10000) goto redo;
            }
        } 

        env->p[i].x[0] = x;
        env->p[i].x[1] = y;
        env->p[i].x[2] = z;
        env->p[i].v[0] = 0;
        env->p[i].v[1] = 0;
        env->p[i].v[2] = 0;

        DBG(2) printf("%ld] %e %e %e\n", i, x, y, z);
    }
}

void ic_uniform_random_shell(env_t *env)
{
    env->N      = 10;
    env->M      = 1e22 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));

    shell(env, 0,0,0, env->radius, 0, env->N);
}

void ic_two_shells(env_t *env)
{
    env->N      = 1000;
    env->M      = 1e22 / 4;
    env->radius = 1e9;
    env->eps    = 1e8;
    env->p      = malloc(env->N * sizeof(*(env->p)));

    shell(env, -env->radius,0,0, env->radius/5, 0,env->N/2);
    shell(env,  env->radius,0,0, env->radius/5, env->N/2,env->N);

    fflush(stdout);
}

void ic_figure8(env_t *env)
{
    double kv = 1e11;
    double kx = 1e9;

    env->N = 3;
    env->M = 1 * pow(kv,2);
    env->p = malloc(env->N * sizeof(*(env->p)));
    env->radius = 1 << 30;
    env->eps=1e8;

    printf(MASST"\n", env->M);

    double vfac= kv/sqrt(kx);
    double rfac= kx;
    env->p[0].x[0] = -0.995492 * rfac;
    env->p[0].x[1] = -0.0      * rfac;
    env->p[0].x[2] = -0.0      * rfac;
    env->p[0].v[0] = -0.347902 * vfac;
    env->p[0].v[1] = -0.53393  * vfac;
    env->p[0].v[2] =  0.0      * vfac;

    env->p[1].x[0] =  0.995492 * rfac;
    env->p[1].x[1] =  0.0      * rfac;
    env->p[1].x[2] =  0.0      * rfac;
    env->p[1].v[0] = -0.347902 * vfac;
    env->p[1].v[1] = -0.53393  * vfac;
    env->p[1].v[2] =  0.0      * vfac;

    env->p[2].x[0] = 0.0       * rfac;
    env->p[2].x[1] = 0.0       * rfac;
    env->p[2].x[2] = 0.0       * rfac;
    env->p[2].v[0] = 0.695804  * vfac;
    env->p[2].v[1] = 1.06786   * vfac;
    env->p[2].v[2] = 0.0       * vfac;
}

void generate_ics(env_t *env)
{
    //ic_random(env);
    //ic_uniform_random_shell(env);
    ic_figure8(env);
    //ic_two_shells(env);
    //ic_sphere(env);
}

int save_image_png(char *filename, unsigned char *image, uint32_t nrows, uint32_t ncols)
{
    int i;

    char *fname = malloc(strlen(filename)+4+1);
    sprintf(fname, "%s.png", filename);

    FILE *fp = fopen(fname, "wb");

    if (fp == NULL)
    {
        fprintf(stderr, "Can't open %s\n", fname);
        free(fname);
        return 1;
    }
    free(fname);

    png_structp png_ptr = png_create_write_struct
       (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) return 1;

    png_init_io(png_ptr, fp);

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
       png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
       return 1;
    }


    png_set_IHDR(png_ptr, info_ptr, ncols, nrows,
           8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
           PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);

    int row_stride = ncols * 3;

    for (i=0; i < nrows; i++)
    {
        png_bytep row_pointer = & (image[i * row_stride]);
        png_write_row(png_ptr, row_pointer);
    }

    png_write_end(png_ptr, info_ptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(fp);

    return 0;
}

char save_color[3] = {255,0,0};
void save(env_t *env)
{
    size_t i;
    //memset(env->image, 0, env->img_rows*env->img_cols*3*sizeof(*env->image));
    for (i=0; i < env->N; i++)
    {
        int32_t c = (env->p[i].x[0] + env->radius) / (2*env->radius / env->img_cols);
        int32_t r = (-env->p[i].x[1] + env->radius) / (2*env->radius / env->img_rows);
        if (r >= env->img_rows) r = env->img_rows-1; else if (r < 0) r = 0;
        if (c >= env->img_cols) c = env->img_cols-1; else if (c < 0) c = 0;
#if 0
        env->image[3*(r*env->img_cols + c) + 0] = save_color[0];
        env->image[3*(r*env->img_cols + c) + 1] = save_color[1];
        env->image[3*(r*env->img_cols + c) + 2] = save_color[2];
#else
        env->image[3*(r*env->img_cols + c) + 0] = env->p[i].color[0];
        env->image[3*(r*env->img_cols + c) + 1] = env->p[i].color[1];
        env->image[3*(r*env->img_cols + c) + 2] = env->p[i].color[2];
#endif
    }

    char *fname = malloc(strlen("mss3d.")+5+1);
    sprintf(fname, "mss3d.%05ld", env->step);

    //save_image_png(fname, env->image, env->img_rows, env->img_cols);
}

void finalize(env_t *env)
{
#if 1
    save_image_png("mss3d", env->image, env->img_rows, env->img_cols);
#endif
}

#if 0
env_t *g_env = NULL;
void periodic_report(int sig)
{
    char buf[100];
    assert(sig == SIGALARM);

    if (g_env == NULL) return;
    sprintf(buf, "Step %ld/%ld   \r", g_env.step, g_env.Nsteps);
    write(STDOUT_FILENO, buf, sizeof(buf));
    fsync(STDOUT_FILENO);
}
#endif

int main(int argc, char **argv)
{
    size_t i,j;
    env_t env;

    time_t seed = time(NULL);
    printf("Random seed %ld\n", seed);
    srand48(seed);

    fpE = fopen("ENERGY", "w");
    assert(fpE != NULL);

    //--------------------------------------------------------------------------
    // Generate initial conditions and save them for later comparison.
    //--------------------------------------------------------------------------
    env.N            = 100;
    env.M            = 0;
    env.dt           = 2;
    env.radius       = 0;
    env.p0           = NULL;
    env.p            = NULL;
    env.F[0]         = 0;
    env.F[1]         = 0;
    env.F[2]         = 0;

    generate_ics(&env);
    assert(env.N > 0);
    assert(env.p != NULL);

    env.Nsteps       = 10000;
    env.output_every = 10;
    env.p0           = malloc(env.N * sizeof(*env.p0));
    env.F[0]         = malloc(env.N * sizeof(*env.F[0]));
    env.F[1]         = malloc(env.N * sizeof(*env.F[1]));
    env.F[2]         = malloc(env.N * sizeof(*env.F[2]));

    assert(env.Nsteps%2 == 0);

    env.img_rows = ceil(env.Nsteps / env.dt) / env.output_every;
    env.img_cols = 1024;
    env.img_rows = 1024;
    env.image = calloc(env.img_rows, env.img_cols * 3);
    
    memcpy(env.p0, env.p, env.N * sizeof(*env.p0));

    for (i=0; i < env.N/2; i++)
    {
        env.p[i].color[0] = 255; env.p[i].color[1] = 0; env.p[i].color[2] = 0;
    }
    for (i=env.N/2; i < env.N; i++)
    {
        env.p[i].color[0] = 0;   env.p[i].color[1] = 0; env.p[i].color[2] = 255;
    }

    save(&env);

    //--------------------------------------------------------------------------
    // Run the simulation. Half way through, reverse the velocities.
    //--------------------------------------------------------------------------
    for (env.step = 0; env.step < env.Nsteps; env.step++)
    {
#if 0
        if (env.step == env.Nsteps/2)
        {
            save_color[0] = 0;
            save_color[1] = 255;
            save_color[2] = 0;
            for (i=0; i < env.N; i++) 
            {
                switch (0)
                {
                    case 1:
                        if (i == 0) 
                        {
                            printf("Not reversing velocity for particle %ld\n", i);
                            continue;
                        }
                        break;
                    case 2:
                        if (i==0) 
                        {
                            printf("Flipping 1 velocity bit for particle %ld\n", i);
                            env.p[i].v[0] ^= 1<<0;
                            env.p[i].v[1] ^= 1<<0;
                            env.p[i].v[2] ^= 1<<0;
                        }
                        break;
                    case 3:
                        if (i == env.N-1)
                        {
                            printf("Removing particle %ld\n", i);
                            env.N--;
                        }
                        break;
                }

                env.p[i].v[0] *= -1;
                env.p[i].v[1] *= -1;
                env.p[i].v[2] *= -1;
            }
        }
#endif

        if ((env.step % 1) == 0)
        {
            printf("Step %ld/%ld   \r", env.step, env.Nsteps);
            fflush(stdout);
        }

        if ((env.step % env.output_every) == 0)
            save(&env);

        integrate(&env, env.step == env.Nsteps);
    }

    for (i=0; i < env.N/2; i++)
    {
        env.p[i].color[0] = 0; env.p[i].color[1] = 255; env.p[i].color[2] = 0;
    }
    for (i=env.N/2; i < env.N; i++)
    {
        env.p[i].color[0] = 255;   env.p[i].color[1] = 255; env.p[i].color[2] = 255;
    }

    save(&env);

    finalize(&env);

    //--------------------------------------------------------------------------
    // Check that our final conditions are the same as the initial ones.
    //--------------------------------------------------------------------------
    size_t n_bad = 0;
    for (i=0; i < env.N; i++)
    {
        if ((env.p[i].x[0] != env.p0[i].x[0])
        ||  (env.p[i].x[1] != env.p0[i].x[1])
        ||  (env.p[i].x[2] != env.p0[i].x[2])
        ||  (env.p[i].v[0] != -env.p0[i].v[0])
        ||  (env.p[i].v[1] != -env.p0[i].v[1])
        ||  (env.p[i].v[2] != -env.p0[i].v[2]))
        {
            n_bad++;
//          printf("IC[%ld] x(%ld %ld %ld) ", i, env.p0[i].x[0], env.p0[i].x[1], env.p0[i].x[2]);
//          printf("v(%ld %ld %ld)\n", env.p0[i].v[0], env.p0[i].v[1], env.p0[i].v[2]);
//          printf("  [%ld] x(%ld %ld %ld) ", i, env.p[i].x[0], env.p[i].x[1], env.p[i].x[2]);
//          printf("v(%ld %ld %ld)\n", env.p[i].v[0], env.p[i].v[1], env.p[i].v[2]);
        }
    }
    
    if (n_bad)
        printf("%ld/%ld particles are not in the right place.\n", n_bad, env.N);
    else
        printf("PERFECT REVERSAL!\n");

    fclose(fpE);

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

