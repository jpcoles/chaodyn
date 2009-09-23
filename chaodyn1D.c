#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <png.h>

typedef int64_t pos_t;              /* pc               */
typedef int64_t vel_t;              /* km/s             */
typedef float   acc_t;              /* km/s^2           */
typedef float   rho_t;              /* Msun/pc          */
typedef float  tyme_t;              /* Myr              */
typedef unsigned char class_t;

//#define G 4.3e-3                    /* pc/Msun(km/s)^2  */
//#define G 4.3e-3                    /* pc/Msun(km/s)^2  */


#define Myr  (1e6 * 365 * 24 * 60 * 60)   // [s]
#define kpc  (3.08568025e19)              // [m]
#define Msun (1.98892e30)                 // [kg]

#define G (6.67300e-11)                                     // [m^3 kg^-1 s^-2]
#define L (1e-10 * kpc)
#define T (1     * Myr)
#define M (1e-10 * Msun)

const double _2piG = 2 * M_PI * G / (L*L*L) * M * (T*T);

typedef struct
{
    pos_t x[1];
    vel_t v[1];
    class_t class;
} particle_t;

typedef struct
{
    uint32_t nr;
    uint32_t nc;
    uint64_t size;
    unsigned char *image;
} image_t;

typedef struct
{
    size_t   N;
    particle_t *p0;
    particle_t *p;
    acc_t   *a;
    rho_t    rho;
    tyme_t   dt;
    tyme_t   t;
    tyme_t   start_time;
    tyme_t   end_time;

    pos_t    maxR;
    vel_t    maxV;

    size_t step;
    size_t output_every;

    size_t Nclasses;
    char **class_color;

    image_t snapshot;
    image_t phase_space_diagram;

} env_t;

int p_compare(const void *a0, const void *b0)
{
    const particle_t *a = (const particle_t *)a0;
    const particle_t *b = (const particle_t *)b0;

    assert(a->x[0] != b->x[0]);
//  if (a->x[0] < b->x[0]) return -1;
//  if (a->x[0] > b->x[0]) return +1;
//  assert(0);
    return 2*(a->x[0] > b->x[0]) - 1;
    //return a->x[0] - b->x[0];
}

void sort(env_t *env)
{
    qsort(env->p, env->N, sizeof(*(env->p)), p_compare);
}

void integrate(env_t *env)
{
    size_t i; 
    const size_t  N   = env->N;
    const rho_t   rho = env->rho;
    const tyme_t  dt  = env->dt;
    acc_t      *a = env->a;
    particle_t *p = env->p;

    //printf("%ld %ld %f  %ld %ld %f\n", p[0].x[0], p[0].v[0], a[0], p[1].x[0], p[1].v[0], a[1]);

    /* Drift */ for (i=0; i < N; i++) p[i].x[0] += (pos_t)(p[i].v[0] * dt/2);
    
    /* Force */ 
    sort(env);
    for (i=0; i < N; i++)
    {
        const int N_on_the_right = N - i - 1;
        const int N_on_the_left  = i;
        a[i] = _2piG * rho * (N_on_the_right - N_on_the_left);
    }

    assert(a[0] > 0);
    assert(a[N-1] < 0);

    /* Kick  */ for (i=0; i < N; i++) p[i].v[0] += (vel_t)(a[i]   * dt);
    /* Drift */ for (i=0; i < N; i++) p[i].x[0] += (pos_t)(p[i].v[0] * dt/2);

    //printf("%ld %ld %f  %ld %ld %f\n", p[0].x[0], p[0].v[0], a[0], p[1].x[0], p[1].v[0], a[1]);
}

void ic_random(env_t *env)
{
    size_t i;
    size_t N = env->N;
    for (i=0; i < N; i++)
    {
        env->p[i].x[0] = (pos_t)(env->maxR * (2*drand48()-1));
        env->p[i].v[0] = 0;
    }
}

void ic_uniform_plus_jitter(env_t *env)
{
    size_t i;
    size_t N = env->N;
    pos_t R = env->maxR / 2;
    float  dx = 2.0*R / (N-1);
    //printf("%f %ld\n", dx, N);
    for (i=0; i < N; i++)
    {
        env->p[i].x[0] = (pos_t)(i*dx - R) + (pos_t)(0.2*dx*(drand48()-0.5));
        //env->p[i].x[0] = (pos_t)(i*dx - env->maxR) + (pos_t)(0.5*dx*(drand48()-0.5));
        //env->s[i].x = (2*i*dx + (pos_t)(dx*(drand48()-0.5))) - env->maxR;
        //printf("%i\n", env->s[i].x);
        env->p[i].v[0] = 0;
        env->p[i].class = i / (N/4);
        //env->p[i].class = i / (N/31);
    }
}

void ic_cos(env_t *env)
{
    size_t i;
    size_t N = env->N;
    pos_t R = env->maxR / 2;
    float l = 4*R;
    for (i=0; i < N/2; i++)
    {
        env->p[i].x[0] = -l * (1-cos(i*M_PI/2/N/2)) - 10001;
        env->p[i].v[0] = 0;
        env->p[i].class = i / (N/4);
    }

    for (i=N/2; i < N; i++)
    {
        env->p[i].x[0] = l * (1-cos((i-N/2)*M_PI/2/N/2)) + 10001;
        env->p[i].v[0] = 0;
        env->p[i].class = i / (N/4);
    }
}

void generate_ics(env_t *env)
{
    //ic_random(env);
    //ic_uniform_plus_jitter(env);
    ic_cos(env);
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

void capture(env_t *env)
{
    size_t i;
    uint32_t r = env->step / env->output_every;

    uint32_t nc = env->snapshot.nc;
    uint32_t nr = env->snapshot.nr;

    for (i=0; i < env->N; i++)
    {
        int32_t c = (env->p[i].x[0] + env->maxR) / (2.0*env->maxR) * nc;
        if (!(0 <= c && c < env->snapshot.nc)) continue;
        env->snapshot.image[3*(r*nc + c) + 0] = env->class_color[env->p[i].class][0];
        env->snapshot.image[3*(r*nc + c) + 1] = env->class_color[env->p[i].class][1];
        env->snapshot.image[3*(r*nc + c) + 2] = env->class_color[env->p[i].class][2];
    }
}

void capture_phase_space(env_t *env, int clear)
{
    size_t i;
    uint32_t nc = env->phase_space_diagram.nc;
    uint32_t nr = env->phase_space_diagram.nr;

    if (clear)
        memset(env->phase_space_diagram.image, 0, env->phase_space_diagram.size);

    for (i=0; i < env->N; i++)
    {
        int32_t r = (env->p[i].v[0] + env->maxV) / (2.0*env->maxV) * nr;
        int32_t c = (env->p[i].x[0] + env->maxR) / (2.0*env->maxR) * nc;
        if (!(0 <= r && r < nr)) continue;
        if (!(0 <= c && c < nc)) continue;
        env->phase_space_diagram.image[3*(r*nc + c) + 0] = env->class_color[env->p[i].class][0];
        env->phase_space_diagram.image[3*(r*nc + c) + 1] = env->class_color[env->p[i].class][1];
        env->phase_space_diagram.image[3*(r*nc + c) + 2] = env->class_color[env->p[i].class][2];
    }
}

void save_sim(env_t *env)
{
    save_image_png("chaodyn", env->snapshot.image, env->snapshot.nr, env->snapshot.nc);
}

void save_phase_space_diagram(env_t *env)
{
    char *fname = malloc(3+1+5+1);
    sprintf(fname, "chaodyn.%05ld", env->step);
    save_image_png(fname, env->phase_space_diagram.image, 
                          env->phase_space_diagram.nr, 
                          env->phase_space_diagram.nc);
    free(fname);
}

int main(int argc, char **argv)
{
    size_t i;
    env_t env;

    env.N            = 500;
    env.p            = malloc(env.N * sizeof(*env.p));
    env.a            = malloc(env.N * sizeof(*env.a));
    env.rho          = 5e5 * Msun/pow(kpc,2) / (M/pow(L,2)) / env.N; //1e2;
    env.dt           = 2;

    // Smaller values of maxR leads to few clumps forming
    env.maxR         = 2*kpc / L; //1L << 32;

    env.maxV         = .008*(kpc/Myr) / (L/T); //1L << 26;
    env.start_time   = 0;
    env.end_time     = 4000;
    env.t            = env.start_time;
    env.output_every = 4;

    env.Nclasses       = 32;
    env.class_color    = malloc(env.Nclasses * sizeof(*env.class_color));

    srand48(0);

    for (i=0; i < env.Nclasses; i++)
    {
        env.class_color[i]    = malloc(3 * sizeof(**env.class_color));
        env.class_color[i][0] = 255;
        env.class_color[i][1] = 255;
        env.class_color[i][2] = 255;
    }

#define CC(i,r,g,b) do { \
    env.class_color[i][0] = r; \
    env.class_color[i][1] = g; \
    env.class_color[i][2] = b; \
} while (0)

    CC( 0, 255,   0,   0);
    CC( 1, 255, 255,   0);
    CC( 2,   0, 255,   0);
    CC( 3, 255, 255, 255);
    CC( 4,   0, 255,   0);
    CC( 5,   0, 255, 255);
    CC( 6,   0,   0, 255);
    CC( 7, 192,   0,   0);
    CC( 8, 255, 192,   0);
    CC( 9, 255,   0, 192);
    CC(10, 255, 192, 192);
    CC(11,   0, 192,   0);
    CC(12,   0, 192, 192);
    CC(13,   0,   0, 192);
    CC(14, 128,   0,   0);
    CC(15, 192, 128,   0);
    CC(16, 192,   0, 128);
    CC(17, 192, 128, 128);
    CC(18,   0, 128,   0);
    CC(19,   0, 128, 128);
    CC(20,   0,   0, 128);
    CC(21,  96,   0,   0);
    CC(22, 128,  96,   0);
    CC(23, 128,   0,  96);
    CC(24, 128,  96,  96);
    CC(25,   0,  96,   0);
    CC(26,   0,  96,  96);
    CC(27,   0,   0,  96);
    CC(28, 96,  64,  64);
    CC(29,   0,  64,   0);
    CC(30,   0,  64,  64);
    CC(31,   0,   0,  64);

    env.snapshot.nr    = ceil(abs(env.end_time - env.start_time)) / env.output_every;
    env.snapshot.nc    = 1024;
    env.snapshot.size  = 3 * env.snapshot.nr * env.snapshot.nc;
    env.snapshot.image = malloc(env.snapshot.size); memset(env.snapshot.image,0,env.snapshot.size);

    env.phase_space_diagram.nr    = 1024;
    env.phase_space_diagram.nc    = 1024;
    env.phase_space_diagram.size  = 3 * env.phase_space_diagram.nr * env.phase_space_diagram.nc;
    env.phase_space_diagram.image = malloc(env.phase_space_diagram.size); 

    generate_ics(&env);

    for (env.step  = 1,
         env.t     = env.start_time;
         env.t    <= env.end_time;
         env.t++,
         env.step++)
    {
        if ((env.step % env.output_every) == 0) capture(&env);
#if 1
        //if ((env.step % env.output_every) == 0) 
        if ((env.step % 10) == 0) 
        {
            capture_phase_space(&env, 1);
            save_phase_space_diagram(&env);
        }
#endif

#if 1
        if (env.step == 2000)
        {
#if 1
            for (i=0; i < env.N; i++) 
            {
#if 1
                if (i == 0) continue;
//              if (i==0) env.p[i].v[0] *= 1+1e-3;
//              env.p[i].v[0] *= 1 + (drand48() - 0.5) * 0.1;
#endif
                env.p[i].v[0] *= -1;
            }

//          env.N--;
#endif
        }
#endif

        printf("Step %ld\n", env.step);
        integrate(&env);
    }

    save_sim(&env);

    return 0;
}

