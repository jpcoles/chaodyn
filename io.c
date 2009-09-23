#include <png.h>
#include <stdlib.h>
#include "mss3domp.h"
#include "io.h"

//==============================================================================
//                                  save_sim
//==============================================================================
void save_sim(env_t *env, char *fname)
{
    header_t header = 
    {
        "mss3domp v0.1   ",
        WITH_INTEGERS,
        env->N,
        env->M,
        env->eps,
        env->dt,
        env->opt.Nclasses,
        env->step,
        0,
        env->radius,
        0,
        env->seed,
        env->opt.img_rows,
        env->opt.img_cols,
        env->units
    };

    VL(1) LOG(CL "Writing to %s\n", fname);
    FILE *fp = fopen(fname, "w");
    fwrite(&header, sizeof(header), 1, fp);
    fwrite(env->p, sizeof(*(env->p)), env->N, fp);
    fclose(fp);
}

//==============================================================================
//                                  save_ic
//==============================================================================
void save_ic(env_t *env)
{
    char *fname = malloc(strlen(env->opt.tag)+1+2+1+3+1);
    sprintf(fname, "%s.ic.sim", env->opt.tag);
    save_sim(env, fname);
    free(fname);
}

//==============================================================================
//                                    save
//==============================================================================
void save(env_t *env)
{
    char *fname = malloc(strlen(env->opt.tag)+1+5+1+3+1);
    sprintf(fname, "%s.%05ld.sim", env->opt.tag, env->step);
    save_sim(env, fname);
    free(fname);
}

//==============================================================================
//                                    load
//==============================================================================
void load(env_t *env, char *fname)
{
    header_t h;
    
    VL(1) LOG(CL "Reading from %s\n", fname);

    FILE *fp = fopen(fname, "r"); assert(fp != NULL);
    fread(&h, sizeof(h), 1, fp);

    if (h.with_integers != WITH_INTEGERS)
    {
        if (WITH_INTEGERS)
            eprintf("This version of mss3domp was compiled with integers enabled, but\n"
                    "the given simulation input was generated using floats. The two\n"
                    "are not compatible.\n");
        else
            eprintf("This version of mss3domp was compiled with floats enabled, but\n"
                    "the given simulation input was generated using integers. The two\n"
                    "are not compatible.\n");
        exit(EXIT_FAILURE);
    }

    env->N                = h.N;
    env->M                = h.M;
    env->eps              = h.eps;
    env->dt               = h.dt;
    env->opt.Nclasses     = h.Nclasses;
    env->step             = h.step;
    env->radius           = h.radius;
    env->seed             = h.seed;
    env->opt.img_rows     = h.img_rows;
    env->opt.img_cols     = h.img_cols;
    env->units            = h.units;

    env->p = malloc(env->N * sizeof(*env->p)); assert(env->p != NULL);
    fread(env->p, sizeof(*(env->p)), env->N, fp);
    fclose(fp);

    switch (env->opt.dump_sim)
    {
        default: assert(0);
        case 0:  break;

        case 1:
            if (h.with_integers)
            {
                LOG("%s Compiled with integers\n", h.version);
                LOG("Number of Particles : %ld\n", h.N);
                LOG("Mass                : "I_MASST"\n", h.M);
                LOG("Softening           : "I_SOFTT"\n", h.eps);
                LOG("dt                  : "I_TIMET"\n", h.dt);
                LOG("Current Timestep    : %ld\n", h.step);
                LOG("Start Radius        : "I_POST"\n", h.radius);
                LOG("Random Number Seed  : %ld\n", h.seed);
                LOG("Image size          : %ix%i\n", h.img_rows, h.img_cols);
            }
            else
            {
                LOG("%s Compiled with floats\n", h.version);
                LOG("Number of Particles : %ld\n", h.N);
                LOG("Mass                : "F_MASST"\n", h.M);
                LOG("Softening           : "F_SOFTT"\n", h.eps);
                LOG("dt                  : "F_TIMET"\n", h.dt);
                LOG("Current Timestep    : %ld\n", h.step);
                LOG("Start Radius        : "F_POST"\n", h.radius);
                LOG("Random Number Seed  : %ld\n", h.seed);
                LOG("Image size          : %ix%i\n", h.img_rows, h.img_cols);
            }
            break;
        case 2:
            break;
    }

}

//==============================================================================
//                               save_snapshot
//==============================================================================
void save_snapshot(env_t *env)
{
    char *fname = malloc(strlen(env->opt.tag)+1+5+1+3+1);
    sprintf(fname, "%s.%05ld.png", env->opt.tag, env->step);
    VL(1) LOG(CL "Writing image to %s\n", fname);
    save_image_png(fname, env->snapshot.image, env->snapshot.nr, env->snapshot.nc);
    free(fname);
}

//==============================================================================
//                              save_path_image
//==============================================================================
void save_path_image(env_t *env)
{
    char *fname = malloc(strlen(env->opt.tag)+1+4+1+3+1);
    sprintf(fname, "%s.path.png", env->opt.tag);
    VL(1) LOG(CL "Writing image to %s\n", fname);
    save_image_png(fname, env->path.image, env->path.nr, env->path.nc);
    free(fname);
}

//==============================================================================
//                           save_comparison_image
//==============================================================================
void save_comparison_image(env_t *env)
{
    char *fname = malloc(strlen(env->opt.tag)+1+7+1+3+1);
    sprintf(fname, "%s.compare.png", env->opt.tag);
    VL(1) LOG(CL "Writing image to %s\n", fname);
    save_image_png(fname, env->snapshot.image, env->snapshot.nr, env->snapshot.nc);
    free(fname);
}

//==============================================================================
//                               save_image_png
//==============================================================================
int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols)
{
    int i;
    FILE *fp = fopen(fname, "wb");

    if (fp == NULL)
    {
        eprintf("Can't open %s\n", fname);
        return 1;
    }

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
