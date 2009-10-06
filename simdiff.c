#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "chaodyn.h"

void usage()
{
    eprintf("usage: simdiff FILES...\n");
    exit(2);
}

int main(int argc, char **argv)
{
    FILE *fp0, *fp1;
    header_t h0;
    header_t h1;

    size_t i;


    if (argc != 3) usage();

    fp0 = fopen(argv[1], "r"); if (fp0 == NULL) { eprintf("%s: %s\n", argv[1], strerror(errno)); exit(EXIT_FAILURE); }
    fp1 = fopen(argv[2], "r"); if (fp1 == NULL) { eprintf("%s: %s\n", argv[2], strerror(errno)); exit(EXIT_FAILURE); }

    fread(&h0, sizeof(h0), 1, fp0);
    fread(&h1, sizeof(h1), 1, fp1);

#define CHECK(c, msg, ...) do { \
    if (!(c)) { eprintf(msg, ##__VA_ARGS__); exit(EXIT_FAILURE); } \
} while (0)

#define CHKPROP(s0,s1,p, t) do { \
    if ((((s0).p) != ((s1).p))) { eprintf("%-10s: "t" "t"\n", #p, ((s0).p), ((s1).p)); } \
} while (0)

    CHECK(h0.with_integers == h1.with_integers,
        "Data type mismatch. One uses integers the other uses floats.\n");

    CHECK(h0.with_integers == WITH_INTEGERS,
        "Data type mismatch. This version of simdiff was compiled with WITH_INTEGERS=%i.\n",
        WITH_INTEGERS);

    CHECK(h0.N == h1.N, "Number of particles mismatch %"PRIu64" != %"PRIu64"\n", h0.N, h1.N);

    CHKPROP(h0, h1, M,        "%e");
    CHKPROP(h0, h1, eps,      "%e");
    CHKPROP(h0, h1, dt,       "%ld");
    CHKPROP(h0, h1, Nclasses, "%ld");
    CHKPROP(h0, h1, step,     "%ld");
    CHKPROP(h0, h1, radius,   "%ld");
    CHKPROP(h0, h1, seed,     "%ld");

    double diff0 = 0;
    double diff1 = 0;
    double diff2 = 0;
    double diff3 = 0;
#if 1
    double dv0[3] = {0,0,0};
    double dv1[3] = {0,0,0};
    double dv2[3] = {0,0,0};
#endif
    
    particle_t p0, p1;
    size_t nbad = 0;
    for (i=0; i < h0.N; i++)
    {
        fread(&p0, sizeof(p0), 1, fp0);
        fread(&p1, sizeof(p1), 1, fp1);

#if 1
        dv0[0] += pow(p0.x[0] - p1.x[0], 2);
        dv0[1] += pow(p0.x[1] - p1.x[1], 2);
        dv0[2] += pow(p0.x[2] - p1.x[2], 2);
 
        dv1[0] += pow(p0.v[0] - p1.v[0], 2);
        dv1[1] += pow(p0.v[1] - p1.v[1], 2);
        dv1[2] += pow(p0.v[2] - p1.v[2], 2);
#else

#if 0
        dv0[0] += pow(p0.x[0],2);
        dv0[1] += pow(p0.x[1],2);
        dv0[2] += pow(p0.x[2],2);

        dv1[0] += pow(p0.v[0],2);
        dv1[1] += pow(p0.v[1],2);
        dv1[2] += pow(p0.v[2],2);
#endif

        diff0 += pow(p0.x[0] - p1.x[0],2) + pow(p0.x[1] - p1.x[1],2) + pow(p0.x[2] - p1.x[2],2);
        diff1 += pow(p0.v[0] - p1.v[0],2) + pow(p0.v[1] - p1.v[1],2) + pow(p0.v[2] - p1.v[2],2);

        if ((p0.x[0] != p1.x[0]) ||  (p0.x[1] != p1.x[1]) ||  (p0.x[2] != p1.x[2])
//      ||  (p0.v[0] != p1.v[0]) ||  (p0.v[1] != p1.v[1]) ||  (p0.v[2] != p1.v[2])
        ||  (p0.v[0] != -p1.v[0]) ||  (p0.v[1] != -p1.v[1]) ||  (p0.v[2] != -p1.v[2])
            )
        {
            nbad++;
#if 0
            printf("Particle %"PRIu64":\n"
                   "    x("POST" "POST" "POST")\n    v("VELT" "VELT" "VELT")\n"
                   "    x("POST" "POST" "POST")\n    v("VELT" "VELT" "VELT")\n",
                   i,
                   p0.x[0], p0.x[1], p0.x[2], p0.v[0], p0.v[1], p0.v[2],
                   p1.x[0], p1.x[1], p1.x[2], p1.v[0], p1.v[1], p1.v[2]);
#endif
        }
    }

    diff0 = sqrt(diff0);
    diff1 = sqrt(diff1);

    double r = sqrt(pow(dv0[0],2) + pow(dv0[1],2) + pow(dv0[2],2)),
           t = acos(((double)dv0[2]) / r),
           p = atan2(dv0[1], dv0[0]);

    r = log10(r);

    double x = r * sin(t) * cos(p),
           y = r * sin(t) * sin(p),
           z = r * cos(t);

    r = sqrt(pow(dv1[0],2) + pow(dv1[1],2) + pow(dv1[2],2)),
    t = acos(dv1[2] / r),
    p = atan2(dv1[1], dv1[0]);

    r = log10(r);

    double vx = r * sin(t) * cos(p),
           vy = r * sin(t) * sin(p),
           vz = r * cos(t);

    printf("%ld %ld %e %e %e %e %e %e %e %e\n", h0.step, h1.step, diff0, diff1, x,y,z, vx,vy,vz);
    //printf("%ld %ld %e %e %e %e %e %e %e %e %e %e\n", h0.step, h1.step, diff0, diff1, diff2, diff3, dv0[0], dv0[1], dv0[2], dv1[0], dv1[1], dv1[2]);

    if (nbad == 0)
        eprintf("PERFECT PHASE-SPACE MATCH!\n");
    else
        eprintf("%ld/%ld particles are not in the right place.\n", nbad, h0.N);

    fclose(fp0);
    fclose(fp1);

    return EXIT_SUCCESS;
}

