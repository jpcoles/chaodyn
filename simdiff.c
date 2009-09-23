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

    fp0 = fopen(argv[1], "r"); assert(fp0 != NULL);
    fp1 = fopen(argv[2], "r"); assert(fp1 != NULL);

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

    double diff[3] = {0,0,0};
    double diff2 = 0;
    
    particle_t p0, p1;
    size_t nbad = 0;
    for (i=0; i < h0.N; i++)
    {
        fread(&p0, sizeof(p0), 1, fp0);
        fread(&p1, sizeof(p1), 1, fp1);

#if 0
        diff[0] += pow(p0.x[0] - p1.x[0], 2);
        diff[1] += pow(p0.x[1] - p1.x[1], 2);
        diff[2] += pow(p0.x[2] - p1.x[2], 2);
#endif

#if 0
        diff[0] += abs(p0.x[0] - p1.x[0]);
        diff[1] += abs(p0.x[1] - p1.x[1]);
        diff[2] += abs(p0.x[2] - p1.x[2]);
#endif

#if 1
        diff2 += sqrt(pow(p0.x[0] - p1.x[0],2) + pow(p0.x[1] - p1.x[1], 2) + pow(p0.x[2] - p1.x[2],2));
#endif

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

    double d = sqrt(pow(diff[0],2) + pow(diff[1],2) + pow(diff[2],2)) / h0.N;
    printf("diff: %e %e\n", d, diff2);

    if (nbad == 0)
        eprintf("PERFECT PHASE-SPACE MATCH!\n");
    else
        eprintf("%ld/%ld particles are not in the right place.\n", nbad, h0.N);

    fclose(fp0);
    fclose(fp1);

    return EXIT_SUCCESS;
}

