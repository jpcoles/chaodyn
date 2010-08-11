#define _GNU_SOURCE
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "chaodyn.h"
#include "io.h"

void usage()
{
    eprintf("usage: sim2png [--geom=WxH] SIMFILES...\n");
    exit(2);
}

int main(int argc, char **argv)
{
    uint32_t width  = 1024;
    uint32_t height = 1024;

    char *select_fname = NULL;
    size_t *select_ids = NULL;
    size_t nids = 0;
    char color_split = 0;

    uint32_t xaxis=0, yaxis=1;

    FILE *select_fp = NULL;

    char *fname = NULL;

    char class_color[2][3] = { {255,0,0}, {0,255,0} };

    env_t *env = calloc(1, sizeof(env_t));

    if (argc < 2) usage();

    while(1)
    {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"geom",            1, 0, 0},
            {"select",          1, 0, 0},
            {"xy",              0, 0, 0},
            {"yz",              0, 0, 0},
            {"xz",              0, 0, 0},
            {0, 0, 0, 0}
        };

        int c = getopt_long(argc, argv, "X:o:", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) 
        {
            #define OPTSTR(s) (!strcmp(s, long_options[option_index].name))
            case 0:
                if OPTSTR("geom")
                {
                    char *endptr;
                    width = strtoll(optarg, &endptr, 10);
                    if (*endptr != 'x' || *endptr != 'X') usage();
                    endptr++;
                    if (*endptr == 0) usage();
                    height = strtoll(endptr, &endptr, 10);
                    if (*endptr != 0) usage();
                }
                else if OPTSTR("select")
                {
                    //select_fname = malloc(strlen(optarg)+1);
                    //strcpy(select_fname, optarg);
                    select_fname = optarg;
                }
                else if OPTSTR("xy") { xaxis=0; yaxis=1; break; }
                else if OPTSTR("yz") { xaxis=1; yaxis=2; break; }
                else if OPTSTR("xz") { xaxis=0; yaxis=2; break; }
                break;

            case 'o':
                fname = optarg;
                break;

            case 'X':
                if (!strcmp("color-split", optarg))
                    color_split = 1;
                else
                    usage();
                break;

            default:
                usage();
                break;
        }
    }

    if (fname == NULL)
        fname = "sim2png.png";

    if (select_fname)
    {
        if (!strcmp("-", select_fname)) 
            select_fp = stdin;
        else
            select_fp = fopen(select_fname, "r"); assert(select_fp != NULL);

        nids = 0;
        size_t allocated_ids = 64;

        select_ids = malloc(allocated_ids * sizeof(*select_ids));

        while (!feof(select_fp))
        {
            size_t id;
            if (fscanf(select_fp, "%ld", &id) == 0) break;
            eprintf("%ld\n", id);
            if (nids == allocated_ids)
            {
                allocated_ids *= 2;
                select_ids = realloc(select_ids, allocated_ids);
            }

            select_ids[nids] = id;
            nids++;
        }

        if (select_fp != stdin)
            fclose(select_fp);
    }

    unsigned char *image = calloc(3 * width * height, sizeof(unsigned char));

    while (optind < argc)
    {
        eprintf("%s\n", argv[optind]);
        load(env, argv[optind++]);

        env->radius *= 2;

        size_t i;
//      for (i=0; i < env->N; i++)
//          env->p[i].class = i / (env->N/2);

        #define PLOT(i) do { \
            int32_t c = ( env->p[i].x[xaxis] + env->radius) / (2*env->radius / width);\
            int32_t r = (-env->p[i].x[yaxis] + env->radius) / (2*env->radius / height);\
            if (!(0 <= r && r < height)) break; \
            if (!(0 <= c && c < width))  break; \
            image[3*(r*width + c) + 0] = class_color[env->p[i].class][0]; \
            image[3*(r*width + c) + 1] = class_color[env->p[i].class][1]; \
            image[3*(r*width + c) + 2] = class_color[env->p[i].class][2]; \
        } while (0)
            //image[3*(r*width + c) + 0] = 0;\
            //image[3*(r*width + c) + 1] = 255;\
            //image[3*(r*width + c) + 2] = 0;\

        if (select_ids)
        {
            for (i=0; i < nids; i++)
                PLOT(select_ids[i]);
        }
        else
        {
            for (i=0; i < env->N; i++)
                PLOT(i);
        }

        free(env->p);
    }

    //save_image_png("sim2png.png", image, height, width);
    save_image_png(fname, image, height, width);

    free(env);

    return EXIT_SUCCESS;
}

