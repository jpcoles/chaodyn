#ifndef IO_H
#define IO_H

#include "mss3domp.h"

void save_sim(env_t *env, char *fname);
void save_ic(env_t *env);
void save(env_t *env);
void load(env_t *env, char *fname);
void save_snapshot(env_t *env);
void save_path_image(env_t *env);
void save_comparison_image(env_t *env);
int save_image_png(char *fname, unsigned char *image, uint32_t nrows, uint32_t ncols);

#endif

