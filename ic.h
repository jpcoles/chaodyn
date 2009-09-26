#ifndef IC_H
#define IC_H

#include "chaodyn.h"

struct iclist_s
{
    char *name;
    void (*f)(env_t *env);
};

struct iclist_s *iclist();

void ic_random(env_t *env);
void sphere(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1, class_t class);
void shell(env_t *env, pos_t x0, pos_t y0, pos_t z0, pos_t R, size_t N0, size_t N1, class_t class);
void ic_uniform_random_shell(env_t *env);
void ic_two_shells(env_t *env);
void ic_one_shell(env_t *env);
void ic_colliding_halos(env_t *env);
void ic_sphere(env_t *env);
void ic_two_points(env_t *env);
void ic_figure8(env_t *env);

#endif
