#ifndef _physics_h
#define _physics_h

#include "coordinate.h"

#define STEP_SIZE 1.0 /* the step size use in the integration */

void push_lst(plst_t* head, pcord_t new_part);

void multi_push(plst_t* head, pcord_t new_part[], int n_to_append);

void remove_from_list(plst_t* head, int n);

pcord_t get_part(plst_t* head, int n);

void convert_to_plst_t(pcord_t from[], int n_particles,
		       plst_t* to);

int feuler(pcord_t *a,
	   float time_step) ;

float wall_collide(pcord_t *p,
		   cord_t wall) ;

float mpi_wall_collide(pcord_t *p, cord_t wall,
		       int rank, int num_p);

void update_locparticels(pcord_t *p, pcord_t *recv_left,
			 int left_tot_count, pcord_t *recv_right,
			 int right_tot_count, int* totsize);

float collide(pcord_t *p1,
	      pcord_t *p2) ;

void interact(pcord_t *p1,
	      pcord_t *p2,
	      float t) ;


#endif
