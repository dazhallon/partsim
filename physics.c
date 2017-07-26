#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "physics.h"


#ifndef sqr
#define sqr(a) ((a)*(a))
#endif

#ifndef sign
#define sign(a) ((a) > 0 ? 1 : -1)
#endif

void push_lst(plst_t** head, pcord_t new_part) {
  plst_t* new;
  new = (plst_t*) malloc(sizeof(plst_t));
  new->val = new_part;
  new->next = *head;

  *head = new;
}

void multi_push(plst_t** head, pcord_t new_part[],
		int n_to_append){
  int i;

  for(i = 0; 0 < n_to_append; ++i){
    push_lst(head, new_part[i]);
  }
}

void remove_from_list(plst_t** head, int n){
  plst_t *current = *head, *temp = NULL;
  int i;

  // iterate to the particle before the one to remove
  for (i = 0; i < n-1; ++i){
    current = current->next;
  }

  // temporally save the particle to remove 
  temp = current->next;
  current->next = temp->next;

  free(temp);
}

plst_t* get_part(plst_t* head, int n){
  plst_t *current = head;
  int i;

  
  for (i = 0; i<n; ++i){
    if (current->next != NULL)
      current = current->next;
    else{
      fprintf(stderr, "The function get_part had a to large int n = %d compared to the size of plst_t* head\n",
	      n);
      exit(1);
    }
  }
  
  return current;
}


void alter_part(plst_t* head, int where, pcord_t part){
  plst_t* current = NULL;
  current = (plst_t*) malloc(sizeof(plst_t));
  current = head;
  int i;
  for (i = 0; i< where; ++i){
    if (current->next != NULL)
      current = current->next;
    else{
      fprintf(stderr, "Use to large index in a plst_t *\n");
      exit(1);
    }
  }
  current->val = part;
}
void convert_to_plst_t(pcord_t from[], int n_particles,
		       plst_t* to){
  int i;

  for (i = 0; i < n_particles; ++i){
    push_lst(&to, from[i]);
  }
}

int feuler(plst_t *a, float time_step){
    a->val.x = a->val.x + time_step* a->val.vx ;
    a->val.y = a->val.y + time_step* a->val.vy ;	    
    return 0 ;
}

float wall_collide(plst_t *p, cord_t wall){
    float gPreassure = 0.0 ;
    
    if(p->val.x < wall.x0){
	p->val.vx = -p->val.vx ;
	p->val.x  = wall.x0 + (wall.x0-p->val.x);
	gPreassure += 2.0*fabs(p->val.vx) ;
    }
    if(p->val.x > wall.x1){
	p->val.vx = -p->val.vx ;
	p->val.x  = wall.x1 - (p->val.x-wall.x1);
	gPreassure += 2.0*fabs(p->val.vx) ;
    }
    if(p->val.y < wall.y0){
	p->val.vy = -p->val.vy ;
	p->val.y  = wall.y0 + (wall.y0-p->val.y);	
	gPreassure += 2.0*fabs(p->val.vy) ;
    }
    if(p->val.y > wall.y1){
	p->val.vy = -p->val.vy ;
	p->val.y  = wall.y1 - (p->val.y-wall.y1);
	gPreassure += 2.0*fabs(p->val.vy);
    }
    return gPreassure ;
}

//Since all walls would not add to the pressure in a box: we use this function.
float mpi_wall_collide(pcord_t *p, cord_t wall, int rank, int num_p){
    float gPreassure = 0.0;
    
    if(p->x < wall.x0 && rank == 0){
	p->vx = -p->vx ;
	p->x  = wall.x0 + (wall.x0-p->x);
	gPreassure += 2.0*fabs(p->vx) ;
    }
    if(p->x > wall.x1 && rank == num_p-1){
	p->vx = -p->vx ;
	p->x  = wall.x1 - (p->x-wall.x1);
	gPreassure += 2.0*fabs(p->vx) ;
    }
    if(p->y < wall.y0){
	p->vy = -p->vy ;
	p->y  = wall.y0 + (wall.y0-p->y);	
	gPreassure += 2.0*fabs(p->vy) ;
    }
    if(p->y > wall.y1){
	p->vy = -p->vy ;
	p->y  = wall.y1 - (p->y-wall.y1);
	gPreassure += 2.0*fabs(p->vy) ;
    }
    
    return gPreassure ;
}



float collide(pcord_t *p1, pcord_t *p2){
    double a,b,c;
    double temp,t1,t2;

    a=sqr(p1->vx-p2->vx)+sqr(p1->vy-p2->vy);
    b=2*((p1->x - p2->x)*(p1->vx - p2->vx)+(p1->y - p2->y)*(p1->vy - p2->vy));
    c=sqr(p1->x-p2->x)+sqr(p1->y-p2->y)-4*1*1;
    
    if (a!=0.0){	
	temp=sqr(b)-4*a*c;
	if (temp>=0){
	    temp=sqrt(temp);
	    t1=(-b+temp)/(2*a);
	    t2=(-b-temp)/(2*a);
      
	    if (t1>t2){
		temp=t1;
		t1=t2;
		t2=temp;
	    }
	    
	    if ((t1>=0)&(t1<=1))
		return t1;
	    else if ((t2>=0)&(t2<=1))
		return t2;    
	}
    }	
    return -1;
}



void interact(plst_t *p1,plst_t *p2, float t){
    float c,s,a,b,tao;
    pcord_t p1temp,p2temp;
  
    if (t>=0){

	/* Move to impact point */
	(void)feuler(p1,t);
	(void)feuler(p2,t);
    
	/* Rotate the coordinate system around p1*/
	p2temp.x=p2->val.x-p1->val.x;
	p2temp.y=p2->val.y-p1->val.y;
    
	/* Givens plane rotation, Golub, van Loan p. 216 */
	a=p2temp.x;
	b=p2temp.y;
	if (p2->val.y==0){
	    c=1;s=0;
	}
	else{
	    if (fabs(b)>fabs(a)){
		tao=-a/b;
		s=1/(sqrt(1+sqr(tao)));
		c=s*tao;
	    }
	    else{
		tao=-b/a;
		c=1/(sqrt(1+sqr(tao)));
		s=c*tao;
	    }
	}
    
	p2temp.x=c * p2temp.x+s * p2temp.y; /* This should be equal to 2r */
	p2temp.y=0.0;
    
	p2temp.vx= c* p2->val.vx + s* p2->val.vy;
	p2temp.vy=-s* p2->val.vx + c* p2->val.vy;
	p1temp.vx= c* p1->val.vx + s* p1->val.vy;
	p1temp.vy=-s* p1->val.vx + c* p1->val.vy;
    
	/* Assume the balls has the same mass... */
	p1temp.vx=-p1temp.vx;
	p2temp.vx=-p2temp.vx;
    
	p1->val.vx = c * p1temp.vx - s * p1temp.vy;
	p1->val.vy = s * p1temp.vx + c * p1temp.vy;
	p2->val.vx = c * p2temp.vx - s * p2temp.vy;
	p2->val.vy = s * p2temp.vx + c * p2temp.vy;

	/* Move the balls the remaining time. */
	c=1.0-t;
	(void)feuler(p1,c);
	(void)feuler(p2,c);
    }

}


void update_locparticels(pcord_t *p, pcord_t *recv_left, 
			 int left_tot_count,  pcord_t *recv_right,
			 int right_tot_count, int* totsize){
  int i, right_count, left_count;
  // Fill in the empty spaces in p
  //printf(stderr,"Updatating locparticles\n");
  //fprintf(stderr, "left_tot_count = %d right_tot_count = %d\n", left_tot_count,
  //	  right_tot_count);
  //printf(stderr,"scount[rank] = %d\n", *totsize);
  if(*totsize < 0 || left_tot_count > 500 || right_tot_count > 500 ||
     *totsize > 500 || left_tot_count < 0 || right_tot_count < 0)
    exit(2);

  
  for (i=right_count=left_count=0; i < *totsize; ++i){
    if(right_count < right_tot_count){
      p[i] = recv_right[right_count];
      ++right_count;
    }
    else if( left_count < left_tot_count) {
      p[i] = recv_left[left_count];
      ++left_count;
    }
  }
  
  // if there are more particles that have come in than gone out
  if (i < *totsize-1) { 
    i = *totsize-1;
    while(1) {
      if (right_count < right_tot_count){
	p[i] = recv_right[right_count];
	++i;
	++right_count;
      }
      else if (left_count < left_tot_count){
	p[i] = recv_left[left_count];
	++i;
	++left_count;
      }
      else {
	break;
      }
    }
  }

  // find the totsize
  //  *totsize = sizeof(p)/sizeof(p[0]); // will be zero.
  *totsize += left_tot_count + right_tot_count;
  //rintf(stderr,"scount[rank] = %d\n", *totsize);
  if(*totsize < 0 || left_tot_count > 500 || right_tot_count > 500 ||
     *totsize > 500 || left_tot_count < 0 || right_tot_count < 0)
    exit(3);
}
