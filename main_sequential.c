#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

// to do:
/* fix the code when we check for collsiion and so forth, here we use the old way. */
//Feel free to change this program to facilitate parallelization.
void print_pcord(pcord_t* p){
  pcord_t ptemp;
  int i;
  for(i = 0; i < 100; ++i){
    ptemp = p[i];
    printf("p[%d]: x=%f y=%f vx=%f vy=%f\n",
	   i, ptemp.x, ptemp.y, ptemp.vx, ptemp.vy);
  }
}
void w_from_to_pcord(pcord_t to, pcord_t from){
  to.x = from.x; to.y = from.y;
  to.vx = from.vx; to.vy = from.vy;
}

float rand1(){
  return (float)( rand()/(float) RAND_MAX );
}

void init_collisions(bool *collisions, unsigned int max){
  unsigned int i;
  for(i=0;i<max;++i)
    collisions[i]=0;
}

bool in_local(pcord_t p, cord_t wall) // checks if particle p is in wall wall.
{
  return  (p.x > wall.x0 && p.x < wall.x1);
}

bool exited_left(pcord_t p, cord_t wall, int rank){
  if(rank == 0)
    return false;
  else 
    return p.x < wall.x0;
}

bool exited_right(pcord_t p, cord_t wall, int rank, int num_p){
  if(rank == num_p -1)
    return false;
  else
    return p.x > wall.x1;
}

int main(int argc, char** argv){
  
  unsigned int time_stamp = 0, time_max, rnk;
  int send_left_count, send_right_count;
  int recv_right_count, recv_left_count=0, num_flags= 0;
  float pressure = 0, sum = 0;

  //MPI stuff
  int rank, num_p;
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_p);
  MPI_Status status;
  //init_mpi_pcord();

  
  /* create a mpi type for pcord_t */
  const int nitems=4;
  int          blocklengths[4] = {1,1,1,1};
  MPI_Datatype ptypes[5] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
  MPI_Datatype MPI_PCORD;
  MPI_Aint     poffsets[5];

  poffsets[0] = offsetof(pcord_t, x);
  poffsets[1] = offsetof(pcord_t, y);
  poffsets[2] = offsetof(pcord_t, vx);
  poffsets[3] = offsetof(pcord_t, vy);
  //  poffsets[4] = offsetof(pcord_t, to_remove);
  
  MPI_Type_create_struct(nitems, blocklengths, poffsets, ptypes, &MPI_PCORD);
  MPI_Type_commit(&MPI_PCORD);

  /* create a type for cord_t */
  MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
  MPI_Datatype MPI_CORD;
  MPI_Aint offsets[4];

  offsets[0] = offsetof(cord_t, x0);
  offsets[1] = offsetof(cord_t, x1);
  offsets[2] = offsetof(cord_t, y0);
  offsets[3] = offsetof(cord_t, y1);

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_CORD);
  MPI_Type_commit(&MPI_CORD);

  		
  // Declare send arrays used in mpi functions
  int scounts[num_p], displs[num_p];
  //pcord_t *sendparticles, *locparticles;
  plst_t  *locparticles = NULL;
  pcord_t **sendparticles; // first dimension stans for the rank of the process, the other for index of particle.
  pcord_t *recvparticles;
  // parse arguments
  if(argc != 2) {
    fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
    fprintf(stderr, "For example: %s 10\n", argv[0]);
    exit(-1);
  }

  time_max = atoi(argv[1]);

  /* Initialize */
  // 1. set the walls
  cord_t wall;
  wall.y0 = wall.x0 = 0;
  wall.x1 = BOX_HORIZ_SIZE;
  wall.y1 = BOX_VERT_SIZE;

  cord_t locwall;
  cord_t walls[num_p]; //needed in other to distr. the particels

  // 2. allocate particle buffer and initialize the particles
  pcord_t *particles = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  bool *collisions=(bool*) malloc(INIT_NO_PARTICLES*sizeof(bool));
  sendparticles = (pcord_t**) malloc(num_p*INIT_NO_PARTICLES*sizeof(pcord_t));
  recvparticles = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  //sendparticles = (plst_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  // locparticles = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  locparticles = (plst_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));									
  pcord_t *send_left = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  pcord_t *send_right = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  pcord_t *recv_left = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));
  pcord_t *recv_right = (pcord_t*) malloc(INIT_NO_PARTICLES*sizeof(pcord_t));

  srand( time(NULL) + 1234 );

  
  // set up walls 
  walls[rank].y0 = 0;
  walls[rank].y1 = BOX_VERT_SIZE;
  walls[rank].x0 = rank*BOX_HORIZ_SIZE/num_p;
  if (rank == num_p-1)
    walls[rank].x1 = BOX_HORIZ_SIZE;
  else
    walls[rank].x1 = walls[rank].x0 + BOX_HORIZ_SIZE/num_p;
  
  for(rnk = 0; rnk < num_p; ++rnk){
    scounts[rnk] = displs[rnk] = 0;
  }
      
  float r, a;
  int i;
  // send init particles and designate to processes. 
  for(i=0; i<INIT_NO_PARTICLES/num_p; i++) {
    // initialize random position
    recvparticles[i].x = locwall.x0 + rand1()*BOX_HORIZ_SIZE/4;
    recvparticles[i].y = locwall.y0 + rand1()*BOX_VERT_SIZE/4;
	
    // initialize random velocity
    r = rand1()*MAX_INITIAL_VELOCITY;
    a = rand1()*2*PI;
    recvparticles[i].vx = r*cos(a);
    recvparticles[i].vy = r*sin(a);
    //    recvparticles[i].to_remove=0;

    push_lst(locparticles, recvparticles[i]);
	      
  }
  
  
  // For each proces, do:
  unsigned int p, pp;
  int pressure_count =0;

  pcord_t current_part;
  pcord_t sub_current_part;

  // convert recvparticles into locarticles, so thath we get a linked list to work with later
  // convert_to_plst_t(recvparticles, INIT_NO_PARTICLES/4, locparticles); // bug here

  
  /* Main loop */
  for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp
    init_collisions(collisions, scounts[rank] + num_flags);
    for(p=0; p<scounts[rank]; p++) { // for all particles
      current_part = get_part(locparticles, p);
      if(collisions[p])// || locparticles[p].to_remove != 0) 
	/* check for collisions */
	// If a local boundary is hit, check if collide with any particle
	// Check if the particle will remain or will go to an other local box
	for(pp=p+1; pp<scounts[rank] + num_flags; pp++) {
	  if(collisions[pp]) continue; //|| locparticles[p].to_remove != 0) continue;
	  sub_current_part = get_part(locparticles, p);
	  float t=collide(&current_part, &sub_current_part);
	  if(t!=-1){ // collision
	    collisions[p]=collisions[pp]=1;
	    interact(&current_part, &sub_current_part, t);
	    break; // only check collision of two particles
	  }
	}
    }

    // move particles that has not collided with another
    for(p=send_left_count=send_right_count=0,
	  recv_left_count=recv_right_count=0;
	p < scounts[rank] + num_flags; p++) {
      current_part = get_part(locparticles, p);
      if(!collisions[p]) { //&& locparticles[p].to_remove == 0){
	feuler(&current_part, 1);
	pressure += mpi_wall_collide(&current_part, locwall,
				     rank, num_p);
	++pressure_count;
	fprintf(stderr, "Pressure was added from rank %d\n", rank);
      }
      else if (exited_left(current_part, locwall, rank)) {// &&
	//      locparticles[p].to_remove==0){
	send_left[send_left_count] = current_part;
	//locparticles[p].to_remove = 1;
	++send_left_count;
	if(num_p == 1 || rank == 0)
	  exit(5);
	printf("sending %d to left\trank %d\n", send_left_count, rank);
      }
      else if (exited_right(current_part, locwall, rank, num_p)) {// &&
	  //locparticles[p].to_remove==0){
	  send_right[send_right_count] = current_part;
	  //locparticles[p].to_remove = 1;
	  ++send_right_count;
	  if(num_p == 1 || rank == num_p - 1)
	    exit(6);
	  printf("sending %d to right\trank %d\n", send_right_count, rank);
      }
    }

    
    /* Send scheme(for the tags used in messages): 
       0 = send to left, 1 = send to right */
    //Send particles to other boxes

    if (rank-1 >= 0) {
      MPI_Send(send_left, send_left_count, MPI_PCORD, rank-1, 0, comm);
          printf("rank %d sent to left\n", rank);
    }

    if (rank+1 < num_p) {
      MPI_Send(send_right, send_right_count, MPI_PCORD, rank+1, 1, comm);
        printf("rank %d sent to right\n", rank);
    }

    // check if something was sent. Using probe to find recv count.
    if(rank-1 >= 0){
      MPI_Probe(rank-1, 1, comm, &status);
      MPI_Get_count(&status, MPI_PCORD, &recv_left_count);
      //fprintf(stderr, "recv_left_count = %d after MPI_probe and MPI_Get_count\n", recv_left_count);
      MPI_Recv(recv_left, recv_left_count, MPI_CORD, rank-1, 1, comm,
	       &status);
    }

    if(rank+1 < num_p){
      MPI_Probe(rank+1, 0, comm, &status);
      MPI_Get_count(&status, MPI_PCORD, &recv_right_count);
      MPI_Recv(recv_right, recv_right_count, MPI_CORD, rank+1, 0, comm,
	       &status);
    }
    
    // update locparticels
    // think about what happens when we tranfer particles between boxes
    multi_push(locparticles, recv_right, recv_right_count);
    multi_push(locparticles, recv_left, recv_left_count); //prob. need a special case hanfle for when recv_right/left == NULL
  }
  
  // Gather pressure
  int tot_pressure_count;
  MPI_Barrier(comm);
  MPI_Reduce(&pressure, &sum, 1, MPI_FLOAT, MPI_SUM, 0, comm);
  MPI_Reduce(&pressure_count, &tot_pressure_count, 1, MPI_INT, MPI_SUM, 0, comm);
  if (rank == 0)
    fprintf(stderr, "Average pressure = %f pressure count = %d\n",
	    sum/(WALL_LENGTH*time_max), tot_pressure_count);

  // free stuff
  free(particles);
  free(collisions);
  free(sendparticles);
  free(locparticles);
  free(send_left);
  free(send_right);
  free(recv_left);
  free(recv_right);
  MPI_Type_free(&MPI_CORD);
  MPI_Type_free(&MPI_PCORD);
  
  MPI_Finalize();
  return 0;
}

