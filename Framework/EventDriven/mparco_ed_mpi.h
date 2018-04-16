* Header containing the event-driven primitives of the MParCO framework 
* Author: Allberson/Corrêa

#ifndef GLO_ED_H
#define GLO_ED_H
#include <mpi.h>

#ifdef _OPENMPI_
#define ED_NULL NULL
#endif
#ifdef _MPICH_
#define ED_NULL 0
#endif

/*
 * Configuration constants: in each case, options DEFAULT and OFF are
 * the smallest and greatest values, respectively
 */
#define ED_TERMINATION                 0
#define ED_TERMINATION_DEFAULT         0
#define ED_TERMINATION_DIFFUSING       1
#define ED_TERMINATION_OFF             2

#define ED_DIFFUSING_ACK               100
#define ED_DIFFUSING_TERM              200

#define ED_MSG_ORDER                   2
#define ED_MSG_ORDER_DEFAULT           0
#define ED_MSG_ORDER_FIFO              1
#define ED_MSG_ORDER_CAUSAL            2
#define ED_MSG_ORDER_OFF               3

#define ED_SMALL_MEM_BUF               1024
#define ED_MEDIUM_MEM_BUF              8192
#define ED_LARGE_MEM_BUF               32768

/*
 * Setup functions
 */

/*
 * Sets the event function of the following computation
 */
int MPI_Event_set_model(int(*event)(MPI_Status *), int(*spontaneous_event)(MPI_Status *));
// collective properties are only assured for the processes in the virtual topology and messages exchanged with
// event functions. The following function allows inspecting the virtual topology
int MPI_Event_neighbors_test(int dest);

/*
 * Sets a specified configuration parameter with a specified value.
 * Returns MPI_SUCCESS in case of success and error code otherwise
 */
int MPI_Event_set_config_param(int param, int value);
int MPI_Event_run(MPI_Comm comm, MPI_Status *comp_status);

/*
 * Communication functions callable from inside model functions only
 */

// callable from event and spontaneous event

#define MPI_Event_Send _ED_MPI_Send_
extern int (*MPI_Event_Send)( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
#define MPI_Event_Isend _ED_MPI_Isend_
extern int (*MPI_Event_Isend)( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );

// callable from event 

#define MPI_Event_Recv _ED_MPI_Recv_
extern int (*MPI_Event_Recv)( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );

#endif

