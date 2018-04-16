* Header containing the pulse-driven primitives of the MParCO framework 
* Author: Allberson/Correa

#ifndef GLO_PD_H
#define GLO_PD_H
#include <mpi.h>

#ifdef _OPENMPI_
#define PD_NULL NULL
#endif
#ifdef _MPICH_
#define PD_NULL 0
#endif

/*
 * Configuration constants: in each case, options DEFAULT and OFF are
 * the smallest and greatest values, respectively
 */
#define PD_TERMINATION                 0  // Termination mode
#define PD_TERMINATION_DEFAULT         0  // Default: PD_TERMINATION_OFF
#define PD_TERMINATION_DIFFUSING       1  // Diffusing computations: Dijkstra's algorithm
#define PD_TERMINATION_OFF             2  // No termination implemented

#define PD_PULSE_GENERATION            1  // Pulse generation mode
#define PD_PULSE_GENERATION_DEFAULT    0  // Default: PD_PULSE_GENERATION_OFF
#define PD_PULSE_GENERATION_AUTOMATIC  1  // Automatically generate pulses up to the first one that can receive a given message
#define PD_PULSE_GENERATION_EMPTY      2  // Generate empty pulses up to the first one that can receive a given message
#define PD_PULSE_GENERATION_POSTPONE   3  // Postpone a given message to be considered in the first pulse in the future that can receive it
#define PD_PULSE_GENERATION_OFF        4  // No pulses are generated when a message arrives too early

#define PD_MSG_ORDER                   2  // Message ordering mode
#define PD_MSG_ORDER_DEFAULT           0  // Default: PD_MSG_ORDER_OFF
#define PD_MSG_ORDER_FIFO              1
#define PD_MSG_ORDER_CAUSAL            2
#define PD_MSG_ORDER_OFF               3

#define PD_SMALL_MEM_BUF               1024
#define PD_MEDIUM_MEM_BUF              8192
#define PD_LARGE_MEM_BUF               32768

/*
 * Setup functions
 */

/*
 * Sets the event and pulse functions of the following computation
 */
int MPI_Pulse_set_model(int(*event)(MPI_Status *), int(*pulse)(long long));
// collective properties are only assured for the processes in the virtual topology and messages exchanged with
// pulse functions. The following function allows inspecting the virtual topology
int MPI_Pulse_neighbors_test(int dest);

/*
 * Sets a specified configuration parameter with a specified value.
 * Returns MPI_SUCCESS in case of success and error code otherwise
 */
int MPI_Pulse_set_config_param(int param, int value);
int MPI_Pulse_run(MPI_Comm comm, MPI_Status *comp_status);

/*
 * Communication functions callable from inside model functions only
 */

// callable from event

#define MPI_Pulse_Recv _PD_MPI_Recv_
extern int (*MPI_Pulse_Recv)( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );

// postpone message originating current event after it is received.
// The count and datatype (size of elements in buf) arguments should match the arguments provided to the receive call of the current event.
#define MPI_Pulse_postpone_msg _PD_MPI_postpone_msg_
extern int (*MPI_Pulse_postpone_msg)(long long targetpulse);

#define MPI_Pulse_resume_pulse_counter _PD_MPI_resume_pulse_counter_
extern int (*MPI_Pulse_resume_pulse_counter)();

// callable from pulse

#define MPI_Pulse_Send _PD_MPI_Send_
extern int (*MPI_Pulse_Send)( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
#define MPI_Pulse_Isend _PD_MPI_Isend_
extern int (*MPI_Pulse_Isend)( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );

#define MPI_Pulse_advance_pulse_counter _PD_MPI_advance_pulse_counter_
extern int (*MPI_Pulse_advance_pulse_counter)(int newpulse);

#define MPI_Pulse_suspend_pulse_counter _PD_MPI_suspend_pulse_counter_
extern int (*MPI_Pulse_suspend_pulse_counter)();

#endif
