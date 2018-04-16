* Implementation of the event-driven primitives of the MParCO framework 
* Author: Allberson/Corrêa

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <heap.h>
#include <repository.h>
#include <mparco_ed_mpi.h>

/*
 * Types
 */

typedef struct P_msg {
	MPI_Status *    status;
	void *          buf;
	MPI_Datatype    datatype;
	long long       event;
} Postponed_msg;

typedef struct Ps_msg {
	void *          buf;
	int             count;
	MPI_Datatype    datatype;
	int             dest;
	int             tag;
	struct Ps_msg * next;
} Event_msg;

typedef struct M_h {
	long long    event;             // event
	int          datatypesz;
	int          count;
} Msg_header;

/*
 * Datatypes
 */
static MPI_Datatype MPI_MSG_HEADER    = ED_NULL;
static MPI_Datatype MPI_HEADED_PACKED = ED_NULL;

/*
 * Configuration parameters
 */

static int config_term      = ED_TERMINATION_DEFAULT;
static int config_msg_order = ED_MSG_ORDER_OFF;

/*
 * Topology
 */

static MPI_Comm topology_comm        = ED_NULL;
static int      topology_type        = MPI_UNDEFINED;
static int      topology_rank        = -1;
static int      topology_size        = -1;
static int      topology_neigh_size  = 0;
static int *    topology_neigh_array = ED_NULL;   // is sorted in ascending order

/*
 * Model functions
 */
static int (*model_event)(MPI_Status *);
static int (*init)();
/*
 * Computation
 */
static long long    comp_event;
static long long    comp_advance;
static MPI_Status * comp_event_msg;
static int          comp_postponed_msg;
static int          comp_msg_was_received;
static void *       comp_buf;
static MPI_Datatype comp_datatype;
static int          comp_is_running           = 0;
static int          comp_msg_header_sz;

/*
 * Termination Detection
 */
static int  termination_parent        =    -1;    
static int  termination_deficit       =     0;
static int  termination_term          =     0;
static int  termination_sent_term     =     0;
static int  termination_lterm         =     0;
static int *termination_ack_array     =     ED_NULL;

/*
 * Memory management
 */

// heap of postponed messages
static Postponed_msg ** postmsg;
static size_t           npostmsg;
static size_t           postmsgsz;
// repository of postponed messages
static void *           postmsgrepo;
// repository of event messages
static void *           postmsgrepo;
// repository of MPI_Status
static void *           statrepo;
// repository of small size messages
static void *           smallmembufrepo;
// repository of small size messages
static void *           mediummembufrepo;
// repository of small size messages
static void *           largemembufrepo;

/*
 * Communication functions
 */

// callable from event and spontaneous event

static int _ED_Send_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
static int _ED_Send_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
static int _ED_Send_with_event_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
int (*MPI_Event_Send)( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) = _ED_Send_no_impl_;

static int _ED_Isend_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );
static int _ED_Isend_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );
int (*MPI_Event_Isend)( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) = _ED_Isend_no_impl_;

// callable from event

static int _ED_Recv_no_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );
static int _ED_Recv_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );
int (*MPI_Event_Recv)( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) = _ED_Recv_no_impl_;

/*
 * Local functions
 */

static inline void * getMemoryBuffer(int l) {
	if (0 <= l && l <= ED_SMALL_MEM_BUF)
		return getFromRepo(smallmembufrepo);
	if (l <= ED_MEDIUM_MEM_BUF)
		return getFromRepo(mediummembufrepo);
	if (l <= ED_LARGE_MEM_BUF)
		return getFromRepo(largemembufrepo);
	return (void *) malloc(l);
}

static inline void releaseMemoryBuffer(void * m, int l) {
	if (0 <= l && l <= ED_SMALL_MEM_BUF)
		release(smallmembufrepo, m);
	else if (l <= ED_MEDIUM_MEM_BUF)
		release(mediummembufrepo, m);
	else if (l <= ED_LARGE_MEM_BUF)
		release(largemembufrepo, m);
	else
		free(m);
}

static inline MPI_Datatype create_msg_header_datatype() {
	// this implementation follows instructions at http://www.cc.gatech.edu/projects/ihpcl/mpichdoc/www3/MPI_Type_struct.html
	MPI_Datatype newtype;
	int blen[4];
	MPI_Aint indices[4];
	MPI_Datatype oldtypes[4];
	Msg_header h;

    blen[0] = 1; indices[0] = 0; oldtypes[0] = MPI_LONG_LONG;
    blen[1] = 1; indices[1] = (char *) &h.datatypesz - (char *) &h; oldtypes[1] = MPI_INT;
    blen[2] = 1; indices[2] = (char *) &h.count - (char *) &h; oldtypes[2] = MPI_INT;
    blen[3] = 1; indices[3] = sizeof(Msg_header); oldtypes[3] = MPI_UB;
    MPI_Type_struct( 4, blen, indices, oldtypes, &newtype );
    MPI_Type_commit(&newtype);

    return newtype;
}

static inline MPI_Datatype create_headed_packed_datatype() {
	// this implementation follows instructions at http://www.cc.gatech.edu/projects/ihpcl/mpichdoc/www3/MPI_Type_struct.html
	MPI_Datatype newtype;
    MPI_Type_contiguous(sizeof(MPI_PACKED), MPI_BYTE, &newtype);
    MPI_Type_commit(&newtype);

    return newtype;
}

static int neigcmp(const void *n1, const void *n2) {
	return *(int *) n1 - *(int *) n2;
}

/*
 * Sets the graph of the following computation
 * Error codes returned are originated at MPI functions
 */
static int get_topology(MPI_Comm comm) {
	int i, j;
	int rank;
	int size;
	int ngsize;
	int *ngarray;
	int retcode;

	MPI_Topo_test ( comm, &topology_type );
	switch (topology_type) {
		case MPI_UNDEFINED:
			if ((retcode = MPI_Comm_rank (comm, &rank)) != MPI_SUCCESS)
				return retcode;
			if ((retcode = MPI_Comm_size (comm, &size)) != MPI_SUCCESS)
				return retcode;
			ngsize = size - 1;
			ngarray = (int *) calloc(ngsize, sizeof(int));	
			for (i = 0, j = 0; i < size; i++)
				if (i != rank)
				ngarray[j++] = i;
			break;
		case MPI_CART: {
			int ndims;
			if ((retcode = MPI_Cartdim_get ( comm, &ndims )) != MPI_SUCCESS)
				return retcode;
			int dims[ndims];
			int periods[ndims];
			int coords[ndims];
			if ((retcode = MPI_Cart_get ( comm, ndims, dims, periods, coords )) != MPI_SUCCESS)
				return retcode;
			if ((retcode = MPI_Cart_rank ( comm, coords, &rank )) != MPI_SUCCESS)
				return retcode;
			size = 1;
			int d;
			for (d = 0; d < ndims; d++) {
				size *= dims[d];
				if (coords[d] < dims[d] - 1)
					ngsize++;
				if (coords[d] > 0)
					ngsize++;
			}
			ngarray = (int *) calloc(ngsize, sizeof(int));
			int i;
			for (d = 0, i = 0; d < ndims; d++) {
				if (coords[d] < dims[d] - 1) {
					coords[d]++;
					if ((retcode = MPI_Cart_rank ( comm, coords, &ngarray[i++] )) != MPI_SUCCESS)
						return retcode;
					coords[d]--;
				}
				if (coords[d] > 0) {
					coords[d]--;
					if ((retcode = MPI_Cart_rank ( comm, coords, &ngarray[i++] )) != MPI_SUCCESS)
						return retcode;
					coords[d]++;
				}
			}
			break;
		}
		case MPI_GRAPH:
			if ((retcode = MPI_Comm_rank (comm, &rank)) != MPI_SUCCESS)
				return retcode;
			int nedges;
			if ((retcode = MPI_Graphdims_get ( comm, &size, &nedges )) != MPI_SUCCESS)
				return retcode;
			if ((retcode = MPI_Graph_neighbors_count ( comm, rank, &ngsize )) != MPI_SUCCESS)
				return retcode;
			ngarray = (int *) calloc(ngsize, sizeof(int));
			if ((retcode = MPI_Graph_neighbors ( comm, rank, ngsize, ngarray )) != MPI_SUCCESS)
				return retcode;
			break;
		default:
			return MPI_ERR_OP;
	}

	topology_comm = comm;
	topology_rank = rank;
	topology_size = size;
	topology_neigh_size = ngsize;
	topology_neigh_array = ngarray;
	

	if (topology_neigh_array != ED_NULL)
		qsort(topology_neigh_array, topology_neigh_size, sizeof(int), neigcmp);

	return MPI_SUCCESS;
}

int MPI_Event_set_model(int(*event)(MPI_Status *), int(*init_fun)()) {
	if (event == ED_NULL || init_fun == ED_NULL)
		return MPI_ERR_ARG;
	model_event = event;
	init=init_fun;	
	
	return MPI_SUCCESS;
}

int MPI_Event_set_config_param(int param, int value) {
	switch (param) {
		case ED_TERMINATION:
			if (value < ED_TERMINATION_DEFAULT || value > ED_TERMINATION_OFF)
				return MPI_ERR_ARG;
			config_term = value;
			break;
		case ED_MSG_ORDER:
			if (value < ED_MSG_ORDER_DEFAULT || value > ED_MSG_ORDER_OFF)
				return MPI_ERR_ARG;
			config_msg_order = value;
			break;
		default:
			return MPI_ERR_ARG;
			break;
	}

	return MPI_SUCCESS;
}

int MPI_Event_run(MPI_Comm comm, MPI_Status *comp_status) {
	int mpi_return;
	int i;
	
	if ((mpi_return = get_topology(comm)) != MPI_SUCCESS)
		return mpi_return;

	comp_is_running = 1;
	
	termination_ack_array = (int *) calloc(topology_size, sizeof(int));
	
	for(i=0; i<topology_size;i++){
		termination_ack_array[i]=0;
	}
	

	smallmembufrepo = newRepo(ED_SMALL_MEM_BUF, 32);
	mediummembufrepo = newRepo(ED_MEDIUM_MEM_BUF, 8);
	largemembufrepo = newRepo(ED_LARGE_MEM_BUF, 2);

	postmsgrepo = newRepo(sizeof(Postponed_msg), 16);
	statrepo = newRepo(sizeof(MPI_Status), 16);
	postmsg = (Postponed_msg **) getMemoryBuffer(sizeof(Postponed_msg *) << 4);
	postmsgsz = 16;
	npostmsg = 0;

	int (*comp_send_impl)( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
	
	comp_send_impl = _ED_Send_impl_;
	comp_msg_header_sz = 0;
	
	init();

	if(topology_rank == 0){

	
		comp_event_msg = ED_NULL;
		MPI_Event_Send = comp_send_impl;
		termination_lterm = model_event(comp_event_msg);
		MPI_Event_Send = _ED_Send_no_impl_;

						
	}
	
	while ((config_term == ED_TERMINATION_DEFAULT && !termination_lterm) || (config_term == ED_TERMINATION_DIFFUSING && termination_term < topology_neigh_size)) {
		int tryevent;

		
			
				comp_event_msg = (MPI_Status *) getFromRepo(statrepo);
			
				if ((mpi_return = MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, topology_comm, comp_event_msg)) != MPI_SUCCESS) {
					
					*comp_status = *comp_event_msg;
					return mpi_return;
				}	
							
				if (comp_event_msg->MPI_TAG == ED_DIFFUSING_ACK || comp_event_msg->MPI_TAG == ED_DIFFUSING_TERM){ 
				
					if (comp_event_msg->MPI_TAG == ED_DIFFUSING_ACK) {
						MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, ED_DIFFUSING_ACK, topology_comm, comp_event_msg);
						//printf("i am %d and i've received ack %d from %d and my deficit BEFORE is %d \n", topology_rank, i, comp_event_msg->MPI_SOURCE, termination_deficit );
						termination_deficit -= i;
						//printf("i am %d and i've received ack %d from %d and my deficit AFTER is %d \n", topology_rank, i, comp_event_msg->MPI_SOURCE, termination_deficit );
					}
					else if (comp_event_msg->MPI_TAG == ED_DIFFUSING_TERM) {
						MPI_Recv(&i, 1, MPI_INT, MPI_ANY_SOURCE, ED_DIFFUSING_TERM, topology_comm, comp_event_msg);
						if (termination_term == 0 && topology_rank !=0 ) {
							for (i = 0; i < topology_neigh_size; i++) {
								//printf("i am %d sending difusing term to %d\n", topology_rank, topology_neigh_array[i]);
								MPI_Event_Send = comp_send_impl;
		        	                                MPI_Event_Send(&i, 1, MPI_INT, topology_neigh_array[i], ED_DIFFUSING_TERM);	
								MPI_Event_Send = _ED_Send_no_impl_;
							} 
	 					}
					        termination_term++;
					}
				}		
				
				else {
					
					comp_msg_was_received = 0;
					MPI_Event_Recv = _ED_Recv_impl_;
					MPI_Event_Send = comp_send_impl;
							
					termination_lterm = model_event(comp_event_msg);
		
					//*comp_status = *comp_event_msg;
					MPI_Event_Recv = _ED_Recv_no_impl_;
					MPI_Event_Send = _ED_Send_no_impl_;
									
					if (!comp_msg_was_received) {
						*comp_status = *comp_event_msg;
						comp_status->MPI_ERROR = MPI_ERR_OTHER;
						return MPI_ERR_IN_STATUS;
					}
	
					if (comp_event_msg != ED_NULL){
						if (config_term == ED_TERMINATION_DIFFUSING){
						//printf("rank %d received message---------------------source %d acks pending %d to him. my parent %d\n",topology_rank, comp_event_msg->MPI_SOURCE, termination_ack_array[comp_event_msg->MPI_SOURCE], termination_parent);
						}
						if (config_term == ED_TERMINATION_DIFFUSING && termination_ack_array[comp_event_msg->MPI_SOURCE] >= 16 && comp_event_msg->MPI_SOURCE != termination_parent) {
							//printf("i am %d sending acks %d to %d tagout %d\n",topology_rank, termination_ack_array[comp_event_msg->MPI_SOURCE], comp_event_msg->MPI_SOURCE, ED_DIFFUSING_ACK);	
							MPI_Event_Send = comp_send_impl;					
							MPI_Event_Send(&termination_ack_array[comp_event_msg->MPI_SOURCE], 1, MPI_INT, comp_event_msg->MPI_SOURCE, ED_DIFFUSING_ACK);
							MPI_Event_Send = _ED_Send_no_impl_;

							termination_ack_array[comp_event_msg->MPI_SOURCE] = 0;
			    			}
						release(statrepo, comp_event_msg);
					}
				}

		if (config_term == ED_TERMINATION_DIFFUSING && termination_lterm == 1) {
			//printf("i am %d sai do while term deficit %d\n", topology_rank, termination_deficit);
			for (i = 0; i < topology_size; i++) {
				if(i != topology_rank)
				//printf("i am %d and i have acks %d to %d \n",topology_rank, termination_ack_array[i], i);

				if (termination_ack_array[i] > 0 && i != topology_rank && i != termination_parent) {
					//printf("i am %d and i am sending flush acks %d to %d \n",topology_rank, termination_ack_array[i], i);
					MPI_Event_Send = comp_send_impl;			
		   			MPI_Event_Send(&termination_ack_array[i], 1, MPI_INT, i, ED_DIFFUSING_ACK);
					MPI_Event_Send = _ED_Send_no_impl_;
			        	termination_ack_array[i] = 0;
		        	}
	    		}
	
    			if(termination_deficit == 0){
		       		if (termination_parent != -1) {
					//printf("i am %d and i am sending term with parent to %d \n",topology_rank, termination_parent);
					MPI_Event_Send = comp_send_impl;	
					MPI_Event_Send(&termination_ack_array[termination_parent], 1, MPI_INT, termination_parent, ED_DIFFUSING_ACK);
					MPI_Event_Send = _ED_Send_no_impl_;
					termination_ack_array[termination_parent] = 0;
					termination_parent = -1;
		       		}
	       			else  {
					if(termination_sent_term == 0){
							for (i = 0; i < topology_neigh_size; i++) {
							    //printf("i am %d and i am sending term to %d \n",topology_rank, topology_neigh_array[i]);
							    MPI_Event_Send = comp_send_impl;	
							    MPI_Event_Send(&i, 1, MPI_INT, topology_neigh_array[i], ED_DIFFUSING_TERM);
							    MPI_Event_Send = _ED_Send_no_impl_;
							}
							termination_sent_term = 1;	
						}
		                	}
				}
			}
	}
	

	

	freeRepo(postmsgrepo);
	freeRepo(statrepo);
	
	for (i = 0; i < npostmsg; i++) {
		int l;
		MPI_Get_count(postmsg[i]->status, postmsg[i]->datatype, &l);
		if (l > ED_LARGE_MEM_BUF)
			free(postmsg[i]->buf);
	}
	if (postmsgsz*sizeof(Postponed_msg *) > ED_LARGE_MEM_BUF)
		free(postmsg);
	freeRepo(smallmembufrepo);
	freeRepo(mediummembufrepo);
	freeRepo(largemembufrepo);
	//freeBaseRepo();

	if (MPI_MSG_HEADER != ED_NULL)
		MPI_Type_free(&MPI_MSG_HEADER);
	if (MPI_HEADED_PACKED != ED_NULL)
		MPI_Type_free(&MPI_HEADED_PACKED);

	if (topology_neigh_array != ED_NULL)
		free(topology_neigh_array);

	comp_is_running = 0;

	return MPI_SUCCESS;
}


static int _ED_Recv_no_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _ED_Recv_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) {
	int i=0;
	MPI_Event_Recv = _ED_Recv_no_impl_;
	comp_msg_was_received = 1;
	
	comp_buf = buf;
	comp_datatype = datatype;
	int ret = MPI_Recv(buf, count, datatype, comp_event_msg->MPI_SOURCE, comp_event_msg->MPI_TAG, topology_comm, status);
	*comp_event_msg = *status;
	if(config_term == ED_TERMINATION_DIFFUSING){	
		if (termination_parent  == -1 && topology_rank != 0 && termination_lterm == 1) {
			termination_parent = comp_event_msg->MPI_SOURCE;
	       		//printf("i am %d my parent is %d\n",topology_rank, termination_parent); 	
	        } 
		termination_ack_array[comp_event_msg->MPI_SOURCE]++;
		//printf("i am %d and adding an ack to %d. %d acks to him\n",topology_rank, comp_event_msg->MPI_SOURCE, termination_ack_array[comp_event_msg->MPI_SOURCE]); 
	}
	return ret;
	
	//  this test should be made only in case postpone message function returns an MPI_Request
//	int cancelled;
//	MPI_Test_cancelled(comp_event_msg, &cancelled);
//	if (cancelled) {
//		return MPI_ERR_OP;
//	}

}

int MPI_Event_neighbors_test(int dest) {
	return comp_is_running && (topology_neigh_size == topology_size - 1 || bsearch(&dest, topology_neigh_array, topology_neigh_size, sizeof(int), neigcmp) == ED_NULL);
}

static int _ED_Send_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _ED_Send_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) {
	//printf("%d is sending\n",topology_rank);	
	if(config_term == ED_TERMINATION_DIFFUSING && tag != ED_DIFFUSING_TERM && tag != ED_DIFFUSING_ACK){
		termination_deficit++;
		//printf("i am %d and adding 1 to termination: %d.sending message to %d\n",topology_rank, termination_deficit, dest); 	
	}
	return MPI_Send( buf, count, datatype, dest, tag, topology_comm );
}

static int _ED_Isend_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _ED_Isend_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) {
	
	if(config_term == ED_TERMINATION_DIFFUSING){
	termination_deficit++;
	}
	return MPI_Isend( buf, count, datatype, dest, tag, topology_comm, request );
}

