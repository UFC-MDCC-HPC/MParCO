* Implementation of the pulse-driven primitives of the MParCO framework 
* Author: Allberson/Correa

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <heap.h>
#include <repository.h>
#include <mparco_pd_mpi.h>

/*
 * Types
 */

typedef struct P_msg {
	MPI_Status *    status;
	void *          buf;
	MPI_Datatype    datatype;
	long long       pulse;
} Postponed_msg;

typedef struct Ps_msg {
	void *          buf;
	int             count;
	MPI_Datatype    datatype;
	int             dest;
	int             tag;
	struct Ps_msg * next;
} Pulse_msg;

typedef struct M_h {
	long long    pulse;             // pulse
	int          datatypesz;
	int          count;
} Msg_header;

/*
 * Datatypes
 */
static MPI_Datatype MPI_MSG_HEADER    = PD_NULL;
static MPI_Datatype MPI_HEADED_PACKED = PD_NULL;

/*
 * Configuration parameters
 */

static int config_term      = PD_TERMINATION_OFF;
static int config_pulse_gen = PD_PULSE_GENERATION_AUTOMATIC;
static int config_msg_order = PD_MSG_ORDER_OFF;

/*
 * Topology
 */

static MPI_Comm topology_comm        = PD_NULL;
static int      topology_type        = MPI_UNDEFINED;
static int      topology_rank        = -1;
static int      topology_size        = -1;
static int      topology_neigh_size  = 0;
static int *    topology_neigh_array = NULL;   // is sorted in ascending order

/*
 * Model functions
 */
static int (*model_event)(MPI_Status *);
static int (*model_pulse)(long long);

/*
 * Computation
 */
static long long    comp_pulse;
static long long    comp_advance;
static int          comp_suspended;
static MPI_Status * comp_event_msg;
static int          comp_postponed_msg;
static int          comp_msg_was_received;
static void *       comp_buf;
static MPI_Datatype comp_datatype;
static int          comp_is_running           = 0;
static Pulse_msg *  comp_pulse_msg;
static int          comp_msg_header_sz;

/*
 * Memory management
 */

// heap of postponed messages
static Postponed_msg ** postmsg;
static size_t           npostmsg;
static size_t           postmsgsz;
// repository of postponed messages
static void *           postmsgrepo;
// repository of pulse messages
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

// callable from event

static int _PD_Recv_no_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );
static int _PD_Recv_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status );
int (*MPI_Pulse_Recv)( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) = _PD_Recv_no_impl_;

static int _PD_Postpone_msg_no_impl_(long long targetpulse);
static int _PD_Postpone_msg_impl_(long long targetpulse);
int (*MPI_Pulse_postpone_msg)(long long targetpulse) = _PD_Postpone_msg_no_impl_;
static int postpone_msg(void * buf, long long targetpulse);

static int _PD_Resume_pulse_counter_no_impl_();
static int _PD_Resume_pulse_counter_impl_();
int (*MPI_Pulse_resume_pulse_counter)() = _PD_Resume_pulse_counter_no_impl_;

// callable from pulse

static int _PD_Send_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
static int _PD_Send_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
static int _PD_Send_with_pulse_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
int (*MPI_Pulse_Send)( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) = _PD_Send_no_impl_;

static int _PD_Isend_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );
static int _PD_Isend_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request );
int (*MPI_Pulse_Isend)( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) = _PD_Isend_no_impl_;

static int _PD_Advance_pulse_counter_no_impl_(int newpulse);
static int _PD_Advance_pulse_counter_impl_(int newpulse);
int (*MPI_Pulse_advance_pulse_counter)(int newpulse) = _PD_Advance_pulse_counter_no_impl_;

static int _PD_Suspend_pulse_counter_no_impl_();
static int _PD_Suspend_pulse_counter_impl_();
int (*MPI_Pulse_suspend_pulse_counter)() = _PD_Suspend_pulse_counter_no_impl_;

/*
 * Local functions
 */

static inline void * getMemoryBuffer(int l) {
	if (0 <= l && l <= PD_SMALL_MEM_BUF)
		return getFromRepo(smallmembufrepo);
	if (l <= PD_MEDIUM_MEM_BUF)
		return getFromRepo(mediummembufrepo);
	if (l <= PD_LARGE_MEM_BUF)
		return getFromRepo(largemembufrepo);
	return (void *) malloc(l);
}

static inline void releaseMemoryBuffer(void * m, int l) {
	if (0 <= l && l <= PD_SMALL_MEM_BUF)
		release(smallmembufrepo, m);
	else if (l <= PD_MEDIUM_MEM_BUF)
		release(mediummembufrepo, m);
	else if (l <= PD_LARGE_MEM_BUF)
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
			ngarray = NULL;
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

	if (topology_neigh_array != NULL)
		qsort(topology_neigh_array, topology_neigh_size, sizeof(int), neigcmp);

	return MPI_SUCCESS;
}

int MPI_Pulse_set_model(int(*event)(MPI_Status *), int(*targetpulse)(long long)) {
	if (event == NULL || targetpulse == NULL)
		return MPI_ERR_ARG;
	model_event = event;
	model_pulse = targetpulse;

	return MPI_SUCCESS;
}

int MPI_Pulse_set_config_param(int param, int value) {
	switch (param) {
		case PD_TERMINATION:
			if (value < PD_TERMINATION_DEFAULT || value > PD_TERMINATION_OFF)
				return MPI_ERR_ARG;
			config_term = value;
			break;
		case PD_PULSE_GENERATION:
			if (value < PD_PULSE_GENERATION_DEFAULT || value > PD_PULSE_GENERATION_OFF)
				return MPI_ERR_ARG;
			config_pulse_gen = value;
			break;
		case PD_MSG_ORDER:
			if (value < PD_MSG_ORDER_DEFAULT || value > PD_MSG_ORDER_OFF)
				return MPI_ERR_ARG;
			config_msg_order = value;
			break;
		default:
			return MPI_ERR_ARG;
			break;
	}

	return MPI_SUCCESS;
}

int MPI_Pulse_run(MPI_Comm comm, MPI_Status *comp_status) {
	int mpi_return;
	int comp_flag;

	if ((mpi_return = get_topology(comm)) != MPI_SUCCESS)
		return mpi_return;

	comp_is_running = 1;

	smallmembufrepo = newRepo(PD_SMALL_MEM_BUF, 32);
	mediummembufrepo = newRepo(PD_MEDIUM_MEM_BUF, 8);
	largemembufrepo = newRepo(PD_LARGE_MEM_BUF, 2);

	postmsgrepo = newRepo(sizeof(Postponed_msg), 16);
	statrepo = newRepo(sizeof(MPI_Status), 16);
	postmsg = (Postponed_msg **) getMemoryBuffer(sizeof(Postponed_msg *) << 4);
	postmsgsz = 16;
	npostmsg = 0;

	int (*comp_send_impl)( void *buf, int count, MPI_Datatype datatype, int dest, int tag );
	if (config_pulse_gen == PD_PULSE_GENERATION_AUTOMATIC) {
		comp_send_impl = _PD_Send_with_pulse_impl_;
		MPI_MSG_HEADER = create_msg_header_datatype();
		MPI_HEADED_PACKED = create_headed_packed_datatype();
		MPI_Pack_size(1, MPI_MSG_HEADER, topology_comm, &comp_msg_header_sz);
	}
	else {
		comp_send_impl = _PD_Send_impl_;
		comp_msg_header_sz = 0;
	}
	MPI_Pulse_Send = comp_send_impl;
	MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_impl_;
	MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_impl_;

	comp_pulse = 0;
	comp_suspended = 0;
	int lterm = model_pulse(comp_pulse);
	comp_advance = 0;

	MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_no_impl_;
	MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_no_impl_;
	MPI_Pulse_Send = _PD_Send_no_impl_;

	while (!lterm) {
		int tryevent;
		do {
			tryevent = 0;

			while (comp_advance > 0 && (npostmsg == 0 || comp_pulse < postmsg[0]->pulse - 1)) {
				long long newpulse = comp_pulse+comp_advance;
				if (npostmsg > 0 && newpulse > postmsg[0]->pulse - 1)
					newpulse = postmsg[0]->pulse - 1;
				comp_advance -= newpulse-comp_pulse;
				switch (config_pulse_gen) {
				case PD_PULSE_GENERATION_AUTOMATIC:
					if (comp_pulse < newpulse) {
						MPI_Pulse_Send = comp_send_impl;
						MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_impl_;
						MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_impl_;
						while (comp_pulse < newpulse && !lterm)
							lterm = model_pulse(++comp_pulse);
						MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_no_impl_;
						MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_no_impl_;
						MPI_Pulse_Send = _PD_Send_no_impl_;

						if (lterm)
							goto TERMINATION;
					}
					break;
				case PD_PULSE_GENERATION_EMPTY:
					comp_pulse = newpulse;
					break;
				default:
					return MPI_ERR_OP;;
					break;
				}
			}

			comp_flag = 0;

			if (npostmsg > 0 && postmsg[0]->pulse == comp_pulse + 1) {
				comp_flag = 1;
				comp_event_msg = postmsg[0]->status;
				comp_postponed_msg = 1;
			}
			else {
				comp_event_msg = (MPI_Status *) getFromRepo(statrepo);
				if (comp_suspended) {
					if ((mpi_return = MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, topology_comm, comp_event_msg)) != MPI_SUCCESS) {
						*comp_status = *comp_event_msg;
						return mpi_return;
					}
				}
				else
					if ((mpi_return = MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, topology_comm, &comp_flag, comp_event_msg)) != MPI_SUCCESS) {
						*comp_status = *comp_event_msg;
						return mpi_return;
					}
				if (comp_suspended || comp_flag) {
					if (config_pulse_gen == PD_PULSE_GENERATION_AUTOMATIC) {
						int count;
						MPI_Get_count(comp_event_msg, MPI_PACKED, &count);
						void * buf = getMemoryBuffer(count);
						int ret = MPI_Recv(buf, count, MPI_PACKED, comp_event_msg->MPI_SOURCE, comp_event_msg->MPI_TAG, topology_comm, comp_event_msg);

						int position = 0;
						Msg_header h;
						MPI_Unpack(buf, count, &position, &h, 1, MPI_MSG_HEADER, topology_comm);
						comp_datatype = MPI_HEADED_PACKED;

						MPI_Status_set_elements(comp_event_msg, MPI_BYTE, count-position);
						int p = h.pulse >= comp_pulse ? h.pulse + 1 : comp_pulse + 1;
						postpone_msg(buf, p);
						comp_flag = 0;
						tryevent = 1;
					}
					else
						comp_postponed_msg = 0;
				}
			}

			if (comp_suspended || comp_flag) {
				comp_msg_was_received = 0;
				MPI_Pulse_Recv = _PD_Recv_impl_;
				MPI_Pulse_postpone_msg = _PD_Postpone_msg_impl_;
				MPI_Pulse_resume_pulse_counter = _PD_Resume_pulse_counter_impl_;
				if ((mpi_return = model_event(comp_event_msg)) != MPI_SUCCESS) {
					*comp_status = *comp_event_msg;
					MPI_Pulse_Recv = _PD_Recv_no_impl_;
					MPI_Pulse_postpone_msg = _PD_Postpone_msg_no_impl_;
					MPI_Pulse_resume_pulse_counter = _PD_Resume_pulse_counter_no_impl_;
					return mpi_return;
				}
				MPI_Pulse_Recv = _PD_Recv_no_impl_;
				MPI_Pulse_postpone_msg = _PD_Postpone_msg_no_impl_;
				MPI_Pulse_resume_pulse_counter = _PD_Resume_pulse_counter_no_impl_;
				if (!comp_msg_was_received) {
					*comp_status = *comp_event_msg;
					comp_status->MPI_ERROR = MPI_ERR_OTHER;
					return MPI_ERR_IN_STATUS;
				}
				if (comp_event_msg != NULL)
					release(statrepo, comp_event_msg);
				tryevent = 1;
			}
		} while (tryevent);

		MPI_Pulse_Send = comp_send_impl;
		MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_impl_;
		MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_impl_;

		lterm = model_pulse(++comp_pulse);

		MPI_Pulse_suspend_pulse_counter = _PD_Suspend_pulse_counter_no_impl_;
		MPI_Pulse_advance_pulse_counter = _PD_Advance_pulse_counter_no_impl_;
		MPI_Pulse_Send = _PD_Send_no_impl_;
	}

	TERMINATION:

	freeRepo(postmsgrepo);
	freeRepo(statrepo);
	int i;
	for (i = 0; i < npostmsg; i++) {
		int l;
		MPI_Get_count(postmsg[i]->status, postmsg[i]->datatype, &l);
		if (l > PD_LARGE_MEM_BUF)
			free(postmsg[i]->buf);
	}
	if (postmsgsz*sizeof(Postponed_msg *) > PD_LARGE_MEM_BUF)
		free(postmsg);
	freeRepo(smallmembufrepo);
	freeRepo(mediummembufrepo);
	freeRepo(largemembufrepo);
	freeBaseRepo();

	if (MPI_MSG_HEADER != PD_NULL)
		MPI_Type_free(&MPI_MSG_HEADER);
	if (MPI_HEADED_PACKED != PD_NULL)
		MPI_Type_free(&MPI_HEADED_PACKED);

	if (topology_neigh_array != NULL)
		free(topology_neigh_array);

	comp_is_running = 0;

	return MPI_SUCCESS;
}

static int pulsecmp(const void *m1, const void *m2) {
	// since field pulse is long long
	int ret = (*(Postponed_msg **) m1)->pulse > (*(Postponed_msg **) m2)->pulse;
	return ret ? 1 : -1;
}

// this should return an MPI_Request to allow further inspections (like MPI_Cancel or other)?
static int postpone_msg(void * buf, long long targetpulse) {
	Postponed_msg * newmsg = (Postponed_msg *) getFromRepo(postmsgrepo);
	newmsg->buf = buf;
	newmsg->datatype = comp_datatype;
	newmsg->pulse = targetpulse;
	newmsg->status = comp_event_msg;
	comp_event_msg = NULL;

	if (npostmsg == postmsgsz) {
		int length = postmsgsz * sizeof(Postponed_msg *);
		Postponed_msg ** newpm = (Postponed_msg **) getMemoryBuffer(length + (sizeof(Postponed_msg *) << 4));
		memcpy(newpm, postmsg, length);
		releaseMemoryBuffer(postmsg, length);
		postmsg = newpm;
		postmsgsz += 16;
	}

	heapoffer(&newmsg, postmsg, npostmsg++, sizeof(Postponed_msg *), pulsecmp);

	return MPI_SUCCESS;
}

static int _PD_Postpone_msg_no_impl_(long long targetpulse) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

// this should return an MPI_Request to allow further inspections (like MPI_Cancel or other)?
static int _PD_Postpone_msg_impl_(long long targetpulse) {
	if (targetpulse <= comp_pulse)
		return MPI_ERR_OP;

//  is that possible this message has been canceled?
//	int canceled;
//	MPI_Test_cancelled(comp_event_msg, &canceled);
//	if (canceled) {
//		comp_has_event_msg = 0;
//		return MPI_ERR_OP;
//	}

	int l;
	MPI_Get_count(comp_event_msg, comp_datatype, &l);
	void * buf = getMemoryBuffer(l);
	memcpy(buf, comp_buf, l);

	return postpone_msg(buf, targetpulse);
}

static int _PD_Resume_pulse_counter_no_impl_() {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Resume_pulse_counter_impl_() {
	comp_suspended = 0;

	return MPI_SUCCESS;
}

static int _PD_Advance_pulse_counter_no_impl_(int newpulse) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Advance_pulse_counter_impl_(int newpulse) {
	comp_advance = newpulse-comp_pulse-1;

	return MPI_SUCCESS;
}

static int _PD_Suspend_pulse_counter_no_impl_() {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Suspend_pulse_counter_impl_() {
	comp_suspended = 1;

	return MPI_SUCCESS;
}

static int _PD_Recv_no_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Recv_impl_( void *buf, int count, MPI_Datatype datatype, MPI_Status *status ) {
	MPI_Pulse_Recv = _PD_Recv_no_impl_;
	comp_msg_was_received = 1;

	if (!comp_postponed_msg) {
		comp_buf = buf;
		comp_datatype = datatype;
		int ret = MPI_Recv(buf, count, datatype, comp_event_msg->MPI_SOURCE, comp_event_msg->MPI_TAG, topology_comm, status);
		*comp_event_msg = *status;
		return ret;
	}

//  this test should be made only in case postpone message function returns an MPI_Request
//	int cancelled;
//	MPI_Test_cancelled(comp_event_msg, &cancelled);
//	if (cancelled) {
//		return MPI_ERR_OP;
//	}

	int sz;
	MPI_Type_size( datatype, &sz );
	size_t length = count*sz;
	size_t postlength;

	if (postmsg[0]->datatype == MPI_HEADED_PACKED) {
		MPI_Get_count(comp_event_msg, MPI_PACKED, &postlength);
		postlength += comp_msg_header_sz;

		int position = 0;
		Msg_header h;
		MPI_Unpack(postmsg[0]->buf, postlength, &position, &h, 1, MPI_MSG_HEADER, topology_comm);
		MPI_Unpack(postmsg[0]->buf, postlength, &position, buf, h.count < count ? h.count : count, datatype, topology_comm);

		comp_buf = buf;
		comp_datatype = datatype;
	}
	else {
		MPI_Get_count( comp_event_msg, postmsg[0]->datatype, &postlength );
		memcpy(buf, postmsg[0]->buf, length < postlength ? length : postlength);
	}

	comp_postponed_msg = 0;
	releaseMemoryBuffer(postmsg[0]->buf, postlength);
	release(postmsgrepo, postmsg[0]);
	heappoll(postmsg, npostmsg--, sizeof(Postponed_msg *), pulsecmp);

	*status = *comp_event_msg;

	return length < postlength - comp_msg_header_sz ? MPI_ERR_TRUNCATE : MPI_SUCCESS;
}

int MPI_Pulse_neighbors_test(int dest) {
	return comp_is_running && (topology_neigh_size == topology_size - 1 || bsearch(&dest, topology_neigh_array, topology_neigh_size, sizeof(int), neigcmp) == NULL);
}

static int _PD_Send_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Send_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) {
	return MPI_Send( buf, count, datatype, dest, tag, topology_comm );
}

static int _PD_Send_with_pulse_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag ) {
	int insize;
	MPI_Pack_size(count, datatype, topology_comm, &insize);
	int outsize = comp_msg_header_sz+insize;
	void * outbuf = getMemoryBuffer(outsize);

	Msg_header h = {comp_pulse, 0, count};
	MPI_Type_size(datatype, &h.datatypesz);

	int position = 0;
	MPI_Pack(&h, 1, MPI_MSG_HEADER, outbuf,	outsize, &position, topology_comm);
	MPI_Pack(buf, count, datatype, outbuf, outsize, &position, topology_comm);

	int ret = MPI_Send( outbuf, position, MPI_PACKED, dest, tag, topology_comm );

	releaseMemoryBuffer(outbuf, outsize);

	return ret;
}

static int _PD_Isend_no_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) {
	// use MPI_ERRORS_RETURN
	return MPI_ERR_OP;
}

static int _PD_Isend_impl_( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Request *request ) {
	return MPI_Isend( buf, count, datatype, dest, tag, topology_comm, request );
}
