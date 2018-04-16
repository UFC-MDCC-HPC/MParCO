/*******************************************************************************************************************************
 * Solving the maximum stable(independent) set problem using the algorithm described in
 *
 * E. Tomita and T. Kameda
 * An efficient branch-and-bound algorithm for finding a maximum clique with computational experiments
 * Journal of Global Optimization, 2007, 37:1, pages 95--111
 * http://dx.doi.org/10.1007/s10898-006-9039-7
 * 
 * Author: Allberson/Corrêa
 *
 * Parallel
 * --------
 *
 * A parallel (pthread based) version is obtained compiling it with _PT_PARALLEL_ defined
 *
 * Parallel parameters:
 *
 * 1) number of threads
 *
 * Distributed
 * -----------
 *
 * A distributed (MPI based) version is obtained compiling it with _PD_DISTRIBUTED_ (pulse-driven model)
 * or _ED_DISTRIBUTED_ (event-driven model) defined
 *
 * Distributed parameters:
 *
 * 1) number of donation attempts before locally terminating (#define ATTEMPTS)
 *
 * Termination:
 *
 * currently, a simplified termination detection mechanism is used. A fixed binary termination tree is determined and used
 * to propagate termination waves. In the event driven model, this same tree is used to start the computation as a
 * diffusing computation (in the sense of Dijkstra's termination algorithm).
 *
 */

#include <time.h>

#include <mcr.h>
#include <stab.h>
#include <repository.h>

#ifdef _PT_PARALLEL_
#include <pthread.h>
#define _PARALLEL_
#endif

#ifdef _PD_DISTRIBUTED_
#include <mparco_pd_mpi.h>
#define _DISTRIBUTED_
#define RECV_MSG 					MPI_Pulse_Recv
#define SEND_MSG 					MPI_Pulse_Send
#define SET_CONFIG_PARAM 			MPI_Pulse_set_config_param
#define MODEL_RUN 					MPI_Pulse_run
#define DIST_TERMINATION 			PD_TERMINATION
#define DIST_TERMINATION_DEFAULT	PD_TERMINATION_DEFAULT
#endif

#ifdef _ED_DISTRIBUTED_
#include <mparco_ed_mpi.h>
#define _DISTRIBUTED_
#define RECV_MSG 					MPI_Event_Recv
#define SEND_MSG 					MPI_Event_Send
#define SET_CONFIG_PARAM 			MPI_Event_set_config_param
#define MODEL_RUN 					MPI_Event_run
#define DIST_TERMINATION 			ED_TERMINATION
#define DIST_TERMINATION_DEFAULT 	ED_TERMINATION_DEFAULT
#endif

#ifdef _DISTRIBUTED_
#define ATTEMPTS (n >> 1)

// message tags

#define EVENT_GENER_TAG		0
#define NO_EVENT_GENER_TAG	1
#define DONATION_TAG		2   // donates subproblems, represented by an ordered subvetor of vertices
#define NO_DONATION_TAG		3   // denies donation, which occurs in response a pairing message when there are no subproblems to donate
#define PAIRING_TAG			4   // requests donation
#define SOLUTION_TAG    	5   // contains a lower bound for the maximum stable set
#define SOLUTION_ACK_TAG	6   // acknowledges a lower bound received
#define START_COMP_TAG      7
#define LOCAL_TERM_TAG		8
#define GLOBAL_TERM_TAG		9
#define PAIR_GLOBAL_TERM_TAG 10
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef _WEIGHTED_
#define VTYPE               double
#define VTYPEREPO 			doublevecrepo
#define MPI_VTYPE 			MPI_DOUBLE
#define WEIGHT(i) 			weight[i]
#define COVERSORT(a, b, c) 	coverWeightedSort(a, b, weight, c)
#define	MCRSORT(a, b, c, d)	mcrWeightedSort(a, b, weight, d)
#define MCR(a, b, c)		wmcr(a, b, weight, c)
#else
#define VTYPE               int
#define VTYPEREPO 			intvecrepo
#define MPI_VTYPE 			MPI_INT
#define WEIGHT(i) 			1
#define COVERSORT(a, b, c) 	coverSort(a, b, c)
#define	MCRSORT(a, b, c, d)	mcrSort(a, b, c, d)
#define MCR(a, b, c)		mcr(a, b, c)
#endif

typedef struct SP {
	int 		q;		// index in the search queue
	int * 		R;		// vertex set
	VTYPE * 	cover;	// (weighted) clique cover
	VTYPE 		l;  	// (weighted) level in the search tree
	int 		r;		// vertex set size
	int 		p;		// index in the vertex set of the next branch
	int 		f;		// index in the vertex set of the last branch
	int 		t;		// thread
} subproblem;

// parameters
static int nthreads;

// repositories
static void ** intvecrepo;
static void ** subprepo;

// graph's adjacency matrix
static adjMatrix g;

// the number of vertices and edges
static int 			n = 0;
static long long 	m = 0;
// parallel tasks
static int rank;
static int size;
// largest (size or weight) of a stable set
static VTYPE 	qmax;
static int 		szqmax;

#ifdef _WEIGHTED_
static void ** doublevecrepo;
// vertex weights
static double * weight;
#endif

// current subproblem
static subproblem * current;
// number of explored nodes
static long long * nnodes;
// number of threads at current expansion
static int nt;
static subproblem ** 	newsubp;
static int        *  	newq;

#ifdef _DISTRIBUTED_
static double start, end;
#else
struct timespec start, end;
#endif
static double elapsed;

#ifdef _PARALLEL_
static pthread_t *		threads;
static pthread_mutex_t	subp_mutex;
static pthread_cond_t	subp_cond;
static int				waiting_threads;
static int				all_activated;
static int				global_term;
#endif

// generates a new subproblem of the curretn subproblem of a specified thread
// the vertex vertex of the new subproblem is generated and sorted
static inline void generate(int curt) {

	int * curR;
	int Rp;
	int * R;
	int i;
	int r;
	int curp = current->p;

	newq[curt] = 1;

	curR = current->R;
	Rp = curR[curp - curt];
	R = newsubp[curt]->R;

	for (i = 0, r = 0; i < curp - curt; i++)
		if (!hasEdge(curR[i], Rp))
			R[r++] = curR[i];

	if (r > 0) { // new subproblem
		newq[curt] = 0;
		subproblem * new = newsubp[curt];

		new->l = current->l + WEIGHT(Rp);
		new->r = r;
		new->p = r - 1;
		new->f = 0;
		new->t = curt;
		COVERSORT(r, R, new->cover);

		nnodes[curt]++;
	}

}

#ifdef _DISTRIBUTED_

// backtracking control variables are global when parallelism is active
static subproblem ** subpqueue;

static int it_per_pulse;
static int it;

#ifdef _ED_DISTRIBUTED_
static int 		pair_event_generation;
static int 		i_am_in_computation;
static int		pair_knows_my_gterm;
static int		pair_terminated;
#endif
static int 		event_generated;
static int		npending_event_gen;

static int 		pending_donation;
static int 		donations_attempted;
static void * 	donation_buf;
static int 		donation_bufsz;
static int 		donation_is_available;

static int 		npairing_requests;
static int * 	pairing_inquirer;

static int 		pending_sol_sending;
static int 		npending_ack_sol;
static int * 	pending_ack_sol;
static VTYPE * 	qout;
static VTYPE 	qoutgenerated;
static int * 	solution_founder;
static int 		nsol_founder;

static int 		parent;
static int 		nchildren;
static int 		nterminated_children;
static int 		locally_terminated;
static int 		globally_terminated;

static long long my_don_request;
static long long my_don_deny;
static long long my_don;
static long long my_event_gen;
static long long don_request;
static long long don_deny;
static long long don;
static long long event_gen;

static
#ifdef _PD_DISTRIBUTED_
inline
#endif
int init() {
	pairing_inquirer = (int *) calloc(n + (n << 1),sizeof(int));
	int start = n;
	pending_ack_sol = &pairing_inquirer[start];
	start += n;
	solution_founder = &pairing_inquirer[start];

	qout = (VTYPE *) calloc(n, sizeof(VTYPE));
	qoutgenerated = 0;

	subpqueue = (subproblem **) calloc(nthreads*n, sizeof(subproblem *));
	newsubp = (subproblem **) calloc(nthreads, sizeof(subproblem *));
	newq = (int *) calloc(nthreads, sizeof(int));

	nnodes = (long long *) calloc(nthreads, sizeof(long long));

	it_per_pulse = 3;
	it = 0;
	elapsed = 0;

#ifdef _ED_DISTRIBUTED_
	if (!(rank & 0x01))
		pair_event_generation = (rank + 1);
	else
		pair_event_generation = (rank - 1);
	i_am_in_computation = 0;
	pair_knows_my_gterm = 0;
	pair_terminated = 0;
#endif
	event_generated = 0;
	npending_event_gen = 0;

	pending_donation = 0;
	npairing_requests = 0;
	donations_attempted = 0;

	// int: f, r
	// int or double: l, qmax
	// int *: R
	// int * or double *: cover
	int bufsz;
	MPI_Pack_size(2+n, MPI_INT, MPI_COMM_WORLD, &bufsz);
	donation_bufsz = bufsz;
	MPI_Pack_size(2+n, MPI_VTYPE, MPI_COMM_WORLD, &bufsz);
	donation_bufsz += bufsz;

	donation_buf = malloc(donation_bufsz);
	donation_is_available = 0;

	pending_sol_sending = 0;
	pending_ack_sol[rank] = 1;
	npending_ack_sol = 0;
	nsol_founder = 0;

	parent = ((rank + 1) >> 1) - 1;

	if (rank + 1 > size >> 1)
		nchildren = 0;
	else if (rank + 1 == size >> 1)
		nchildren = 1 + (size & 1);
	else
		nchildren = 2;

	nterminated_children = 0;
	locally_terminated = 0;
	globally_terminated = 0;

	my_don_request = 0;
	my_don_deny = 0;
	my_don = 0;
	my_event_gen = 0;
	don_request = 0;
	don_deny = 0;
	don = 0;
	event_gen = 0;

	return 0;
}

int event(MPI_Status *st) {

#ifdef _ED_DISTRIBUTED_
	if (!st && rank)
		return 0;

	if (!st) {
		i_am_in_computation = 1;
		event_generated++;
		if (nchildren) {
			SEND_MSG(&qmax, 1, MPI_VTYPE, 1 + (rank << 1), START_COMP_TAG);
			if (nchildren & 0x02)
				SEND_MSG(&qmax, 1, MPI_VTYPE, 2 + (rank << 1), START_COMP_TAG);
		}
	}
	else
#else
	if (!st)
		return 0;
#endif
		switch (st->MPI_TAG) {
		case DONATION_TAG:
			RECV_MSG(donation_buf, donation_bufsz, MPI_PACKED, st);
			donation_is_available = 1;
			pending_donation = 0;
			donations_attempted = 0;
			my_don++;
			break;
		case NO_DONATION_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			pending_donation = 0;
			donations_attempted++;
			my_don_deny++;
			break;
		case PAIRING_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			pairing_inquirer[npairing_requests++] = st->MPI_SOURCE;
			don_request++;
			break;
		case SOLUTION_TAG:
			RECV_MSG(&qout[nsol_founder], 1, MPI_VTYPE, st);
			solution_founder[nsol_founder++] = st->MPI_SOURCE;
			break;
#ifdef _ED_DISTRIBUTED_
		case START_COMP_TAG:
			RECV_MSG(&qoutgenerated, 1, MPI_VTYPE, st);
			event_generated++;
			i_am_in_computation = 1;
			if (nchildren) {
				SEND_MSG(NULL, 0, MPI_CHAR, 1 + (rank << 1), START_COMP_TAG);
				if (nchildren & 0x02)
					SEND_MSG(NULL, 0, MPI_CHAR, 2 + (rank << 1), START_COMP_TAG);
			}
			break;
		case NO_EVENT_GENER_TAG:
			event_generated++;
		case EVENT_GENER_TAG:
			event_generated++;
			npending_event_gen--;
			RECV_MSG(&qoutgenerated, 1, MPI_VTYPE, st);
			event_gen++;
			break;
		case PAIR_GLOBAL_TERM_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			pair_terminated = 1;
			break;
#endif
		case SOLUTION_ACK_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			npending_ack_sol--;
			pending_ack_sol[st->MPI_SOURCE] = 0;
			break;
		case LOCAL_TERM_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			nterminated_children++;
			break;
		case GLOBAL_TERM_TAG:
			RECV_MSG(NULL, 0, MPI_CHAR, st);
			globally_terminated = 1;
			break;
		default:
			break;
	}

#ifdef _PD_DISTRIBUTED_
	return 0;
}

int pulse(long long pulsecnt) {
	if (pulsecnt == 0)
		init();
	else {
#endif
		// incorporate donation
		if (donation_is_available) {
			int position = 0;
			MPI_Unpack(donation_buf, donation_bufsz, &position, &current->f, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(donation_buf, donation_bufsz, &position, &current->r, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(donation_buf, donation_bufsz, &position, &current->l, 1, MPI_VTYPE, MPI_COMM_WORLD);
			VTYPE qout;
			MPI_Unpack(donation_buf, donation_bufsz, &position, &qout, 1, MPI_VTYPE, MPI_COMM_WORLD);
			if (qout > qmax)
				qmax = qout;
			current->p = current->r - 1;
//			current->R = (int *) getFromRepo(intvecrepo[current->t]);
			MPI_Unpack(donation_buf, donation_bufsz, &position, current->R, current->r, MPI_INT, MPI_COMM_WORLD);
//			current->cover = (VTYPE *) getFromRepo(VTYPEREPO[current->t]);
			MPI_Unpack(donation_buf, donation_bufsz, &position, current->cover, current->r, MPI_VTYPE, MPI_COMM_WORLD);
			current->q = 0;

			donation_is_available = 0;
		}

		// donate
		while (npairing_requests > 0) {
			int dest = pairing_inquirer[--npairing_requests];
			int d = 0;
			subproblem * subp;
			while (d < current->q) {
				subp = subpqueue[d];
				if (subp->p > subp->f && subp->l + subp->cover[subp->p - 1] > qmax)
					break;
				d++;
			}

			if (d < current->q) {
				int position = 0;
				MPI_Pack(&subp->f, 1, MPI_INT, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				MPI_Pack(&subp->p, 1, MPI_INT, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				MPI_Pack(&subp->l, 1, MPI_VTYPE, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				MPI_Pack(&qmax, 1, MPI_VTYPE, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				MPI_Pack(subp->R, subp->p, MPI_INT, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				MPI_Pack(subp->cover, subp->p, MPI_VTYPE, donation_buf, donation_bufsz, &position, MPI_COMM_WORLD);
				SEND_MSG(donation_buf, position, MPI_PACKED, dest, DONATION_TAG);

				subp->f = subp->p;
				don++;
			}
			else {
				SEND_MSG(NULL, 0, MPI_CHAR, dest, NO_DONATION_TAG);
				don_deny++;
			}
		}

		// get new solution
		VTYPE maxqout = 0;
		if (qoutgenerated > maxqout)
			maxqout = qoutgenerated;
		while (nsol_founder) {
			SEND_MSG(NULL, 0, MPI_CHAR, solution_founder[--nsol_founder], SOLUTION_ACK_TAG);
			if (qout[nsol_founder] > maxqout)
				maxqout = qout[nsol_founder];
		}
		int newqmax = 0;
		if (maxqout > qmax) {
			qmax = maxqout;
			newqmax = 1;
		}

		// propagate new solution
		if ((!locally_terminated && newqmax) || pending_sol_sending) { // cannot generate new chain when locally terminated
			if (!locally_terminated && newqmax)
				pending_sol_sending = 2;

			if (npending_ack_sol < size - 1) {
				int pair = lrand48() % size;

#ifdef _ED_DISTRIBUTED_
				while (pair == rank || pending_ack_sol[pair] || (npending_ack_sol < size - 2 && pair == pair_event_generation)) {
#else
				while (pair == rank || pending_ack_sol[pair]) {
#endif
					if (++pair == rank || pair == size)
							if (++pair >= size)
								pair = 0;
				}

				npending_ack_sol++;
				pending_ack_sol[pair] = 1;
				SEND_MSG(&qmax, 1, MPI_VTYPE, pair, SOLUTION_TAG);

				if (pending_sol_sending)
					pending_sol_sending--;
			}
		}
#ifdef _PD_DISTRIBUTED_
	}
#endif
#else
static int expand() {
	printf("EXPAND\n");

	subproblem * subpqueue[nthreads*current->r];
	newsubp = (subproblem **) calloc(nthreads, sizeof(subproblem *));
	newq = (int *) calloc(nthreads, sizeof(int));
	nnodes = (long long *) calloc(nthreads, sizeof(long long));
#endif
	int q = current->q;
	int p = current->p;
	while (q >= 0) {

		if (!(nnodes[0] & 0x07FF)) {
#ifndef _EXECUTABLE_
			elapsed = 0;
#else
#ifdef _DISTRIBUTED_
			end = MPI_Wtime();
			elapsed = end - start;
#else
			clock_gettime(CLOCK_MONOTONIC, &end);
			elapsed = (end.tv_sec - start.tv_sec);
			elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
#endif
#endif
		}

		while (elapsed < 1800 && p >= current->f && current->l + current->cover[p] > qmax) {
//			if (q > 0)
//				printf("RANK=%d qmax=%d f=%d p=%d l=%d cover=%d r=%d ub=%d\n", rank, qmax, current->f, p, current->l, current->cover[p], current->r, subpqueue[0]->cover[subpqueue[0]->p]);
//			else
//				printf("RANK=%d qmax=%d f=%d p=%d l=%d cover=%d r=%d ub=%d\n", rank, qmax, current->f, p, current->l, current->cover[p], current->r, current->cover[current->p]);

			int t;

#ifdef _PARALLEL_
			pthread_mutex_lock(&subp_mutex);
			all_activated = 0;
			for (nt = 0; nt < nthreads && p - nt >= current->f && current->l + current->cover[p - nt] > qmax; nt++)
				if (newsubp[nt] == NULL) {
					newsubp[nt] = (subproblem *) getFromRepo(subprepo[nt]);
					newsubp[nt]->R = (int *) getFromRepo(intvecrepo[nt]);
					newsubp[nt]->cover = (VTYPE *) getFromRepo(VTYPEREPO[nt]);
				}
			pthread_cond_broadcast(&subp_cond);
			pthread_cond_wait(&subp_cond, &subp_mutex);
			pthread_mutex_unlock(&subp_mutex);
#else

			for (nt = 0; nt < nthreads && p - nt >= current->f && current->l + current->cover[p - nt] > qmax; nt++)
				if (newsubp[nt] == NULL) {
					newsubp[nt] = (subproblem *) getFromRepo(subprepo[nt]);
					newsubp[nt]->R = (int *) getFromRepo(intvecrepo[nt]);
					newsubp[nt]->cover = (VTYPE *) getFromRepo(VTYPEREPO[nt]);
				}

			for (t = 0; t < nt; t++) {
				generate(t);
			}
#endif

			// push current subproblem
			current->p = p-nt;
			subpqueue[q++] = current;
			// push new subproblems
			int tt, firstt = -1, maxtt = -1;
			for (tt = 0; tt < nt; tt++)
				if (newq[tt]) {
					newq[0]++;
					if (maxtt == -1 || WEIGHT(current->R[p - tt]) > WEIGHT(current->R[p - maxtt]))
						maxtt = tt;
				}
				else {
					if (firstt < 0)
						firstt = tt;
					else {
						subpqueue[q++] = newsubp[tt];
						newsubp[tt] = NULL;
					}
				}

			// update best sol
			if (newq[0] && current->l >= qmax) {
				qmax = current->l + WEIGHT(current->R[p - maxtt]);

#ifdef _DISTRIBUTED_
				if (npending_ack_sol < size - 1) {
					int pair = lrand48() % size;

#ifdef _ED_DISTRIBUTED_
					while (pair == rank || pending_ack_sol[pair] || (npending_ack_sol < size - 2 && pair == pair_event_generation)) {
#else
					while (pair == rank || pending_ack_sol[pair]) {
#endif
						if (++pair == rank || pair == size)
							if (++pair >= size)
								pair = 0;
					}

					npending_ack_sol++;
					pending_ack_sol[pair] = 1;
					SEND_MSG(&qmax, 1, MPI_VTYPE, pair, SOLUTION_TAG);

					if (pending_sol_sending)
						pending_sol_sending--;
				}
				else
					pending_sol_sending = 2;
#endif
			}
			// set next subproblem
			if (firstt >= 0) {
				current = newsubp[firstt];
				newsubp[firstt] = NULL;
				current->q = q;
			}
			else
				// remove current from queue
				q--;
			p = current->p;

#ifdef _DISTRIBUTED_
			if (++it == it_per_pulse) {
				it = 0;
#ifdef _ED_DISTRIBUTED_
				if (event_generated) {
					event_generated = 0;
					npending_event_gen++;
					SEND_MSG(&qmax, 1, MPI_VTYPE, pair_event_generation, EVENT_GENER_TAG);
					my_event_gen++;
				}
#endif
				return 0;
			}
#endif
		}
		q--;  // backtracking
		if (q >= 0) {
			release(intvecrepo[current->t], current->R);
			release(VTYPEREPO[current->t], current->cover);
			release(subprepo[current->t], current);
			current = subpqueue[q];
			current->q = q;
			p = current->p;
		}
		else
			current->q = q = -1;
	}

#ifdef _DISTRIBUTED_
	if (!pending_donation && donations_attempted < ATTEMPTS) {
		pending_donation = 1;
		int pair = lrand48() % size;
		if (pair == rank)
			if (++pair == size)
				pair = 0;
		SEND_MSG(NULL, 0, MPI_CHAR, pair, PAIRING_TAG);
		my_don_request++;
	}
#ifdef _ED_DISTRIBUTED_
	else if (i_am_in_computation && !npending_ack_sol && !pending_sol_sending && donations_attempted >= ATTEMPTS && !locally_terminated && nterminated_children == nchildren) {
#else
	else if (!npending_ack_sol && !pending_sol_sending && donations_attempted >= ATTEMPTS && !locally_terminated && nterminated_children == nchildren) {
#endif
		if (parent >= 0)
			SEND_MSG(NULL, 0, MPI_CHAR, parent, LOCAL_TERM_TAG);
		else
			globally_terminated = 1;
		locally_terminated = 1;
	}

#ifdef _ED_DISTRIBUTED_
	if (!locally_terminated) {
		if (event_generated) {
			event_generated = 0;
			npending_event_gen++;
			SEND_MSG(&qmax, 1, MPI_VTYPE, pair_event_generation, EVENT_GENER_TAG);
			my_event_gen++;
		}
	}
	else
		if (event_generated == 1) {
			event_generated = 0;
			SEND_MSG(&qmax, 1, MPI_VTYPE, pair_event_generation, NO_EVENT_GENER_TAG);
			my_event_gen++;
		}
#endif

	if (!globally_terminated)
		return 0;

	if (nchildren) {
		SEND_MSG(NULL, 0, MPI_CHAR, 1 + (rank << 1), GLOBAL_TERM_TAG);
		if (nchildren > 1)
			SEND_MSG(NULL, 0, MPI_CHAR, 2 + (rank << 1), GLOBAL_TERM_TAG);
		nchildren = 0;
	}

#ifdef _ED_DISTRIBUTED_
	if (!pair_knows_my_gterm && npending_event_gen <= 0) {
		pair_knows_my_gterm = 1;
		SEND_MSG(NULL, 0, MPI_CHAR, pair_event_generation, PAIR_GLOBAL_TERM_TAG);
	}

	if (!pair_terminated || !pair_knows_my_gterm)
		return 0;
#endif

	free(donation_buf);
	free(subpqueue);
	free(pairing_inquirer);
	free(qout);
#endif

	int l;
	for (l = 0; l < nthreads; l++)
		if (newsubp[l] != NULL) {
			release(intvecrepo[l], newsubp[l]->R);
			release(VTYPEREPO[l], newsubp[l]->cover);
			release(subprepo[l], newsubp[l]);
		}

	free(newsubp);

	return 1;
}

#ifdef _PARALLEL_
void *subpGeneration(void *threadid) {
	int tid = (int) threadid;
	int t;

	pthread_mutex_lock(&subp_mutex);
	waiting_threads++;
	if (waiting_threads < nthreads)
		pthread_cond_wait(&subp_cond, &subp_mutex);
	else
		pthread_cond_broadcast(&subp_cond);

	while (!global_term) {
		while (all_activated == nt)
			pthread_cond_wait(&subp_cond, &subp_mutex);

		if (global_term)
			break;

		waiting_threads--;
		t = all_activated++;

		pthread_mutex_unlock(&subp_mutex);

		generate(t);

		pthread_mutex_lock(&subp_mutex);

		waiting_threads++;
		if (all_activated == nt && waiting_threads == nthreads)
			pthread_cond_broadcast(&subp_cond);
		else
			pthread_cond_wait(&subp_cond, &subp_mutex);
	}

	pthread_mutex_unlock(&subp_mutex);
	pthread_exit(NULL);
}
#endif

// external call must set the graph in stab.h
#ifndef _EXECUTABLE_
#ifdef _WEIGHTED_
static VTYPE wmcr(int l, int * vset, double * w, VTYPE * ccover) {
#else
static VTYPE mcr(int l, int * vset, VTYPE * ccover) {
#endif
#else
static VTYPE mcr(int l, int * vset, VTYPE * ccover) {
#endif

	int t;

#ifndef _EXECUTABLE_
	clock_gettime(CLOCK_MONOTONIC, &start);

	nthreads = 1;
	rank = 0;
	size = 1;
	intvecrepo = calloc(nthreads, sizeof(void *));
#ifdef _WEIGHTED_
	weight = w;
	doublevecrepo = calloc(nthreads, sizeof(void *));
#endif
#endif

	subprepo = calloc(nthreads, sizeof(void *));
	for (t = 0; t < nthreads; t++) {
		subprepo[t] = newRepo(sizeof(subproblem), l < 100 ? l : 100);
#ifndef _EXECUTABLE_
		intvecrepo[t] = newRepo(l*sizeof(int), l < 100 ? l : 100);
#ifdef _WEIGHTED_
		doublevecrepo[t] = newRepo(l*sizeof(double), l < 100 ? l : 100);
#endif
#endif
	}

	current = (subproblem *) getFromRepo(subprepo[0]);
	current->f = rank == size - 1 ? 0 : l-1;
	current->r = l;
	current->R = vset;
	current->cover = ccover;
	current->q = 0;
	current->l = 0;
	current->p = current->r - 1;
	current->t = 0;

#ifdef _DISTRIBUTED_
	SET_CONFIG_PARAM(DIST_TERMINATION, DIST_TERMINATION_DEFAULT);
#ifdef _PD_DISTRIBUTED_
	MPI_Pulse_set_config_param(PD_PULSE_GENERATION, PD_PULSE_GENERATION_OFF);
	MPI_Pulse_set_model(event, pulse);
#else
	MPI_Event_set_model(event, init);
#endif

	MPI_Status status;
	MODEL_RUN(MPI_COMM_WORLD, &status);
#else
	expand();
#endif

	release(subprepo[0], current);

	for (t = 0; t < nthreads; t++) {
#ifndef _EXECUTABLE_
		freeRepo(intvecrepo[t]);
#ifdef _WEIGHTED_
		freeRepo(doublevecrepo[t]);
#endif
#endif
		freeRepo(subprepo[t]);
	}

#ifndef _EXECUTABLE_
	freeBaseRepo();
	free(intvecrepo);
#ifdef _WEIGHTED_
	free(doublevecrepo);
#endif
#endif
	free(subprepo);

	return qmax;
}

#ifdef _EXECUTABLE_

int main (int argc, char *argv[]) {

	int t;

#ifdef _DISTRIBUTED_
	MPI_Init(&argc, &argv);

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	srand48(rank+time(NULL));
#else
	rank = 0;
	size = 1;
#endif

	if (rank == 0) {
		int argi = 1;
		if (argc <= argi) {
			printf("Please, specify the graph file in the DIMACS format.\n");
			return 0;
		}

		FILE * graphFile = fopen(argv[argi], "r");
		if (graphFile == NULL) {
			printf("Could not open %s.\n", argv[argi]);
			return 0;
		}

		if (readGraphFile(graphFile, complement, &n, &m)) {
			fclose(graphFile);
			printf("Could not read graph file %s.", argv[argi]);
			return 0;
		}
		fclose(graphFile);
		getStabGraph(&g);
		setAdjMatrix((&g));

#ifdef _WEIGHTED_
		argi++;
		if (argc <= argi) {
			printf("Please, specify the weight file.\n");
			return 0;
		}

		FILE * weightFile = fopen(argv[argi], "r");
		if (weightFile == NULL) {
			printf("Could not open %s.\n", argv[argi]);
			return 0;
		}

		weight = (double *) calloc(n, sizeof(double));
		if (readWeightFile(weightFile, n, weight)) {
			fclose(weightFile);
			printf("Could not read weight file %s.", argv[argi]);
			return 0;
		}
		fclose(weightFile);
#endif

		argi++;
		if (argc > argi)
			nthreads = atoi(argv[argi]);
		else
			nthreads = 1;
	}

#ifdef _DISTRIBUTED_
	int gsz[3] = {n,m,nthreads};
	MPI_Bcast(gsz, 3, MPI_INT, 0, MPI_COMM_WORLD);
	n = gsz[0];
	m = gsz[1];
	nthreads = gsz[2];
	if (n > 0 && g.matrixsize == 0) {
		newAdjMatrix(n);
		getAdjMatrix((&g));
		setStabGraph(&g);
	}
	MPI_Bcast(g.g, g.matrixsize, MPI_CHAR, 0, MPI_COMM_WORLD);
#ifdef _WEIGHTED_
	if (weight == NULL)
		weight = (double *) calloc(n, sizeof(double));
	MPI_Bcast(weight, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#endif

#ifdef _PARALLEL_
	threads = (pthread_t *) calloc(nthreads, sizeof(pthread_t));
	pthread_mutex_init(&subp_mutex, NULL);
	pthread_cond_init (&subp_cond, NULL);
	waiting_threads = 0;
	global_term = 0;
	all_activated = 0;
	nt = 0;

	int rc;
	for(t=0; t<nthreads; t++){
		printf("In main: creating thread %d\n", t);
		rc = pthread_create(&threads[t], NULL, subpGeneration, (void *) t);
		if (rc){
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	pthread_mutex_lock(&subp_mutex);
	if (waiting_threads < nthreads)
		pthread_cond_wait(&subp_cond, &subp_mutex);
	pthread_mutex_unlock(&subp_mutex);
#endif

#ifdef _DISTRIBUTED_
	start = MPI_Wtime();
#else
	clock_gettime(CLOCK_MONOTONIC, &start);
#endif

	intvecrepo = calloc(nthreads, sizeof(void *));
#ifdef _WEIGHTED_
	doublevecrepo = calloc(nthreads, sizeof(void *));
#endif
	for (t = 0; t < nthreads; t++) {
		intvecrepo[t] = newRepo(n*sizeof(int), n);
#ifdef _WEIGHTED_
		doublevecrepo[t] = newRepo(n*sizeof(double), n);
#endif
	}

	int i, ad;
	VTYPE maxantideg;
	// vertex set
	int * vset = (int *) getFromRepo(intvecrepo[0]);
	// clique cover
	VTYPE * ccover = (VTYPE *) getFromRepo(VTYPEREPO[0]);

	ad = -1;  // only matters for the weighted case
	if (rank == 0) {
		for (i = 0; i < n; i++)
			vset[i] = i;
		int antideg[n];
		i = MCRSORT(n, vset, antideg, &maxantideg);
		ad = antideg[i];
	}

#ifdef _DISTRIBUTED_
	int csz[3] = {i, ad, maxantideg};
	MPI_Bcast(csz, 3, MPI_INT, 0, MPI_COMM_WORLD);
	i = csz[0];
	ad = csz[1];
	maxantideg = csz[2];
	MPI_Bcast(vset, n, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	// size of the root subproblem
	int l;

	l = n - rank;
	if (l < i)
		i = l - 1;

	if (ad == i) { // first i vertices constitute a stab set
		qmax = ++i;
		int j;
		for (j = 0; j < i; j++)
			ccover[j] = 1;
	}
	else {
		COVERSORT(++i, vset, ccover);
		qmax = 1;
	}

	for (; i < l; i++)
		ccover[i] = ccover[i-1] + WEIGHT(vset[i]) < maxantideg ? ccover[i-1] + WEIGHT(vset[i]) : maxantideg;

	mcr(l, vset, ccover);

	release(intvecrepo[0], vset);
	release(VTYPEREPO[0], ccover);

	for (t = 0; t < nthreads; t++) {
		freeRepo(intvecrepo[t]);
#ifdef _WEIGHTED_
		freeRepo(doublevecrepo[t]);
#endif
	}
	freeBaseRepo();

	free(intvecrepo);
#ifdef _WEIGHTED_
	free(doublevecrepo);
#endif

#ifdef _DISTRIBUTED_
	end = MPI_Wtime();
	elapsed = end - start;
#else
	clock_gettime(CLOCK_MONOTONIC, &end);
	elapsed = (end.tv_sec - start.tv_sec);
	elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;

	struct rusage relapsed;
	getrusage(RUSAGE_SELF,&relapsed);
#endif

#ifdef _PARALLEL_
	pthread_mutex_lock(&subp_mutex);
	global_term = 1;
	all_activated = 0;
	nt = nthreads;
	pthread_cond_broadcast(&subp_cond);
	pthread_mutex_unlock(&subp_mutex);
#endif

	freeAdjMatrix();
	printf("Number of threads: %d\n", nthreads);
#ifdef _WEIGHTED_
	free(weight);
	printf("Largest size found(%d): %8.5g\n", rank, qmax);
#else
	printf("Largest size found(%d): %d\n", rank, qmax);
#endif
	//printf("Number of explored nodes(%d): ", rank);
	int tnodes = 0;
	for (t = 0; t < nthreads; t++) {
		//printf("[%d]=%lld ", t, nnodes[t]);
		tnodes += nnodes[t];
	}
	printf("total=%d\n", tnodes);
	printf("time=%8.5g\n", elapsed);

	free(nnodes);

#ifdef _DISTRIBUTED_
	printf("requests=%lld\ndeny=%lld\ndonations=%lld\ngens=%lld\n", don_request, my_don_deny, my_don, event_gen);
	//printf("Sent (%d): %lld pairing requests, %lld don deny, %lld donations, and %lld event gens\n", rank, my_don_request, don_deny, don, my_event_gen);
	MPI_Finalize();
#else
	printf("Elapsed time in detail: user time = %ld microsecond(s); system time = %ld microsecond(s).\n", relapsed.ru_utime.tv_sec*1000000+relapsed.ru_utime.tv_usec, relapsed.ru_stime.tv_sec*1000000+relapsed.ru_stime.tv_usec);
#endif

#ifdef _PARALLEL_
	free(threads);
	pthread_mutex_destroy(&subp_mutex);
	pthread_cond_destroy(&subp_cond);
	pthread_exit(NULL);
#endif

	return 0;
}
#else
#ifdef _WEIGHTED_
double getWeightedStabSet(adjMatrix * g, int nvset, int * vset, double * w, VTYPE qq, double ldens, double udens) {
	weight = w;
#else
int getStabSet(adjMatrix * g, int nvset, int * vset, VTYPE qq, double ldens, double udens) {
#endif
	int i, ad;
	int maxantideg;

	printf("mcr\n");

	adjMatrix gg;
	getStabGraph(&gg);
	setStabGraph(g);

	setAdjMatrix(g);

	printf("matrices set\n");

	ad = -1;  // only matters for the weighted case
	VTYPE antideg[nvset];
	i = MCRSORT(nvset, vset, antideg, &maxantideg);
	ad = antideg[i];

	int s, ss;
	for (s = 0; s < nvset; s++) {
		antideg[s] = 0;
		for (ss = 0; ss < s; ss++)
			if (!hasEdge(vset[s], vset[ss]))
				antideg[s] += WEIGHT(vset[ss]);
	}

	double antidens;
	int nenum = 1;
	VTYPE antidegsum, nvsum, wsum;
	antidegsum=nvsum=0;
	wsum=WEIGHT(vset[0]);
	for (s = 1; s < nvset; s++) {
		antidegsum += antideg[s];
		nvsum += wsum;
		wsum += WEIGHT(vset[s]);
		antidens = ((double) antidegsum)/nvsum;
		if (antidens < ldens || antidens > udens)
			nenum = s + 1;
	}

	for (s = 0; nvset >= nenum; s++) {
		int z = vset[nvset-1];
		int a, b;
		for (a = 0, b = 0; b < nvset-1; b++)
			if (!hasEdge(z, vset[b]))
				vset[a++] = vset[b];
		nvset = a;
	}

	if (!nvset)
		return s;

	if (s) {
		ad = -1;
		i = MCRSORT(nvset, vset, antideg, &maxantideg);
		ad = antideg[i];
	}

	// clique cover
	VTYPE * ccover = (VTYPE *) calloc(nvset, sizeof(VTYPE));
	if (ad == i) { // first i vertices constitute a stab set
		int l;
		for (l = 0; l < i; l++)
			ccover[l] = 1;
	}
	else
		COVERSORT(++i, vset, ccover);

	for (; i < nvset; i++)
		ccover[i] = ccover[i-1] + WEIGHT(vset[i]) < maxantideg ? ccover[i-1] + WEIGHT(vset[i]) : maxantideg;

	printf("greedy done s=%d nvset=%d\n", s, nvset);

	qmax = qq-s;
	VTYPE ret = s + MCR(nvset, vset, ccover);

	setStabGraph(&gg);
	free(ccover);

	return ret;
}
#endif

