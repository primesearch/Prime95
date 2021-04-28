/**************************************************************
 *
 *	ecm.c
 *
 *	ECM, P-1, and P+1 factoring routines
 *
 *	Original author:  Richard Crandall - www.perfsci.com
 *	Adapted to Mersenne numbers and optimized by George Woltman
 *	Further optimizations from Paul Zimmerman's GMP-ECM program
 *	Other important ideas courtesy of Peter Montgomery, Mihai Preda, Pavel Atnashev.
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	c. 1998-2021 Mersenne Research, Inc.
 *	All Rights Reserved.
 *
 *************************************************************/

/* Global variables */

int	QA_IN_PROGRESS = FALSE;
int	QA_TYPE = 0;
giant	QA_FACTOR = NULL;
int	PRAC_SEARCH = 7;

/* Macros for readability */

#define primes_less_than(x)		((double)(x) / (log ((double)(x)) - 1.0))
#define isMersenne(k,b,n,c)		((k) == 1.0 && (b) == 2 && (c) == -1)
#define isGeneralizedFermat(k,b,n,c)	((k) == 1.0 && isPowerOf2 (n) && (c) == 1)
#define isPowerOf2(n)			(((n) & ((n)-1)) == 0)

/**********************************************************************************************************************/
/*                                ECM, P-1, and P+1 best stage 2 implementation routines                              */
/**********************************************************************************************************************/

/* Various D values that we will consider in creating an ECM, P-1 or P+1 plan */

#define NUM_D		21
struct D_data {
	int	D;
	int	numrels;			// Number of values less than D that are relatively prime to D
	int	first_missing_prime;		// First prime that does not divide D
	int	second_missing_prime;		// Second prime that does not divide D
} D_data[NUM_D] = {
	{ 2 * 3 * 5, 2 * 4 / 2, 7, 11 },				// 30, 4
	{ 2 * 3 * 7, 2 * 6 / 2, 5, 11 },				// 42, 6
	{ 2 * 3 * 5 * 2, 2 * 4 / 2 * 2, 7, 11 },			// 60, 8
//	{ 2 * 5 * 7, 4 * 6 / 2, 3, 11 },				// 70, 12
	{ 2 * 3 * 5 * 3, 2 * 4 / 2 * 3, 7, 11 },			// 90, 12
	{ 2 * 3 * 5 * 4, 2 * 4 / 2 * 4, 7, 11 },			// 120, 16
	{ 2 * 3 * 5 * 5, 2 * 4 / 2 * 5, 7, 11 },			// 150, 20
	{ 2 * 3 * 5 * 7, 2 * 4 * 6 / 2, 11, 13 },			// 210, 24
	{ 2 * 3 * 5 * 11, 2 * 4 * 10 / 2, 7, 13 },			// 330, 40
	{ 2 * 3 * 5 * 7 * 2, 2 * 4 * 6 / 2 * 2, 11, 13 },		// 420, 48
	{ 2 * 3 * 7 * 11, 2 * 6 * 10 / 2, 5, 13 },			// 462, 60
	{ 2 * 3 * 5 * 7 * 3, 2 * 4 * 6 / 2 * 3, 11, 13 },		// 630, 72
	{ 2 * 3 * 5 * 7 * 4, 2 * 4 * 6 / 2 * 4, 11, 13 },		// 840, 96
//	{ 2 * 5 * 7 * 11, 4 * 6 * 10 / 2, 3, 13 },			// 770, 120
	{ 2 * 3 * 5 * 7 * 5, 2 * 4 * 6 / 2 * 5, 11, 13 },		// 1050, 120
	{ 2 * 3 * 5 * 7 * 6, 2 * 4 * 6 / 2 * 6, 11, 13 },		// 1260, 144
	{ 2 * 3 * 5 * 7 * 7, 2 * 4 * 6 / 2 * 7, 11, 13 },		// 1470, 168
	{ 2 * 3 * 5 * 7 * 8, 2 * 4 * 6 / 2 * 8, 11, 13 },		// 1680, 192
	{ 2 * 3 * 5 * 7 * 9, 2 * 4 * 6 / 2 * 9, 11, 13 },		// 1890, 216
	{ 2 * 3 * 5 * 7 * 11, 2 * 4 * 6 * 10 / 2, 13, 17 },		// 2310, 240
	{ 2 * 3 * 5 * 7 * 13, 2 * 4 * 6 * 12 / 2, 11, 17 },		// 2730, 288
//	{ 2 * 3 * 5 * 11 * 13, 2 * 4 * 10 * 12 / 2, 7, 17 },		// 4290, 480
	{ 2 * 3 * 5 * 7 * 11 * 2, 2 * 4 * 6 * 10 / 2 * 2, 13, 17 },	// 4620, 480
//	{ 2 * 3 * 7 * 11 * 13, 2 * 6 * 10 * 12 / 2, 5, 17 },		// 6006, 720
	{ 2 * 3 * 5 * 7 * 11 * 3, 2 * 4 * 6 * 10 / 2 * 3, 13, 17 }	// 6930, 720
};
#define MAX_D		6930
#define MAX_RELPRIMES	720

/* Select the best D value for the given the number of temporary gwnums that can be allocated.  We trade off more D steps vs. better */
/* prime pairing vs. different B2 start points using the ECM, P-1, or P+1 costing function. */
/* Returns the cost.  Cost function can return more information, such as best D value, B2_start, B2_end. */

double best_stage2_impl_internal (
	uint64_t C_start,		/* Starting point for bound #2 */
	uint64_t C,			/* Bound #2 */
	int	totrels,		/* Number of gwnum temporaries used for storing relative prime data */
	double	numprimes,		/* Estimated number of primes to be processed in stage 2 */
	int	max_bitarray_size,	/* Maximum size of the work bit array in MB */
	double	(*cost_func)(int, uint64_t, uint64_t, uint64_t, int, double, double, double, void *), /* ECM, P-1 or P+1 costing function */
	void	*cost_func_data)	/* User-supplied data to pass to the costing function */
{
	uint64_t B2_start, B2_end, numDsections, bitarraymaxDsections;
	double	density, multiplier, pairing_percentage, numpairs, numsingles;
	double	cost, best_cost;
	int	D, i, best_i, j;

/* Kludge to make one pass finding the best D value and a second pass to re-call the cost function using the best D. */
/* Re-calling the cost function allows it to pass back any data that P-1, P+1, or ECM may need to save. */

	best_cost = 1.0e99;
	for (j = 0; j < 2; j++) {

/* Try various values of D until we find the best one. */

	    for (i = 0; i < NUM_D; i++) {

/* On second pass only cost the best D from the first pass */

		if (j == 1) {
			if (best_cost > 1.0e98) break;    
			i = best_i;
		}

/* Check if this D value would require using too many gwnum temporaries */

		if (D_data[i].numrels > totrels) break;

/* We move the smaller primes to a composite value higher up in the B1 to B2 range to improve our pairing chances and */
/* reduce the number of D sections to process.  Calculate B2_start - the first prime that cannot be moved higher up. */

		D = D_data[i].D;
		B2_end = round_up_to_multiple_of (C - D / 2, D) + D / 2;
		B2_start = (uint64_t) floor ((double) (B2_end / (double) D_data[i].first_missing_prime));
		if (B2_start < C_start) B2_start = C_start;
		B2_start = round_down_to_multiple_of (B2_start + D / 2, D) - D / 2;
		numDsections = (B2_end - B2_start) / D;
		bitarraymaxDsections = divide_rounding_down ((uint64_t) max_bitarray_size * 8000000, (uint64_t) totrels);

/* Estimate our prime pairing percentage from the prime density and excess number of relative primes. */
/* This formula is by no means perfect, but seems to be reasonably accurate over the range of 1x to 7x relative primes. */

		density = (numprimes / (B2_end - B2_start)) * ((double) D / ((double) D_data[i].numrels * 2));
		multiplier = (double) totrels / (double) D_data[i].numrels;
		pairing_percentage = 1.0 - pow (1.0 - density, (5.0 - multiplier) * 0.05 + pow (multiplier, 0.85));
		numpairs = (pairing_percentage * numprimes) / 2.0;
		numsingles = numprimes - numpairs * 2.0;

/* Calculate the cost of this stage 2 plan */

		cost = (*cost_func) (D, B2_start, numDsections, bitarraymaxDsections, D_data[i].numrels, multiplier, numpairs, numsingles, cost_func_data);

/* On second pass, break out of loop after costing the best D from the first pass */

		if (j == 1) break;

/* Remember best cost and best D */

		if (cost < best_cost) {
			best_cost = cost;
			best_i = i;
		}
	    }
	}

/* Return best cost */

	return  (cost);
}

/* Binary search for the best number of gwnums to allocate and the best D value to use.  Caller specifies the maximum number of gwnums */
/* that can be allocated.  We trade off more D steps vs. better prime pairing vs. different B2 start points using the ECM, P-1, or P+1 costing function. */
/* Returns the cost.  Cost function can return more information, such as best D value, B2_start, B2_end. */

double best_stage2_impl (
	uint64_t C_start,		/* Starting point for bound #2 */
	uint64_t C,			/* Bound #2 */
	int	maxtotrels,		/* Maximum number of gwnum temporaries that can be used for storing relative prime data */
	double	(*cost_func)(int, uint64_t, uint64_t, uint64_t, int, double, double, double, void *), /* ECM, P-1, or P+1 costing function */
	void	*cost_func_data)	/* User-supplied data to pass to the costing function */
{
	int	max_bitarray_size;	/* Maximum size of the work bit array in MB */
	double	numprimes;		/* Estimated number of primes */
	struct best_stage2_cost {
		int	totrels;
		double	cost;
	} best[3], midpoint;

/* Return infinite cost if maxtotrels is less than D=30 needs */

	if (maxtotrels < 4) return (1.0e99);

/* Determine the maximum bit arraysize.  The costing functions must double (or more) the stage 2 setup costs if the bit array must be split in chunks. */

	max_bitarray_size = IniGetInt (INI_FILE, "MaximumBitArraySize", 250);
	if (max_bitarray_size > 2000) max_bitarray_size = 2000;
	if (max_bitarray_size < 1) max_bitarray_size = 1;

/* Estimate the number of primes between B1 and B2 */

	numprimes = primes_less_than (C) - primes_less_than (C_start);

/* Prepare for a possible binary search looking for the lowest cost stage 2 implementation varying the number of gwnums */
/* available for relative primes data. */

	if (maxtotrels > 4) {
		best[1].totrels = maxtotrels - 1;
		best[1].cost = best_stage2_impl_internal (C_start, C, best[1].totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);
	}
	best[2].totrels = maxtotrels;
	best[2].cost = best_stage2_impl_internal (C_start, C, best[2].totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);

/* We hope that most of the time a binary search is not necessary.  That is, using all available memory is best */

	if (maxtotrels == 4 || best[1].cost > best[2].cost) return (best[2].cost);

/* A binary search is required */

	best[0].totrels = 4;
	best[0].cost = best_stage2_impl_internal (C_start, C, best[0].totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);

/* Work until midpoint is better than the start point. */
/* The search code requires best[1] is better than best[0] and best[2]. */

	while (best[0].cost < best[1].cost) {
		best[2] = best[1];
		best[1].totrels = (best[0].totrels + best[2].totrels) / 2;
		best[1].cost = best_stage2_impl_internal (C_start, C, best[1].totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);
	}

/* Now we can do a binary search */

	while (best[2].totrels - best[0].totrels > 2) {
		// Work on the bigger of the lower section and upper section
		if (best[1].totrels - best[0].totrels > best[2].totrels - best[1].totrels) {		// Work on lower section
			midpoint.totrels = (best[0].totrels + best[1].totrels) / 2;
			midpoint.cost = best_stage2_impl_internal (C_start, C, midpoint.totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);
			if (midpoint.cost < best[1].cost) {			// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			midpoint.totrels = (best[1].totrels + best[2].totrels) / 2;
			midpoint.cost = best_stage2_impl_internal (C_start, C, midpoint.totrels, numprimes, max_bitarray_size, cost_func, cost_func_data);
			if (midpoint.cost < best[1].cost) {			// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Redo the best stage 2 implementation.  This lets the costing function return the best data.  Return best cost. */

	return (best_stage2_impl_internal (C_start, C, best[1].totrels, numprimes, max_bitarray_size, cost_func, cost_func_data));
}

/**********************************************************************************************************************/
/*                                        ECM, P-1, and P+1 prime pairing routines                                    */
/**********************************************************************************************************************/

#define bitset_prime(array,p)	bitset (array, ((p) - C_start) / D * numrels * 2 + map_relprime_to_index ((p) % D, D, &reldata))
#define bitclr_prime(array,p)	bitclr (array, ((p) - C_start) / D * numrels * 2 + map_relprime_to_index ((p) % D, D, &reldata))
#define bittst_prime(array,p)	bittst (array, ((p) - C_start) / D * numrels * 2 + map_relprime_to_index ((p) % D, D, &reldata))
#define bit_prime(k,j)		(C_start + (k) * D + D / 2 + map_index_to_relprime (j, D, &reldata))
#define bitset_work(array,k,j)	bitset (array, (k) * totrels + ((j) < totrels ? totrels - 1 - (j) : (j) - totrels))

// Map a D value to the D_data index
int map_D_to_index (int D)
{
	int	i;
	for (i = 0; i < NUM_D; i++) if (D_data[i].D == D) return (i);
	ASSERTG (0);
	return (0);
}

// Data for mapping relprimes
struct relmap_data {
	int	relmap1[MAX_RELPRIMES];
	int	relmap2[MAX_D];
	int	map_relprimes;
};

// Routines to manage relative primes to d (we only calculate half of the relative primes)
void initmap (int D, struct relmap_data *data)
{
	int	relp;
	data->map_relprimes = 0;
	memset (data->relmap1, 0, sizeof (data->relmap1));
	memset (data->relmap2, 0, sizeof (data->relmap2));
	for (relp = 1; relp < D / 2; relp += 2) {
		if (_intgcd (relp, D) != 1) continue;
		data->relmap1[data->map_relprimes] = relp;
		data->relmap2[relp] = data->map_relprimes;
		data->map_relprimes++;
	}
}
// Take number between 1 and d - 1, return the index for the prime bitmap
int map_relprime_to_index (int n, int D, struct relmap_data *data)
{
	ASSERTG (n < D);
	if (n < D / 2) return (data->map_relprimes + data->relmap2[n]);
	else return (data->map_relprimes - 1 - data->relmap2[D - n]);
}
// Take index between 0 and relprimes * 2 - 1 from the prime bitmap, return the relative prime
int map_index_to_relprime (int n, int D, struct relmap_data *data)
{
	ASSERTG (n < D);
	if (n < data->map_relprimes) return (- data->relmap1[data->map_relprimes - 1 - n]);
	else return (data->relmap1[n - data->map_relprimes]);
}
// Return nth relative prime to d
int map_nth_relprime (int n, int D, struct relmap_data *data)
{
	ASSERTG (n < data->map_relprimes);
	if (n < data->map_relprimes) return (data->relmap1[n]);
	else return (D - data->relmap1[data->map_relprimes - 1 - n]);
}
// Return true if n is relative prime to d
int map_is_relprime (int n, int D, struct relmap_data *data)
{
	if (n > D) n = n % D;
	if (n > D / 2) n = D - n;
	return (n == 1 || data->relmap2[n]);
}

// Count the number of primes a particular prime can pair with
// This was found to be slow for high multipliers (an N^2 algorithm) and worse pairing than just matching the first available prime!!
#ifdef NOT_USED_ANYMORE
int count_possible_pairs (char *primes, char *relocated, uint64_t numDsections, int numrels, int totrels, int multiplier, uint64_t k, int j)
{
	int	i, count, j_pair, work_j, work_j_pair;
	int64_t	bit, bit_pair;

	j_pair = numrels * 2 - 1 - j;
	bit = k * numrels * 2 + j;
	bit_pair = k * numrels * 2 + j_pair;
	work_j = j + (totrels - numrels);
	work_j_pair = j_pair + (totrels - numrels);

	count = 0;
	for (i = -(multiplier / 2 + 1); i <= (multiplier / 2 + 1); i++) {
		if (work_j + i * numrels * 2 < 0 || work_j + i * numrels * 2 >= totrels * 2) continue;
		if (bit + i * numrels * 4 < 0 || bit + i * numrels * 4 >= (int64_t) numDsections * numrels * 2) continue;
		if (bittst (primes, bit_pair + i * numrels * 4)) count++;
		else if (relocated != NULL && bittst (relocated, bit_pair + i * numrels * 4)) count++;
	}
	return (count);
}
#endif

/* Allocate and fill a bit array in such a way that it maximizes stage 2 prime pairings in ECM, P-1, and P+1.  Large B2 values will result */
/* in large bit arrays -- and if B2 is really large then the bit array must be created in chunks. */

int fill_work_bitarray (
	int	thread_num,		/* For outputting informative messages */
	void	**sieve_info,		/* Prime number sieve to use / initialize */
	int	D,			/* Calculated by best_stage2_impl, best D value ("big step") */
	int	totrels,		/* Calculated by best_stage2_impl, number of relative primes to use for pairing */
	uint64_t first_relocatable,	/* First relocatable prime (same as B1 unless bit arrays must be split or mem change caused a replan) */
	uint64_t last_relocatable,	/* End of relocatable primes (same as B2_start unless a change in available memory caused a new stage 2 plan) */
	uint64_t C_start,		/* First D section to place in the bit array (set to B2_start for first bit array) */
	uint64_t C,			/* Bound #2, end sieving for primes here */
	uint64_t maxDsections,		/* Maximum number of D sections to put in the bit array */
	char	**bitarray)		/* Pointer to work bit array pointer -- may be allocated here */
{
	int	D_index, numrels;
	uint64_t B2_end;		/* D section just after the last D section */
	uint64_t numDsections;		/* Number of D sections in the bit array */
	char	*primes, *relocated;
	uint64_t len;
	uint64_t prime, numprimes;
	char	*work;			/* Work bit array to fill in */
	uint64_t bits_set;		/* Bits set in the work array */
	int	multiplier, loop_count;
	struct relmap_data reldata;

/* Clear pointers to make error cleanup easier */

	primes = relocated = NULL;

/* Look up D in D_data so that we can access number of relative primes to D that are less than D */

	D_index = map_D_to_index (D);
	numrels = D_data[D_index].numrels;

/* Figure out how many D sections will be in this bit array, adjust B2_end and C if we're only processing some of the primes (a split bit array) */

	C_start = round_down_to_multiple_of (C_start - D / 2, D) + D / 2;
	B2_end = round_up_to_multiple_of (C + D / 2, D) - D / 2;
	numDsections = (B2_end - C_start) / D;
	if (numDsections > maxDsections) {
		numDsections = maxDsections;
		B2_end = C = C_start + numDsections * D;
	}

/* Allocate the work bit array if necessary */

	if (*bitarray == NULL) {
		len = divide_rounding_up (totrels * numDsections, 8);
		*bitarray = (char *) malloc ((size_t) len);
		if (*bitarray == NULL) goto oom;
	}

/* Allocate two bit arrays for the primes and relocated primes */

	len = divide_rounding_up (numDsections * numrels * 2, 8);
	primes = (char *) malloc ((size_t) len);
	if (primes == NULL) goto oom;
	relocated = (char *) malloc ((size_t) len);
	if (relocated == NULL) goto oom;

/* Set one bit (or more for relocatables) for each prime between C_start and C. */

	memset (primes, 0, (size_t) len);
	memset (relocated, 0, (size_t) len);
	start_sieve_with_limit (thread_num, first_relocatable, (uint32_t) sqrt((double) C), sieve_info);
	initmap (D, &reldata);
	numprimes = 0;
	while ((prime = sieve (*sieve_info)) <= C) {
		if (prime >= C_start) {
			numprimes++;
			bitset_prime (primes, prime);
		}
		else if (prime >= last_relocatable ||						// Prime is not one we care about
			 D_data[D_index].first_missing_prime * prime > B2_end) {		// No place to relocate prime - do in future bit array
			start_sieve_with_limit (thread_num, C_start, (uint32_t) sqrt((double) C), sieve_info);
			continue;
		}
		else if (D_data[D_index].second_missing_prime * prime > B2_end) {		// Only one place to relocate prime
			do {
				prime = prime * D_data[D_index].first_missing_prime;
			} while (prime < C_start);
			if (prime < B2_end) {
				numprimes++;
				bitset_prime (primes, prime);
			}
		} else {									// Mark all possible relocations
			int	i;
			int	set_a_bit = 0;
			for (i = D_data[D_index].first_missing_prime; i * prime < B2_end; i += 2) {
				if (map_is_relprime (i, D, &reldata) && i * prime >= C_start) {
					bitset_prime (relocated, i * prime);
					set_a_bit = 1;
				}
			}
			numprimes += set_a_bit;
		}
	}

/* Clear the work array */

	work = *bitarray;
	memset (work, 0, (size_t) divide_rounding_up (numDsections * totrels, 8));
	bits_set = 0;

/* Fill a work bit array in such a way that it maximizes stage 2 prime pairings in ECM, P-1, and P+1. */
/* Input is a bit array of non-relocatable primes and a bit array of relocatable primes. */

/* Pair primes making multiple passes over the input bit arrays.  We've tried several algorithms and may try several more! */
/* When loop_count is 1000 we've paired all the bits in the primes array leaving just relocatables -- pair any that can be paired. */
/* When loop_count is 1001 we're left with only unpairable relocatables.  Just output the first occurrence and we're done. */

// An example to help visualize what the code below is doing.  For D=30, numrels=4, multiplier=3:
// The 8 relprimes used for a "j" index are: -13,-11,-7,-1,1,7,11,13.  The 24 relprimes used for a "work_j" index that a particular
// multiple of D can access are: -43,-41,-37,-31,-29,-23,-19,-17,-13,-11,-7,-1,1,7,11,13,17,19,23,29,31,37,41,43
// Prime == 49 has three representations:
//	30 accesses -15 to 75,	relprime = 19 (work_j #17)		pair is 11	-16 bits in the prime bitmap
//	60 accesses 15 to 105,	relprime = -11 (work_j #9)		pair is 71	(the i=0 entry in the prime bitmap)
//	90 accesses 45 to 135,  relprime = -41 (work_j #1)		pair is 131	+16 bits in the prime bitmap

#define NO_MATCH	9999
	multiplier = totrels / numrels;
	for (loop_count = 0; loop_count <= 1001; loop_count++) {
		int	i, j, j_pair, work_j, work_j_pair, match, match_relocatable;
		uint64_t k, prime;
		int64_t	bit, bit_pair;

//GW - we could start at k = multiplier-1 and end at numDsections - multiplier (for a little cost in reduced pairs) but gain in reduced D multiplications
		for (k = 0; k < numDsections; k++) {				// loop through all sections of size D
			for (j = 0; j < numrels * 2; j++) {			// loop through each relative prime to find a possible pairing
				bit = k * numrels * 2 + j;

				// See if this entry in the prime bit array (or relocated array) needs to be paired up
				if (loop_count < 1000) {
					if (!bittst (primes, bit)) continue;
				} else {
					if (!bittst (relocated, bit)) continue;
				}

				// j and j_pair are zero-based indices for accessing the input prime bit arrays
				// for D=30 values are -13,-11,-7,-1,1,7,11,13
				j_pair = numrels * 2 - 1 - j;
				// work_j and work_j_pair are zero based indices to aid in generating output bit array
				// for D=30, multiplier=1, totrels=6 values are -19,-17,-13,-11,-7,-1,1,7,11,13,17,19
				work_j = j + (totrels - numrels);
				work_j_pair = j_pair + (totrels - numrels);
				// bit indexes into the input prime bit arrays
				bit_pair = k * numrels * 2 + j_pair;

// Loop over all possible ways this prime can be represented by multiple of D +/- a relative prime.  If any of these possible representations
// pair with a prime that is wonderful.  If multiple representations pair with a prime simply select the first.  We used to select the pairing that
// gives us the most flexibility for other primes to pair up, but that actually led to a slightly lower pairing percentage.

				match = NO_MATCH;	// no match yet

				// Find all the possible representations for this entry in the prime bit array
				// We only need to do this when looking for pairs from the primes array and we only need to look ahead.
				// Any possible matches against a smaller prime have already been matched with a different prime.
				if (loop_count < 1000)
			        for (i = 0; i <= (multiplier / 2 + 1); i++) {
					// Make sure this representation of the pairing prime is possible
					if (work_j - i * numrels * 2 < 0) continue;
					if (bit_pair + i * numrels * 4 >= (int64_t) numDsections * numrels * 2) continue;

					// See if the pair for this representation is also required in the output work bit array
					if (bittst (primes, bit_pair + i * numrels * 4)) {
						match = i;
						match_relocatable = FALSE;
						break;
					}
				}

				// Find all the other possible representations for this entry in the relocatable bit array
				if (match == NO_MATCH)
			        for (i = -(multiplier / 2 + 1); i <= (multiplier / 2 + 1); i++) {
					// Make sure this representation of the pairing prime is possible
					if (work_j - i * numrels * 2 < 0 || work_j - i * numrels * 2 >= totrels * 2) continue;
					if (bit_pair + i * numrels * 4 < 0 || bit_pair + i * numrels * 4 >= (int64_t) numDsections * numrels * 2) continue;

					// See if there is a pairing using the relocatable primes
					if (bittst (relocated, bit_pair + i * numrels * 4)) {
						match = i;
						match_relocatable = TRUE;
						break;
					}
				}

				// On first pass doing relocatables, only process matches.  A different multiple of the relocatable
				// prime may have a pairing.
				if (loop_count == 1000 && match == NO_MATCH) continue;

				// Clear bit from the primes array
				if (loop_count < 1000) {
//					Bit clear not necessary with our current algorithm that handles all primes in one-pass
//					bitclr (primes, bit);
				}
				// Clear bit from the relocated array for all possible relocations of the small prime
				else {
					ASSERTG (bittst (relocated, bit));
					prime = bit_prime (k, j);
					for (i = D_data[D_index].first_missing_prime; prime >= i * first_relocatable; i += 2) {
						if (map_is_relprime (i, D, &reldata)) while (prime % i == 0) prime = prime / i;
					}
					for (i = D_data[D_index].first_missing_prime; i * prime < B2_end; i += 2) {
						if (map_is_relprime (i, D, &reldata) && i * prime >= C_start) {
							ASSERTG (bittst_prime (relocated, i * prime));
							bitclr_prime (relocated, i * prime);
						}
					}
				}

				// Set the work output bit for this unpaired prime
				if (match == NO_MATCH) {
					bitset_work (work, k, work_j);
					bits_set++;
				}

				// Clear bit from the primes array for the matching prime of this pair
				else if (!match_relocatable) {
					bitset_work (work, k + match, work_j - match * numrels * 2);
					bits_set++;
					ASSERTG (bittst (primes, bit_pair + match * numrels * 4));
					bitclr (primes, bit_pair + match * numrels * 4);
				}

				// Clear bit from the relocated array for all possible relocations of the matching small prime
				else {						//matching an optional
					bitset_work (work, k + match, work_j - match * numrels * 2);
					bits_set++;
					ASSERTG (bittst (relocated, bit_pair + match * numrels * 4));
					prime = bit_prime (k, j_pair) + match * D * 2;
					for (i = D_data[D_index].first_missing_prime; prime >= i * first_relocatable; i += 2) {	// Reduce optional composite to underlying prime
						if (map_is_relprime (i, D, &reldata)) while (prime % i == 0) prime = prime / i;
					}
					ASSERTG (prime >= first_relocatable);
					for (i = D_data[D_index].first_missing_prime; i * prime < B2_end; i += 2) {
						if (map_is_relprime (i, D, &reldata) && i * prime >= C_start) {
							ASSERTG (bittst_prime (relocated, i * prime));
							bitclr_prime (relocated, i * prime);
						}
					}
				}
			}
		}
		// This rather odd loop_count test is a holdover from tests that did multiple passes looking for better pairing
		if (loop_count < 999) loop_count = 999;
	}

/* Free two temporary bit arrays */

	free (primes);
	free (relocated);

/* Output an informational message on our pairing efficiency */

	{
		char buf[120];
		sprintf (buf, "D: %d, relative primes: %d, stage 2 primes: %" PRIu64 ", pair%%=%5.2f\n",
			 D, totrels, numprimes, (double) (numprimes - bits_set) * 2.0 / (double) numprimes * 100.0);
		OutputStr (thread_num, buf);
	}

/* All done */

	return (0);

/* Out of memory exit */

oom:	free (primes);
	free (relocated);
	return (OutOfMemory (thread_num));
}

/* When a bit array completes, we know that all relocatable primes that could be relocated to that bit array are processed. */
/* Calculate the new first prime that needs relocating. */

uint64_t calc_new_first_relocatable (
	int	D,			/* Calculated by best_stage2_impl, best D value ("big step") */
	uint64_t C_done)		/* Bit arrays have been completed to this point */
{
	int	D_index = map_D_to_index (D);
	return (C_done / D_data[D_index].first_missing_prime);
}

/********************/
/* Utility routines */
/********************/

/* Test if N is a probable prime */
/* Compute i^(N-1) mod N for i = 3,5,7 */

int isProbablePrime (
	gwhandle *gwdata,
	giant	N)
{
	int	i, j, len, retval;
	gwnum	t1, t2;
	giant	x;

	if (isone (N)) return (TRUE);

	retval = TRUE;		/* Assume it is a probable prime */
	t1 = gwalloc (gwdata);
	len = bitlen (N);
	for (i = 3; retval && i <= 7; i += 2) {
		t2 = gwalloc (gwdata);
		dbltogw (gwdata, (double) 1.0, t1);
		dbltogw (gwdata, (double) i, t2);
		gwfft (gwdata, t2, t2);
		for (j = 1; j <= len; j++) {
			gwsquare (gwdata, t1);
			if (bitval (N, len-j)) gwfftmul (gwdata, t2, t1);
		}
		gwfree (gwdata, t2);
		x = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 10);
		if (gwtogiant (gwdata, t1, x)) retval = FALSE;	/* Technically, prime status is unknown on an unexpected error */
		else {
			modgi (&gwdata->gdata, N, x);
			iaddg (-i, x);
			if (!isZero (x)) retval = FALSE;	/* Not a prime */
		}
		pushg (&gwdata->gdata, 1);
	}
	gwfree (gwdata, t1);
	return (retval);
}

// From Alex Kruppa, master of all things ECM, the following formula computes the value of a curve when using B2 values that are not 100 * B1.
// curve_worth = 0.11 + 0.89 * (log10(B2 / B1) / 2) ^ 1.5
// B2 = B1 * 10 ^ ((((curve_worth - 0.11) / 0.89) ^ (1 / 1.5)) * 2)

#define kruppa_adjust_ratio(B2B1ratio)	(0.11 + 0.89 * pow (_log10(B2B1ratio) / 2.0, 1.5))
#define kruppa_adjust(B2,B1)		kruppa_adjust_ratio ((double)(B2) / (double)(B1))
#define kruppa_unadjust(worth,B1)	(uint64_t) round((B1) * pow (10.0, (pow (((worth) - 0.11) / 0.89, 1.0 / 1.5) * 2)))

/*************************************************/
/* ECM structures and setup/termination routines */
/*************************************************/

/* We used to compute ell_dbl and ell_add by first computing (x+z) and (x-z).  This structure gave as a place to keep track if addsub4 had been called. */

struct xz {
	gwnum	x;		/* x or FFT(x) */
	gwnum	z;		/* z or FFT(z) */
};

/* Data maintained during ECM process */

#define POOL_3MULT		1	/* Modinv algorithm that takes 3 multiplies (9 FFTs) using 100% more gwnums */
#define POOL_3POINT44MULT	2	/* Modinv algorithm that takes 3.444 multiplies (10.333 FFTs) using 33% more gwnums */
#define POOL_3POINT57MULT	3	/* Modinv algorithm that takes 3.573 multiplies (10.714 FFTs) using 14% more gwnums */
#define POOL_N_SQUARED		4	/* Use O(N^2) multiplies modinv algorithm using no extra gwnums */

#define ECM_STATE_STAGE1_INIT		0	/* Selecting sigma for curve */
#define ECM_STATE_STAGE1		1	/* In middle of stage 1 */
#define ECM_STATE_MIDSTAGE		2	/* Stage 2 initialization for the first time */
#define ECM_STATE_STAGE2		3	/* In middle of stage 2 (processing a bit array) */
#define ECM_STATE_GCD			4	/* Stage 2 GCD */

typedef struct {
	gwhandle gwdata;	/* GWNUM handle */
	int	thread_num;	/* Worker thread number */
	struct work_unit *w;	/* Worktodo.txt entry */
	unsigned long curve;	/* Curve # starting with 1 */
	int	state;		/* Curve state defined above */
	double	sigma;		/* Sigma for the current curve */
	uint64_t B;		/* Bound #1 (a.k.a. B1) */
	uint64_t C;		/* Bound #2 (a.k.a. B2) */
	int	optimal_B2;	/* TRUE if we calculate optimal bound #2 given currently available memory.  FALSE for a fixed bound #2. */
	uint64_t average_B2;	/* Average Kruppa-adjusted bound #2 work done on ECM curves thusfar */

	gwnum	Ad4;		/* Pre-computed value used for doubling */
	readSaveFileState read_save_file_state;	/* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	void	*sieve_info;	/* Prime number sieve */
	uint64_t stage1_prime;	/* Prime number being processed */

	struct xz xz;		/* The stage 1 value being computed */
	gwnum	gg;		/* The stage 2 accumulated value */

	int	pool_type;	/* Modinv algorithm type to use */
	int	pool_count;	/* Count values in the modinv pool */
	int	poolz_count;	/* Count values in the modinv poolz */
	gwnum	pool_modinv_value;/* Value we will eventually do a modinv on */
	gwnum	*pool_values;	/* Array of values to normalize */
	gwnum	*poolz_values;	/* Array of z values we are normalizing */
	unsigned long modinv_count; /* Stats - count of modinv calls */

	int	D;		/* Stage 2 loop increment */
	int	E;		/* Number of mQx values to pool together into one stage 2 modular inverse */
	int	totrels;	/* Number relatively prime nQx values used */
	uint64_t B2_start;	/* Starting point of first D section to be processed in stage 2 (an odd multiple of D/2) */
	uint64_t numDsections;	/* Number of D sections to process in stage 2 */
	uint64_t Dsection;	/* Current D section being processed in stage 2 */

	char	*bitarray;	/* Bit array for prime pairings in each D section */
	uint64_t bitarraymaxDsections;	/* Maximum number of D sections per bit array */
	uint64_t bitarrayfirstDsection; /* First D section in the bit array */
	uint64_t first_relocatable; /* First relocatable prime (same as B1 unless bit arrays must be split or mem change caused a replan) */
	uint64_t last_relocatable; /* Last relocatable prime for filling bit arrays (unless mem change causes a replan) */
	uint64_t C_done;	/* Stage 2 completed thusfar (updated every D section that is completed) */

	int	stage2_numvals;	/* Number of gwnums used in stage 2 */
	int	TWO_FFT_STAGE2;	/* Type of ECM stage 2 to execute */
	gwnum	*nQx;		/* Array of relative primes data used in stage 2 */
	struct xz Qm, Qprevm, QD; /* Values used to calculate successive D values in stage 2 */
	int	Qm_state;	/* For 4-FFT continuation, TRUE if Qm has been returned by mQ_next */
	gwnum	*mQx;		/* Array of calculated D values when modular inverse pooling is active */
	int	mQx_count;	/* Number of values remaining in the mQx array */
	int	mQx_retcount;	/* Next mQx array entry to return */

	double	pct_mem_to_use;	/* If we get memory allocation errors, we progressively try using less and less. */
} ecmhandle;

/* Forward declarations */

void normalize_pool_term (ecmhandle *ecmdata);
void mQ_term (ecmhandle *ecmdata);

/* Perform cleanup functions */

void ecm_cleanup (
	ecmhandle *ecmdata)
{
	normalize_pool_term (ecmdata);
	mQ_term (ecmdata);
	free (ecmdata->nQx), ecmdata->nQx = NULL;
	free (ecmdata->bitarray), ecmdata->bitarray = NULL;
	gwdone (&ecmdata->gwdata);
	end_sieve (ecmdata->sieve_info), ecmdata->sieve_info = NULL;
}

/**************************************************************
 *	ECM Functions
 **************************************************************/

/* This routine initializes an xz pair with two allocated gwnums */

__inline int alloc_xz (			/* Returns TRUE if successful */
	ecmhandle *ecmdata,
	struct xz *arg)
{
	arg->x = gwalloc (&ecmdata->gwdata);
	if (arg->x == NULL) return (FALSE);
	arg->z = gwalloc (&ecmdata->gwdata);
	if (arg->z == NULL) return (FALSE);
	return (TRUE);
}

/* This routine cleans up an xz pair with two allocated gwnums */

__inline void free_xz (
	ecmhandle *ecmdata,
	struct xz *arg)
{
	gwfree (&ecmdata->gwdata, arg->x); arg->x = NULL;
	gwfree (&ecmdata->gwdata, arg->z); arg->z = NULL;
}

/* Macro to swap to xz structs */

#define xzswap(a,b)	{ struct xz t; t = a; a = b; b = t; }

/* Computes 2P=(out.x:out.z) from P=(in.x:in.z), uses the global variable Ad4. */
/* Input arguments may be in FFTed state.  Out argument can be same as input argument. */
/* Scratch xz argument can equal in but cannot equal out. */

void ell_dbl_xz_scr (
	ecmhandle *ecmdata,
	struct xz *in,		/* Input value to double */
	struct xz *out,		/* Output value */
	struct xz *scr)		// Scratch registers (only the .x gwnum is used)
{				/* 10 or 11 FFTs, 4 adds */
	gwnum	t2, t3, t4;

	ASSERTG (scr != out);

	/* If we have extra_bits to square the output of a gwadd, then we can use an algorithm that minimizes FFTs. */
	if (square_safe (&ecmdata->gwdata, 1)) {
		t3 = scr->x;
		t2 = out->x;
		t4 = out->z;

		gwsetmulbyconst (&ecmdata->gwdata, 4);
		gwmul3 (&ecmdata->gwdata, in->x, in->z, t3, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); /* t3 = 4*x*z */
		gwsub3o (&ecmdata->gwdata, in->x, in->z, t2, GWADD_DELAY_NORMALIZE);			/* Compute x - z */
		gwsquare2 (&ecmdata->gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = (x - z)^2 */
		gwmul3 (&ecmdata->gwdata, t2, ecmdata->Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	/* t4 = t2 * Ad4 */
		gwaddmul4 (&ecmdata->gwdata, t2, t3, t4, out->x, GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT); /* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&ecmdata->gwdata, t4, t3, t3, out->z, GWMUL_STARTNEXTFFT);			/* outz = (t4 + t3) * t3 */
	}

	/* This algorithm uses an extra FFT to assure that it will work even if there are no extra bits.  We've made a conscious decision to */
	/* spend an extra FFT to make ell_add as clean as possible.  Ell_adds are ten times more common than ell_dbls.  This version lets */
	/* ell_add FFT x and z without any worries. */
	else {
		int	t1_safe = addmul_safe (&ecmdata->gwdata, 2,0,0);

		t2 = out->x;
		t3 = out->z;
		t4 = scr->x;

		gwsquare2 (&ecmdata->gwdata, in->x, scr->x, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(t1_safe));	/* x^2 */
		gwsquare2 (&ecmdata->gwdata, in->z, scr->z, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(t1_safe));	/* z^2 */
		gwsetmulbyconst (&ecmdata->gwdata, 2);
		gwmul3 (&ecmdata->gwdata, in->x, in->z, t3, GWMUL_MULBYCONST);					/* t3 = 2*x*z */

		gwadd3o (&ecmdata->gwdata, scr->x, scr->z, t2, GWADD_DELAY_NORMALIZE);				/* x^2 + z^2 */
		gwsub3o (&ecmdata->gwdata, t2, t3, t2, GWADD_DELAYNORM_IF(t1_safe));				/* t2 = x^2 - 2xz + z^2 */
		gwadd3o (&ecmdata->gwdata, t3, t3, t3, GWADD_FORCE_NORMALIZE);					/* t3 = 4*x*z */

		gwmul3 (&ecmdata->gwdata, t2, ecmdata->Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);		/* t4 = t2 * Ad4 */
		gwaddmul4 (&ecmdata->gwdata, t2, t3, t4, out->x, GWMUL_FFT_S2 | GWMUL_FFT_S3 | GWMUL_STARTNEXTFFT); /* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&ecmdata->gwdata, t4, t3, t3, out->z, GWMUL_STARTNEXTFFT);				/* outz = (t4 + t3) * t3 */
	}
}

/* Like ell_dbl_xz_scr, but the output arguments are not partially FFTed.  The input argument is assumed to never be used again. */

void ell_dbl_xz_scr_last (
	ecmhandle *ecmdata,
	struct xz *in,		/* Input value to double */
	struct xz *out,		/* Output value */
	struct xz *scr)		// Scratch registers (only the .x gwnum is used)
{				/* 10 or 11 FFTs, 4 adds */
	gwnum	t2, t3, t4;

	ASSERTG (scr != out);

	/* If we have extra_bits to square the output of a gwadd, then we can use an algorithm that minimizes FFTs. */
	if (square_safe (&ecmdata->gwdata, 1)) {
		t3 = scr->x;
		t2 = out->x;
		t4 = out->z;

		gwsetmulbyconst (&ecmdata->gwdata, 4);
		gwmul3 (&ecmdata->gwdata, in->x, in->z, t3, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_MULBYCONST | GWMUL_STARTNEXTFFT); /* t3 = 4*x*z */
		gwsub3o (&ecmdata->gwdata, in->x, in->z, t2, GWADD_DELAY_NORMALIZE);			/* Compute x - z */
		gwsquare2 (&ecmdata->gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = (x - z)^2 */
		gwmul3 (&ecmdata->gwdata, t2, ecmdata->Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	/* t4 = t2 * Ad4 */
		gwaddmul4 (&ecmdata->gwdata, t2, t3, t4, out->x, GWMUL_FFT_S2 | GWMUL_FFT_S3);		/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&ecmdata->gwdata, t4, t3, t3, out->z, 0);					/* outz = (t4 + t3) * t3 */
	}

	/* This algorithm uses an extra FFT to assure that it will work even if there are no extra bits.  We've made a conscious decision to */
	/* spend an extra FFT to make ell_add as clean as possible.  Ell_adds are ten times more common than ell_dbls.  This version lets */
	/* ell_add FFT x and z without any worries. */
	else {
		int	t1_safe = addmul_safe (&ecmdata->gwdata, 2,0,0);

		t2 = out->x;
		t3 = out->z;
		t4 = scr->x;

		gwsquare2 (&ecmdata->gwdata, in->x, scr->x, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(t1_safe));	/* x^2 */
		gwsquare2 (&ecmdata->gwdata, in->z, scr->z, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT_IF(t1_safe));	/* z^2 */
		gwsetmulbyconst (&ecmdata->gwdata, 2);
		gwmul3 (&ecmdata->gwdata, in->x, in->z, t3, GWMUL_MULBYCONST);					/* t3 = 2*x*z */

		gwadd3o (&ecmdata->gwdata, scr->x, scr->z, t2, GWADD_DELAY_NORMALIZE);				/* x^2 + z^2 */
		gwsub3o (&ecmdata->gwdata, t2, t3, t2, GWADD_DELAYNORM_IF(t1_safe));				/* t2 = x^2 - 2xz + z^2 */
		gwadd3o (&ecmdata->gwdata, t3, t3, t3, GWADD_FORCE_NORMALIZE);					/* t3 = 4*x*z */

		gwmul3 (&ecmdata->gwdata, t2, ecmdata->Ad4, t4, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);		/* t4 = t2 * Ad4 */
		gwaddmul4 (&ecmdata->gwdata, t2, t3, t4, out->x, GWMUL_FFT_S2 | GWMUL_FFT_S3);			/* outx = (t1 = t2 + t3) * t4 */
		gwaddmul4 (&ecmdata->gwdata, t4, t3, t3, out->z, 0);						/* outz = (t4 + t3) * t3 */
	}
}

/* Adds Q=(in2.x:in2.z) and R=(in1.x:in1.z) and puts the result in (out.x:out.z).  Assumes that Q-R=P or R-Q=P where P=(diff.x:diff.z). */
/* Input arguments may be in FFTed format.  Out argument can be same as any of the 3 input arguments. */
/* Scratch xz argument cannot equal in1, in2, or diff. */

void ell_add_xz_scr (
	ecmhandle *ecmdata,
	struct xz *in1,
	struct xz *in2,
	struct xz *diff,
	struct xz *out,
	struct xz *scr)		// Scratch registers
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;
	int	options;

	ASSERTG (scr != in1 && scr != in2 && scr != diff);
	t1 = scr->z;
	t2 = scr->x;
	gwmulmulsub5 (&ecmdata->gwdata, in1->x, in2->z, in1->z, in2->x, t1, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	/* t1 = x1z2 - z1x2 */
	gwmulmulsub5 (&ecmdata->gwdata, in1->x, in2->x, in1->z, in2->z, t2, GWMUL_STARTNEXTFFT);				/* t2 = x1x2 - z1z2 */
	gwsquare2 (&ecmdata->gwdata, t1, t1, GWMUL_STARTNEXTFFT);				/* t1 = t1^2 */
	gwsquare2 (&ecmdata->gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = t2^2 */
	options = (diff == out) ? GWMUL_STARTNEXTFFT : GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT;
	gwmul3 (&ecmdata->gwdata, t2, diff->z, t2, options);					/* t2 = t2 * zdiff (will become outx) */
	gwmul3 (&ecmdata->gwdata, t1, diff->x, t1, options);					/* t1 = t1 * xdiff (will become outz) */
	if (out != scr) xzswap (*scr, *out);
}

// Like ell_add_xz_scr except that out is not equal any input (including diff) so that scratch register can be inferred.
// This is simply a shortcut for better readability.

#define ell_add_xz(h,i1,i2,dif,o)	ell_add_xz_scr(h,i1,i2,dif,o,o);

// Like ell_add_xz_scr except that a scratch register is allocated.

int ell_add_xz_noscr (
	ecmhandle *ecmdata,
	struct xz *in1,
	struct xz *in2,
	struct xz *diff,
	struct xz *out)
{
	struct xz scr;
	if (!alloc_xz (ecmdata, &scr)) return (OutOfMemory (ecmdata->thread_num));
	ell_add_xz_scr (ecmdata, in1, in2, diff, out, &scr);
	free_xz (ecmdata, &scr);
	return (0);
}

/* Like ell_add_xz_scr but in1, in2 and diff are assumed to not be used again.  Out argument never has its forward FFT begun. */

void ell_add_xz_last (
	ecmhandle *ecmdata,
	struct xz *in1,
	struct xz *in2,
	struct xz *diff,
	struct xz *out,
	struct xz *scr)
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;

	ASSERTG (scr != in1 && scr != in2 && scr != diff);
	t1 = scr->z;
	t2 = scr->x;
	gwmulmulsub5 (&ecmdata->gwdata, in1->x, in2->z, in1->z, in2->x, t1, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	/* t1 = x1z2 - z1x2 */
	gwmulmulsub5 (&ecmdata->gwdata, in1->x, in2->x, in1->z, in2->z, t2, GWMUL_STARTNEXTFFT);				/* t2 = x1x2 - z1z2 */
	gwsquare2 (&ecmdata->gwdata, t1, t1, GWMUL_STARTNEXTFFT);				/* t1 = t1^2 */
	gwsquare2 (&ecmdata->gwdata, t2, t2, GWMUL_STARTNEXTFFT);				/* t2 = t2^2 */
	gwmul3 (&ecmdata->gwdata, t2, diff->z, t2, 0);						/* t2 = t2 * zdiff (will become outx) */
	gwmul3 (&ecmdata->gwdata, t1, diff->x, t1, 0);						/* t1 = t1 * xdiff (will become outz) */
	if (out != scr) xzswap (*scr, *out);
}


/* Perform an elliptic multiply using an algorithm developed by Peter Montgomery.  Basically, we try to find a near optimal Lucas */
/* chain of additions that generates the number we are multiplying by.  This minimizes the number of calls to ell_dbl and ell_add. */

/* The costing function assigns an ell_dbl call a cost of 10 and an ell_add call a cost of 12. */
/* This cost estimates the number of forward and inverse transforms performed. */

#define swap(a,b)	{uint64_t t=a;a=b;b=t;}

int lucas_cost (
	uint64_t n,
	uint64_t d)
{
	uint64_t e;//, dmod3, emod3;
	unsigned long c;

	if (d >= n || d <= n/2) return (999999999);		/* Catch invalid costings */

	c = 0;
	while (n != 1) {
	    e = n - d;
	    d = d - e;

	    c += 12;

	    while (d != e) {
		if (d < e) {
			swap (d,e);
		}
//		if (d <= e + (e >> 2)) {
//			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
//				t = d;
//				d = (d+d-e)/3;
//				e = (e+e-t)/3;
//				c += 36;
//				continue;
//			}
//			if (dmod3 == emod3 && (d&1) == (e&1)) {
//				d = (d-e) >> 1;
//				c += 22;
//				continue;
//			}
//		}
		if (100 * d <= 296 * e) {
			d = d-e;
			c += 12;
		} else if ((d&1) == (e&1)) {
			d = (d-e) >> 1;
			c += 22;
		} else if ((d&1) == 0) {
			d = d >> 1;
			c += 22;
//		} else if ((dmod3 = d%3) == 0) {
//			d = d/3-e;
//			c += 46;
//		} else if (dmod3 == 3 - (emod3 = e%3)) {
//			d = (d-e-e)/3;
//			c += 46;
//		} else if (dmod3 == emod3) {
//			d = (d-e)/3;
//			c += 46;
		} else {
			e = e >> 1;
			c += 22;
		}
	    }
	    c += 10;
	    if (d == 1) break;
	    n = d;
	    d = (uint64_t) ((double) n * 0.6180339887498948);
	}

	return (c);
}

int lucas_mul (
	ecmhandle *ecmdata,
	struct xz *A,
	uint64_t n,
	uint64_t d,
	int	last_mul)	// TRUE if the last multiply should not use the GWMUL_STARTNEXTFFT option
{
	uint64_t e;//, dmod3, emod3;
	struct xz B, C, T;//, S

	if (!alloc_xz (ecmdata, &B)) goto oom;
	if (!alloc_xz (ecmdata, &C)) goto oom;
//      if (!alloc_xz (ecmdata, &S)) goto oom;
	if (!alloc_xz (ecmdata, &T)) goto oom;

	while (n != 1) {
	    ell_dbl_xz_scr (ecmdata, A, &B, &C);				/* B = 2*A, scratch reg = C */
										/* C = A (but we delay setting that up) */

	    e = n - d;
	    d = d - e;

	    // To save two gwcopies setting C=A, we handle the most common case for the first iteration of the following loop.
	    // I've only seen three cases that end up doing the gwcopies, n=3, n=11 and n=17.  With change to 2.96 there are a couple more.
	    if (e > d && 100 * e <= 296 * d) {
		swap (d, e);
		xzswap (*A, B);							/* swap A & B, thus diff C = B */
		ell_add_xz (ecmdata, A, &B, &B, &C);				/* B = A+B */
		xzswap (B, C);							/* C = B */
		d = d-e;
	    }
	    else if (d > e && 100 * d <= 296 * e) {
		ell_add_xz (ecmdata, A, &B, A, &C);				/* B = A+B */
		xzswap (B, C);							/* C = B */
		d = d-e;
	    } else {
		gwcopy (&ecmdata->gwdata, A->x, C.x);
		gwcopy (&ecmdata->gwdata, A->z, C.z);				/* C = A */
	    }

	    while (d != e) {
		if (d < e) {
			swap (d, e);
			xzswap (*A, B);
		}
		// These cases were in Peter Montgomery's original PRAC implementation.  I've found that we do fewer FFTs with
		// these cases removed.  Should they be added back in the if statements above for gwcopies and in lucas_cost need to be restored.
		// SHOULD RE-EVALUATE with addition of d > 4 * e check (where 4 may not be optimal)
//		if (d <= e + (e >> 2)) {
//			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
//				ell_add_xz (ecmdata, A, &B, &C, &S);		/* S = A+B */
//				ell_add_xz (ecmdata, A, &S, &B, &T);		/* T = A+S */
//				ell_add_xz_scr (ecmdata, &S, &B, A, &B, &T);	/* B = B+S, scratch reg = T */
//				xzswap (T, *A);					/* A = T */
//				t = d;
//				d = (d+d-e)/3;
//				e = (e+e-t)/3;
//				continue;
//			}
//			if (dmod3 == emod3 && (d&1) == (e&1)) {
//				ell_add_xz_scr (ecmdata, A, &B, &C, &B, &T);	/* B = A+B, scratch reg = T */
//				ell_dbl_xz_scr (ecmdata, A, A, &T);		/* A = 2*A, scratch reg = T */
//				d = (d-e) >> 1;
//				continue;
//			}
//		}
		if (100 * d <= 296 * e) {					/* d <= 2.96 * e (Montgomery used 4.00) */
			ell_add_xz_scr (ecmdata, A, &B, &C, &C, &T);		/* B = A+B, scratch reg = T */
			xzswap (B, C);						/* C = B */
			d = d-e;
		} else if ((d&1) == (e&1)) {
			ell_add_xz_scr (ecmdata, A, &B, &C, &B, &T);		/* B = A+B, scratch reg = T */
			ell_dbl_xz_scr (ecmdata, A, A, &T);			/* A = 2*A, scratch reg = T */
			d = (d-e) >> 1;
		} else if ((d&1) == 0) {
			ell_add_xz_scr (ecmdata, A, &C, &B, &C, &T);		/* C = A+C, scratch reg = T */
			ell_dbl_xz_scr (ecmdata, A, A, &T);			/* A = 2*A, scratch reg = T */
			d = d >> 1;
		}
		// These cases were in Peter Montgomery's original PRAC implementation.  I've found that optimal addition chains
		// rarely (never?) use these rules.  Should they be added back lucas_cost needs to be restored too.
		// SHOULD RE-EVALUATE with addition of d > 4 * e check (where 4 may not be optimal)
//		else if ((dmod3 = d%3) == 0) {
//			ell_dbl_xz_scr (ecmdata, A, &S, &T);			/* S = 2*A, scratch reg = T */
//			ell_add_xz (ecmdata, A, &B, &C, &T);			/* T = A+B */
//			ell_add_xz_scr (ecmdata, &S, &T, &C, &C, &T);		/* B = S+T, scratch reg = T */
//			ell_add_xz_scr (ecmdata, &S, A, A, A, &S);		/* A = S+A, scratch reg = S */
//			xzswap (B, C);						/* C = B */
//			d = d/3-e;
//		} else if (dmod3 == 3 - (emod3 = e%3)) {
//			ell_add_xz (ecmdata, A, &B, &C, &S);			/* S = A+B */
//			ell_add_xz_scr (ecmdata, A, &S, &B, &B, &T);		/* B = A+S, scratch reg = T */
//			ell_dbl_xz_scr (ecmdata, A, &S, &T);			/* S = 2*A, scratch reg = T */
//			ell_add_xz_scr (ecmdata, &S, A, A, A, &T);		/* A = S+A, scratch reg = T */
//			d = (d-e-e)/3;
//		} else if (dmod3 == emod3) {
//			ell_add_xz (ecmdata, A, &B, &C, &T);			/* T = A+B */
//			ell_add_xz_scr (ecmdata, A, &C, &B, &C, &S);		/* C = A+C, scratch reg = S */
//			xzswap (T, B);						/* B = T */
//			ell_dbl_xz_scr (ecmdata, A, &S, &T);			/* S = 2*A, scratch reg = T */
//			ell_add_xz_scr (ecmdata, &S, A, A, A, &T);		/* A = S+A, scratch reg = T */
//			d = (d-e)/3;
//		}
		else {
			ell_add_xz_scr (ecmdata, &B, &C, A, &C, &T);		/* C = C-B, scratch reg = T */
			ell_dbl_xz_scr (ecmdata, &B, &B, &T);			/* B = 2*B, scratch reg = T */
			e = e >> 1;
		}
	    }

	    if (d == 1 && last_mul) {
		ell_add_xz_last (ecmdata, &B, A, &C, A, &T);			/* A = A+B, scratch reg = T */
		break;
	    } else {
		ell_add_xz_scr (ecmdata, &B, A, &C, A, &T);			/* A = A+B, scratch reg = T */
	    }

	    n = d;
	    d = (uint64_t) ((double) n * 0.6180339887498948);
	}
	free_xz (ecmdata, &B);
	free_xz (ecmdata, &C);
//      free_xz (ecmdata, &S);
	free_xz (ecmdata, &T);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
#undef swap

/* Try a series of Lucas chains to find the cheapest. */
/* First try v = (1+sqrt(5))/2, then (2+v)/(1+v), then (3+2*v)/(2+v), */
/* then (5+3*v)/(3+2*v), etc.  Finally, execute the cheapest. */

__inline int lucas_cost_several (uint64_t n, uint64_t *d) {
	int	i, c, min;
	uint64_t testd;
	for (i = 0, testd = *d - PRAC_SEARCH / 2; i < PRAC_SEARCH; i++, testd++) {
		c = lucas_cost (n, testd);
		if (i == 0 || c < min) min = c, *d = testd;
	}
	return (min);
}

int ell_mul (
	ecmhandle *ecmdata,
	struct xz *arg,
	uint64_t n,
	int	last_mul)	// TRUE if the this is the last mul in a series and we do not want to apply the GWMUL_STARTNEXTFFT option to the result
{
	unsigned long zeros;
	int	stop_reason;

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		int	c, min;
		uint64_t d, mind;

		mind = (uint64_t) ceil((double) 0.6180339887498948 * n);		/*v=(1+sqrt(5))/2*/
		min = lucas_cost_several (n, &mind);

		d = (uint64_t) ceil ((double) 0.7236067977499790 * n);			/*(2+v)/(1+v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.5801787282954641 * n);			/*(3+2*v)/(2+v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6328398060887063 * n);			/*(5+3*v)/(3+2*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6124299495094950 * n);			/*(8+5*v)/(5+3*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6201819808074158 * n);			/*(13+8*v)/(8+5*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6172146165344039 * n);			/*(21+13*v)/(13+8*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6183471196562281 * n);			/*(34+21*v)/(21+13*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6179144065288179 * n);			/*(55+34*v)/(34+21*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6180796684698958 * n);			/*(89+55*v)/(55+34*v)*/
		c = lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		stop_reason = lucas_mul (ecmdata, arg, n, mind, zeros == 0 && last_mul);
		if (stop_reason) return (stop_reason);
	}
	while (zeros--) {
		struct xz scr;
		if (!alloc_xz (ecmdata, &scr)) return (OutOfMemory (ecmdata->thread_num));
		if (zeros || !last_mul) ell_dbl_xz_scr (ecmdata, arg, arg, &scr);
		else ell_dbl_xz_scr_last (ecmdata, arg, arg, &scr);
		free_xz (ecmdata, &scr);
	}
	return (0);
}

/* Test if factor divides N, return TRUE if it does */

int testFactor (
	gwhandle *gwdata,
	struct work_unit *w,
	giant	f)		/* Factor to test */
{
	giant	tmp;
	int	divides_ok;

/* See if this is a valid factor */

	tmp = popg (&gwdata->gdata, f->sign + 5);	/* Allow room for mul by KARG */
	itog (w->b, tmp);
	powermod (tmp, w->n, f);
	dblmulg (w->k, tmp);
	iaddg (w->c, tmp);
	modgi (&gwdata->gdata, f, tmp);
	divides_ok = isZero (tmp);
	pushg (&gwdata->gdata, 1);
	if (!divides_ok) return (FALSE);

/* If QAing, see if we found the expected factor */

	if (QA_IN_PROGRESS) {
		tmp = popg (&gwdata->gdata, f->sign + 5);
		gtog (f, tmp);
		modg (QA_FACTOR, tmp);
		divides_ok = isZero (tmp);
		pushg (&gwdata->gdata, 1);
		if (!divides_ok) {
			char	buf[200];
			strcpy (buf, "ERROR: Factor not found. Expected ");
			gtoc (QA_FACTOR, buf+strlen(buf), 150);
			strcat (buf, "\n");
			OutputBoth (MAIN_THREAD_NUM, buf);
		}
	}

/* All done, return success */

	return (TRUE);
}

/* Set N, the number we are trying to factor */

int setN (
	int	thread_num,
	struct work_unit *w,
	giant	*N)		/* k*b^n+c as a giant */
{
	unsigned long bits, p;
	char	buf[2500];

/* Create the binary representation of the number we are factoring */
/* Allocate 5 extra words to handle any possible k value. */

	bits = (unsigned long) (w->n * _log2 (w->b));
	*N = allocgiant ((bits >> 5) + 5);
	if (*N == NULL) return (OutOfMemory (thread_num));

/* This special code comes from Serge Batalov */

	if (IniGetInt (INI_FILE, "PhiExtensions", 0) &&
	    w->k == 1.0 && w->b == 2 && w->c == -1) {		/*=== this input means Phi(n,2) with n semiprime ===*/
		unsigned int i,k,q,knownSmallMers[] = {3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
						       4253, 4423, 9689, 9941, 11213, 19937, 999999999}; /* for now, just the cases where w->n = p * q, and 2^q-1 is prime */
		for (i=0; (q=knownSmallMers[i]) < w->n || q*q <= w->n; i++) if ((w->n%q) == 0) {
			giant tmp = allocgiant ((bits >> 5) + 5);
			if (!tmp) return (OutOfMemory (thread_num));
			ultog (1, tmp);
			ultog (1, *N);
			gshiftleft (w->n-w->n/q, *N);
			for (k=2; k < q; k++) {
				gshiftleft (w->n/q, tmp);
				addg (tmp, *N);
			}
			iaddg (1, *N);
			if (q != w->n/q) {
				ultog (w->b, tmp);
				power (tmp, w->n/q);
				iaddg (w->c, tmp);
				divg (tmp, *N);
			}
			free (tmp);
			/* w->minimum_fftlen = w->n; */ /*=== too late to do this here. Moved before gwsetup() --SB. */
			if (!w->known_factors || !strcmp (w->known_factors, "1")) {
				p = sprintf (buf, "M%lu", w->n/q); if(q != w->n/q) p += sprintf (buf+p, "/M%d", q);
				w->known_factors = (char *) malloc (p+1);
				memcpy (w->known_factors, buf, p+1);
			}
			return (0);
		}
        }

	if (IniGetInt (INI_FILE, "PhiExtensions", 0) &&
	    w->k == 1.0 && labs(w->c) == 1 && (w->n%3) == 0) {		/*=== this input means Phi(3,-b^(n/3)) ===*/
		giant	tmp = allocgiant ((bits >> 5) + 5);
		if (tmp == NULL) return (OutOfMemory (thread_num));
		ultog (w->b, tmp);
		power (tmp, w->n/3);
		gtog (tmp, *N);
		squareg (*N);
		if (w->c == 1) subg (tmp, *N); else addg (tmp, *N);
		iaddg (1, *N);
		free (tmp);
		/* w->minimum_fftlen = w->n; */ /*=== too late to do this here. Moved before gwsetup() --SB. */
		if (!w->known_factors) {
			p = sprintf (buf, "(%lu^%lu%+ld)", w->b, w->n/3, w->c);
			w->known_factors = (char *) malloc (p+1);
			memcpy (w->known_factors, buf, p+1);
		}
		return (0);
	}

/* Standard code for working on k*b^n+c */

	ultog (w->b, *N);
	power (*N, w->n);
	dblmulg (w->k, *N);
	iaddg (w->c, *N);

/* If we have a list of known factors then process it */

	if (w->known_factors != NULL) {
		char	*comma, *p;
		giant	tmp, f;

		tmp = allocgiant ((bits >> 5) + 5);
		if (tmp == NULL) return (OutOfMemory (thread_num));
		f = allocgiant ((bits >> 5) + 5);
		if (f == NULL) return (OutOfMemory (thread_num));

/* Process each factor */

		for (p = w->known_factors; ; ) {

/* Get the factor - raise error is it is less than or equal to one */

			ctog (p, f);
			if (gsign (f) < 1 || isone (f)) {
				char	msg[1000], kbnc[80];
				strcpy (buf, p);
				buf[10] = 0;
				gw_as_string (kbnc, w->k, w->b, w->n, w->c);
				sprintf (msg, "Error parsing known factors of %s near: '%s'\n", kbnc, buf);
				OutputStr (thread_num, msg);
				free (f);
				free (tmp);
				deleteWorkToDoLine (thread_num, w, FALSE);
				return (STOP_ABORT);
			}

/* Divide N by factor - then verify the factor */

			gtog (*N, tmp);
			divg (f, tmp);
			mulg (tmp, f);
			if (gcompg (f, *N)) {
				char	kbnc[80];
				strcpy (buf, p);
				comma = strchr (buf, ',');
				if (comma != NULL) *comma = 0;
				gw_as_string (kbnc, w->k, w->b, w->n, w->c);
				sprintf (buf+strlen(buf), " does not divide %s\n", kbnc);
				OutputBoth (thread_num, buf);
				free (f);
				free (tmp);
				deleteWorkToDoLine (thread_num, w, FALSE);
				return (STOP_ABORT);
			}
			gtog (tmp, *N);

/* Skip to next factor in list */

			comma = strchr (p, ',');
			if (comma == NULL) break;
			p = comma + 1;
		}
		free (f);
		free (tmp);
		return (0);
	}

/* Return success */

	return (0);
}

/* Do a GCD of the input value and N to see if a factor was found. */
/* The GCD is returned in factor iff a factor is found. */
/* This routine used to be interruptible and thus returns a stop_reason. */
/* Since switching to GMP's mpz code to implement the GCD this routine is no longer interruptible. */

int gcd (
	gwhandle *gwdata,
	int	thread_num,
	gwnum	gg,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found if any */
{
	giant	v;
	mpz_t	a, b;

/* Assume a factor will not be found */

	*factor = NULL;

/* Convert input number to binary */

	v = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 10);
	if (v == NULL) goto oom;
	gwunfft (gwdata, gg, gg);		// Just in case caller partially FFTed gg
	if (gwtogiant (gwdata, gg, v)) {	// On unexpected error, return no factor found
		pushg (&gwdata->gdata, 1);
		return (0);
	}

/* Do the GCD */

	mpz_init (a);
	mpz_init (b);
	gtompz (v, a);
	gtompz (N, b);
	pushg (&gwdata->gdata, 1);
	mpz_gcd (a, a, b);

/* If a factor was found, save it in FAC */

	if (mpz_cmp_ui (a, 1) && mpz_cmp (a, b)) {
		*factor = allocgiant ((int) divide_rounding_up (mpz_sizeinbase (a, 2), 32));
		if (*factor == NULL) goto oom;
		mpztog (a, *factor);
	}

/* Cleanup and return */

	mpz_clear (a);
	mpz_clear (b);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (thread_num));
}

/* Computes the modular inverse of a number.  This is done using the */
/* extended GCD algorithm.  If a factor is accidentally found, it is */
/* returned in factor.  Function returns stop_reason if it was */
/* interrupted by an escape. */

int ecm_modinv (
	ecmhandle *ecmdata,
	gwnum	b,
	giant	N,			/* Number we are factoring */
	giant	*factor)		/* Factor found, if any */
{
	giant	v;

/* Convert input number to binary */

	v = popg (&ecmdata->gwdata.gdata, ((int) ecmdata->gwdata.bit_length >> 5) + 10);
	if (v == NULL) goto oom;
	if (gwtogiant (&ecmdata->gwdata, b, v)) {
		// On unexpected, should-never-happen error, return out-of-memory for lack of a better error message
		goto oom;
	}

#ifdef MODINV_USING_GIANTS

	int	stop_reason;

/* Let the invg code use gwnum b's memory.  This code is slower, but at least it is interruptible. */
/* Compute 1/v mod N */

	gwfree_temporarily (&ecmdata->gwdata, b);
	stop_reason = invgi (&ecmdata->gwdata.gdata, ecmdata->thread_num, N, v);
	if (stop_reason == GIANT_OUT_OF_MEMORY)
		stop_reason = OutOfMemory (ecmdata->thread_num);
	gwrealloc_temporarily (&ecmdata->gwdata, b);
	if (stop_reason) {
		pushg (&ecmdata->gwdata.gdata, 1);
		return (stop_reason);
	}

/* If a factor was found, save it in FAC */

	if (v->sign < 0) {
		negg (v);
		*factor = allocgiant (v->sign);
		if (*factor == NULL) goto oom;
		gtog (v, *factor);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		*factor = NULL;
		gianttogw (&ecmdata->gwdata, v, b);
	}

/* Use the faster GMP library to do an extended GCD which gives us 1/v mod N */

#else
	{
	mpz_t	__v, __N, __gcd, __inv;

/* Do the extended GCD */

	mpz_init (__v);
	mpz_init (__N);
	mpz_init (__gcd);
	mpz_init (__inv);
	gtompz (v, __v);
	gtompz (N, __N);
	mpz_gcdext (__gcd, __inv, NULL, __v, __N);
	mpz_clear (__v);

/* If a factor was found (gcd != 1 && gcd != N), save it in FAC */

	if (mpz_cmp_ui (__gcd, 1) && mpz_cmp (__gcd, __N)) {
		*factor = allocgiant ((int) divide_rounding_up (mpz_sizeinbase (__gcd, 2), 32));
		if (*factor == NULL) goto oom;
		mpztog (__gcd, *factor);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		*factor = NULL;
		if (mpz_sgn (__inv) < 0) mpz_add (__inv, __inv, __N);
		mpztog (__inv, v);
		gianttogw (&ecmdata->gwdata, v, b);
	}

/* Cleanup and return */

	mpz_clear (__gcd);
	mpz_clear (__inv);
	mpz_clear (__N);
	}
#endif

/* Clean up */

	pushg (&ecmdata->gwdata.gdata, 1);

/* Increment count and return */

	ecmdata->modinv_count++;
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}


/* Computes the modular inverse of an array of numbers */
/* Uses extra multiplications to make only one real modinv call */
/* Uses the simple formula 1/a = b * 1/ab, 1/b = a * 1/ab */
/* If we accidentally find a factor it is returned in factor. */
/* Return stop_reason if there is a user interrupt. */

int grouped_modinv (
	ecmhandle *ecmdata,
	gwnum	*b,
	unsigned int size,
	gwnum	*tmp,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	unsigned int i;
	gwnum	*orig_tmp;
	int	stop_reason;

/* Handle group of 1 as a special case */

	if (size == 1) return (ecm_modinv (ecmdata, *b, N, factor));

/* Handle an odd size */

	orig_tmp = tmp;
	if (size & 1) {
		gwswap (b[0], *tmp);
		tmp++;
	}

/* Multiply each pair of numbers */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (&ecmdata->gwdata, b[i], b[i]);
		gwfft (&ecmdata->gwdata, b[i+1], b[i+1]);
		gwfftfftmul (&ecmdata->gwdata, b[i], b[i+1], *tmp);
		tmp++;
	}

/* Recurse */

	stop_reason = grouped_modinv (ecmdata, orig_tmp, (size+1) / 2, tmp, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) return (0);

/* Handle an odd size */

	if (size & 1) {
		gwswap (b[0], *orig_tmp);
		orig_tmp++;
	}

/* Now perform multiplications on each pair to get the modular inverse */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (&ecmdata->gwdata, *orig_tmp, *orig_tmp);
		gwfftfftmul (&ecmdata->gwdata, *orig_tmp, b[i], b[i]);
		gwfftfftmul (&ecmdata->gwdata, *orig_tmp, b[i+1], b[i+1]);
		gwswap (b[i], b[i+1]);
		orig_tmp++;
	}

/* All done, return */

	return (0);
}

/* Takes a point (a,b) and multiplies it by a value such that b will be one */
/* If we accidentally find a factor it is returned in factor.  Function */
/* returns stop_reason if it was interrupted by an escape. */

int normalize (
	ecmhandle *ecmdata,
	gwnum	a,
	gwnum	b,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	int	stop_reason;

/* Compute the modular inverse and scale up the first input value */

	stop_reason = ecm_modinv (ecmdata, b, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) return (0);
	gwmul (&ecmdata->gwdata, b, a);
	return (0);
}

/* Initialize the normalize (modular inverse) pool. */

int normalize_pool_init (
	ecmhandle *ecmdata,
	int	max_pool_size)
{
	ecmdata->pool_count = 0;
	ecmdata->poolz_count = 0;
	ecmdata->pool_values = NULL;
	ecmdata->poolz_values = NULL;
	ecmdata->pool_modinv_value = NULL;

/* Allocate memory for modular inverse pooling arrays of pointers to gwnums */

	ecmdata->pool_values = (gwnum *) malloc (max_pool_size * sizeof (gwnum));
	if (ecmdata->pool_values == NULL) goto oom;
	if (ecmdata->pool_type == POOL_3MULT) {		// 3MULT pooling requires 100% more gwnums
		ecmdata->poolz_values = (gwnum *) malloc (max_pool_size * sizeof (gwnum));
		if (ecmdata->poolz_values == NULL) goto oom;
	}
	if (ecmdata->pool_type == POOL_3POINT44MULT) {	// 3POINT44MULT pooling requires 33% plus a couple more gwnums
		ecmdata->poolz_values = (gwnum *) malloc ((divide_rounding_up (max_pool_size, 3) + 2) * sizeof (gwnum));
		if (ecmdata->poolz_values == NULL) goto oom;
	}
	if (ecmdata->pool_type == POOL_3POINT57MULT) {	// 3POINT57MULT pooling requires 14% plus 5 more gwnums
		ecmdata->poolz_values = (gwnum *) malloc ((divide_rounding_up (max_pool_size, 7) + 5)* sizeof (gwnum));
		if (ecmdata->poolz_values == NULL) goto oom;
	}

/* All done */

	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Free data associated with the normalize (modular inverse) pool. */

void normalize_pool_term (
	ecmhandle *ecmdata)
{
	int	i;

	if (ecmdata->poolz_values != NULL) {
		for (i = 0; i < ecmdata->poolz_count; i++) gwfree (&ecmdata->gwdata, ecmdata->poolz_values[i]);
		free (ecmdata->poolz_values); ecmdata->poolz_values = NULL;
	}

	free (ecmdata->pool_values); ecmdata->pool_values = NULL;
	gwfree (&ecmdata->gwdata, ecmdata->pool_modinv_value); ecmdata->pool_modinv_value = NULL;
	ecmdata->pool_count = 0;
	ecmdata->poolz_count = 0;
}

/* Adds a point (a,b) to the list of numbers that need normalizing. */
/* This is done in such a way as to minimize the amount of memory used. */

/* This is an interesting bit of code with a variety of algorithms available.  Assuming there are N pairs to normalize, then you can: */
/* 1) Use 2*N memory and 3 multiplies per point. */
/* 2) Use 1.33*N memory and use 3.444 multiplies per point. */
/* 3) Use 1.14*N memory and use 3.573 multiplies per point. */
/* 4) Use N memory and use O(N^2) multiplies per point. */

int add_to_normalize_pool_internal (
	ecmhandle *ecmdata,
	gwnum	a,		/* Number to multiply by 1/b.  This input can be in a full or partial FFTed state.  Not preserved! */
	gwnum	b,		/* Preserved, will not be FFTed in place */
	giant	*factor,	/* Factor found, if any */
	int	relinquish_b)	/* TRUE if caller is relinquishing ownership of b, saves a gwcopy */
{
#define allocAndFFTb(d)		if(relinquish_b){gwfft(&ecmdata->gwdata,b,b);d=b;}\
				else{gwnum t=gwalloc(&ecmdata->gwdata);if(t==NULL)goto oom;gwfft(&ecmdata->gwdata,b,t);d=t;}
#define allocAndCopyb(d)	if(relinquish_b)d=b;\
				else{gwnum t=gwalloc(&ecmdata->gwdata);if(t==NULL)goto oom;gwcopy(&ecmdata->gwdata,b,t);d=t;}
#define copyAndFFTb(d)		if(relinquish_b){gwswap(b,d);gwfree(&ecmdata->gwdata,b);}else gwfft(&ecmdata->gwdata,b,d);

/* Switch off the type of pooling we are going to do */

	switch (ecmdata->pool_type) {

/* Implement 3-multiply algorithm with 100% memory overhead */

	case POOL_3MULT:

/* If this is the first call allocate memory for the gwnum we use */

		if (ecmdata->pool_count == 0) {
			allocAndFFTb (ecmdata->pool_modinv_value);
		}

/* Otherwise, multiply a by the accumulated b values */

		else {
			// Accumulate previous b value (delayed so that we can use GWMUL_STARTNEXTFFT on all but last b value)
			if (ecmdata->pool_count > 1) {
				gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			}
			// Multiply a by accumulated b values
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, a, a, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// Remember b
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;

/* Implement 3.44-multiply algorithm with 33% memory overhead */
// The basic algorithm here takes 4 (a,b) points and converts them as follows:
//	a1							b1			one copy
//	a1*b2	    a2*b1					b1*b2			7 FFTs (FFT b1, FFT b2, FFTINVFFT a1*b2, FFTINVFFT a2*b1, INVFFT b1*b2)
//	a1*b2	    a2*b1       a3,b3				b1*b2
//	a1*b2       a2*b1       a3*b4,b3*b4 a4*b3		b1*b2			7 FFTs (FFT b3, FFT b4, FFTINVFFT a3*b4, FFTINVFFT a4*b3, INVFFT b3*b4)
//	a1*b2*b3*b4 a2*b1*b3*b4 a3*b4*b1*b2 a4*b3*b1*b2		b1*b2*b3*b4		11 FFTs (FFT b1*b2, FFT b3*b4, FFTINVFFT a1* a2* a3* a4*, INVFFT b1*b2*b3*b4)
//	inv = 1/(b1*b2*b3*b4)
//	a4 *= inv									9 FFTs (FFT inv, FFTINVFFT a1* a2* a3* a4*)
//	a3 *= inv								
//	a2 *= inv
//	a1 *= inv
// The trick is that subsequent groups of 4 operate on (1,product_of_b) and 3 more (a,b) points.  For example, the steps needed for the next 3 (a,b) points:
// On 5th input:
//	let b0 = pooled_modinv
//	fft b5 to alloc poolz[0] - will become 1/(b1*b2*b3*b4)
//	pool[4] = a5*b0 will become a5/b5
//	pooled_modinv *= b5
//	b5	    a5*b0					b0*b5			5 FFTs (FFT b0, FFT b5, FFTINVFFT a5*b0, INVFFT b0*b5)
// On 6th input:
//	fft b6 to tmp1, will free later
//	pool[5] = a6 will become a6/b6
//	b5	    a5*b0       a6,b6				b0*b5
// On 7th input:
//	tmp2 = fft b7
//	pool[5] = a6*b7
//	pool[6] = a7*b6 will become a6/b6
//	tmp1 = b6*b7
//	b5       a5*b0       a6*b7,b6*b7 a7*b6			b0*b5			7 FFTs (FFT b6, FFT b7, FFTINVFFT a6*b7, FFTINVFFT a7*b6, INVFFT b6*b7)
//	poolz[0] *= b6*b7
//	pool[4] *= b6*b7
//	pool[5] *= b0*b5
//	pool[6] *= b0*b5
//	pooled_modinv *= b6*b7
//	free tmp1,tmp2
//	b5*b6*b7 a5*b0*b6*b7 a6*b7*b0*b5 a7*b6*b0*b5		b0*b5*b6*b7		10 FFTs (FFT b0*b5, FFT b6*b7, INVFFT b5*, FFTINVFFT a5* a6* a7*, INVFFT b0*b5*b6*b7)
//	inv = 1/(b0*b5*b6*b7)
//	a7 *= inv									9 FFTs (FFT inv, FFTINVFFT a5* a5* a6* a7*)
//	a6 *= inv								
//	a5 *= inv
//	poolz[0] *= inv
// Careful counting shows that each 3 additional points costs 31 FFTs

	case POOL_3POINT44MULT:

/* Handle the first call for the first set of four points.  Allocate and set pool_modinv_value (the product_of_b's) */

		if (ecmdata->pool_count == 0) {
			// FFT b1
			allocAndFFTb (ecmdata->pool_modinv_value);
		}

/* Handle the second call for the first set of four points */

		else if (ecmdata->pool_count == 1) {
			// FFT b2 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
			// a1 *= b2
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_values[0], GWMUL_STARTNEXTFFT);
			// a2 *= b1
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool_modinv_value *= b2 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
		}

/* Handle the third call for the first set of four points */

		else if (ecmdata->pool_count == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b3 (temporarily store in poolz)
			copyAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count-1]);
		}

/* Handle the fourth call for the first set of four points */

		else if (ecmdata->pool_count == 3) {
			gwnum	tmp1, tmp2;
			// FFT b4
			allocAndFFTb (tmp2);
			// a3 *= b4
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], tmp2, ecmdata->pool_values[2], GWMUL_STARTNEXTFFT);
			// a4 *= b3
			tmp1 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			gwmul3 (&ecmdata->gwdata, a, tmp1, a, GWMUL_STARTNEXTFFT);
			// tmp1 = b3*b4
			gwmul3 (&ecmdata->gwdata, tmp1, tmp2, tmp1, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, tmp2);
			// pool[0] *= b3*b4
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], tmp1, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool[1] *= b3*b4
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], tmp1, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// pool[2] *= b1*b2
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], ecmdata->pool_modinv_value, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// a4 *= b1*b2
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pooled_modinv *= b3*b4 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, tmp1, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
		}

/* Handle the first call for a new set of three points */

		else if (ecmdata->pool_count % 3 == 1) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b5
			copyAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count-1]);
			// a5 *= pool_modinv_value
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool_modinv_value *= b5 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
		}

/* Handle the second call for a set of three points */

		else if (ecmdata->pool_count % 3 == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b6 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
		}

/* Handle the third call for a set of three points */

		else {
			gwnum	tmp1, tmp2;
			// FFT b7
			allocAndFFTb (tmp2);
			// a6 *= b7
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], tmp2, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_STARTNEXTFFT);
			// a7 *= b6
			tmp1 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			gwmul3 (&ecmdata->gwdata, a, tmp1, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// tmp1 = b6*b7
			gwmul3 (&ecmdata->gwdata, tmp1, tmp2, tmp1, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, tmp2);
			// poolz[0] *= b6*b7
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-2], tmp1, ecmdata->poolz_values[ecmdata->poolz_count-2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool[4] *= b6*b7
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], tmp1, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// pool[5] *= b0*b5
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// a7 *= b0*b5
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pooled_modinv *= b6*b7 (delay this so we can use startnextfft)
			// gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, tmp1, ecmdata->pool_modinv_value, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;

/* Implement 3.573-multiply algorithm with 14% memory overhead */
// The basic algorithm here takes 8 (a,b) points and converts them as follows:
//	a1							b1		one copy
//	a1*b2	    a2*b1					b1*b2		7 FFTs (FFT b1, FFT b2, FFTINVFFT a2*b1, INVFFT b1*b2)
//	a1*b2	    a2*b1       a3,b3				b1*b2
//	a1*b2       a2*b1       a3*b4,b3*b4 a4*b3		b1*b2		7 FFTs (FFT b3, FFT b4, FFTINVFFT a3*b4, FFTINVFFT a4*b3, INVFFT b3*b4)
//	a5,b5
//	a5*b6	    a6*b5					b5*b6		7 FFTs (FFT b5, FFT b6, FFTINVFFT a5*b6, FFTINVFFT a6*b5, INVFFT b5*b6)
//	a5*b6	    a6*b5       a7,b7				b5*b6
//	a5*b6       a6*b5       a7*b8,b7*b8 a8*b7		b5*b6		7 FFTs (FFT b7, FFT b8, FFTINVFFT a7*b8, FFTINVFFT a8*b7, INVFFT b7*b8)
//		b1*b2*b3*b4	b5*b6*b7*b8					6 FFTs (FFT b1*b2, FFT b3*b4, FFT b5*b6, FFT b7*b8, 2 INVFFT)
//		b12*b5678  b34*b5678  b56*b1234  b78*b1234			6 FFTs (FFT b1*b2*b3*b4, FFT b5*b6*b7*b8, 4 INVFFT)
//	a1*b2*b345678 a2*b1*b345678 a3*b4*b125678 a4*b3*b125678			10 FFTs (FFT b345678, FFT b125678, FFTINVFFT a1* a2* a3* a4*)
//	a5*b6*b123478 a6*b5*b123478 a7*b8*b123456 a8*b7*b123456			10 FFTs (FFT b123478, FFT b123456, FFTINVFFT a5* a6* a7* a8*)
//							b1*b2*b3*b4*b5*b6*b7*b8	1 FFT (INVFFT)
//	8 multiplies by inv							17 FFTs
// The trick is that subsequent groups of 8 operate on (1,product_of_b) and 7 more (a,b) points.  For example, the steps needed for the next 7 (a,b) points:
// On 9th input:
//	let b0 = pooled_modinv
//	fft b9 to alloc poolz[0] - will become 1/(b1*b2*b3*b4*b5*b6*b7*b8)
//	b9	    a9*b0					b0*b9		5 FFTs (FFT b0, FFT b9, FFTINVFFT a9*b0, INVFFT b0*b9)
// Then 3 pairs do this:
//	a10,b10
//	a10*b11,b10*b11 a11*b10							7 FFTs (FFT b10, FFT b11, FFTINVFFT a10*b11, FFTINVFFT a11*b10, INVFFT b10*b11)
//	a12,b12
//	a12*b13,b12*b13   a13*b12						7 FFTs (FFT b12, FFT b13, FFTINVFFT a12*b13, FFTINVFFT a13*b12, INVFFT b12*b13)
//	a14,b14			
//	a14*b15,b14*b15   a15*b14						7 FFTs (FFT b14, FFT b15, FFTINVFFT a14*b15, FFTINVFFT a15*b14, INVFFT b14*b15)
// After the last pair, we do a lot of combining:
//		b0*b9*b10*b11	b12*b13*b14*b15					6 FFTs (FFT b0*b9, FFT b10*b11, FFT b12*b13, FFT b14*b15, 2 INVFFT)
//		b02*b12131415  b1011*b12131415  b1213*b091011  b1415*b091011	6 FFTs (FFT b0*b9*b10*b11, FFT b12*b13*b14*b15, 4 INVFFT)
//	b9*b10-15  a9*b0*b10-15 a10*b11*b0212-15 a11*b10*b0212-15		9 FFTs (FFT b10-15, FFT b0212-15, INVFFT b9*, FFTINVFFT a9* a10* a11*)
//	a12*b13*b0910111415 a13*b0910111415 a14*b15*b09-13 a15*b14*b09-13	10 FFTs (FFT b0910111415, FFT b09-13, FFTINVFFT a12* a13* a14* a15*)
//						b0*b9*b10*b11*b12*b13*b14*b15	1 FFT (INVFFT)
//	8 multiplies by inv							17 FFTs
// Careful counting shows that each 7 additional points costs 75 FFTs

	case POOL_3POINT57MULT:

/* Handle the first call for the first set of eight points.  Allocate and set pool_modinv_value (the product_of_b's) */

		if (ecmdata->pool_count == 0) {
			// FFT b1
			allocAndFFTb (ecmdata->pool_modinv_value);
		}

/* Handle the second call for the first set of eight points */

		else if (ecmdata->pool_count == 1) {
			// FFT b2 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
			// a1 *= b2
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_values[0], GWMUL_STARTNEXTFFT);
			// a2 *= b1
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool_modinv_value *= b2 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
		}

/* Handle the third call for the first set of eight points */

		else if (ecmdata->pool_count == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b3 (temporarily store in poolz)
			copyAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count-1]);
		}

/* Handle the fourth and sixth call for the first set of eight points */

		else if (ecmdata->pool_count == 3 || ecmdata->pool_count == 5) {
			// FFT b4 or b6
			gwnum	tmp;
			allocAndFFTb (tmp);
			// a3 *= b4 or a5 *= b6
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], tmp, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// a4 *= b3 or a6 *= b5
			gwmul3 (&ecmdata->gwdata, a, ecmdata->poolz_values[ecmdata->poolz_count-1], a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b3 *= b4 or b5 *= b6
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], tmp, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, tmp);
		}

/* Handle the fifth and seventh call for first set of eight points */

		else if (ecmdata->pool_count == 4 || ecmdata->pool_count == 6) {
			// FFT b5 or b7 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
		}

/* Handle the eighth call for first set of eight points */

		else if (ecmdata->pool_count == 7) {
			gwnum	b7, b8, b34, b56, b78, b5678, b6vals;
			// Allocate two more temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b7 = ecmdata->poolz_values[2];
			b56 = ecmdata->poolz_values[1];
			b34 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// FFT b8
			allocAndFFTb (b8);
			// a7 *= b8
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], b8, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b78 = b7*b8
			gwmul3 (&ecmdata->gwdata, b7, b8, b8, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			b78 = b8;
			// a8 *= b7
			gwmul3 (&ecmdata->gwdata, a, b7, a, GWMUL_STARTNEXTFFT);
			// b56 * b78
			b5678 = b7;
			gwmul3 (&ecmdata->gwdata, b56, b78, b5678, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b34 * b5678
			gwmul3 (&ecmdata->gwdata, b34, b5678, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a1*b2)*b345678
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b6vals, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b345678
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b6vals, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// b12 * b5678
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b5678, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a3*b4)*b125678
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], b6vals, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a4*b3)*b125678
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[3], b6vals, ecmdata->pool_values[3], GWMUL_STARTNEXTFFT);
			// b12 * b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b34, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// b78 * b1234
			gwmul3 (&ecmdata->gwdata, b78, ecmdata->pool_modinv_value, b6vals, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a5*b6)*b123478
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[4], b6vals, ecmdata->pool_values[4], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a6*b5)*b123478
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[5], b6vals, ecmdata->pool_values[5], GWMUL_STARTNEXTFFT);
			// b56 * b1234
			gwmul3 (&ecmdata->gwdata, b56, ecmdata->pool_modinv_value, b6vals, GWMUL_STARTNEXTFFT);
			// (a7*b8)*b123456
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[6], b6vals, ecmdata->pool_values[6], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a8*b7)*b123456
			gwmul3 (&ecmdata->gwdata, a, b6vals, a, GWMUL_STARTNEXTFFT);
			// b1*b2*b3*b4*b5*b6*b7*b8 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b5678, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			ecmdata->poolz_values[0] = b5678;
			ecmdata->poolz_count = 1;
			// Free temps
			gwfree (&ecmdata->gwdata, b34);
			gwfree (&ecmdata->gwdata, b56);
			gwfree (&ecmdata->gwdata, b78);
			gwfree (&ecmdata->gwdata, b6vals);
		}

/* Handle the first call for a new set of seven points */

		else if (ecmdata->pool_count % 7 == 1) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b9
			copyAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count-1]);
			// a9 *= pool_modinv_value
			gwmul3 (&ecmdata->gwdata, a, ecmdata->pool_modinv_value, a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool_modinv_value *= b9 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
		}

/* Handle the second call for a set of seven points */

		else if (ecmdata->pool_count % 7 == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// FFT b10 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
		}

/* Handle the third and fifth call for a set of seven points */

		else if (ecmdata->pool_count % 7 == 3 || ecmdata->pool_count % 7 == 5) {
			// FFT b11 or b13
			gwnum	tmp;
			allocAndFFTb (tmp);
			// a10 *= b11 or a12 *= b13
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], tmp, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// a11 *= b10 or a13 *= b12
			gwmul3 (&ecmdata->gwdata, a, ecmdata->poolz_values[ecmdata->poolz_count-1], a, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b10 *= b11 or b12 *= b13
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], tmp, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, tmp);
		}

/* Handle the fourth and sixth call for a set of seven points */

		else if (ecmdata->pool_count % 7 == 4 || ecmdata->pool_count % 7 == 6) {
			// FFT b12 or b14 (temporarily store in poolz)
			allocAndFFTb (ecmdata->poolz_values[ecmdata->poolz_count++]);
		}

/* Handle the last call for a set of seven points */

		else {
			gwnum	b14, b15, b1011, b1213, b1415, b12131415, b6vals;
			// Allocate two more temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b14 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			b1213 = ecmdata->poolz_values[ecmdata->poolz_count-2];
			b1011 = ecmdata->poolz_values[ecmdata->poolz_count-3];
			ecmdata->poolz_count -= 3;
			// FFT b15
			allocAndFFTb (b15);
			// a14 *= b15
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], b15, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1415 = b14*b15
			gwmul3 (&ecmdata->gwdata, b14, b15, b15, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			b1415 = b15;
			// a15 *= b14
			gwmul3 (&ecmdata->gwdata, a, b14, a, GWMUL_STARTNEXTFFT);
			// b1213 * b1415
			b12131415 = b14;
			gwmul3 (&ecmdata->gwdata, b1213, b1415, b12131415, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1011 * b12131415
			gwmul3 (&ecmdata->gwdata, b1011, b12131415, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b9*b101112131415
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b6vals, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b101112131415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-6], b6vals, ecmdata->pool_values[ecmdata->pool_count-6], GWMUL_STARTNEXTFFT);
			// b09 * b12131415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b12131415, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a10*b11)*b0912131415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-5], b6vals, ecmdata->pool_values[ecmdata->pool_count-5], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a11*b10)*b0912131415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-4], b6vals, ecmdata->pool_values[ecmdata->pool_count-4], GWMUL_STARTNEXTFFT);
			// b09 * b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1011, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// b1415 * b091011
			gwmul3 (&ecmdata->gwdata, b1415, ecmdata->pool_modinv_value, b6vals, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a12*b13)*b0910111415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-3], b6vals, ecmdata->pool_values[ecmdata->pool_count-3], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a13*b12)*b0910111415
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], b6vals, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// b1213 * b091011
			gwmul3 (&ecmdata->gwdata, b1213, ecmdata->pool_modinv_value, b6vals, GWMUL_STARTNEXTFFT);
			// (a14*b15)*b0910111213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], b6vals, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a15*b14)*b0910111213
			gwmul3 (&ecmdata->gwdata, a, b6vals, a, GWMUL_STARTNEXTFFT);
			// b0*b9*b10*b11*b12*b13*b14*b15 (delay this multiply so that we can use startnextfft)
			//gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b12131415, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			ecmdata->poolz_values[ecmdata->poolz_count++] = b12131415;
			// Free temps
			gwfree (&ecmdata->gwdata, b1011);
			gwfree (&ecmdata->gwdata, b1213);
			gwfree (&ecmdata->gwdata, b1415);
			gwfree (&ecmdata->gwdata, b6vals);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;

/* Implement N^2 algorithm using only 1 extra gwnum */
/* The algorithm costs (in FFTs):	(N+3)^2 - 13 + modinv_cost	*/
/* To find the break even point:					*/
/*	(N+3)^2 + modinv_cost = 2 * ((N/2+3)^2 + modinv_cost)		*/
/*	(N+3)^2 + modinv_cost = (N/2+3)^2 * 2 + 2 * modinv_cost		*/
/*	(N+3)^2 = (N/2+3)^2 * 2 + modinv_cost				*/
/*	N^2 + 6*N + 9 = (N^2/4 + 3*N + 9) * 2 + modinv_cost		*/
/*	N^2 + 6*N + 9 - (N^2/4 + 3*N + 9) * 2 = modinv_cost		*/
/*	N^2/2 - 9 = modinv_cost						*/
/*	N = sqrt(2 * (modinv_cost + 9))					*/
/* If modinv_cost = 4000, breakeven is 89.5				*/

	case POOL_N_SQUARED:

/* If this is the first call allocate memory for the gwnum we use */

		if (ecmdata->pool_count == 0) {
			allocAndCopyb (ecmdata->pool_modinv_value);
		}

/* Otherwise, multiply a by the accumulated b values and multiply all previous a's by this b */

		else {
			int	i;
			gwnum	tmp;
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, a, a, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			allocAndFFTb (tmp);
			for (i = 0; i < ecmdata->pool_count; i++)
				gwmul3 (&ecmdata->gwdata, tmp, ecmdata->pool_values[i], ecmdata->pool_values[i], GWMUL_STARTNEXTFFT);
			gwmul3 (&ecmdata->gwdata, tmp, ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, tmp);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;
	}
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

// The original add_to_normalize_pool interface
int add_to_normalize_pool (
	ecmhandle *ecmdata,
	gwnum	a,		/* Number to multiply by 1/b.  This input can be in a full or partial FFTed state.  Not preserved! */
	gwnum	b,		/* Preserved, will not be FFTed in place */
	giant	*factor)	/* Factor found, if any */
{
	return (add_to_normalize_pool_internal (ecmdata, a, b, factor, FALSE));
}

// The new add_to_normalize_pool interface. "_ro" stands for relinquish ownership.  Caller relinquishes ownership of b.  Saves a gwcopy.
int add_to_normalize_pool_ro (
	ecmhandle *ecmdata,
	gwnum	a,		/* Number to multiply by 1/b.  This input can be in a full or partial FFTed state.  Not preserved! */
	gwnum	b,		/* Caller loses ownership of this gwnum! */
	giant	*factor)	/* Factor found, if any */
{
	return (add_to_normalize_pool_internal (ecmdata, a, b, factor, TRUE));
}

/* Takes each point from add_to_normalize_pool and normalizes it. */
/* If we accidentally find a factor, it is returned in factor. */
/* Return stop_reason if there is a user interrupt. */

int normalize_pool (
	ecmhandle *ecmdata,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	int	i, stop_reason;

/* For 3-mult pooling, accumulate last b value (was delayed so that we can use GWMUL_STARTNEXTFFT on all but the last pooled b value) */

	if (ecmdata->pool_type == POOL_3MULT && ecmdata->pool_count > 1) {
		gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, 0);
	}

/* For 3.44-mult pooling finish off a partial group of four.  Accumulate last b value (was delayed so that we can use */
/* GWMUL_STARTNEXTFFT on all but the last pooled b value). */

	if (ecmdata->pool_type == POOL_3POINT44MULT) {

		/* Handle having only two of first set of four points. */
		if (ecmdata->pool_count == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[0], ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[0]);
			ecmdata->poolz_count = 0;
		}

		/* Handle having only three of first set of four points. */
		/* Do much of the work that would have happened if the fourth point had been added to the pool. */
		else if (ecmdata->pool_count == 3) {
			gwnum	b3 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// pool[0] *= b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b3, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool[1] *= b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b3, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// pool[2] *= b1*b2
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], ecmdata->pool_modinv_value, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pooled_modinv *= b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b3, ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, b3);
		}

		/* Handle having only the first point in a set of three additional points */
		else if (ecmdata->pool_count > 4 && ecmdata->pool_count % 3 == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, 0);
		}

		/* Handle having only two of an additional set of three points */
		/* Do much of the work that would have happened if the last point had been added to the pool */
		else if (ecmdata->pool_count > 4 && ecmdata->pool_count % 3 == 0) {
			gwnum	b6 = ecmdata->poolz_values[--ecmdata->poolz_count];
			// poolz[0] *= b6
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b6, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pool[4] *= b6
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], b6, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// pool[5] *= b0*b5
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// pooled_modinv *= b6
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b6, ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, b6);
		}

		/* Handle having a complete set of points */
		else if (ecmdata->pool_count >= 4 && ecmdata->pool_count % 3 == 1) {
			// Do the delayed multiply for pool_modinv
			ecmdata->poolz_count--;
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count], ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count]);
		}
	}

/* For 3.57-mult pooling finish off a partial group of eight.  Accumulate last b value (was delayed so that we can use */
/* GWMUL_STARTNEXTFFT on all but the last pooled b value). */

	if (ecmdata->pool_type == POOL_3POINT57MULT) {

		/* Handle having only two of first set of eight points */
		if (ecmdata->pool_count == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[0], ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[0]);
			ecmdata->poolz_count = 0;
		}

		/* Handle having only three of first set of eight points */
		else if (ecmdata->pool_count == 3) {
			gwnum	b3;
			// Load values temporarily stored in poolz
			b3 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// (a1*b2)*b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b3, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b3, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// a3*b12
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], ecmdata->pool_modinv_value, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1*b2*b3
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b3, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b3);
		}

		/* Handle having only four of first set of eight points */
		else if (ecmdata->pool_count == 4) {
			gwnum	b34;
			// Load values temporarily stored in poolz
			b34 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// (a1*b2)*b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b34, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b34, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// (a3*b4)*b12
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], ecmdata->pool_modinv_value, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a4*b3)*b12
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[3], ecmdata->pool_modinv_value, ecmdata->pool_values[3], GWMUL_STARTNEXTFFT);
			// b1*b2*b3*b4
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b34, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b34);
		}

		/* Handle having only five of first set of eight points */
		else if (ecmdata->pool_count == 5) {
			gwnum	b34, b5, b6vals;
			// Allocate more temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b5 = ecmdata->poolz_values[1];
			b34 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// b34 * b5
			gwmul3 (&ecmdata->gwdata, b34, b5, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a1*b2)*b345
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b6vals, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b345
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b6vals, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// b12 * b5
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b5, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a3*b4)*b125
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], b6vals, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a4*b3)*b125
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[3], b6vals, ecmdata->pool_values[3], GWMUL_STARTNEXTFFT);
			// b12 * b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b34, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// a5*b1234
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[4], ecmdata->pool_modinv_value, ecmdata->pool_values[4], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1*b2*b3*b4*b5
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b5, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b34);
			gwfree (&ecmdata->gwdata, b5);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having only six of first set of eight points */
		else if (ecmdata->pool_count == 6) {
			gwnum	b34, b56, b6vals;
			// Allocate more temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b56 = ecmdata->poolz_values[1];
			b34 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// b34 * b56
			gwmul3 (&ecmdata->gwdata, b34, b56, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a1*b2)*b3456
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b6vals, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b3456
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b6vals, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// b12 * b56
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b56, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a3*b4)*b1256
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], b6vals, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a4*b3)*b1256
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[3], b6vals, ecmdata->pool_values[3], GWMUL_STARTNEXTFFT);
			// b12 * b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b34, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// (a5*b6)*b1234
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[4], ecmdata->pool_modinv_value, ecmdata->pool_values[4], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a6*b5)*b1234
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[5], ecmdata->pool_modinv_value, ecmdata->pool_values[5], GWMUL_STARTNEXTFFT);
			// b1*b2*b3*b4*b5*b6
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b56, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b34);
			gwfree (&ecmdata->gwdata, b56);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having seven of first set of eight points. */
		/* Do much of the work that would have happened if the fourth point had been added to the pool. */
		else if (ecmdata->pool_count == 7) {
			gwnum	b7, b34, b56, b567, b6vals;
			// Allocate more temporaries
			b567 = gwalloc (&ecmdata->gwdata);
			if (b567 == NULL) goto oom;
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b7 = ecmdata->poolz_values[2];
			b56 = ecmdata->poolz_values[1];
			b34 = ecmdata->poolz_values[0];
			ecmdata->poolz_count = 0;
			// b56 * b7
			gwmul3 (&ecmdata->gwdata, b56, b7, b567, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b34 * b567
			gwmul3 (&ecmdata->gwdata, b34, b567, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a1*b2)*b34567
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[0], b6vals, ecmdata->pool_values[0], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a2*b1)*b34567
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[1], b6vals, ecmdata->pool_values[1], GWMUL_STARTNEXTFFT);
			// b12 * b567
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b567, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a3*b4)*b12567
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[2], b6vals, ecmdata->pool_values[2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a4*b3)*b12567
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[3], b6vals, ecmdata->pool_values[3], GWMUL_STARTNEXTFFT);
			// b12 * b34
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b34, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// b7 * b1234
			gwmul3 (&ecmdata->gwdata, b7, ecmdata->pool_modinv_value, b6vals, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a5*b6)*b12347
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[4], b6vals, ecmdata->pool_values[4], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a6*b5)*b12347
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[5], b6vals, ecmdata->pool_values[5], GWMUL_STARTNEXTFFT);
			// b56 * b1234
			gwmul3 (&ecmdata->gwdata, b56, ecmdata->pool_modinv_value, b6vals, GWMUL_STARTNEXTFFT);
			// a7*b123456
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[6], b6vals, ecmdata->pool_values[6], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1*b2*b3*b4*b5*b6*b7
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b567, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b34);
			gwfree (&ecmdata->gwdata, b56);
			gwfree (&ecmdata->gwdata, b7);
			gwfree (&ecmdata->gwdata, b567);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having only the first point in a set of seven additional points */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 2) {
			// Do the delayed multiply for pool_modinv
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, 0);
		}

		/* Handle having only two of an additional set of seven points */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 3) {
			gwnum	b10;
			// Load values temporarily stored in poolz
			b10 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			ecmdata->poolz_count -= 1;
			// b9*b10
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b10, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b10
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], b10, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// a10*b09
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b0*b9*b10
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b10, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b10);
		}

		/* Handle having only three of an additional set of seven points */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 4) {
			gwnum	b1011;
			// Load values temporarily stored in poolz
			b1011 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			ecmdata->poolz_count -= 1;
			// b9*b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b1011, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-3], b1011, ecmdata->pool_values[ecmdata->pool_count-3], GWMUL_STARTNEXTFFT);
			// (a10*b11)*b09
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a11*b10)*b09
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_STARTNEXTFFT);
			// b0*b9*b10*b11
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1011, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b1011);
		}

		/* Handle having only four of an additional set of seven points */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 5) {
			gwnum	b1011, b12, b6vals;
			// Allocate temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b12 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			b1011 = ecmdata->poolz_values[ecmdata->poolz_count-2];
			ecmdata->poolz_count -= 2;
			// b1011 * b12
			gwmul3 (&ecmdata->gwdata, b1011, b12, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b9*b101112
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b6vals, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b101112
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-4], b6vals, ecmdata->pool_values[ecmdata->pool_count-4], GWMUL_STARTNEXTFFT);
			// b09 * b12
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b12, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a10*b11)*b0912
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-3], b6vals, ecmdata->pool_values[ecmdata->pool_count-3], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a11*b10)*b0912
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], b6vals, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// b09 * b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1011, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// a12*b091011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b0*b9*b10*b11*b12
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b12, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b1011);
			gwfree (&ecmdata->gwdata, b12);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having only five of an additional set of seven points */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 6) {
			gwnum	b1011, b1213, b6vals;
			// Allocate temporaries
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b1213 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			b1011 = ecmdata->poolz_values[ecmdata->poolz_count-2];
			ecmdata->poolz_count -= 2;
			// b1011 * b1213
			gwmul3 (&ecmdata->gwdata, b1011, b1213, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b9*b10111213
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b6vals, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b10111213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-5], b6vals, ecmdata->pool_values[ecmdata->pool_count-5], GWMUL_STARTNEXTFFT);
			// b09 * b1213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1213, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a10*b11)*b091213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-4], b6vals, ecmdata->pool_values[ecmdata->pool_count-4], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a11*b10)*b091213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-3], b6vals, ecmdata->pool_values[ecmdata->pool_count-3], GWMUL_STARTNEXTFFT);
			// b09 * b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1011, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// (a12*b13)*b091011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a13*b12)*b091011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], ecmdata->pool_modinv_value, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_STARTNEXTFFT);
			// b0*b9*b10*b11*b12*b13
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1213, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b1011);
			gwfree (&ecmdata->gwdata, b1213);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having only six of an additional set of seven points */
		/* Do much of the work that would have happened if the seventh point had been added to the pool. */
		else if (ecmdata->pool_count > 8 && ecmdata->pool_count % 7 == 0) {
			gwnum	b14, b1011, b1213, b121314, b6vals;
			// Allocate temporaries
			b121314 = gwalloc (&ecmdata->gwdata);
			if (b121314 == NULL) goto oom;
			b6vals = gwalloc (&ecmdata->gwdata);
			if (b6vals == NULL) goto oom;
			// Load values temporarily stored in poolz
			b14 = ecmdata->poolz_values[ecmdata->poolz_count-1];
			b1213 = ecmdata->poolz_values[ecmdata->poolz_count-2];
			b1011 = ecmdata->poolz_values[ecmdata->poolz_count-3];
			ecmdata->poolz_count -= 3;
			// b1213 * b14
			gwmul3 (&ecmdata->gwdata, b1213, b14, b121314, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b1011 * b121314
			gwmul3 (&ecmdata->gwdata, b1011, b121314, b6vals, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b9*b1011121314
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count-1], b6vals, ecmdata->poolz_values[ecmdata->poolz_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a9*b0)*b1011121314
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-6], b6vals, ecmdata->pool_values[ecmdata->pool_count-6], GWMUL_STARTNEXTFFT);
			// b09 * b121314
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b121314, b6vals, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			// (a10*b11)*b09121314
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-5], b6vals, ecmdata->pool_values[ecmdata->pool_count-5], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a11*b10)*b09121314
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-4], b6vals, ecmdata->pool_values[ecmdata->pool_count-4], GWMUL_STARTNEXTFFT);
			// b09 * b1011
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b1011, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			// b14 * b091011
			gwmul3 (&ecmdata->gwdata, b14, ecmdata->pool_modinv_value, b6vals, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a12*b13)*b09101114
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-3], b6vals, ecmdata->pool_values[ecmdata->pool_count-3], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// (a13*b12)*b09101114
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-2], b6vals, ecmdata->pool_values[ecmdata->pool_count-2], GWMUL_STARTNEXTFFT);
			// b1213 * b091011
			gwmul3 (&ecmdata->gwdata, b1213, ecmdata->pool_modinv_value, b6vals, GWMUL_STARTNEXTFFT);
			// a14*b0910111213
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_values[ecmdata->pool_count-1], b6vals, ecmdata->pool_values[ecmdata->pool_count-1], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			// b0*b9*b10*b11*b12*b13*b14
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, b121314, ecmdata->pool_modinv_value, 0);
			// Free temps
			gwfree (&ecmdata->gwdata, b14);
			gwfree (&ecmdata->gwdata, b1011);
			gwfree (&ecmdata->gwdata, b1213);
			gwfree (&ecmdata->gwdata, b121314);
			gwfree (&ecmdata->gwdata, b6vals);
		}

		/* Handle having a complete set of points */
		else if (ecmdata->pool_count >= 8 && ecmdata->pool_count % 7 == 1) {
			// Do the delayed multiply for pool_modinv
			ecmdata->poolz_count--;
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count], ecmdata->pool_modinv_value, 0);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count]);
		}
	}

/* Compute the modular inverse */

	gwunfft (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->pool_modinv_value);
	stop_reason = ecm_modinv (ecmdata, ecmdata->pool_modinv_value, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) goto exit;

/* Now invert each value.  Switch off the type of pooling we are doing. */

	switch (ecmdata->pool_type) {

/* Implement 3 multiply algorithm */

	case POOL_3MULT:
		for (i = ecmdata->pool_count-1; ; i--) {
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->pool_values[i], ecmdata->pool_values[i], GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			if (i == 0) break;
			ecmdata->poolz_count--;
			gwmul3 (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count], ecmdata->pool_modinv_value, ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count]);
		}
		break;

/* Implement 3.44 multiply algorithm */

	case POOL_3POINT44MULT:
		for (i = ecmdata->pool_count-1; i >= 0; i--) {
			// Invert one pooled value
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->pool_values[i], ecmdata->pool_values[i], GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			if (i < 4 || i % 3 != 1) continue;
			// Compute inverse for next group of 4
			ecmdata->poolz_count--;
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count]);
		}
		break;

/* Implement 3.57 multiply algorithm */

	case POOL_3POINT57MULT:
		for (i = ecmdata->pool_count-1; i >= 0; i--) {
			// Invert one pooled value
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->pool_values[i], ecmdata->pool_values[i], GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
			if (i < 8 || i % 7 != 1) continue;
			// Compute inverse for next group of 8
			ecmdata->poolz_count--;
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->poolz_values[ecmdata->poolz_count], ecmdata->pool_modinv_value, GWMUL_STARTNEXTFFT);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[ecmdata->poolz_count]);
		}
		break;

/* Implement N^2 algorithm */

	case POOL_N_SQUARED:
		for (i = 0; i < ecmdata->pool_count; i++)
			gwmul3 (&ecmdata->gwdata, ecmdata->pool_modinv_value, ecmdata->pool_values[i], ecmdata->pool_values[i], GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);
		break;
	}

/* Cleanup and reinitialize */

exit:	ecmdata->pool_count = 0;
	ecmdata->poolz_count = 0;
	gwfree (&ecmdata->gwdata, ecmdata->pool_modinv_value); ecmdata->pool_modinv_value = NULL;
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

 
/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed (called $\sigma$ in the Tech. Report),

	u/v = (s^2 - 5)/(4s)

Then starting point is (x_1, y_1) where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/
int choose12 (
	ecmhandle *ecmdata,
	struct xz *xz,		/* Return the curve starting value here */
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	gwnum	xs, zs, t2, t3;
	int	stop_reason;

	xs = gwalloc (&ecmdata->gwdata);
	if (xs == NULL) goto oom;
	zs = gwalloc (&ecmdata->gwdata);
	if (zs == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	t3 = ecmdata->Ad4 = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Ad4 == NULL) goto oom;

	dbltogw (&ecmdata->gwdata, ecmdata->sigma, zs);
	gwsquare2 (&ecmdata->gwdata, zs, xs, 0);	/* s^2 */
	gwsmalladd (&ecmdata->gwdata, -5.0, xs);	/* u = s^2 - 5 */
	gwsmallmul (&ecmdata->gwdata, 4.0, zs);		/* v = 4*s */
	if (xz != NULL) {
		gwsquare2 (&ecmdata->gwdata, xs, xz->x, GWMUL_STARTNEXTFFT);
		gwsafemul (&ecmdata->gwdata, xs, xz->x);	/* x = u^3 */
		gwsquare2 (&ecmdata->gwdata, zs, xz->z, GWMUL_STARTNEXTFFT);
		gwsafemul (&ecmdata->gwdata, zs, xz->z);	/* z = v^3 */
	}

	/* Now for A. */
	gwsub3o (&ecmdata->gwdata, zs, xs, t2, GWADD_GUARANTEED_OK | GWADD_SQUARE_INPUT);
	gwsquare2 (&ecmdata->gwdata, t2, t3, GWMUL_STARTNEXTFFT);
	gwmul (&ecmdata->gwdata, t3, t2);						/* (v-u)^3 */
	gwadd3o (&ecmdata->gwdata, xs, xs, t3, GWADD_GUARANTEED_OK | GWADD_MUL_INPUT);
	gwadd3o (&ecmdata->gwdata, xs, t3, t3, GWADD_GUARANTEED_OK | GWADD_MUL_INPUT);
	gwadd3o (&ecmdata->gwdata, zs, t3, t3, GWADD_GUARANTEED_OK | GWADD_MUL_INPUT);	/* 3u+v */
	gwmul3 (&ecmdata->gwdata, t2, t3, t2, 0);					/* An = (v-u)^3 (3u+v) */
	gwmul3 (&ecmdata->gwdata, xs, zs, ecmdata->Ad4, GWMUL_STARTNEXTFFT);
	gwsquare2 (&ecmdata->gwdata, xs, xs, GWMUL_STARTNEXTFFT);
	gwmul3 (&ecmdata->gwdata, xs, ecmdata->Ad4, ecmdata->Ad4, 0);			/* u^3 * v */
	gwsmallmul (&ecmdata->gwdata, 4.0, ecmdata->Ad4);				/* An/Ad is now A + 2 */

	/* Normalize so that An is one */
	stop_reason = normalize (ecmdata, ecmdata->Ad4, t2, N, factor);
	if (stop_reason) return (stop_reason);

	/* For extra speed, precompute Ad * 4 */
	gwsmallmul (&ecmdata->gwdata, 4.0, ecmdata->Ad4);

	/* Even more speed, save FFT of Ad4 */
	gwfft (&ecmdata->gwdata, ecmdata->Ad4, ecmdata->Ad4);

/* Clean up temporaries */

	gwfree (&ecmdata->gwdata, xs);
	gwfree (&ecmdata->gwdata, zs);
	gwfree (&ecmdata->gwdata, t2);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Print message announcing the start of this curve */

void curve_start_msg (
	ecmhandle *ecmdata)
{
	char	buf[120];

	if (ecmdata->curve != 1) OutputStr (ecmdata->thread_num, "\n");
	sprintf (buf, "%s ECM curve #%ld", gwmodulo_as_string (&ecmdata->gwdata), ecmdata->curve);
	title (ecmdata->thread_num, buf);

	if (!ecmdata->optimal_B2 || ecmdata->state >= ECM_STATE_STAGE2)
		sprintf (buf,
			 "ECM on %s: curve #%ld with s=%.0f, B1=%" PRIu64 ", B2=%" PRIu64 "\n",
			 gwmodulo_as_string (&ecmdata->gwdata), ecmdata->curve, ecmdata->sigma, ecmdata->B, ecmdata->C);
	else
		sprintf (buf,
			 "ECM on %s: curve #%ld with s=%.0f, B1=%" PRIu64 ", B2=TBD\n",
			 gwmodulo_as_string (&ecmdata->gwdata), ecmdata->curve, ecmdata->sigma, ecmdata->B);
	OutputStr (ecmdata->thread_num, buf);
}

/* These routines manage the computing of Q^m in stage 2 */

int mQ_init (
	ecmhandle *ecmdata,
	uint64_t m)
{
	uint64_t mask;
	struct xz *V, *Vn, *Vn1, T;

/* Allocate the two multiples of D needed for a string of ell_adds */

	if (!alloc_xz (ecmdata, &ecmdata->Qm)) goto oom;
	if (!alloc_xz (ecmdata, &ecmdata->Qprevm)) goto oom;
	if (!alloc_xz (ecmdata, &T)) goto oom;

/* Compute Q * mD and Q * (m+1)D) */
/* We know that m is divisible by D.  Find top bit in multiplier. */

	m = m / ecmdata->D;
	for (mask = 1ULL << 63; (m & mask) == 0; mask >>= 1);

/* Init V values processing second mask bit.  Caller has computed Q * D and stored it in QD. */

	V = &ecmdata->QD;
	Vn = &ecmdata->Qprevm;
	Vn1 = &ecmdata->Qm;
	mask >>= 1;
	if (m & mask) {
		ell_dbl_xz_scr (ecmdata, V, Vn1, &T);				// V2, scratch T
		ell_add_xz_scr (ecmdata, V, Vn1, V, Vn, Vn);			// next n   = V3 = V1 + V2, diff 1, scratch Vn
		ell_dbl_xz_scr (ecmdata, Vn1, Vn1, &T);				// next n+1 = V4 = V2 + V2, scratch T
	}
	else {
		ell_dbl_xz_scr (ecmdata, V, Vn, Vn1);				// next n   = V2, scratch Vn1
		ell_add_xz_scr (ecmdata, V, Vn, V, Vn1, Vn1);			// next n+1 = V3 = V1 + V2, diff 1, scratch Vn1
	}

/* Loop processing one multiplier bit at a time.  Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1} */

	for (mask >>= 1; mask; mask >>= 1) {
	        if (m & mask) {
			ell_add_xz_scr (ecmdata, Vn, Vn1, V, Vn, &T);		// next n   = n + (n+1), diff 1, scratch T
			ell_dbl_xz_scr (ecmdata, Vn1, Vn1, &T);			// next n+1 = (n+1) + (n+1), scratch T
		}
		else {
			ell_add_xz_scr (ecmdata, Vn, Vn1, V, Vn1, &T);		// next n+1 = n + (n+1), diff 1, scratch T
			ell_dbl_xz_scr (ecmdata, Vn, Vn, &T);			// next n   = n + n, scratch T
		}
	}
	free_xz (ecmdata, &T);

/* There will be no more ell_dbl calls */

	gwfree (&ecmdata->gwdata, ecmdata->Ad4);
	ecmdata->Ad4 = NULL;

/* Initialize Qm_state for mQ_next to use as it sees fit */

	ecmdata->Qm_state = 0;

/* Init the arrays used in pooled normalizes of mQx values */

	if (ecmdata->TWO_FFT_STAGE2) {
		int	i;
		ecmdata->mQx = (gwnum *) malloc (ecmdata->E * sizeof (gwnum));
		if (ecmdata->mQx == NULL) goto oom;
		for (i = 0; i < ecmdata->E; i++) ecmdata->mQx[i] = NULL;
		ecmdata->mQx_count = 0;
	}
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
int mQ_next (
	ecmhandle *ecmdata,
	gwnum	*retx,
	gwnum	*retz,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	int	stop_reason;

/* The non-normalized case - multiply the last Q^m value by Q^D to get the next Q^m value */

	if (!ecmdata->TWO_FFT_STAGE2) {
		// On first call return the Q^mD value calculated by mQ_init and stored in Qprevm
		if (ecmdata->Qm_state == 0) {
			ecmdata->Qm_state = 1;
			*retx = ecmdata->Qprevm.x;
			*retz = ecmdata->Qprevm.z;
		}
		// On second call return the Q^(m+1)D value calculated by mQ_init and stored in Qm
		else if (ecmdata->Qm_state == 1) {
			ecmdata->Qm_state = 2;
			*retx = ecmdata->Qm.x;
			*retz = ecmdata->Qm.z;
		}
		// Compute next Q^(m+i)D value and return it from Qm
		else {
			stop_reason = ell_add_xz_noscr (ecmdata, &ecmdata->Qm, &ecmdata->QD, &ecmdata->Qprevm, &ecmdata->Qprevm);
			if (stop_reason) return (stop_reason);
			xzswap (ecmdata->Qm, ecmdata->Qprevm);
			*retx = ecmdata->Qm.x;
			*retz = ecmdata->Qm.z;
		}
		// Return success
		return (0);
	}

/* The normalized case - batch up a bunch of Q^m values and normalize them.  Then return them one at a time. */
/* Obviously retz need not be returned since it is always one. */

	if (ecmdata->mQx_count == 0) {
		for ( ; ecmdata->mQx_count < ecmdata->E; ecmdata->mQx_count++) {
			struct xz nextQm;
			nextQm.x = NULL;
			nextQm.z = NULL;

			// If no more Qm's need to be calculated, break out of loop if we've used up all the Qm values else use the NULL nextQm
			if (ecmdata->Dsection + ecmdata->mQx_count + 2 >= ecmdata->numDsections) {
				if (ecmdata->Qprevm.x == NULL) break;
			}

			// Allocate memory for nextQm, compute next multiple of D
			else {
				// Convoluted allocation scheme for nextQm (reduces peak mem usage a hair).  Allocate nextQm from the mQx array.
				// This re-uses all the mQx gwnums, then we allocate as necessary.  Doing so means there are no unused mQx entries
				// when we compute the very last Qm value.
				if (ecmdata->mQx_count*2 < ecmdata->E) {
					nextQm.x = ecmdata->mQx[ecmdata->mQx_count*2];
					ecmdata->mQx[ecmdata->mQx_count*2] = NULL;
				}
				if (ecmdata->mQx_count*2+1 < ecmdata->E) {
					nextQm.z = ecmdata->mQx[ecmdata->mQx_count*2+1];
					ecmdata->mQx[ecmdata->mQx_count*2+1] = NULL;
				}
				if (nextQm.x == NULL) {
					nextQm.x = gwalloc (&ecmdata->gwdata);
					if (nextQm.x == NULL) return (OutOfMemory (ecmdata->thread_num));
				}
				if (nextQm.z == NULL) {
					nextQm.z = gwalloc (&ecmdata->gwdata);
					if (nextQm.z == NULL) return (OutOfMemory (ecmdata->thread_num));
				}

				// Now compute next multiple of D
				ell_add_xz (ecmdata, &ecmdata->Qm, &ecmdata->QD, &ecmdata->Qprevm, &nextQm);
			}

			// Before discarding Qprevm add it to the normalize pool
			ecmdata->mQx[ecmdata->mQx_count] = ecmdata->Qprevm.x;
			stop_reason = add_to_normalize_pool_ro (ecmdata, ecmdata->mQx[ecmdata->mQx_count], ecmdata->Qprevm.z, factor);
			if (stop_reason) return (stop_reason);
			if (*factor != NULL) return (0);

			// Shuffle data along
			ecmdata->Qprevm = ecmdata->Qm;
			ecmdata->Qm = nextQm;
		}

		// Normalize the entire pool, freeing some gwnums (amount depends on pool type)
		stop_reason = normalize_pool (ecmdata, N, factor);
		if (stop_reason) return (stop_reason);
		if (*factor != NULL) return (0);
		mallocFreeForOS ();

		// Return first entry in the mQx array
		ecmdata->mQx_retcount = 0;
	}

	/// Return next entry from the pool
	*retx = ecmdata->mQx[ecmdata->mQx_retcount++];
	ecmdata->mQx_count--;
	return (0);
}
void mQ_term (
	ecmhandle *ecmdata)
{
	int	i;
	free_xz (ecmdata, &ecmdata->Qm);
	free_xz (ecmdata, &ecmdata->Qprevm);
	free_xz (ecmdata, &ecmdata->QD);
	if (ecmdata->mQx != NULL) {
		for (i = 0; i < ecmdata->E; i++) gwfree (&ecmdata->gwdata, ecmdata->mQx[i]);
		free (ecmdata->mQx); ecmdata->mQx = NULL;
	}
}

/* Record the amount of memory being used by this thread.  Until we get to stage 2, ECM uses 9 gwnums (x, z, AD4, 6 for ell_mul). */

void ecm_stage1_memory_usage (
	int	thread_num,
	ecmhandle *ecmdata)
{
	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&ecmdata->gwdata, 9));
}

/* Cost out a stage 2 plan with the given D, normalize_pool algorithm, and 2 vs. 4 FFT stage 2 setting. */

struct ecm_stage2_cost_data {
	int	impl;		/* 4 possible implementations.  2 vs. 4 FFT, N^2 vs. 3-mult pooling. */
	int	numvals;	/* Total number of gwnum temporaries available */
	double	gcd_cost;	/* Cost (in FFTs) of a GCD */
	double	modinv_cost;	/* Cost (in FFTs) of a modular inverse */
				/* Data returned from cost function follows */
	int	D;		/* D value for big steps */
	int	totrels;	/* Total number of relative primes used for pairing */
	int	stage2_numvals;	/* Total number of gwnum temporaries to be used in stage 2 */
	uint64_t B2_start;	/* Stage 2 start */
	uint64_t numDsections;	/* Number of D sections to process */
	uint64_t bitarraymaxDsections; /* Maximum number of D sections per bit array (in case work bit array must be split) */
	int	pool_type;	/* Modular inverse pooling implementation */
	int	TWO_FFT_STAGE2;	/* 2 vs. 4 FFT implementation */
	int	E;		/* In 2-FFT stage 2, number of D sections pooled for a single modular inverse */
};

double ecm_stage2_cost (
	int	D,		/* Stage 2 big step */
	uint64_t B2_start,	/* Stage 2 start point */
	uint64_t numDsections,	/* Number of stage 2 D sections */
	uint64_t bitarraymaxDsections, /* Maximum number of D sections per bit array (in case work bit array must be split) */
	int	numrels,	/* Number of relative primes less than D */
	double	multiplier,	/* Multiplier * numrels relative primes will be used in stage 2 */
	double	numpairs,	/* Estimated number of paired primes in stage 2 */
	double	numsingles,	/* Estimated number of prime singles in stage 2 */
	void	*data)		/* ECM specific costing data */
{
	int	totrels;	/* Total number of relative primes in use */
	int	e;
	double	numgcdsections;
	double	cost;
	struct ecm_stage2_cost_data *cost_data = (struct ecm_stage2_cost_data *) data;

/* MQ_init requires B2_start is at least 4 times D */

	if (B2_start < 4 * D) return (1.0e99);

/* Define constants for the number of transforms for various operations. */

#define ELL_ADD_COST		12
#define N_SQUARED_POOL_COST(n)	(((double)((n)+3)) * ((double)((n)+3)) - 13.0 + cost_data->modinv_cost)
#define MULT3_POOL_COST(n)	((double)(n) * 9 + cost_data->modinv_cost)
#define MULT3POINT44_POOL_COST(n) ((double)(n) * 10.3333333333333 + cost_data->modinv_cost)
#define MULT3POINT57_POOL_COST(n) ((double)(n) * 10.7142857142857 + cost_data->modinv_cost)
#define NQX_SETUP_COST(fft4)	((((fft4) || D % 3 != 0) ? 1.0 : 0.666666666667) * multiplier * ((double) D / 4) * ELL_ADD_COST)

/* Calculate the number of temporaries used for relative primes.  This will leave some available for modinv pooling. */

	totrels = (int) (numrels * multiplier);

/* Cost (in FFTs) of a 4FFT stage 2 using this D		*/
/* The cost will be:						*/
/*	multiplier * D / 4 ell_add_xzs + pool_cost (totrels) +	*/
/*	totrels (each nQx must be FFTed) +			*/
/*	(C-B)/D ell_add_xzs +					*/
/*	(C-B)/D * 2 (each mQx mQz must be FFTed) +		*/
/*	#primes*4						*/
/* The memory consumed will be:					*/
/*	9 + totrels gwnums in main stage 2 loop if N^2 pooling	*/
/* or	8 + totrels*2 gwnums in nQx setup if 3N pooling		*/
/* or	8 + totrels*1.33+2 gwnums in nQx setup if 3.44N pooling */
/* or	8 + totrels*1.14+5 gwnums in nQx setup if 3.57N pooling */

//GW:  what about doing the pool in chunks O(n^2) grows quickly!  There is a best "e" value
	if (cost_data->impl == 0) {
		cost = NQX_SETUP_COST (TRUE) + N_SQUARED_POOL_COST (totrels) + totrels +
			(double) numDsections * ELL_ADD_COST + numDsections * 2 +
			(numpairs + numsingles) * 4;
		cost_data->TWO_FFT_STAGE2 = FALSE;
		cost_data->pool_type = POOL_N_SQUARED;
		cost_data->E = 0;
		cost_data->stage2_numvals = 9 + totrels;
	}

	if (cost_data->impl == 1) {
		// We need 8 + totrels * 2 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (8 + totrels * 2 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations
		cost = NQX_SETUP_COST (TRUE) + MULT3_POOL_COST (totrels) + totrels +
			(double) numDsections * ELL_ADD_COST + numDsections * 2 +
			(numpairs + numsingles) * 4;
		cost_data->TWO_FFT_STAGE2 = FALSE;
		cost_data->pool_type = POOL_3MULT;
		cost_data->E = 0;
		cost_data->stage2_numvals = 8 + totrels * 2;
	}

	if (cost_data->impl == 2) {
		// We need 8 + totrels * 1.33 + 2 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (8 + totrels * 4 / 3 + 2 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations
		cost = NQX_SETUP_COST (TRUE) + MULT3POINT44_POOL_COST (totrels) + totrels +
			(double) numDsections * ELL_ADD_COST + numDsections * 2 +
			(numpairs + numsingles) * 4;
		cost_data->TWO_FFT_STAGE2 = FALSE;
		cost_data->pool_type = POOL_3POINT44MULT;
		cost_data->E = 0;
		cost_data->stage2_numvals = 8 + totrels * 4 / 3 + 2;
	}

	if (cost_data->impl == 3) {
		// We need 8 + totrels * 1.14 + 5 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (8 + totrels * 8 / 7 + 5 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations
		cost = NQX_SETUP_COST (TRUE) + MULT3POINT57_POOL_COST (totrels) + totrels +
			(double) numDsections * ELL_ADD_COST + numDsections * 2 +
			(numpairs + numsingles) * 4;
		cost_data->TWO_FFT_STAGE2 = FALSE;
		cost_data->pool_type = POOL_3POINT57MULT;
		cost_data->E = 0;
		cost_data->stage2_numvals = 8 + totrels * 8 / 7 + 5;
	}

/* Cost (in FFTs) of a 2FFT stage 2 using this D		*/
/* The cost will be:						*/
/*	multiplier * D / 6 ell_add_xzs + pool_cost (totrels) +	*/
/*	totrels (each nQx must be FFTed) +			*/
/*	(C-B)/D ell_add_xzs +					*/
/*	(C-B)/D (each mQx must be FFTed) +			*/
/*	(C-B)/D/E * (pool_cost (e) + modinv_cost) +		*/
/*	#primes*2						*/
/* The memory consumed will be:					*/
/*	max (13 + totrels + 1, 9 + totrels + e + 1) gwnums if N^2 pooling */
/* or	max (13 + totrels * 2, 9 + totrels + e * 2) gwnums if 3N pooling */
/* or	max (13 + totrels * 2, 9 + totrels + e * 4 / 3 + 2) gwnums if 3.44N pooling */
/* or	max (13 + totrels * 2, 9 + totrels + e * 8 / 7 + 5) gwnums if 3.57N pooling */

//GW: feedback loop that might make optimizing faster.  choose e such that every gcd section is the same size
	if (cost_data->impl == 4) {
		int	beste = (int) sqrt (2 * (cost_data->modinv_cost + 9)); // Best value for E using O(N^2) pooling
		// We need 13 + totrels + 1 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (13 + totrels + 1 > cost_data->numvals) return (1.0e99);
		e = cost_data->numvals - 9 - totrels - 1;
		if (e <= 0) return (1.0e99);
		if (e > beste) e = beste;
		numgcdsections = ceil ((double) numDsections / e);
		cost = NQX_SETUP_COST (FALSE) + N_SQUARED_POOL_COST (totrels) + totrels +
		       (double) numDsections * (ELL_ADD_COST + 1) +
		       numgcdsections * N_SQUARED_POOL_COST (e) +
		       (numpairs + numsingles) * 2.0;
		cost_data->TWO_FFT_STAGE2 = TRUE;
		cost_data->pool_type = POOL_N_SQUARED;
		cost_data->E = e;
		cost_data->stage2_numvals = _intmax (13 + totrels + 1, 9 + totrels + e + 1);
	}

	if (cost_data->impl == 5) {
		// We need 13 + totrels * 2 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (13 + totrels * 2 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations when totrels > E
		e = (cost_data->numvals - 9 - totrels) / 2;
		numgcdsections = ceil ((double) numDsections / e);
		e = (unsigned long) ceil ((double) numDsections / numgcdsections);
		cost = NQX_SETUP_COST (FALSE) + MULT3_POOL_COST (totrels) + totrels +
		       (double) numDsections * (ELL_ADD_COST + 1) +
		       numgcdsections * MULT3_POOL_COST (e) +
		       (numpairs + numsingles) * 2.0;
		cost_data->TWO_FFT_STAGE2 = TRUE;
		cost_data->pool_type = POOL_3MULT;
		cost_data->E = e;
		cost_data->stage2_numvals = _intmax (13 + totrels * 2, 9 + totrels + e * 2);
	}

	if (cost_data->impl == 6) {
		// We need 13 + totrels * 4 / 3 + 2 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (13 + totrels * 4 / 3 + 2 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations when totrels > E
// GW: This e calculation may need some refinement
		e = (cost_data->numvals - 9 - totrels - 2) * 3 / 4;
		numgcdsections = ceil ((double) numDsections / e);
		e = (unsigned long) ceil ((double) numDsections / numgcdsections);
		cost = NQX_SETUP_COST (FALSE) + MULT3POINT44_POOL_COST (totrels) + totrels +
		       (double) numDsections * (ELL_ADD_COST + 1) +
		       numgcdsections * MULT3POINT44_POOL_COST (e) +
		       (numpairs + numsingles) * 2.0;
		cost_data->TWO_FFT_STAGE2 = TRUE;
		cost_data->pool_type = POOL_3POINT44MULT;
		cost_data->E = e;
		cost_data->stage2_numvals = _intmax (13 + totrels * 4 / 3 + 2, 9 + totrels + e * 4 / 3 + 2);
	}

	if (cost_data->impl == 7) {
		// We need 13 + totrels * 8 / 7 + 5 available gwnums for the nQx setup.  Return infinite cost if there aren't enough temporaries available.
		if (13 + totrels * 8 / 7 + 5 > cost_data->numvals) return (1.0e99);
//GW: seems shameful to dedicate many temporaries to use-once nQx setup computations when totrels > E
// GW: This e calculation may need some refinement
		e = (cost_data->numvals - 9 - totrels - 5) * 7 / 8;
		numgcdsections = ceil ((double) numDsections / e);
		e = (unsigned long) ceil ((double) numDsections / numgcdsections);
		cost = NQX_SETUP_COST (FALSE) + MULT3POINT57_POOL_COST (totrels) + totrels +
		       (double) numDsections * (ELL_ADD_COST + 1) +
		       numgcdsections * MULT3POINT57_POOL_COST (e) +
		       (numpairs + numsingles) * 2.0;
		cost_data->TWO_FFT_STAGE2 = TRUE;
		cost_data->pool_type = POOL_3POINT57MULT;
		cost_data->E = e;
		cost_data->stage2_numvals = _intmax (13 + totrels * 8 / 7 + 5, 9 + totrels + e * 8 / 7 + 5);
	}

/* Return data ECM implementation will need */

	cost_data->D = D;
	cost_data->totrels = totrels;
	cost_data->B2_start = B2_start;
	cost_data->numDsections = numDsections;
	cost_data->bitarraymaxDsections = bitarraymaxDsections;

/* Return the cost of executing this stage 2 */

	return (cost);
}

/* Given a number of temporaries derived from a memory limit, choose best value for D, 2 or 4 FFT stage 2, and best algorithm for normalize_pool. */
/* Returns the cost (in FFTs) of implementing stage 2 as well as other cost data needed to implement stage 2 */

double ecm_stage2_impl_given_numvals (
	ecmhandle *ecmdata,
	int	numvals,				/* Number of gwnum temporaries available */
	struct ecm_stage2_cost_data *return_cost_data)	/* Returned extra data from ECM costing function */
{
	int	impl;				/* Four possible stage 2 implementations */
	double	cost, best_cost;		/* Best cost for each of the 4 possible stage 2 implementations */
	struct ecm_stage2_cost_data cost_data;	/* Extra data passed to and returned from ECM costing function */

/* If memory is really tight, then the 4 FFT - O(n^2) pooling is the most memory efficient ECM implementation.  This will be our default */
/* plan.  Note: D=30 (4 nQx values) requires 17 gwnums.  The next D value (60) requires 21 gwnums. */

//GW: Test for at least minimum numvals?

/* Calculate the GCD and modular inverse cost in terms of number of transforms. */
/* The costs come from the timing code running ECM on M604 and spreadsheeting. */
/* Since GCDs are single-threaded we double the costs for multi-threaded ECM runs. */

	cost_data.gcd_cost = 320.53 * log (ecmdata->gwdata.bit_length) - 3302.0;
	if (cost_data.gcd_cost < 100.0) cost_data.gcd_cost = 100.0;
	cost_data.modinv_cost = 570.16 * log (ecmdata->gwdata.bit_length) - 6188.4;
	if (cost_data.modinv_cost < 1.25 * cost_data.gcd_cost) cost_data.modinv_cost = 1.25 * cost_data.gcd_cost;
	if (gwget_num_threads (&ecmdata->gwdata)) cost_data.gcd_cost *= 2.0, cost_data.modinv_cost *= 2.0;

/* Find the least costly stage 2 plan looking at the four combinations of 2 vs. 4 FFT and pool type N^2 vs. 3-MULT vs. 3.44-MULT vs. 3.57-MULT. */

	best_cost = 1.0e99;
	for (impl = 0; impl <= 7; impl++) {

/* Check for QA'ing a specific ECM implementation type */

		if (QA_TYPE != 0 && QA_TYPE != impl + 1) continue;

/* Cost out an ECM stage 2 implementation.  2 vs. 4 FFT, N^2 vs. 3-MULT vs. 3.44-MULT vs. 3.57-MULT pooling.  All implementations must account for 9 gwnum */
/* temporaries required by the main stage 2 loop (6 for computing mQx, gg, 2 for ell_add_xz_noscr temps).  Keep track of the best implementation. */

		cost_data.impl = impl;
		cost_data.numvals = numvals;
		cost = best_stage2_impl (ecmdata->C_done, ecmdata->C, numvals - 9, &ecm_stage2_cost, &cost_data);
		if (cost < best_cost) {
			best_cost = cost;
			*return_cost_data = cost_data;
		}

//GW:  play with e = num modinvs.  That is, break the N_SQUARED pool in half or thirds, etc.
//			do we need a binary search on breaking up the 3N pooling into multiple segments?
	}

/* Return our best implementation */

	return (best_cost);
}

/* Choose most effective B2 for an ECM curve given number of gwnums we are allowed to allocate */

void ecm_choose_B2 (
	ecmhandle *ecmdata,
	unsigned long numvals)
{
	struct ecm_stage2_cost_data cost_data;	/* Extra data passed to ECM costing function */
	double	B1_cost, B2_cost;
	struct ecm_stage2_efficiency {
		int	i;
		double	efficiency;
	} best[3];
	char	buf[100];

// From Alex Kruppa, master of all things ECM, the following formula compensates for using B2 values that are not 100 * B1.
// curve_worth = 0.11 + 0.89 * (log10(B2 / B1) / 2) ^ 1.5

#define kruppa(x,B2mult)	x.i = B2mult; \
				ecmdata->C = x.i * ecmdata->B; \
				B2_cost = ecm_stage2_impl_given_numvals (ecmdata, numvals, &cost_data); \
				x.efficiency = kruppa_adjust_ratio (x.i) / (B1_cost + B2_cost + cost_data.gcd_cost);

// The cost of B1 (in FFTs) is about 25.48 * B1 (measured at 25.42 for B1=50000, 25.53 for B1=250000).

	B1_cost = 25.48 * (double) ecmdata->B;

/* Look for the best B2 which is likely between 50 * B1 and 150 * B1.  If optimal is not between these bounds, don't worry */
/* we'll locate the optimal spot anyway. */

	kruppa (best[0], 50);
	kruppa (best[1], 100);
	kruppa (best[2], 150);

/* Handle case where midpoint has worse efficiency than the start point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[0].efficiency > best[1].efficiency) {
		best[2] = best[1];
		kruppa (best[1], (best[0].i + best[2].i) / 2);
	}

/* Handle case where midpoint has worse efficiency than the end point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[2].efficiency > best[1].efficiency) {
		best[0] = best[1];
		best[1] = best[2];
		kruppa (best[2], best[1].i * 2);
	}

/* Find the best B2.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].i + 2 != best[2].i) {
		struct ecm_stage2_efficiency midpoint;

		ASSERTG (best[1].efficiency >= best[0].efficiency);
		ASSERTG (best[1].efficiency >= best[2].efficiency);

		// Work on the bigger of the lower section and upper section
		if (best[1].i - best[0].i > best[2].i - best[1].i) {		// Work on lower section
			kruppa (midpoint, (best[0].i + best[1].i) / 2);
			if (midpoint.efficiency > best[1].efficiency) {		// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			kruppa (midpoint, (best[1].i + best[2].i) / 2);
			if (midpoint.efficiency > best[1].efficiency) {		// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Return the best B2 */

	ecmdata->C = best[1].i * ecmdata->B;
	sprintf (buf, "Optimal B2 is %d*B1 = %" PRIu64 ".\n", best[1].i, ecmdata->C);
	OutputStr (ecmdata->thread_num, buf);
}

/* Using the current memory limit, choose best value for D, 2 or 4 FFT stage 2, and best algorithm for normalize_pool. */

int ecm_stage2_impl (
	ecmhandle *ecmdata)
{
	unsigned int memory, min_memory, desired_memory;	/* Memory is in MB */
	unsigned int bitarray_size, max_bitarray_size;		/* Bitarray memory in MB */
	int	numvals;			/* Number of gwnums we can allocate */
	struct ecm_stage2_cost_data cost_data;	/* Extra data passed returned from ECM costing function */
	int	stop_reason;

/* Calculate (roughly) the memory required for the bit-array.  A typical stage 2 implementation is D=210, B2_start=B2/11, numrels=24, totrels=7*24. */
/* A typical optimal B2 is 130*B1.  This means bitarray memory = B2 * 10/11 / 210 * 24 * 7 / 8 bytes. */

	bitarray_size = divide_rounding_up ((unsigned int) ((ecmdata->optimal_B2 ? 130 * ecmdata->B : ecmdata->C) * 1680 / 18480), 1 << 20);
	max_bitarray_size = IniGetInt (INI_FILE, "MaximumBitArraySize", 250);
	if (bitarray_size > max_bitarray_size) bitarray_size = max_bitarray_size;

/* Get available memory.  We need 13 gwnums to do the smallest stage 2. */
/* If continuing from a stage 2 save file then we desire as many temporaries as the save file used.  Otherwise, assume 120 temporaries will */
/* allow us to do a reasonable efficient stage 2 implementation. */
//gw: is 120 reasonable?

	min_memory = bitarray_size + cvt_gwnums_to_mem (&ecmdata->gwdata, 13);
	desired_memory = bitarray_size + cvt_gwnums_to_mem (&ecmdata->gwdata, ecmdata->state >= ECM_STATE_STAGE2 ? ecmdata->stage2_numvals : 120);
	stop_reason = avail_mem (ecmdata->thread_num, min_memory, desired_memory, &memory);
	if (stop_reason) return (stop_reason);

/* Factor in the multiplier that we set to less than 1.0 when we get unexpected memory allocation errors. */
/* Make sure we can still allocate 13 temporaries. */

	memory = (unsigned int) (ecmdata->pct_mem_to_use * (double) memory);
	if (memory < min_memory) return (avail_mem_not_sufficient (ecmdata->thread_num, min_memory, desired_memory));
	if (memory < 8) memory = 8;

/* Output a message telling us how much memory is available */

	if (NUM_WORKER_THREADS > 1) {
		char	buf[100];
		sprintf (buf, "Available memory is %dMB.\n", memory);
		OutputStr (ecmdata->thread_num, buf);
	}

/* Figure out how many gwnum values fit in our memory limit.  User nordi had over-allocating memory troubles on Linux testing M1277, presumably */
/* because our estimate genum size was too low.  As a work-around limit numvals to 100,000 by default. */

	numvals = cvt_mem_to_gwnums (&ecmdata->gwdata, memory - bitarray_size);
	if (numvals < 13) numvals = 13;
	if (numvals > 100000) numvals = 100000;

/* Set C_done for future best_stage2_impl calls. */
/* Override B2 with optimal B2 based on amount of memory available. */

	if (ecmdata->state == ECM_STATE_MIDSTAGE) {
		ecmdata->C_done = ecmdata->B;
		if (ecmdata->optimal_B2) ecm_choose_B2 (ecmdata, numvals);
	}

/* If are continuing from a save file that was in stage 2, check to see if we currently have enough memory to continue with the save file's */
/* stage 2 implementation.  Also check if we now have significantly more memory available and stage 2 is not near complete such that a new */
/* stage 2 implementation might give us a faster stage 2 completion. */

//GW: These are rather arbitrary heuristics
	if (ecmdata->state >= ECM_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= ecmdata->stage2_numvals &&				// We have enough memory and
	    (numvals < ecmdata->stage2_numvals * 2 ||				// less than twice as much memory now available or
	     ecmdata->Dsection >= ecmdata->numDsections / 2))			// stage 2 more than half done
		return (0);							// Use old plan

/* Find the least costly stage 2 plan */

	ecm_stage2_impl_given_numvals (ecmdata, numvals, &cost_data);

/* If are continuing from a save file that was in stage 2 and the new plan doesn't look significant better than the old plan, then */
/* we use the old plan and its partially completed bit array. */

	if (ecmdata->state >= ECM_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= ecmdata->stage2_numvals &&				// We have enough memory and
	    cost_data.stage2_numvals < ecmdata->stage2_numvals * 2)		// new plan does not use significantly more memory
		return (0);							// Use old plan

/* If are continuing from a save file that was in stage 2, toss the save file's bit array. */

	if (ecmdata->state >= ECM_STATE_STAGE2) {
		free (ecmdata->bitarray);
		ecmdata->bitarray = NULL;
	}

/* Set all the variables needed for this stage 2 plan */

	ecmdata->stage2_numvals = cost_data.stage2_numvals;
	ecmdata->totrels = cost_data.totrels;
	ecmdata->TWO_FFT_STAGE2 = cost_data.TWO_FFT_STAGE2;
	ecmdata->pool_type = cost_data.pool_type;
	ecmdata->D = cost_data.D;
	ecmdata->E = cost_data.E;
	ecmdata->B2_start = cost_data.B2_start;
	ecmdata->numDsections = cost_data.numDsections;
	ecmdata->bitarraymaxDsections = cost_data.bitarraymaxDsections;
	ecmdata->Dsection = 0;

	if (ecmdata->state < ECM_STATE_STAGE2) {
		ecmdata->first_relocatable = ecmdata->B;
		ecmdata->last_relocatable = ecmdata->B2_start;
		ecmdata->C_done = ecmdata->B2_start;
	}

/* Once the plan is selected, create a bit array maximizing the prime pairings. */

	ecmdata->bitarrayfirstDsection = ecmdata->Dsection;
	stop_reason = fill_work_bitarray (ecmdata->thread_num, &ecmdata->sieve_info, ecmdata->D, ecmdata->totrels,
					  ecmdata->first_relocatable, ecmdata->last_relocatable, ecmdata->C_done, ecmdata->C,
					  ecmdata->bitarraymaxDsections, &ecmdata->bitarray);
	if (stop_reason) return (stop_reason);

	return (0);
}

/* Routines to create and read save files for an ECM factoring job */

#define ECM_MAGICNUM	0x1725bcd9
#define ECM_VERSION	2

void ecm_save (
	ecmhandle *ecmdata)
{
	int	fd;
	struct work_unit *w = ecmdata->w;
	unsigned long sum = 0;

/* Create the intermediate file */

	fd = openWriteSaveFile (&ecmdata->write_save_file_state);
	if (fd < 0) return;

/* Write the file header. */

	if (! write_header (fd, ECM_MAGICNUM, ECM_VERSION, w)) goto writeerr;

/* Write the file data */

	if (! write_long (fd, ecmdata->curve, &sum)) goto writeerr;
	if (! write_longlong (fd, ecmdata->average_B2, NULL)) goto writeerr;
	if (! write_int (fd, ecmdata->state, &sum)) goto writeerr;
	if (! write_double (fd, ecmdata->sigma, NULL)) goto writeerr;
	if (! write_longlong (fd, ecmdata->B, &sum)) goto writeerr;
	if (! write_longlong (fd, ecmdata->C, &sum)) goto writeerr;

	if (ecmdata->state == ECM_STATE_STAGE1) {
		if (! write_longlong (fd, ecmdata->stage1_prime, &sum)) goto writeerr;
	}

	else if (ecmdata->state == ECM_STATE_MIDSTAGE) {
	}

	// Save everything necessary to restart stage 2 without calling ecm_stage2_impl again
	else if (ecmdata->state == ECM_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! write_int (fd, ecmdata->stage2_numvals, &sum)) goto writeerr;
		if (! write_int (fd, ecmdata->totrels, &sum)) goto writeerr;
		if (! write_int (fd, ecmdata->D, &sum)) goto writeerr;
		if (! write_int (fd, ecmdata->E, &sum)) goto writeerr;
		if (! write_int (fd, ecmdata->TWO_FFT_STAGE2, &sum)) goto writeerr;
		if (! write_int (fd, ecmdata->pool_type, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->first_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->last_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->B2_start, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->C_done, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->numDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->Dsection, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->bitarraymaxDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, ecmdata->bitarrayfirstDsection, &sum)) goto writeerr;
		// Output the truncated bit array
		bitarray_numDsections = ecmdata->numDsections - ecmdata->bitarrayfirstDsection;
		if (bitarray_numDsections > ecmdata->bitarraymaxDsections) bitarray_numDsections = ecmdata->bitarraymaxDsections; 		
		bitarray_len = divide_rounding_up (bitarray_numDsections * ecmdata->totrels, 8);
		bitarray_start = divide_rounding_down ((ecmdata->Dsection - ecmdata->bitarrayfirstDsection) * ecmdata->totrels, 8);
		if (! write_array (fd, ecmdata->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto writeerr;
	}

	else if (ecmdata->state == ECM_STATE_GCD) {
	}

/* Write the data values, write_gwnum can handle FFTed or partially FFTed */

	if (ecmdata->state == ECM_STATE_STAGE1) {
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.x, &sum)) goto writeerr;
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.z, &sum)) goto writeerr;
	}
	if (ecmdata->state == ECM_STATE_MIDSTAGE) {
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->QD.x, &sum)) goto writeerr;
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->QD.z, &sum)) goto writeerr;
	}
	if (ecmdata->state == ECM_STATE_STAGE2) {
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->nQx[0], &sum)) goto writeerr;
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->gg, &sum)) goto writeerr;
	}
	if (ecmdata->state == ECM_STATE_GCD) {
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.x, &sum)) goto writeerr;
		if (! write_gwnum (fd, &ecmdata->gwdata, ecmdata->gg, &sum)) goto writeerr;
	}

/* Write the checksum, we're done */

	if (! write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (&ecmdata->write_save_file_state, fd);
	return;

/* An error occurred.  Close and delete the current file. */

writeerr:
	deleteWriteSaveFile (&ecmdata->write_save_file_state, fd);
}

int ecm_old_restore (			/* For version 25 save files */
	ecmhandle *ecmdata,
	int	fd,
	unsigned long filesum)
{
	unsigned long state;
	unsigned long sum = 0;
	uint64_t savefile_B, unused64;

/* Read the file data */

	if (! read_long (fd, &state, &sum)) goto readerr;
	if (! read_long (fd, &ecmdata->curve, &sum)) goto readerr;
	if (! read_double (fd, &ecmdata->sigma, NULL)) goto readerr;
	if (! read_longlong (fd, &savefile_B, &sum)) goto readerr;
	if (! read_longlong (fd, &ecmdata->stage1_prime, &sum)) goto readerr;
	if (! read_longlong (fd, &unused64, &sum)) goto readerr;	// C_processed

/* Handle the should-never-happen case where we have a save file with a smaller bound #1 than the bound #1 we are presently working on. */
/* Restart the curve (and curve counts) from scratch. */

	if (savefile_B < ecmdata->B) {
		OutputStr (ecmdata->thread_num, "ECM save file created with smaller B1 value.  Save file cannot be used.\n");
		goto readerr;
	}

// Convert old state values into new state values.  Restart stage 2 for curves that were in the middle of stage 2. */

	if (state == 0) {
		ecmdata->state = ECM_STATE_STAGE1;
	}
	if (state == 1) {
		ecmdata->state = ECM_STATE_MIDSTAGE;
		OutputStr (ecmdata->thread_num, "Old ECM save file was in stage 2.  Restarting stage 2 from scratch.\n");
	}

/* Old save files did not store B2, assume all the curves were run with the current B2 */

	ecmdata->average_B2 = ecmdata->C;

/* Read the values */

	if (! read_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.x, &sum)) goto readerr;
	if (! read_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.z, &sum)) goto readerr;

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occurred.  Cleanup and return FALSE. */

readerr:
	_close (fd);
	return (FALSE);
}

int ecm_restore (			/* For version 30.4 and later save files */
	ecmhandle *ecmdata)
{
	int	fd;
	struct work_unit *w = ecmdata->w;
	unsigned long version;
	unsigned long sum = 0, filesum;
	uint64_t savefile_B;

/* Open the intermediate file */

	fd = _open (ecmdata->read_save_file_state.current_filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto err;

/* Read the file header */

	if (! read_magicnum (fd, ECM_MAGICNUM)) goto readerr;
	if (! read_header (fd, &version, w, &filesum)) goto readerr;
	if (version == 0 || version > ECM_VERSION) goto readerr;
	if (version == 1) return (ecm_old_restore (ecmdata, fd, filesum));

/* Read the file data */

	if (! read_long (fd, &ecmdata->curve, &sum)) goto readerr;
	if (! read_longlong (fd, &ecmdata->average_B2, NULL)) goto readerr;
	if (! read_int (fd, &ecmdata->state, &sum)) goto readerr;
	if (! read_double (fd, &ecmdata->sigma, NULL)) goto readerr;
	if (! read_longlong (fd, &savefile_B, &sum)) goto readerr;
//GW: should we be overwriting C with savefile_C?
	if (! read_longlong (fd, &ecmdata->C, &sum)) goto readerr;

/* Handle the case where we have a save file with a smaller bound #1 than the bound #1 we are presently working on. */
/* Restart the curve (and curve counts) from scratch. */

	if (savefile_B < ecmdata->B) {
		OutputStr (ecmdata->thread_num, "ECM save file created with smaller B1 value.  Save file cannot be used.\n");
		goto readerr;
	}

/* Read state dependent data */

	if (ecmdata->state == ECM_STATE_STAGE1) {
		if (! read_longlong (fd, &ecmdata->stage1_prime, &sum)) goto readerr;
	}

	else if (ecmdata->state == ECM_STATE_MIDSTAGE) {
	}

	else if (ecmdata->state == ECM_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! read_int (fd, &ecmdata->stage2_numvals, &sum)) goto readerr;
		if (! read_int (fd, &ecmdata->totrels, &sum)) goto readerr;
		if (! read_int (fd, &ecmdata->D, &sum)) goto readerr;
		if (! read_int (fd, &ecmdata->E, &sum)) goto readerr;
		if (! read_int (fd, &ecmdata->TWO_FFT_STAGE2, &sum)) goto readerr;
		if (! read_int (fd, &ecmdata->pool_type, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->first_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->last_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->B2_start, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->C_done, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->numDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->Dsection, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->bitarraymaxDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &ecmdata->bitarrayfirstDsection, &sum)) goto readerr;
		// Read the truncated bit array
		bitarray_numDsections = ecmdata->numDsections - ecmdata->bitarrayfirstDsection;
		if (bitarray_numDsections > ecmdata->bitarraymaxDsections) bitarray_numDsections = ecmdata->bitarraymaxDsections;
		bitarray_len = divide_rounding_up (bitarray_numDsections * ecmdata->totrels, 8);
		bitarray_start = divide_rounding_down ((ecmdata->Dsection - ecmdata->bitarrayfirstDsection) * ecmdata->totrels, 8);
		ecmdata->bitarray = (char *) malloc ((size_t) bitarray_len);
		if (ecmdata->bitarray == NULL) goto readerr;
		if (! read_array (fd, ecmdata->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto readerr;
	}

	else if (ecmdata->state == ECM_STATE_GCD) {
	}

/* Read the gwnum values (there are always two).  Caller will move them to where they need to be for each ECM state. */

	if (! read_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.x, &sum)) goto readerr;
	if (! read_gwnum (fd, &ecmdata->gwdata, ecmdata->xz.z, &sum)) goto readerr;

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occurred.  Cleanup and return FALSE. */

readerr:
	_close (fd);
err:
	return (FALSE);
}


/**************************************************************
 *
 *	Main ECM Function
 *
 **************************************************************/

int ecm (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	ecmhandle ecmdata;
	uint64_t sieve_start, next_prime;
	unsigned long SQRT_B;
	double	last_output, last_output_t, one_over_B;
	double	output_frequency, output_title_frequency;
	double	base_pct_complete, one_bit_pct;
	int	i, min_memory;
	unsigned int memused;
	char	filename[32], buf[255], JSONbuf[4000], fft_desc[200];
	int	res, stop_reason, stage, first_iter_msg;
	giant	N;		/* Number being factored */
	giant	factor;		/* Factor found, if any */
	int	msglen, continueECM, prpAfterEcmFactor;
	char	*str, *msg;
	double	timers[10];

/* Clear pointers to allocated memory */

	N = NULL;
	factor = NULL;
	str = NULL;
	msg = NULL;

/* Begin initializing ECM data structure */
/* MQ_init requires that B is at least 120 (4 times the minimum D) */
/* Choose a default value for the second bound if none was specified */

	memset (&ecmdata, 0, sizeof (ecmhandle));
	ecmdata.thread_num = thread_num;
	ecmdata.w = w;
	ecmdata.B = (uint64_t) w->B1;
	ecmdata.C = (uint64_t) w->B2;
	if (ecmdata.B < 120) {
		OutputStr (thread_num, "Using minimum bound #1 of 120\n");
		ecmdata.B = 120;
	}
	if (ecmdata.C == 0) ecmdata.C = ecmdata.B * 100;
	if (ecmdata.C <= ecmdata.B) ecmdata.C = ecmdata.B;
	ecmdata.pct_mem_to_use = 1.0;				// Use as much memory as we can unless we get allocation errors

/* Decide if we will calculate an optimal B2 when stage 2 begins */

	ecmdata.optimal_B2 = (!QA_IN_PROGRESS && ecmdata.C == 100 * ecmdata.B && IniGetInt (INI_FILE, "ECMBestB2", 1));

/* Little known option to use higher bounds than assigned by PrimeNet */

	if (!QA_IN_PROGRESS) {
		int	mult = IniGetInt (INI_FILE, "ECMBoundsMultiplier", 1);
		if (mult < 1) mult = 1;
		if (mult > 20) mult = 20;
		ecmdata.B *= mult;
		ecmdata.C *= mult;
	}

/* Compute the number we are factoring */

	stop_reason = setN (thread_num, w, &N);
	if (stop_reason) goto exit;

/* Other initialization */

	PRAC_SEARCH = IniGetInt (INI_FILE, "PracSearch", 7);
	if (PRAC_SEARCH < 1) PRAC_SEARCH = 1;
	if (PRAC_SEARCH > 50) PRAC_SEARCH = 50;

/* Clear all timers */

restart:
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Time the giants squaring and multiply code in order to select the */
/* best crossover points.  This should only be done in the release code */
/* (optimized giants library). */

/*#define TIMING1*/
#ifdef TIMING1
if (w->n == 598) {
int i, j;
giant	x, y, z, a, m;
#define TESTSIZE	200
RDTSC_TIMING = 12;
x = allocgiant(TESTSIZE); y = allocgiant (2*TESTSIZE);
z = allocgiant (2*TESTSIZE), a = allocgiant (2*TESTSIZE);
m = allocgiant (TESTSIZE);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE-1] &= 0x00FFFFFF;
m->n[TESTSIZE-1] &= 0x00FFFFFF;
for (i = TESTSIZE; i >= 10; i--) {
	x->sign = i;
	m->sign = i;
	setmulmode (GRAMMAR_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, y);
		start_timer (timers, 0);
		if (B&1) mulg (m, y);
		else squareg (y);
		end_timer (timers, 0);
		if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
		timers[0] = 0;
	}
	setmulmode (KARAT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, z);
		start_timer (timers, 0);
		if (B&1) mulg (m, z);
		else squareg (z);
		end_timer (timers, 0);
		if (timers[2] == 0 || timers[2] > timers[0]) timers[2] = timers[0];
		timers[0] = 0;
	}
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (timers, 0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (timers, 0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  G: ", i);
	print_timer (timers, 1, buf, TIMER_MS | TIMER_CLR);
	strcat (buf, ", K: ");
	print_timer (timers, 2, buf, TIMER_MS | TIMER_CLR);
	strcat (buf, ", F: ");
	print_timer (timers, 3, buf, TIMER_MS | TIMER_NL | TIMER_CLR);
	OutputBoth (thread_num, buf);
	if (gcompg (y, z) != 0)
		i--;
	if (gcompg (y, a) != 0)
		i--;
	Sleep (100);
}
return 0;
}

/* This code lets us time various giants FFT squarings and multiplies */

if (w->n == 601) {
int i, j;
giant	x, a, m;
#define TESTSIZE2	260000
RDTSC_TIMING = 12;
x = allocgiant(TESTSIZE2);
a = allocgiant (2*TESTSIZE2);
m = allocgiant (TESTSIZE2);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE2; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE2-1] &= 0x00FFFFFF;
m->n[TESTSIZE2-1] &= 0x00FFFFFF;
for (i = 30; i < TESTSIZE2/2; i<<=1) {
	x->sign = i;
	m->sign = i;
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (timers, 0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (timers, 0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  , F: ", i);
	print_timer (timers, 3, buf, TIMER_NL | TIMER_CLR | TIMER_MS);
	OutputStr (thread_num, buf);
	Sleep (100);
}
return 0;
}
#endif

/* Include timing code when building the debug version of prime95 */

#ifdef GDEBUG
if (w->n == 600) {
gwhandle gwdata;
void *workbuf;
int j, min_test, max_test, test, cnt, NUM_X87_TESTS, NUM_SSE2_TESTS, NUM_AVX_TESTS, NUM_AVX512_TESTS;
#define timeit(a,n,w) (((void**)a)[0]=w,((uint32_t*)a)[2]=n,gwtimeit(a))

gwinit (&gwdata); gwdata.cpu_flags &= ~CPU_AVX512F;
gwsetup (&gwdata, 1.0, 2, 10000000, -1);
workbuf = (void *) aligned_malloc (400000000, 4096);
memset (workbuf, 0, 400000000);
RDTSC_TIMING = 2;
min_test = IniGetInt (INI_FILE, "MinTest", 0);
max_test = IniGetInt (INI_FILE, "MaxTest", min_test);
NUM_X87_TESTS = timeit (gwdata.asm_data, -1, NULL);
NUM_SSE2_TESTS = timeit (gwdata.asm_data, -2, NULL);
NUM_AVX_TESTS = timeit (gwdata.asm_data, -3, NULL);
NUM_AVX512_TESTS = timeit (gwdata.asm_data, -4, NULL);
//SetThreadPriority (CURRENT_THREAD, THREAD_PRIORITY_TIME_CRITICAL);
for (j = 0; j < NUM_X87_TESTS + NUM_SSE2_TESTS + NUM_AVX_TESTS + NUM_AVX512_TESTS; j++) {
	cnt = 0;
	test = (j < NUM_X87_TESTS ? j :
		j < NUM_X87_TESTS + NUM_SSE2_TESTS ? 1000 + j - NUM_X87_TESTS :
		j < NUM_X87_TESTS + NUM_SSE2_TESTS + NUM_AVX_TESTS ? 2000 + j - NUM_X87_TESTS - NUM_SSE2_TESTS :
			3000 + j - NUM_X87_TESTS - NUM_SSE2_TESTS - NUM_AVX_TESTS);
	if (test == 3000 && CPU_FLAGS & CPU_AVX512F) {
		gwdone (&gwdata);
		gwinit (&gwdata);
		gwsetup (&gwdata, 1.0, 2, 10000000, -1);
	}
	if (min_test && (test < min_test || test > max_test)) continue;
	if (! (CPU_FLAGS & CPU_SSE2) && test >= 1000) break;
	if (! (CPU_FLAGS & CPU_AVX) && test >= 2000) break;
	if (! (CPU_FLAGS & CPU_AVX512F) && test >= 3000) break;
for (i = 1; i <= 50; i++) {
	start_timer (timers, 0);
	timeit (gwdata.asm_data, test, workbuf);
	end_timer (timers, 0);
	if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
	if (i > 1 && timers[0] < 3.0 * timers[1]) {
		if (timers[0] > 1.5 * timers[1])
			i++;
		timers[2] += timers[0];
		cnt++;
	}
	timers[0] = 0;
}
sprintf (buf, "Test %d: ", test);
print_timer (timers, 1, buf, TIMER_CLR);
timers[2] /= cnt;
strcat (buf, ", avg: ");
print_timer (timers, 2, buf, TIMER_NL | TIMER_CLR);
OutputBoth (thread_num, buf);
}
aligned_free (workbuf);
gwdone (&gwdata);
if (min_test) exit (0);
return 0;
}
#endif

#ifdef TIMING604
if (w->n == 604) {
	gmp_randstate_t rstate;
	int	exp;

	gmp_randinit_default (rstate);
	for (exp = 100000; exp < 300000000; exp *= 2) {
		gwhandle gwdata;
		giant	N;
		gwnum	a, b, c;
		int	j;
		mpz_t	r;

		gwinit (&gwdata);
		gwsetup (&gwdata, 1.0, 2, exp, -1);

		a = gwalloc (&gwdata);
		b = gwalloc (&gwdata);
		c = gwalloc (&gwdata);
		N = popg (&gwdata.gdata, ((int) gwdata.bit_length >> 5) + 10);
		mpz_init (r);
		mpz_urandomb (r, rstate, exp);
		mpztog (r, N);
		gianttogw (&gwdata, N, a);
		mpz_urandomb (r, rstate, exp);
		mpztog (r, N);
		gianttogw (&gwdata, N, b);
		mpz_urandomb (r, rstate, exp);
		mpztog (r, N);
		gianttogw (&gwdata, N, c);
		mpz_clear (r);

		start_timer (timers, 0);
		for (j = 0; j < 100; j++) {
			gwsquare2 (&gwdata, c, c, GWMUL_STARTNEXTFFT);
		}
		end_timer (timers, 0);
		sprintf (buf, "Exp: %d, 100 squares: ", exp);
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputBoth (thread_num, buf);

		start_timer (timers, 0);
		for (j = 0; j < 1; j++) {
			mpz_t	__v, __N, __gcd;
			giant	v;
			v = popg (&gwdata.gdata, ((int) gwdata.bit_length >> 5) + 10);
			gwtogiant (&gwdata, b, v);
			/* Do the GCD */
			mpz_init (__v);
			mpz_init (__N);
			mpz_init (__gcd);
			gtompz (v, __v);
			gtompz (N, __N);
			mpz_gcd (__gcd, __v, __N);
			mpz_clear (__gcd);
			mpz_clear (__v);
			mpz_clear (__N);
			pushg (&gwdata.gdata, 1);
		}
		end_timer (timers, 0);
		sprintf (buf, "Exp: %d, GCD: ", exp);
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputBoth (thread_num, buf);

		start_timer (timers, 0);
		for (j = 0; j < 1; j++) {
			giant	v;
			mpz_t	__v, __N, __gcd, __inv;
			v = popg (&gwdata.gdata, ((int) gwdata.bit_length >> 5) + 10);
			gwtogiant (&gwdata, b, v);
			/* Do the extended GCD */
			mpz_init (__v);
			mpz_init (__N);
			mpz_init (__gcd);
			mpz_init (__inv);
			gtompz (v, __v);
			gtompz (N, __N);
			mpz_gcdext (__gcd, __inv, NULL, __v, __N);
			mpz_clear (__v);
			if (mpz_sgn (__inv) < 0) mpz_add (__inv, __inv, __N);
			mpztog (__inv, v);
			gianttogw (&gwdata, v, b);
			mpz_clear (__gcd);
			mpz_clear (__inv);
			mpz_clear (__N);
			pushg (&gwdata.gdata, 1);
		}
		end_timer (timers, 0);
		sprintf (buf, "Exp: %d, modinv: ", exp);
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputBoth (thread_num, buf);
		gwdone (&gwdata);
	}
	return (0);
}
#endif


/* Init filename */

	tempFileName (w, filename);
	uniquifySaveFile (thread_num, filename);

/* Init the random number generator */

	srand ((unsigned) time (NULL));

/* Perform setup functions.  This includes decding how big an FFT to use, allocating memory, calling the FFT setup code, etc. */

/* Setup the gwnum assembly code */

	gwinit (&ecmdata.gwdata);
	gwset_sum_inputs_checking (&ecmdata.gwdata, SUM_INPUTS_ERRCHK);
	if (IniGetInt (LOCALINI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&ecmdata.gwdata);
	if (IniGetInt (INI_FILE, "HyperthreadPrefetch", 0)) gwset_hyperthread_prefetch (&ecmdata.gwdata);
	if (HYPERTHREAD_LL) {
		sp_info->normal_work_hyperthreads = IniGetInt (LOCALINI_FILE, "HyperthreadLLcount", CPU_HYPERTHREADS);
		gwset_will_hyperthread (&ecmdata.gwdata, sp_info->normal_work_hyperthreads);
	}
	gwset_bench_cores (&ecmdata.gwdata, NUM_CPUS);
	gwset_bench_workers (&ecmdata.gwdata, NUM_WORKER_THREADS);
	if (ERRCHK) gwset_will_error_check (&ecmdata.gwdata);
	gwset_num_threads (&ecmdata.gwdata, CORES_PER_TEST[thread_num] * sp_info->normal_work_hyperthreads);
	gwset_thread_callback (&ecmdata.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&ecmdata.gwdata, sp_info);
	gwset_safety_margin (&ecmdata.gwdata, IniGetFloat (INI_FILE, "ExtraSafetyMargin", 0.0));
	gwset_minimum_fftlen (&ecmdata.gwdata, w->minimum_fftlen);
	res = gwsetup (&ecmdata.gwdata, w->k, w->b, w->n, w->c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}

/* A kludge so that the error checking code is not as strict. */

	ecmdata.gwdata.MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);

/* More random initializations */

	gwsetnormroutine (&ecmdata.gwdata, 0, ERRCHK, 0);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);
	first_iter_msg = TRUE;
	calc_output_frequencies (&ecmdata.gwdata, &output_frequency, &output_title_frequency);

/* Optionally do a probable prime test */

	if (IniGetInt (INI_FILE, "ProbablePrimeTest", 0) && isProbablePrime (&ecmdata.gwdata, N)) {
		sprintf (buf, "%s is a probable prime\n", gwmodulo_as_string (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
	}

/* Output a startup message */

	gwfft_description (&ecmdata.gwdata, fft_desc);
	sprintf (buf, "\nUsing %s\n", fft_desc);
	OutputStr (thread_num, buf);
	sprintf (buf, "%5.3f bits-per-word below FFT limit (more than %5.3f allows extra optimizations)\n",
		 ecmdata.gwdata.fft_max_bits_per_word - virtual_bits_per_word (&ecmdata.gwdata), EB_FIRST_ADD);
	OutputStr (thread_num, buf);

/* Check for a continuation file.  Limit number of backup files we try to read in case there is an error deleting bad save files. */

	readSaveFileStateInit (&ecmdata.read_save_file_state, thread_num, filename);
	writeSaveFileStateInit (&ecmdata.write_save_file_state, filename, 0);
	for ( ; ; ) {
		if (! saveFileExists (&ecmdata.read_save_file_state)) {
			/* If there were save files, they are all bad.  Report a message and temporarily abandon the work unit.  We do this in hopes that we */
			/* can successfully read one of the bad save files at a later time.  This sounds crazy, but has happened when OSes get in a funky state. */
			if (ecmdata.read_save_file_state.a_non_bad_save_file_existed) {
				OutputBoth (thread_num, ALLSAVEBAD_MSG);
				return (0);
			}
			/* No save files existed, start from scratch. */
			break;
		}

/* Allocate memory to read gwnums from a save file.  We'll move the read in values to where they need to be depending on the ECM state. */

		if (!alloc_xz (&ecmdata, &ecmdata.xz)) goto oom;
		ecmdata.gg = NULL;

/* Read in the save file.  If the save file is no good ecm_restore will have deleted it.  Loop trying to read a backup save file. */

		if (! ecm_restore (&ecmdata)) {
			free_xz (&ecmdata, &ecmdata.xz);
			/* Close and rename the bad save file */
			saveFileBad (&ecmdata.read_save_file_state);
			continue;
		}

/* Compute Ad4 from sigma while ignoring the x,z starting point */

		curve_start_msg (&ecmdata);
		stop_reason = choose12 (&ecmdata, NULL, N, &factor);
		if (stop_reason) goto exit;

/* Continue in the middle of stage 1 */

		if (ecmdata.state == ECM_STATE_STAGE1) {
			sieve_start = ecmdata.stage1_prime + 1;
			goto restart1;
		}

/* Continue if between stage 1 and stage 2.  Stage 2 init expects Q stored in QD.  Swap xz to QD. */

		if (ecmdata.state == ECM_STATE_MIDSTAGE) {
			xzswap (ecmdata.xz, ecmdata.QD);
			goto restart3;
		}

/* We've finished stage 1, resume stage 2.  Store normalized Q in QD.  Initialize the accumulator. */

		if (ecmdata.state == ECM_STATE_STAGE2) {
			ecmdata.QD.x = ecmdata.xz.x;
			ecmdata.gg = ecmdata.xz.z;
			ecmdata.xz.x = NULL;
			ecmdata.xz.z = NULL;
			ecmdata.QD.z = gwalloc (&ecmdata.gwdata);
			if (ecmdata.QD.z == NULL) goto oom;
			dbltogw (&ecmdata.gwdata, 1.0, ecmdata.QD.z);
			goto restart3;
		}

/* We've finished stage 2, but haven't done the GCD yet.  Place accumulator in gg, leave normalized Q in xz.x. */

		ASSERTG (ecmdata.state == ECM_STATE_GCD);
		gwswap (ecmdata.xz.z, ecmdata.gg);
		goto restart4;
	}

/* Unless a save file indicates otherwise, we are testing our first curve */

	ecmdata.curve = 1;
	ecmdata.average_B2 = 0;

/* Loop processing the requested number of ECM curves */

restart0:
	ecmdata.state = ECM_STATE_STAGE1_INIT;
	ecmdata.pct_mem_to_use = 1.0;				// Use as much memory as we can unless we get allocation errors
	ecm_stage1_memory_usage (thread_num, &ecmdata);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Allocate memory */

	if (!alloc_xz (&ecmdata, &ecmdata.xz)) goto oom;
	ecmdata.gg = NULL;

/* Choose curve with order divisible by 16 and choose a point (x/z) on said curve. */

	do {
		uint32_t hi, lo;
		ecmdata.sigma = (rand () & 0x1F) * 65536.0 * 65536.0 * 65536.0;
		ecmdata.sigma += (rand () & 0xFFFF) * 65536.0 * 65536.0;
		if (CPU_FLAGS & CPU_RDTSC) rdtsc (&hi, &lo);
		ecmdata.sigma += lo ^ hi ^ ((unsigned long) rand () << 16);
	} while (ecmdata.sigma <= 5.0);
	if (w->curve > 5.0 && w->curve < 9007199254740992.0) {
		ecmdata.sigma = w->curve;
		w->curves_to_do = 1;
	}
	curve_start_msg (&ecmdata);
	stop_reason = choose12 (&ecmdata, &ecmdata.xz, N, &factor);
	if (stop_reason) goto exit;
	if (factor != NULL) goto bingo;
	sieve_start = 2;

/* The stage 1 restart point */

restart1:
	ecmdata.state = ECM_STATE_STAGE1;
	ecm_stage1_memory_usage (thread_num, &ecmdata);
	one_over_B = 1.0 / (double) ecmdata.B;
	sprintf (w->stage, "C%luS1", ecmdata.curve);
	w->pct_complete = sieve_start * one_over_B;
	start_timer (timers, 0);
	SQRT_B = (unsigned long) sqrt ((double) ecmdata.B);
	// We guess the max sieve prime for stage 2.  If optimal B2 is less 256 * B, then max sieve prime will be less than 16 * sqrt(B).
	// If our guess is wrong, that's no big deal -- sieve code is smart enough to handle it.
	stop_reason = start_sieve_with_limit (thread_num, sieve_start, 16 * SQRT_B, &ecmdata.sieve_info);
	if (stop_reason) goto exit;
	for (ecmdata.stage1_prime = sieve (ecmdata.sieve_info); ecmdata.stage1_prime <= ecmdata.B; ecmdata.stage1_prime = next_prime) {
		int	saving, count;

		next_prime = sieve (ecmdata.sieve_info);

/* Test for user interrupt, save files, and error checking */

		stop_reason = stopCheck (thread_num);
		saving = testSaveFilesFlag (thread_num);

/* Count the number of prime powers where prime^n <= B */

		if (ecmdata.stage1_prime <= SQRT_B) {
			uint64_t mult, max;
			count = 0;
			max = ecmdata.B / ecmdata.stage1_prime;
			for (mult = ecmdata.stage1_prime; ; mult *= ecmdata.stage1_prime) {
				count++;
				if (mult > max) break;
			}
		} else
			count = 1;

/* Adjust the count of prime powers.  I'm not sure, but I think choose12 means we should include 2 extra twos and 1 extra 3. */

		if (ecmdata.stage1_prime <= 3) count += (ecmdata.stage1_prime == 2) ? 2 : 1;

/* Apply the proper number of primes */

		for ( ; count; count--) {
			int	last_mul = (count == 1 && (saving || stop_reason || next_prime > ecmdata.B));
			int	oom_stop_reason = ell_mul (&ecmdata, &ecmdata.xz, ecmdata.stage1_prime, last_mul);
			if (oom_stop_reason) {	// Out of memory should be the only reason ell_mul fails.  Won't happen, but handle it gracefully anyway
				stop_reason = oom_stop_reason;
				goto exit;
			}
		}

/* Calculate stage 1 percent complete */

		w->pct_complete = ecmdata.stage1_prime * one_over_B;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&ecmdata.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s ECM curve %lu stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve);
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&ecmdata.gwdata);
		}

/* Print a message every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&ecmdata.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s curve %lu stage 1 at prime %" PRIu64 " [%.*f%%].",
				 gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve, ecmdata.stage1_prime, (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			if (ecmdata.stage1_prime != 2)
				OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&ecmdata.gwdata);
			first_iter_msg = FALSE;
		}

/* Check for errors */

		if (gw_test_for_error (&ecmdata.gwdata)) goto err;

/* Write a save file when the user interrupts the calculation and every DISK_WRITE_TIME minutes. */

		if (stop_reason || saving) {
			ecm_save (&ecmdata);
			if (stop_reason) goto exit;
		}
	}

/* Stage 1 complete */

	end_timer (timers, 0);
	sprintf (buf, "Stage 1 complete. %.0f transforms, %lu modular inverses. Time: ", gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&ecmdata.gwdata);
	}

/* If we aren't doing a stage 2, then check to see if we found a factor. */
/* If we are doing a stage 2, then the stage 2 init will do this GCD for us. */

	if (ecmdata.C <= ecmdata.B) {
skip_stage_2:	start_timer (timers, 0);
		stop_reason = gcd (&ecmdata.gwdata, thread_num, ecmdata.xz.z, N, &factor);
		if (stop_reason) {
			ecm_save (&ecmdata);
			goto exit;
		}
		end_timer (timers, 0);
		strcpy (buf, "Stage 1 GCD complete. Time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputStr (thread_num, buf);
		if (factor != NULL) goto bingo;

/* Alexander Kruppa wrote this code to normalize and output the x value along with N and sigma so that it can be used in Paul Zimmermann's */
/* superior GMP-ECM implementation of stage 2. */

		if (IniGetInt (INI_FILE, "GmpEcmHook", 0)) {
			char	*msg, *buf;
			int	msglen, i, leadingzeroes;
			giant	gx;
			char	*hex = "0123456789ABCDEF";

			stop_reason = normalize (&ecmdata, ecmdata.xz.x, ecmdata.xz.z, N, &factor);
			if (stop_reason) goto exit;

			gx = popg (&ecmdata.gwdata.gdata, ((int) ecmdata.gwdata.bit_length >> 5) + 10);
			if (gx == NULL) goto oom;
			if (gwtogiant (&ecmdata.gwdata, ecmdata.xz.x, gx)) goto oom;  // Unexpected error, return oom for lack of a better error message
			modgi (&ecmdata.gwdata.gdata, N, gx);

			msglen = N->sign * 8 + 5;
			buf = (char *) malloc (msglen + msglen + 80);
			if (buf == NULL) goto oom;
			strcpy (buf, "N=0x");
			msg = buf + strlen (buf);
			leadingzeroes = 1; /* Still eat leading zeroes? */
			for (i = 0; i < N->sign * 8; i++) {
				char nibble = ( N->n[N->sign - 1 - i/8] >> ((7-i%8)*4)) & 0xF;
				if (nibble != 0) leadingzeroes = 0;
				if (!leadingzeroes) *msg++ = hex[nibble];
			}
			strcpy (msg, "; QX=0x");
			msg = msg + strlen (msg);
			leadingzeroes = 1;
			for (i = 0; i < gx->sign * 8; i++) {
				char nibble = ( gx->n[gx->sign - 1 - i/8] >> ((7-i%8)*4)) & 0xF;
				if (nibble != 0) leadingzeroes = 0;
				if (!leadingzeroes) *msg++ = hex[nibble];
			}
			strcpy (msg, "; SIGMA=");
			msg = msg + strlen (msg);
			sprintf (msg, "%.0f\n", ecmdata.sigma);
			writeResults (buf);
			free (buf);
			pushg (&ecmdata.gwdata.gdata, 1);
		}

/* Now do the next ECM curve */

		goto more_curves;
	}

/*
   Stage 2:  We support two types of stage 2's here.  One uses less memory and uses fewer extended GCDs, but is slower in accumulating
   each found prime.  Thanks to Richard Crandall and Paul Zimmermann for letting me liberally use their code and ideas here.
   x, z: coordinates of Q at the beginning of stage 2
*/

/* Change state to between stage 1 and 2 */

	ecmdata.state = ECM_STATE_MIDSTAGE;
	sprintf (w->stage, "C%luS2", ecmdata.curve);
	w->pct_complete = 0.0;

/* Thoughout stage 2 init, ecmdata.QD will contain a multiple of Q suitable for saving.  We start QD with Q. */

	xzswap (ecmdata.xz, ecmdata.QD);

/* Entry point for continuing stage 2 from a save file */

restart3:
	start_timer (timers, 0);
	sprintf (buf, "%s ECM curve %lu stage 2 init", gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve);
	title (thread_num, buf);

/* Make sure we will have enough memory to run stage 2 at some time.  We need at least 13 gwnums in stage 2 main loop. */
/* 6 for mQ_next computations (Qm, Qprevm, QD), 4 for nQx, gg, 2 for ell_add_xz_noscr temps. */

replan:	min_memory = cvt_gwnums_to_mem (&ecmdata.gwdata, 13);
	if ((int) max_mem (thread_num) < min_memory) {
		sprintf (buf, "Skipping stage 2 due to insufficient memory -- %dMB needed.\n", min_memory);
		OutputStr (thread_num, buf);
		goto skip_stage_2;
	}

/* Choose the best plan implementation given the currently available memory.  This implementation could be "wait until we have more memory". */

	stop_reason = ecm_stage2_impl (&ecmdata);
	if (stop_reason) goto possible_lowmem;

/* Record the amount of memory this thread will be using in stage 2. */

	memused = cvt_gwnums_to_mem (&ecmdata.gwdata, ecmdata.stage2_numvals);
	memused += (_intmin (ecmdata.numDsections - ecmdata.bitarrayfirstDsection, ecmdata.bitarraymaxDsections) * ecmdata.totrels) >> 23;
	// To dodge possible infinite loop if ecm_stage2_impl allocates too much memory (it shouldn't), decrease the percentage of memory we are allowed to use
	// Beware that replaning may allocate a larger bitarray
	if (set_memory_usage (thread_num, MEM_VARIABLE_USAGE, memused)) {
		ecmdata.pct_mem_to_use *= 0.99;
		free (ecmdata.bitarray); ecmdata.bitarray = NULL;
		goto replan;
	}
	gw_set_max_allocs (&ecmdata.gwdata, ecmdata.stage2_numvals);

/* Output a useful message regarding memory usage and ECM stage 2 implementation */

	if (ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 2 FFTs per prime pair, 3-mult modinv pooling, pool size %d.\n", memused, ecmdata.E);
	else if (ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3POINT44MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 2 FFTs per prime pair, 3.44-mult modinv pooling, pool size %d.\n", memused, ecmdata.E);
	else if (ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3POINT57MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 2 FFTs per prime pair, 3.57-mult modinv pooling, pool size %d.\n", memused, ecmdata.E);
	else if (ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_N_SQUARED)
		sprintf (buf, "Stage 2 uses %uMB of memory, 2 FFTs per prime pair, N^2 modinv pooling, pool size %d.\n", memused, ecmdata.E);
	else if (!ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 4 FFTs per prime pair, 3-mult modinv pooling.\n", memused);
	else if (!ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3POINT44MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 4 FFTs per prime pair, 3.44-mult modinv pooling.\n", memused);
	else if (!ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_3POINT57MULT)
		sprintf (buf, "Stage 2 uses %uMB of memory, 4 FFTs per prime pair, 3.57-mult modinv pooling.\n", memused);
	else if (!ecmdata.TWO_FFT_STAGE2 && ecmdata.pool_type == POOL_N_SQUARED)
		sprintf (buf, "Stage 2 uses %uMB of memory, 4 FFTs per prime pair, N^2 modinv pooling.\n", memused);
	OutputStr (thread_num, buf);

/* Initialize variables for second stage.  Ideally we fill up the nQx array with Q^relprime normalized with only one modular inverse. */

	// Calculate the percent completed by previous bit arrays
	base_pct_complete = (double) (ecmdata.B2_start - ecmdata.last_relocatable) / (double) (ecmdata.C - ecmdata.last_relocatable);
	// Calculate the percent completed by each bit in this bit array
	one_bit_pct = (1.0 - base_pct_complete) / (double) (ecmdata.numDsections * ecmdata.totrels);
	// Calculate the percent completed by previous bit arrays and the current bit array
	w->pct_complete = base_pct_complete + (double) (ecmdata.Dsection * ecmdata.totrels) * one_bit_pct;

/* Allocate nQx array of pointers to relative prime gwnums */

	ecmdata.nQx = (gwnum *) malloc (ecmdata.totrels * sizeof (gwnum));
	if (ecmdata.nQx == NULL) goto lowmem;

/* Allocate some of the memory for modular inverse pooling */

	stop_reason = normalize_pool_init (&ecmdata, ecmdata.stage2_numvals);
	if (stop_reason) goto possible_lowmem;

/* We have two approaches for computing nQx, if memory is really tight we compute every odd multiple of x. */
/* Otherwise, if D is divisible by 3 we compute the 1,5 mod 6 multiples which reduces initialization cost by 33%. */

	if (!ecmdata.TWO_FFT_STAGE2 || ecmdata.D % 3 != 0) {
		struct xz t1, t2, t3, t4;
		struct xz *Q1, *Q2, *Q3, *Qi, *Qiminus2, *Qscratch;
		int	have_QD, totrels;

/* Allocate memory and init values for computing nQx.  We need Q^1, Q^2, Q^3. */

		if (!alloc_xz (&ecmdata, &t1)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t2)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t3)) goto lowmem;
		t4.x = NULL; t4.z = NULL;					// Only needed if more Q2 is needed after QD is calculated

/* Compute Q^2 and place in ecmdata.QD.  This is the value we will save if init is aborted and we create a save file. */
/* Q^2 will be replaced by Q^D later on in this loop.  At all times ecmdata.QD will be a save-able value! */
/* Remember Q^1 as the first nQx value.  After the first stage 2 init, xz.x is already normalized. */

		xzswap (t1, ecmdata.QD);					// Move Q^1 from QD to t1
		Q1 = &t1;							// Q^1
		Q2 = &ecmdata.QD;
		ell_dbl_xz_scr (&ecmdata, Q1, Q2, &t3);				// Q2 = 2 * Q1, scratch = t3
		ecmdata.nQx[0] = Q1->x; totrels = 1;

/* Compute Q^3 */

		Q3 = &t2;
		ell_add_xz_scr (&ecmdata, Q2, Q1, Q1, Q3, &t3);			// Q3 = Q2 + Q1 (diff Q1), scratch = t3

/* Compute the rest of the nQx values (Q^i for i >= 3) */
/* MEMPEAK: 9 + nQx + 1 (AD4, QD.x&z, 6 for t1 through t3, nQx vals, and 1 for modinv_value assuming N^2 pooling) */
/* BUT! If using the minimum number of totrels (e.g. 4, D=30), two of the nQx values are sharing space with Qiminus2.  Reducing MEMPEAK to --- */
/* MEMPEAK: 7 + nQx + 1 (AD4, 6 for t1 through t3, nQx vals, and 1 for modinv_value assuming N^2 pooling) */

		Qi = Q3;
		Qiminus2 = Q1;
		Qscratch = &t3;
		have_QD = FALSE;
		for (i = 3; ; i += 2) {

/* Compute Q^D which we will need in mQ_init.  We need two different ways to do this based on whether D is 0 or 2 mod 4. */

			if (i + i == ecmdata.D) {
				ASSERTG (Q2 == &ecmdata.QD);
				// In tightest memory case we can overwrite Q2 as it will no longer be needed.  Otherwise, move Q2 to fourth temporary.
				if (totrels != ecmdata.totrels) {
					t4 = ecmdata.QD; Q2 = &t4;				// Move Q2 out of ecmdata.QD
					if (!alloc_xz (&ecmdata, &ecmdata.QD)) goto lowmem;	// Allocate QD
				}
				ell_dbl_xz_scr (&ecmdata, Qi, &ecmdata.QD, Qscratch);		// QD = 2 * Qi
				have_QD = TRUE;
			}
			if (i + i-2 == ecmdata.D) {
				// In tightest memory case we can overwrite Q2 as it will no longer be needed.  Otherwise, move Q2 to fourth temporary.
				if (totrels != ecmdata.totrels) {
					t4 = ecmdata.QD; Q2 = &t4;				// Move Q2 out of ecmdata.QD
					if (!alloc_xz (&ecmdata, &ecmdata.QD)) goto lowmem;	// Allocate QD
				}
				ell_add_xz (&ecmdata, Qi, Qiminus2, Q2, &ecmdata.QD);		// QD = i + (i-2), diff Q2
				have_QD = TRUE;
			}

/* Break out of loop when we have all our nQx values */

			if (have_QD && totrels == ecmdata.totrels) break;

/* Get next odd value, but don't destroy Qiminus2 in case it is an nQx value. */

			ell_add_xz (&ecmdata, Qi, Q2, Qiminus2, Qscratch);			// Next odd value,  save in scratch
			if (relatively_prime (i-2, ecmdata.D)) {
				stop_reason = add_to_normalize_pool_ro (&ecmdata, Qiminus2->x, Qiminus2->z, &factor);
				if (stop_reason) goto possible_lowmem;
				if (factor != NULL) goto bingo;
				if (!alloc_xz (&ecmdata, Qiminus2)) goto lowmem;
			}
			xzswap (*Qiminus2, *Qi);
			xzswap (*Qi, *Qscratch);
			if (relatively_prime (i+2, ecmdata.D) && totrels < ecmdata.totrels) ecmdata.nQx[totrels++] = Qi->x;

/* Check for errors and save file creation */

			if (gw_test_for_error (&ecmdata.gwdata)) goto err;

			stop_reason = stopCheck (thread_num);
			if (stop_reason) goto possible_lowmem;
		}

/* Free Q2 if it was moved to t4 */

		free_xz (&ecmdata, &t4);

/* If needed, add any of the last two computed values to the normalization pool */

		for (i = 1; i <= 2; i++) {
			struct xz *needs_work;
			if (t1.x == ecmdata.nQx[totrels-i]) needs_work = &t1;
			else if (t2.x == ecmdata.nQx[totrels-i]) needs_work = &t2;
			else if (t3.x == ecmdata.nQx[totrels-i]) needs_work = &t3;
			else continue;
			stop_reason = add_to_normalize_pool_ro (&ecmdata, needs_work->x, needs_work->z, &factor);
			if (stop_reason) goto possible_lowmem;
			if (factor != NULL) goto bingo;
			needs_work->x = NULL;
			needs_work->z = NULL;
		}

/* Free memory used in computing nQx values */

		free_xz (&ecmdata, &t1);
		free_xz (&ecmdata, &t2);
		free_xz (&ecmdata, &t3);
	}

/* This is the faster 1,5 mod 6 nQx initialization */

	else {
		struct xz t1, t2, t3, t4, t5, t6;
		struct xz *Q1, *Q2, *Q3, *Q5, *Q6, *Q7, *Q11, *Q1mod6, *Q1mod6minus6, *Q5mod6, *Q5mod6minus6, *Qscratch;
		int	have_QD, totrels;

/* Allocate memory and init values for computing nQx.  We need Q^1, Q^5, Q^6, Q^7, Q^11. */
/* We also need Q^4 to compute Q^D in the middle of the nQx init loop. */

		if (!alloc_xz (&ecmdata, &t1)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t2)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t3)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t4)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t5)) goto lowmem;
		if (!alloc_xz (&ecmdata, &t6)) goto lowmem;

/* Compute Q^2 and place in ecmdata.QD.  This is the value we will save if init is aborted and we create a save file. */
/* Q^2 will be replaced by Q^D later on in this loop.  At all times ecmdata.QD will be a save-able value! */
/* Remember Q^1 as the first nQx value.  After the first stage 2 init, xz.x is already normalized. */

		xzswap (t1, ecmdata.QD);					// Move Q^1 from QD to t1
		Q1 = &t1;							// Q^1
		Q2 = &ecmdata.QD;
		ell_dbl_xz_scr (&ecmdata, Q1, Q2, &t6);				// Q2 = 2 * Q1, scratch = t6
		ecmdata.nQx[0] = Q1->x; totrels = 1;

/* Compute Q^5, Q^6, Q^7, Q^11 */

		Q3 = &t2;
		ell_add_xz_scr (&ecmdata, Q2, Q1, Q1, Q3, &t6);			// Q3 = Q2 + Q1 (diff Q1), scratch = t6

		Q5 = &t3;
		ell_add_xz_scr (&ecmdata, Q3, Q2, Q1, Q5, &t6);			// Q5 = Q3 + Q2 (diff Q1), scratch = t6
		if (relatively_prime (5, ecmdata.D)) ecmdata.nQx[totrels++] = Q5->x;

		Q6 = Q3;
		ell_dbl_xz_scr (&ecmdata, Q3, Q6, &t6);				// Q6 = 2 * Q3, scratch = t6, Q3 no longer needed

		Q7 = &t4;
		ell_add_xz_scr (&ecmdata, Q6, Q1, Q5, Q7, &t6);			// Q7 = Q6 + Q1 (diff Q5), scratch = t6
		if (relatively_prime (7, ecmdata.D)) ecmdata.nQx[totrels++] = Q7->x;

		Q11 = &t5;
		ell_add_xz_scr (&ecmdata, Q6, Q5, Q1, Q11, &t6);		// Q11 = Q6 + Q5 (diff Q1), scratch = t6
		if (relatively_prime (11, ecmdata.D)) ecmdata.nQx[totrels++] = Q11->x;

/* Compute the rest of the nQx values (Q^i for i >= 7) */
/* MEMPEAK: 15 + nQx * 2 + (AD4, QD.x&z, 12 for t1 through t6, nQx vals, another nQx for POOL_3MULT poolz values) */
/* But at least 1 nQx and 1 poolz value are still in t1 through t6, so that reduces mempeak by at least two. */

//GW: With 3N pooling, we are not handling case where totrels > 2 * E -- we exceed our memory allocation  (need 2 modinvs)
// with 3.57N pooling this gets a lot better....

		Q1mod6 = Q7;
		Q1mod6minus6 = Q1;
		Q5mod6 = Q11;
		Q5mod6minus6 = Q5;
		Qscratch = &t6;
		have_QD = FALSE;
		for (i = 7; ; i += 6) {

/* Compute Q^D which we will need in mQ_init.  Do this with a single ell_add_xz call when we reach two values that are 4 apart that add to D. */
/* This only works for D = 2 mod 4. */

			if (i + i+4 == ecmdata.D) {
				ASSERTG (Q2 == &ecmdata.QD);
				ell_dbl_xz_scr (&ecmdata, Q2, Q2, Qscratch);				// Q4 = 2 * Q2, Q2 no longer needed
				ell_add_xz_scr (&ecmdata, Q1mod6, Q5mod6, Q2, &ecmdata.QD, Qscratch);	// QD = i + i+4 (diff Q4), Q4 no longer needed
				have_QD = TRUE;
			}

/* Break out of loop when we have all our nQx values */

			if (have_QD && totrels == ecmdata.totrels) break;

/* Get next 1 mod 6 value, but don't destroy Q1mod6minus6 in case it is an nQx value. */

			ell_add_xz (&ecmdata, Q1mod6, Q6, Q1mod6minus6, Qscratch);			// Next 1 mod 6 value, save in scratch
			if (relatively_prime (i-6, ecmdata.D)) {
				stop_reason = add_to_normalize_pool_ro (&ecmdata, Q1mod6minus6->x, Q1mod6minus6->z, &factor);
				if (stop_reason) goto possible_lowmem;
				if (factor != NULL) goto bingo;
				if (!alloc_xz (&ecmdata, Q1mod6minus6)) goto lowmem;
			}
			xzswap (*Q1mod6minus6, *Q1mod6);
			xzswap (*Q1mod6, *Qscratch);
			if (relatively_prime (i+6, ecmdata.D) && totrels < ecmdata.totrels) ecmdata.nQx[totrels++] = Q1mod6->x;

/* Compute Q^D which we will need in mQ_init.  Do this with a single ell_add_xz call when we reach the two values that are 2 apart that add to D. */
/* This only works for D = 0 mod 4. */

			if (i+6 + i+4 == ecmdata.D) {
				ASSERTG (Q2 == &ecmdata.QD);
				ell_add_xz_scr (&ecmdata, Q1mod6, Q5mod6, Q2, &ecmdata.QD, Qscratch);	// QD = i+6 + i+4 (diff Q2), Q2 no longer needed
				have_QD = TRUE;
			}

/* Break out of loop when we have QD and all our nQx values */

			if (have_QD && totrels == ecmdata.totrels) break;

/* Get next 5 mod 6 value, but don't destroy Q5mod6minus6 in case it is an nQx value. */

			ell_add_xz (&ecmdata, Q5mod6, Q6, Q5mod6minus6, Qscratch);			// Next 5 mod 6 value, save in scratch
			if (relatively_prime (i+4-6, ecmdata.D)) {
				stop_reason = add_to_normalize_pool_ro (&ecmdata, Q5mod6minus6->x, Q5mod6minus6->z, &factor);
				if (stop_reason) goto possible_lowmem;
				if (factor != NULL) goto bingo;
				if (!alloc_xz (&ecmdata, Q5mod6minus6)) goto lowmem;
			}
			xzswap (*Q5mod6minus6, *Q5mod6);
			xzswap (*Q5mod6, *Qscratch);
			if (relatively_prime (i+4+6, ecmdata.D) && totrels < ecmdata.totrels) ecmdata.nQx[totrels++] = Q5mod6->x;

/* Check for errors and save file creation */

			if (gw_test_for_error (&ecmdata.gwdata)) goto err;

			stop_reason = stopCheck (thread_num);
			if (stop_reason) goto possible_lowmem;
		}

/* If needed, add any of the last four computed values to the normalization pool */

		for (i = 1; i <= 4; i++) {
			struct xz *needs_work;
			if (t1.x == ecmdata.nQx[totrels-i]) needs_work = &t1;
			else if (t2.x == ecmdata.nQx[totrels-i]) needs_work = &t2;
			else if (t3.x == ecmdata.nQx[totrels-i]) needs_work = &t3;
			else if (t4.x == ecmdata.nQx[totrels-i]) needs_work = &t4;
			else if (t5.x == ecmdata.nQx[totrels-i]) needs_work = &t5;
			else if (t6.x == ecmdata.nQx[totrels-i]) needs_work = &t6;
			else continue;
			stop_reason = add_to_normalize_pool_ro (&ecmdata, needs_work->x, needs_work->z, &factor);
			if (stop_reason) goto possible_lowmem;
			if (factor != NULL) goto bingo;
			needs_work->x = NULL;
			needs_work->z = NULL;
		}

/* Free memory used in computing nQx values */

		free_xz (&ecmdata, &t1);
		free_xz (&ecmdata, &t2);
		free_xz (&ecmdata, &t3);
		free_xz (&ecmdata, &t4);
		free_xz (&ecmdata, &t5);
		free_xz (&ecmdata, &t6);
	}

/* Normalize all the nQx values */
/* MEMPEAK: 3 + nQx gwnums (AD4, QD.x&z, nQx values) + 1 or nQx for pooled_modinv */

	stop_reason = normalize_pool (&ecmdata, N, &factor);
	if (stop_reason) goto possible_lowmem;

/* If we found a factor, we're done */

	if (factor != NULL) goto bingo;

/* Precompute the transforms of nQx */

//	for (i = 0; i < ecmdata.totrels; i++)
//		gwfft (&ecmdata.gwdata, ecmdata.nQx[i], ecmdata.nQx[i]);

/* Init code that computes Q^m, where m is the first D section we are working on */
/* MEMPEAK: 9 + nQx (Ad4, 6 for computing mQx, 2 for mQ_init temporaries, nQx) + another one for gg if resuming a stage 2 */

	stop_reason = mQ_init (&ecmdata, ecmdata.B2_start + ecmdata.Dsection * ecmdata.D + ecmdata.D / 2);
	if (stop_reason) goto possible_lowmem;

/* Now init the accumulator unless this value was read from a continuation file */
/* MEMUSED: 7 + nQx gwnums (6 for computing mQx, gg, nQx values) */

	if (ecmdata.gg == NULL) {
		ecmdata.gg = gwalloc (&ecmdata.gwdata);
		if (ecmdata.gg == NULL) goto lowmem;
		dbltogw (&ecmdata.gwdata, 1.0, ecmdata.gg);
	}

/* Initialization of stage 2 complete */

	sprintf (buf, "%.*f%% of %s ECM curve %lu stage 2 (using %uMB)",
		 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve, memused);
	title (thread_num, buf);

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 init complete. %.0f transforms, %lu modular inverses. Time: ", gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	gw_clear_fft_count (&ecmdata.gwdata);
	ecmdata.modinv_count = 0;

/* Now do stage 2 */

/* We do prime pairing with each loop iteration handling the range m-Q to m+Q where m is a multiple of D and Q is the */
/* Q-th relative prime to D.  Totrels is often much larger than the number of relative primes less than D.  This Preda */
/* optimization (originally from Montgomery) provides us with many more chances to find a prime pairing. */

/* Accumulate (mQx - nQx)(mQz + nQz) - mQx mQz + nQx nQz.		*/
/* Since nQz = 1, we have (the 4 FFT per prime continuation)		*/
/*		== (mQx - nQx)(mQz + 1) - mQx mQz + nQx			*/
/*		== mQx mQz - nQx mQz + mQx - nQx - mQx mQz + nQx	*/
/*		== mQx - nQx mQz					*/
/* If mQz also = 1 (the 2 FFT per prime continuation) then we accumulate*/
/*		== mQx - nQx						*/

//GW:  Should we do our first mQ_next before switching state (to allow going to lowmem where we reduce memory settings and retry)?
//GW:  or can we handle out-of-memory from mQ_next???
	ecmdata.state = ECM_STATE_STAGE2;
	start_timer (timers, 0);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	for ( ; ; ) {
		uint64_t bitnum;
		gwnum	mQx, mQz, t1;

/* Compute this Q^m value */
/* MEMUSED pooling: 7 + nQx + E gwnums (6 for computing mQx, gg, nQx, E normalized D values) */
/* MEMPEAK pooling: 7 + nQx + E + (pooling cost of 1 up to E) + 2 for ell_add temporaries */
/* MEMUSED non-pooling: 7 + nQx gwnums (6 for computing mQx, gg, nQx) */
/* MEMPEAK non-pooling: 9 + nQx (2 for ell_add temporaries) */

		stop_reason = mQ_next (&ecmdata, &mQx, &mQz, N, &factor);
		if (stop_reason) {
			// In case stop_reason is out-of-memory, free some up before calling ecm_save.
			mQ_term (&ecmdata);
			ecm_save (&ecmdata);
			goto exit;
		}
		if (factor != NULL) goto bingo;

/* 4 FFT implementation requires another temporary.  Allocate it here (after mQ_next has freed its 2 ell_add_xz_noscr temporaries) */

		if (!ecmdata.TWO_FFT_STAGE2) {
			t1 = gwalloc (&ecmdata.gwdata);
			if (t1 == NULL) goto oom;
		}

/* Test which relprimes in this D section need to be processed */

		bitnum = (ecmdata.Dsection - ecmdata.bitarrayfirstDsection) * ecmdata.totrels;
		for (i = 0; i < ecmdata.totrels; i++) {

/* Skip this relprime if the pairing routine did not set the corresponding bit in the bitarray. */

			if (! bittst (ecmdata.bitarray, bitnum + i)) continue;

/* 2 FFT per prime continuation - deals with all normalized values */

//GW:dont start next fft if saving,stopping,processing last bit (like p-1)
			if (ecmdata.TWO_FFT_STAGE2) {
				gwsubmul4 (&ecmdata.gwdata, mQx, ecmdata.nQx[i], ecmdata.gg, ecmdata.gg, GWMUL_STARTNEXTFFT);
			}

/* 4 FFT per prime continuation - deals with only nQx values normalized */

			else {
				gwmul3 (&ecmdata.gwdata, ecmdata.nQx[i], mQz, t1, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
				gwsubmul4 (&ecmdata.gwdata, mQx, t1, ecmdata.gg, ecmdata.gg, GWMUL_STARTNEXTFFT);
			}

/* Clear pairing bit from the bitarray in case a save file is written.  Calculate stage 2 percentage. */

			bitclr (ecmdata.bitarray, bitnum + i);
			w->pct_complete = base_pct_complete + (double) (ecmdata.Dsection * ecmdata.totrels + i) * one_bit_pct;
			
/* Check for errors */

			if (gw_test_for_error (&ecmdata.gwdata)) goto err;

/* Output the title every so often */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&ecmdata.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
				sprintf (buf, "%.*f%% of %s ECM curve %lu stage 2 (using %uMB)",
					 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve, memused);
				title (thread_num, buf);
				last_output_t = gw_get_fft_count (&ecmdata.gwdata);
			}

/* Print a message every so often */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&ecmdata.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
				sprintf (buf, "%s curve %lu stage 2 in D-block=%" PRIu64 " [%.*f%%].",
					 gwmodulo_as_string (&ecmdata.gwdata), ecmdata.curve,
					 ecmdata.B2_start + ecmdata.Dsection * ecmdata.D + ecmdata.D / 2,
					 (int) PRECISION, trunc_percent (w->pct_complete));
				end_timer (timers, 0);
				if (first_iter_msg) {
					strcat (buf, "\n");
					clear_timer (timers, 0);
				} else {
					strcat (buf, " Time: ");
					print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
				}
				OutputStr (thread_num, buf);
				start_timer (timers, 0);
				last_output = gw_get_fft_count (&ecmdata.gwdata);
				first_iter_msg = FALSE;
			}

/* Write a save file when the user interrupts the calculation and every DISK_WRITE_TIME minutes. */

			stop_reason = stopCheck (thread_num);
			if (stop_reason || testSaveFilesFlag (thread_num)) {
				ecm_save (&ecmdata);
				if (stop_reason) goto exit;
			}
		}

/* Free 4 FFT temporary so that mQ_next can use it */

		if (!ecmdata.TWO_FFT_STAGE2) gwfree (&ecmdata.gwdata, t1);

/* Move onto the next D section when we are done with all the relprimes */

		ecmdata.Dsection++;
		ecmdata.C_done = ecmdata.B2_start + ecmdata.Dsection * ecmdata.D;
		if (ecmdata.Dsection >= ecmdata.numDsections) break;

/* See if more bitarrays are required to get us to bound #2 */

		if (ecmdata.Dsection < ecmdata.bitarrayfirstDsection + ecmdata.bitarraymaxDsections) continue;

		ecmdata.first_relocatable = calc_new_first_relocatable (ecmdata.D, ecmdata.C_done);
		ecmdata.bitarrayfirstDsection = ecmdata.Dsection;
		stop_reason = fill_work_bitarray (ecmdata.thread_num, &ecmdata.sieve_info, ecmdata.D, ecmdata.totrels,
						  ecmdata.first_relocatable, ecmdata.last_relocatable, ecmdata.C_done, ecmdata.C,
						  ecmdata.bitarraymaxDsections, &ecmdata.bitarray);
		if (stop_reason) {
			ecm_save (&ecmdata);
			goto exit;
		}
	}

/* Move nQx[0] to xz.x in case GCD code calls ecm_save.  Free lots of other stage 2 data. */

	mQ_term (&ecmdata);
	normalize_pool_term (&ecmdata);
	ecmdata.xz.x = ecmdata.nQx[0];
	for (i = 1; i < ecmdata.totrels; i++) gwfree (&ecmdata.gwdata, ecmdata.nQx[i]);
	free (ecmdata.nQx); ecmdata.nQx = NULL;
	free (ecmdata.bitarray); ecmdata.bitarray = NULL;
	mallocFreeForOS ();
	ecm_stage1_memory_usage (thread_num, &ecmdata);		// With the default 10 freed gwnums cached, this should be close to the correct mem usage

/* Compute the new Kruppa-adjusted B2 work completed.  This is needed when curves are run with different optimal B2 values due to */
/* changing available memory.  The Primenet server expects just one B2 value representing the work done. */

	{
		double total_B2 = 0.0;
		if (ecmdata.average_B2 > 0) total_B2 = (ecmdata.curve - 1) * kruppa_adjust (ecmdata.average_B2, ecmdata.B);
		total_B2 += kruppa_adjust (ecmdata.C, ecmdata.B);
		ecmdata.average_B2 = kruppa_unadjust (total_B2 / ecmdata.curve, ecmdata.B);
	}

/* Stage 2 is complete */

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 complete. %.0f transforms, %lu modular inverses. Time: ", gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&ecmdata.gwdata);
	}

/* See if we got lucky! */

restart4:
	ecmdata.state = ECM_STATE_GCD;
	sprintf (w->stage, "C%luS2", ecmdata.curve);
	w->pct_complete = 1.0;
	start_timer (timers, 0);
	stop_reason = gcd (&ecmdata.gwdata, thread_num, ecmdata.gg, N, &factor);
	if (stop_reason) {
		ecm_save (&ecmdata);
		goto exit;
	}
	end_timer (timers, 0);
	strcpy (buf, "Stage 2 GCD complete. Time: ");
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	if (factor != NULL) goto bingo;

/* Check if more curves need to be done */

more_curves:
	gwfreeall (&ecmdata.gwdata);
	ecmdata.xz.x = NULL;
	ecmdata.xz.z = NULL;
	ecmdata.gg = NULL;
	mallocFreeForOS ();
	if (++ecmdata.curve <= w->curves_to_do) goto restart0;

/* Output line to results file indicating the number of curves run */

	sprintf (buf, "%s completed %u ECM %s, B1=%" PRIu64 ",%s B2=%" PRIu64 ", Wi%d: %08lX\n",
		 gwmodulo_as_string (&ecmdata.gwdata), w->curves_to_do, w->curves_to_do == 1 ? "curve" : "curves",
		 ecmdata.B, ecmdata.optimal_B2 ? " average" : "", ecmdata.average_B2, PORT, SEC5 (w->n, ecmdata.B, ecmdata.average_B2));
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"NF", "exponent":45581713, "worktype":"ECM", "b1":50000, "b2":5000000, */
/* "curves":5, "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"NF\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"ECM\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64 ", \"b2\":%" PRIu64, ecmdata.B, ecmdata.average_B2);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"curves\":%u", w->curves_to_do);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", ecmdata.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, ecmdata.B, ecmdata.average_B2));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send ECM completed message to the server.  Although don't do it for puny B1 values. */

	if (!QA_IN_PROGRESS && (ecmdata.B >= 50000 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		truncated_strcpy (pkt.message, sizeof (pkt.message), buf);
		pkt.result_type = PRIMENET_AR_ECM_NOFACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.B1 = (double) ecmdata.B;
		pkt.B2 = (double) ecmdata.average_B2;
		pkt.curves = w->curves_to_do;
		pkt.fftlen = gwfftlen (&ecmdata.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Delete the save file */

	unlinkSaveFiles (&ecmdata.write_save_file_state);

/* Free memory and return */

	stop_reason = STOP_WORK_UNIT_COMPLETE;
exit:	ecm_cleanup (&ecmdata);
	free (N);
	free (factor);
	free (str);
	free (msg);
	return (stop_reason);

/* Low or possibly low on memory in stage 2 init, create save file, reduce memory settings, and try again */

lowmem:	stop_reason = OutOfMemory (thread_num);
possible_lowmem:
	if (ecmdata.state == ECM_STATE_MIDSTAGE) ecm_save (&ecmdata);
	if (stop_reason != STOP_OUT_OF_MEM) goto exit;
	ecm_cleanup (&ecmdata);
	OutputBoth (thread_num, "Memory allocation error.  Trying again using less memory.\n");
	ecmdata.pct_mem_to_use *= 0.8;
	goto restart;

/* We've run out of memory.  Print error message and exit. */

oom:	stop_reason = OutOfMemory (thread_num);
	goto exit;

/* Print a message, we found a factor! */

bingo:	stage = (ecmdata.state > ECM_STATE_MIDSTAGE) ? 2 : (ecmdata.state > ECM_STATE_STAGE1_INIT) ? 1 : 0;
	sprintf (buf, "ECM found a factor in curve #%lu, stage #%d\n", ecmdata.curve, stage);
	writeResults (buf);
	sprintf (buf, "Sigma=%.0f, B1=%" PRIu64 ", B2=%" PRIu64 ".\n", ecmdata.sigma, ecmdata.B, ecmdata.C);
	writeResults (buf);

/* Allocate memory for the string representation of the factor and for */
/* a message.  Convert the factor to a string. */ 

	msglen = factor->sign * 10 + 400;
	str = (char *) malloc (msglen);
	if (str == NULL) goto oom;
	msg = (char *) malloc (msglen);
	if (msg == NULL) goto oom;
	gtoc (factor, str, msglen);

/* Validate the factor we just found */

	if (!testFactor (&ecmdata.gwdata, w, factor)) {
		sprintf (msg, "ERROR: Bad factor for %s found: %s\n", gwmodulo_as_string (&ecmdata.gwdata), str);
		OutputBoth (thread_num, msg);
		OutputStr (thread_num, "Restarting ECM curve from scratch.\n");
		continueECM = TRUE;
		ecmdata.curve--;
		goto bad_factor_recovery;
	}

/* Output the validated factor */

	sprintf (msg, "%s has a factor: %s (ECM curve %lu, B1=%" PRIu64 ", B2=%" PRIu64 ")\n",
		 gwmodulo_as_string (&ecmdata.gwdata), str, ecmdata.curve, ecmdata.B, ecmdata.C);
	OutputStr (thread_num, msg);
	formatMsgForResultsFile (msg, w);
	writeResults (msg);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"F", "exponent":45581713, "worktype":"ECM", "factors":["430639100587696027847"], */
/* "b1":50000, "b2":5000000, "sigma":"123456789123456", "stage":2 */
/* "curves":5, "average-b2":5127843, "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"F\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"ECM\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"factors\":[\"%s\"]", str);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64 ", \"b2\":%" PRIu64, ecmdata.B, ecmdata.C);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"sigma\":%.0f", ecmdata.sigma);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"stage\":%d", stage);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"curves\":%lu", ecmdata.curve);
	if (ecmdata.optimal_B2 && ecmdata.average_B2 != ecmdata.C) sprintf (JSONbuf+strlen(JSONbuf), ", \"average-b2\":%" PRIu64, ecmdata.average_B2);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", ecmdata.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, ecmdata.B, ecmdata.C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* See if the cofactor is prime and set flag if we will be continuing ECM */

	continueECM = IniGetInt (INI_FILE, "ContinueECM", 0);
	prpAfterEcmFactor = IniGetInt (INI_FILE, "PRPAfterECMFactor", bitlen (N) < 100000);
	if (prpAfterEcmFactor || continueECM) divg (factor, N);
	if (prpAfterEcmFactor && isProbablePrime (&ecmdata.gwdata, N)) {
		OutputBoth (thread_num, "Cofactor is a probable prime!\n");
		continueECM = FALSE;
	}

/* Send assignment result to the server.  To avoid flooding the server with small factors from users needlessly redoing */
/* factoring work, make sure the factor is more than 60 bits or so. */

	if (!QA_IN_PROGRESS && (strlen (str) >= 18 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		truncated_strcpy (pkt.message, sizeof (pkt.message), msg);
		pkt.result_type = PRIMENET_AR_ECM_FACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		truncated_strcpy (pkt.factor, sizeof (pkt.factor), str);
		pkt.B1 = (double) ecmdata.B;
		pkt.B2 = (double) ecmdata.average_B2;
		pkt.curves = ecmdata.curve;
		pkt.stage = stage;
		pkt.fftlen = gwfftlen (&ecmdata.gwdata);
		pkt.done = !continueECM;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);

/* If continuing ECM, subtract the curves we just reported from the worktodo count of curves to run.  Otherwise, delete all ECM entries */
/* for this number from the worktodo file. */

		if (continueECM) {
			unlinkSaveFiles (&ecmdata.write_save_file_state);
			w->curves_to_do -= ecmdata.curve;
			stop_reason = updateWorkToDoLine (thread_num, w);
			if (stop_reason) return (stop_reason);
			ecmdata.curve = 0;
		} else {
//bug - how to update worktodo such that all ECM's of this number are deleted???
		}
	}

/* Free memory */

bad_factor_recovery:
	free (str); str = NULL;
	free (msg); msg = NULL;
	free (factor); factor = NULL;

	clear_timer (timers, 0);

/* Since we found a factor, then we likely performed much fewer curves than expected.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	if (!continueECM) {
		unlinkSaveFiles (&ecmdata.write_save_file_state);
		stop_reason = STOP_WORK_UNIT_COMPLETE;
		invalidateNextRollingAverageUpdate ();
		goto exit;
	}

/* Do more curves despite finding a factor */

	goto more_curves;

/* Output a message saying we are restarting */

err:	OutputBoth (thread_num, "SUMOUT error occurred.\n");

/* Sleep five minutes before restarting */

	free (str); str = NULL;
	free (msg); msg = NULL;
	free (factor); factor = NULL;
	ecm_cleanup (&ecmdata);
	stop_reason = SleepFive (thread_num);
	if (stop_reason) return (stop_reason);

/* Restart from last save file */

	goto restart;
}

/* Read a file of ECM tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, sigma, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9991 to run the QA suite */

int ecm_QA (
	int	thread_num,
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	FILE	*fd;

/* Set the title */

	title (thread_num, "QA");

/* Open QA file */

	fd = fopen ("qa_ecm", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa_ecm' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	QA_TYPE = 0;
	for ( ; ; ) {
		struct work_unit w;
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		double	sigma;
		int	stop_reason;

/* Read a line from the file */

		n = 0;
		(void) fscanf (fd, "%lf,%lu,%lu,%ld,%lf,%lu,%lu,%lu,%s\n", &k, &b, &n, &c, &sigma, &B1, &B2_start, &B2_end, fac_str);
		if (n == 0) break;

/* If b is 1, set QA_TYPE */

		if (b == 1) {
			QA_TYPE = c;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		QA_FACTOR = allocgiant ((int) strlen (fac_str));
		ctog (fac_str, QA_FACTOR);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Do the ECM */

		if (B2_start < B1) B2_start = B1;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ECM;
		w.k = k;
		w.b = b;
		w.n = n;
		w.c = c;
		w.B1 = B1;
		w.B2_start = B2_start;
		w.B2 = B2_end;
		w.curves_to_do = 1;
		w.curve = sigma;
		QA_IN_PROGRESS = TRUE;
		stop_reason = ecm (0, sp_info, &w);
		QA_IN_PROGRESS = FALSE;
		free (QA_FACTOR);
		if (stop_reason != STOP_WORK_UNIT_COMPLETE) {
			fclose (fd);
			return (stop_reason);
		}
	}

/* Cleanup */

	fclose (fd);
	return (0);
}

/**************************************************************
 *	P-1 Functions
 **************************************************************/

/* Data maintained during P-1 process */

#define PM1_STATE_STAGE0	0	/* In stage 1, computing 3^exp using a precomputed mpz exp */
#define PM1_STATE_STAGE1	1	/* In stage 1, processing larger primes */
#define PM1_STATE_MIDSTAGE	2	/* Between stage 1 and stage 2 */
#define PM1_STATE_STAGE2	3	/* In middle of stage 2 (processing a bit array) */
#define PM1_STATE_GCD		4	/* Stage 2 GCD */
#define PM1_STATE_DONE		5	/* P-1 job complete */

typedef struct {
	gwhandle gwdata;	/* GWNUM handle */
	int	thread_num;	/* Worker thread number */
	struct work_unit *w;	/* Worktodo.txt entry */
	int	state;		/* One of the states listed above */
	uint64_t B;		/* Bound #1 (a.k.a. B1) */
	uint64_t C;		/* Bound #2 (a.k.a. B2 */
	uint64_t interim_B;	/* B1 we are currently calculating (equals B except when finishing stage 1 from a save file using a different B1). */
	uint64_t B_done;	/* We have completed calculating 3^e to this bound #1 */
	uint64_t first_C_start;	/* First stage 2 starting point (equals B except when worktodo.txt specifies starting point for bound #2). */
	int	optimal_B2;	/* TRUE if we calculate optimal bound #2 given currently available memory.  FALSE for a fixed bound #2. */
	unsigned long max_stage0_prime; /* Maximum small prime that can be used in stage 0 exponent calculation. */
	unsigned long stage0_bitnum; /* Bit number in stage 0 exponent to process next */
	readSaveFileState read_save_file_state;	/* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	void	*sieve_info;	/* Prime number sieve */
	uint64_t stage1_prime;	/* Prime number being processed */

	gwnum	x;		/* The stage 1 value being computed */
	gwnum	gg;		/* The stage 2 accumulated value */

	int	D;		/* Stage 2 loop size */
	int	E;		/* Suyama's power in stage 2 */
	int	totrels;	/* Number relatively prime nQx values used */
	uint64_t B2_start;	/* Starting point of first D section to be processed in stage 2 (an odd multiple of D/2) */
	uint64_t numDsections;	/* Number of D sections to process in stage 2 */
	uint64_t Dsection;	/* Current D section being processed in stage 2 */
	uint64_t interim_C;	/* B2 we are currently calculating (equals C except when finishing stage 2 a save file using different B2) */

	char	*bitarray;	/* Bit array for prime pairings in each D section */
	uint64_t bitarraymaxDsections;	/* Maximum number of D sections per bit array */
	uint64_t bitarrayfirstDsection; /* First D section in the bit array */
	uint64_t first_relocatable; /* First relocatable prime (same as B1 unless bit arrays must be split or mem change caused a replan) */
	uint64_t last_relocatable; /* Last relocatable prime for filling bit arrays (unless mem change causes a replan) */
	uint64_t C_done;	/* Stage 2 completed thusfar (updated every D section that is completed) */

	int	stage2_numvals;	/* Number of gwnums used in stage 2 */
	gwnum	*nQx;		/* Array of relprime data used in stage 2 */
	gwnum	*eQx;		/* Array of data used in stage 2 to compute multiples of D. */

	double	pct_mem_to_use;	/* If we get memory allocation errors, we progressively try using less and less. */
} pm1handle;

/* Perform cleanup functions. */

void pm1_cleanup (
	pm1handle *pm1data)
{

/* Free memory */

	free (pm1data->nQx), pm1data->nQx = NULL;
	free (pm1data->eQx), pm1data->eQx = NULL;
	free (pm1data->bitarray), pm1data->bitarray = NULL;
	gwdone (&pm1data->gwdata);
	end_sieve (pm1data->sieve_info), pm1data->sieve_info = NULL;
}

/* Raises number to the given power */

void pm1_mul (
	pm1handle *pm1data,
	gwnum	xx,
	gwnum	orig_xx_fft,
	uint64_t n)
{
	uint64_t c;

/* Find most significant bit and then ignore it */

	c = 1;
	c <<= 63;
	while ((c&n) == 0) c >>= 1;
	c >>= 1;

/* Handle the second most significant bit */

//GW: use exponentiate??
	gwstartnextfft (&pm1data->gwdata, c > 1);
	gwfft (&pm1data->gwdata, xx, orig_xx_fft);
	gwfftfftmul (&pm1data->gwdata, orig_xx_fft, orig_xx_fft, xx);
	if (c&n) gwfftmul (&pm1data->gwdata, orig_xx_fft, xx);
	c >>= 1;

/* Do the rest of the bits */

	while (c) {
		gwstartnextfft (&pm1data->gwdata, c > 1);
		gwsquare (&pm1data->gwdata, xx);
		if (c&n) gwfftmul (&pm1data->gwdata, orig_xx_fft, xx);
		c >>= 1;
	}
}

/* Code to init "finite differences" for computing successive values of x^(start+i*incr)^E */

int fd_init (
	pm1handle *pm1data,
	uint64_t start,
	unsigned long incr,
	gwnum	x)
{
	int	i, j;
	giant	p;

/* Treat each eQx[i] as a binary value and compute (start+i*incr)^e */

	for (i = 0; i <= (int) pm1data->E; i++) {
		uint64_t val;
		p = allocgiant (pm1data->E * 2);
		if (p == NULL) goto oom;
		val = start + i * incr;
		ulltog (val, p);
		for (j = 2; j <= (int) pm1data->E; j++) ullmulg (val, p);
		pm1data->eQx[i] = (gwnum) p;
	}		

/* Now do the finite differences */

	for (i = 1; i <= (int) pm1data->E; i++) {
		for (j = (int) pm1data->E; j >= i; j--) {
			subg ((giant) pm1data->eQx[j-1], (giant) pm1data->eQx[j]);
		}
	}

/* Now compute each x^difference */

	for (i = 0; i <= (int) pm1data->E; i++) {
		p = (giant) pm1data->eQx[i];
		pm1data->eQx[i] = gwalloc (&pm1data->gwdata);
		if (pm1data->eQx[i] == NULL) goto oom;

/* Test for easy cases */

		ASSERTG (!isZero (p));
		if (isone (p)) {
			gwcopy (&pm1data->gwdata, x, pm1data->eQx[i]);
		}

/* Find most significant bit and then ignore it */

//GW: use exponentiate???  convert exponentiate to gwmul3
//GW: exponentiate uses 4 gwnum temps which may be OK from a working set perspective
		else {
			int	len;

			len = bitlen (p);
			len--;

/* Perform the first squaring using the already FFTed value of x */
/* Then process the second and remaining bits of p */

//GW:  We can reuse this x^2 value (and if we use exponentiate, we could reuse the x^1,3,5,7 values)
			gwsquare2 (&pm1data->gwdata, x, pm1data->eQx[i], GWMUL_STARTNEXTFFT);
			for ( ; ; ) {
				if (bitval (p, len-1))
					gwmul3 (&pm1data->gwdata, x, pm1data->eQx[i], pm1data->eQx[i], GWMUL_STARTNEXTFFT);
				len--;
				if (len == 0) break;
				gwsquare2 (&pm1data->gwdata, pm1data->eQx[i], pm1data->eQx[i], GWMUL_STARTNEXTFFT);
			}
		}
		free (p);
	}
	return (0);

/* Out of memory exit path */

oom:	free (p);
	return (OutOfMemory (pm1data->thread_num));
}

/* Code to compute next x^(start+i*incr)^E value.  Value is returned in eQx[0] partially FFTed. */

void fd_next (
	pm1handle *pm1data)
{
	int	i;

	for (i = 0; i < (int) pm1data->E; i++) {
		gwmul3 (&pm1data->gwdata, pm1data->eQx[i], pm1data->eQx[i+1], pm1data->eQx[i], GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
	}
}

/* Terminate finite differences code */

void fd_term (
	pm1handle *pm1data)
{
	int	i;

/* Free each eQx[i] */

	for (i = 0; i <= (int) pm1data->E; i++)
		gwfree (&pm1data->gwdata, pm1data->eQx[i]);
}

/* Routines to create and read save files for a P-1 factoring job */

#define PM1_MAGICNUM	0x317a394b
//#define PM1_VERSION	2				/* Changed in 29.4 build 7 -- corrected calc_exp bug */
//#define PM1_VERSION	3				/* Changed in 29.8 build 8.  Configurable calc_exp max exponent */
#define PM1_VERSION	4				/* Complete overhaul in 30.4.  New stage 2 code with Preda optimizations */

void pm1_save (
	pm1handle *pm1data)
{
	int	fd;
	struct work_unit *w = pm1data->w;
	unsigned long sum = 0;

/* Create the intermediate file */

	fd = openWriteSaveFile (&pm1data->write_save_file_state);
	if (fd < 0) return;

/* Write the file header */

	if (!write_header (fd, PM1_MAGICNUM, PM1_VERSION, w)) goto writeerr;

/* Write the file data */

	if (! write_int (fd, pm1data->state, &sum)) goto writeerr;

	if (pm1data->state == PM1_STATE_STAGE0) {
		if (! write_longlong (fd, pm1data->interim_B, &sum)) goto writeerr;
		if (! write_long (fd, pm1data->max_stage0_prime, &sum)) goto writeerr;
		if (! write_long (fd, pm1data->stage0_bitnum, &sum)) goto writeerr;
	}

	else if (pm1data->state == PM1_STATE_STAGE1) {
		if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->interim_B, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->stage1_prime, &sum)) goto writeerr;
	}

	else if (pm1data->state == PM1_STATE_MIDSTAGE) {
		if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->C_done, &sum)) goto writeerr;
	}

	// Save everything necessary to restart stage 2 without calling pm1_stage2_impl again
	else if (pm1data->state == PM1_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->C_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->interim_C, &sum)) goto writeerr;
		if (! write_int (fd, pm1data->stage2_numvals, &sum)) goto writeerr;
		if (! write_int (fd, pm1data->totrels, &sum)) goto writeerr;
		if (! write_int (fd, pm1data->D, &sum)) goto writeerr;
		if (! write_int (fd, pm1data->E, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->first_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->last_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->B2_start, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->numDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->Dsection, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->bitarraymaxDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->bitarrayfirstDsection, &sum)) goto writeerr;
		// Output the truncated bit array
		bitarray_numDsections = pm1data->numDsections - pm1data->bitarrayfirstDsection;
		if (bitarray_numDsections > pm1data->bitarraymaxDsections) bitarray_numDsections = pm1data->bitarraymaxDsections; 		
		bitarray_len = divide_rounding_up (bitarray_numDsections * pm1data->totrels, 8);
		bitarray_start = divide_rounding_down ((pm1data->Dsection - pm1data->bitarrayfirstDsection) * pm1data->totrels, 8);
		if (! write_array (fd, pm1data->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto writeerr;
	}

	else if (pm1data->state == PM1_STATE_GCD) {
		if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->C_done, &sum)) goto writeerr;
	}

	else if (pm1data->state == PM1_STATE_DONE) {
		if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pm1data->C_done, &sum)) goto writeerr;
	}

/* Write the gwnum value used in stage 1 and stage 2.  There are occasions where x or gg may be in a partially FFTed state. */

	if (! write_gwnum (fd, &pm1data->gwdata, pm1data->x, &sum)) goto writeerr;
	if (pm1data->state >= PM1_STATE_MIDSTAGE && pm1data->state <= PM1_STATE_GCD) {
		if (! write_gwnum (fd, &pm1data->gwdata, pm1data->gg, &sum)) goto writeerr;
	}

/* Write the checksum, we're done */

	if (! write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (&pm1data->write_save_file_state, fd);
	return;

/* An error occurred.  Close and delete the current file. */

writeerr:
	deleteWriteSaveFile (&pm1data->write_save_file_state, fd);
}

/* Read a version 25 through version 30.3 save file */

int pm1_old_restore (
	pm1handle *pm1data,
	int	fd,
	unsigned long version,
	unsigned long filesum)
{
	unsigned long sum = 0;
	unsigned long state, bitarray_len, unused;
	uint64_t processed, B_done, B, C_done, unused64;

/* Read the first part of the save file, much will be ignored but must be read for backward compatibility */

	if (! read_long (fd, &state, &sum)) goto readerr;
	if (version == 2) pm1data->max_stage0_prime = 13333333;	// The hardwired value prior to version 29.8 build 8
	else if (! read_long (fd, &pm1data->max_stage0_prime, &sum)) goto readerr;
	if (! read_longlong (fd, &B_done, &sum)) goto readerr;
	if (! read_longlong (fd, &B, &sum)) goto readerr;
	if (! read_longlong (fd, &C_done, &sum)) goto readerr;
	if (! read_longlong (fd, &unused64, &sum)) goto readerr;	// C_start
	if (! read_longlong (fd, &unused64, &sum)) goto readerr;	// C_done
	if (! read_longlong (fd, &processed, &sum)) goto readerr;
	if (! read_long (fd, &unused, &sum)) goto readerr;		// D
	if (! read_long (fd, &unused, &sum)) goto readerr;		// E
	if (! read_long (fd, &unused, &sum)) goto readerr;		// rels_done
	if (! read_long (fd, &bitarray_len, &sum)) goto readerr;
	if (bitarray_len) {
		char *bitarray = (char *) malloc (bitarray_len);
		if (bitarray == NULL) goto readerr;
		if (! read_array (fd, bitarray, bitarray_len, &sum)) goto readerr;
		free (bitarray);
	}
	if (! read_longlong (fd, &unused64, &sum)) goto readerr;	// bitarray_first_number
	if (! read_long (fd, &unused, &sum)) goto readerr;		// pairs_set
	if (! read_long (fd, &unused, &sum)) goto readerr;		// pairs_done

/* Depending on the state, some of the values read above are not meaningful. */
/* In stage 0, only B and processed (bit number) are meaningful. */
/* In stage 1, only B_done, B, and processed (prime) are meaningful. */
/* In stage 2, only B_done is useful.  We cannot continue an old stage 2. */
/* When done, only B_done and C_done are meaningful. */

	if (state == 3) {				// PM1_STATE_STAGE0
		pm1data->state = PM1_STATE_STAGE0;
		pm1data->interim_B = B;
		pm1data->stage0_bitnum = (unsigned long) processed;
		// Version 29.4 build 7 changed the calc_exp algorithm which invalidates earlier save files that are in stage 0
		if (version == 1) {
			OutputBoth (pm1data->thread_num, "P-1 save file incompatible with this program version.  Restarting stage 1 from the beginning.\n");
			goto readerr;
		}
	} else if (state == 0) {			// PM1_STATE_STAGE1
//GW: This may not be compatible if we switch to sliding window / big calc_exp approach
		pm1data->state = PM1_STATE_STAGE1;
		pm1data->B_done = B_done;			//GW: set this to processed???
		pm1data->interim_B = B;
		pm1data->stage1_prime = processed;
	} else if (state == 1) {			// PM1_STATE_STAGE2
		pm1data->state = PM1_STATE_MIDSTAGE;
		pm1data->B_done = B_done;
		pm1data->C_done = B_done;
		// Current stage 2 code is incompatible with older save files
		OutputBoth (pm1data->thread_num, "Cannot continue stage 2 from old P-1 save file.  Restarting stage 2 from the beginning.\n");
	} else if (state == 2) {			// PM1_STATE_DONE
		pm1data->state = PM1_STATE_DONE;
		pm1data->B_done = B_done;
		pm1data->C_done = C_done;
	}

/* Read the gwnum values */

	pm1data->x = gwalloc (&pm1data->gwdata);
	if (pm1data->x == NULL) goto readerr;
	if (! read_gwnum (fd, &pm1data->gwdata, pm1data->x, &sum)) goto readerr;

	pm1data->gg = NULL;
	if (state == 1) {
		pm1data->gg = gwalloc (&pm1data->gwdata);
		if (pm1data->gg == NULL) goto readerr;
		if (! read_gwnum (fd, &pm1data->gwdata, pm1data->gg, &sum)) goto readerr;
	}

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);

/* All done */

	return (TRUE);

/* An error occurred.  Cleanup and return. */

readerr:
	_close (fd);
	return (FALSE);
}


/* Read a save file */

int pm1_restore (			/* For version 30.4 and later save files */
	pm1handle *pm1data)
{
	int	fd;
	struct work_unit *w = pm1data->w;
	unsigned long version;
	unsigned long sum = 0, filesum;

/* Open the intermediate file */

	fd = _open (pm1data->read_save_file_state.current_filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto err;

/* Read the file header */

	if (! read_magicnum (fd, PM1_MAGICNUM)) goto readerr;
	if (! read_header (fd, &version, w, &filesum)) goto readerr;
	if (version < 1 || version > PM1_VERSION) goto readerr;
	if (version < 4) return (pm1_old_restore (pm1data, fd, version, filesum));

/* Read the first part of the save file */

	if (! read_int (fd, &pm1data->state, &sum)) goto readerr;

/* Read state dependent data */

	if (pm1data->state == PM1_STATE_STAGE0) {
		if (! read_longlong (fd, &pm1data->interim_B, &sum)) goto readerr;
		if (! read_long (fd, &pm1data->max_stage0_prime, &sum)) goto readerr;
		if (! read_long (fd, &pm1data->stage0_bitnum, &sum)) goto readerr;
	}

	else if (pm1data->state == PM1_STATE_STAGE1) {
		if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->interim_B, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->stage1_prime, &sum)) goto readerr;
	}

	else if (pm1data->state == PM1_STATE_MIDSTAGE) {
		if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
	}

	else if (pm1data->state == PM1_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->interim_C, &sum)) goto readerr;
		if (! read_int (fd, &pm1data->stage2_numvals, &sum)) goto readerr;
		if (! read_int (fd, &pm1data->totrels, &sum)) goto readerr;
		if (! read_int (fd, &pm1data->D, &sum)) goto readerr;
		if (! read_int (fd, &pm1data->E, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->first_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->last_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->B2_start, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->numDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->Dsection, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->bitarraymaxDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->bitarrayfirstDsection, &sum)) goto readerr;
		// Read the truncated bit array
		bitarray_numDsections = pm1data->numDsections - pm1data->bitarrayfirstDsection;
		if (bitarray_numDsections > pm1data->bitarraymaxDsections) bitarray_numDsections = pm1data->bitarraymaxDsections;
		bitarray_len = divide_rounding_up (bitarray_numDsections * pm1data->totrels, 8);
		bitarray_start = divide_rounding_down ((pm1data->Dsection - pm1data->bitarrayfirstDsection) * pm1data->totrels, 8);
		pm1data->bitarray = (char *) malloc ((size_t) bitarray_len);
		if (pm1data->bitarray == NULL) goto readerr;
		if (! read_array (fd, pm1data->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto readerr;
	}

	else if (pm1data->state == PM1_STATE_GCD) {
		if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
	}

	else if (pm1data->state == PM1_STATE_DONE) {
		if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
	}

/* Read the gwnum value used in stage 1 and stage 2 */

	pm1data->x = gwalloc (&pm1data->gwdata);
	if (pm1data->x == NULL) goto readerr;
	if (! read_gwnum (fd, &pm1data->gwdata, pm1data->x, &sum)) goto readerr;

/* Read stage 2 accumulator gwnum */

	pm1data->gg = NULL;
	if (pm1data->state >= PM1_STATE_MIDSTAGE && pm1data->state <= PM1_STATE_GCD) {
		pm1data->gg = gwalloc (&pm1data->gwdata);
		if (pm1data->gg == NULL) goto readerr;
		if (! read_gwnum (fd, &pm1data->gwdata, pm1data->gg, &sum)) goto readerr;
	}

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);

/* All done */

	return (TRUE);

/* An error occurred.  Cleanup and return. */

readerr:
	_close (fd);
err:
	return (FALSE);
}


/* Compute the cost (in squarings) of a particular P-1 stage 2 implementation. */

struct pm1_stage2_cost_data {
	int	E;		/* Brent-Suyama power */
				/* Data returned from cost function follows */
	int	D;		/* D value for big steps */
	int	totrels;	/* Total number of relative primes used for pairing */
	int	stage2_numvals;	/* Total number of gwnum temporaries to be used in stage 2 */
	uint64_t B2_start;	/* Stage 2 start */
	uint64_t numDsections;	/* Number of D sections to process */
	uint64_t bitarraymaxDsections; /* Maximum number of D sections per bit array (in case work bit array must be split) */
};

double pm1_stage2_cost (
	int	D,		/* Stage 2 big step */
	uint64_t B2_start,	/* Stage 2 start point */
	uint64_t numDsections,	/* Number of D sections to process */
	uint64_t bitarraymaxDsections, /* Maximum number of D sections per bit array (if work bit array must be split) */
	int	numrels,	/* Number of relative primes less than D */
	double	multiplier,	/* Multiplier * numrels relative primes will be used in stage 2 */
	double	numpairs,	/* Estimated number of paired primes in stage 2 */
	double	numsingles,	/* Estimated number of prime singles in stage 2 */
	void	*data)		/* P-1 specific costing data */
{
	int	totrels;	/* Total number of relative primes in use */
	struct pm1_stage2_cost_data *cost_data = (struct pm1_stage2_cost_data *) data;
	double	cost;

	ASSERTG (cost_data->E >= 2);

/* Calculate the number of temporaries used for relative primes.  This will leave some available for modinv pooling. */

	totrels = (int) (numrels * multiplier);

/* Compute the nQx setup costs.  They are tiny and if inaccurately estimated can be safely ignored. */

	cost = (cost_data->E + 1) +

/* For the minimum number of relative primes for a D, there are D / 4 fd_next calls at a cost of E multiplies each fd_next call. */
/* If using more than the minimum number of relative primes, increase the cost proportionally. */
/* Also, each nQx value will be FFTed once (a half squaring). */

		multiplier * ((double) D / 4) * (cost_data->E) + totrels * 0.5 +

/* Compute the eQx setup costs.  To calculate eQx values, one fd_init is required with a start point of B. */

		(cost_data->E + 1) * cost_data->E * _log2 (B2_start) * 1.5;

/* Add in the cost of fd_next calls (#sections * E multiplies) */

	cost += (double) (numDsections * cost_data->E);

/* Finally, each prime pair and prime single costs one multiply */

	cost += numpairs + numsingles;

/* Pass 2 FFT multiplications seem to be at least 20% slower than the squarings in pass 1.  This is probably due */
/* to several factors.  These include: better L2 cache usage and no calls to the faster gwsquare routine. */
/* Nov, 2009: On my Macbook Pro, with exponents around 45M and using 800MB memory, pass2 squarings are 40% slower. */
/* Sept. 2020: gwsubmul4 implemented in assembly for AVX and later reducing the overhead to about 10%. */

	cost = cost * 1.10;

/* Return data P-1 implementation will need */
/* Stage 2 memory is totrels gwnums for the nQx array, E+1 gwnums to calculate eQx values, one gwnum for gg. */

	cost_data->stage2_numvals = totrels + cost_data->E + 2;
	cost_data->D = D;
	cost_data->totrels = (int) (numrels * multiplier);
	cost_data->B2_start = B2_start;
	cost_data->numDsections = numDsections;
	cost_data->bitarraymaxDsections = bitarraymaxDsections;

/* Return the resulting cost */

	return (cost);
}

/* Choose the most effective B2 for a P-1 run with a fixed B1 given the number of gwnums we are allowed to allocate. */
/* That is, find the B2 such that investing a fixed cost in either a larger B1 or B2 results in the same increase in chance of finding a factor. */

void pm1_choose_B2 (
	pm1handle *pm1data,
	unsigned long numvals)
{
	struct pm1_stage2_cost_data cost_data;	/* Extra data passed to P-1 costing function */
	struct pm1_stage2_efficiency {
		int	i;
		double	B2_cost;		/* Cost of stage 2 in squarings */
		double	fac_pct;		/* Percentage chance of finding a factor */
	} best[3];
	int	sieve_depth;
	char	buf[160];
	double	takeAwayBits;			/* Bits we get for free in smoothness of P-1 factor */

// Cost out a B2 value
#define p1eval(x,B2mult)	x.i = B2mult; \
				x.B2_cost = best_stage2_impl (pm1data->B, x.i * pm1data->B, numvals - (cost_data.E + 2), &pm1_stage2_cost, &cost_data); \
				x.fac_pct = pm1prob (takeAwayBits, sieve_depth, (double) pm1data->B, (double) (x.i * pm1data->B));

// Return TRUE if x is better than y.  Determined by seeing if taking the increased cost of y's higher B2 and investing it in increasing x's bounds
// results in a higher chance of finding a factor.
#define B1increase(x,y)	(int) ((y.B2_cost - x.B2_cost) / 1.44)			// Each B1 increase costs 1.44 squarings
#define p1compare(x,y)	(pm1prob (takeAwayBits, sieve_depth, (double) (pm1data->B + B1increase(x,y)), (double) (x.i * pm1data->B + B1increase(x,y))) > y.fac_pct)

/* Look for the best B2 which is likely between 5*B1 and 40*B1.  If optimal is not between these bounds, don't worry we'll locate the optimal spot anyway. */

	cost_data.E = IniGetInt (INI_FILE, "BrentSuyama", 2);
	if (cost_data.E < 2) cost_data.E = 2;
	if (cost_data.E > 24) cost_data.E = 24;
	takeAwayBits = isMersenne (pm1data->w->k, pm1data->w->b, pm1data->w->n, pm1data->w->c) ? log2 (pm1data->w->n) + 1.0 :
		       isGeneralizedFermat (pm1data->w->k, pm1data->w->b, pm1data->w->n, pm1data->w->c) ? log2 (pm1data->w->n) : 0.0;
	sieve_depth = (int) pm1data->w->sieve_depth;
	p1eval (best[0], 5);
	p1eval (best[1], 20);
	p1eval (best[2], 40);

/* Handle case where midpoint is worse than the start point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (p1compare (best[0], best[1])) {
		best[2] = best[1];
		p1eval (best[1], (best[0].i + best[2].i) / 2);
	}

/* Handle case where midpoint is worse than the end point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (!p1compare (best[1], best[2])) {
		best[0] = best[1];
		best[1] = best[2];
		p1eval (best[2], (int) (best[1].i + 10));
	}

/* Find the best B2.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].i + 2 != best[2].i) {
		struct pm1_stage2_efficiency midpoint;

		// Work on the bigger of the lower section and upper section
		if (best[1].i - best[0].i > best[2].i - best[1].i) {		// Work on lower section
			p1eval (midpoint, (best[0].i + best[1].i) / 2);
			if (p1compare (midpoint, best[1])) {			// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			p1eval (midpoint, (best[1].i + best[2].i) / 2);
			if (!p1compare (best[1], midpoint)) {			// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Return the best B2 */

	pm1data->C = best[1].i * pm1data->B;
	sprintf (buf, "With trial factoring done to 2^%d, optimal B2 is %d*B1 = %" PRIu64 ".\n", sieve_depth, best[1].i, pm1data->C);
	OutputStr (pm1data->thread_num, buf);
	sprintf (buf, "If no prior P-1, chance of a new factor is %.3g%%\n", best[1].fac_pct * 100.0);
	OutputStr (pm1data->thread_num, buf);
}
#undef p1eval
#undef B1increase
#undef p1compare

/* Choose the best implementation of P-1 stage 2 given the current memory settings.  We may decide there will never be enough memory. */
/* We may decide to wait for more memory to be available. */
/* We choose the best values for D and E that reduce the number of multiplications with the current memory constraints. */

int pm1_stage2_impl (
	pm1handle *pm1data)
{
	unsigned int memory, min_memory, desired_memory;	/* Memory is in MB */
	unsigned int bitarray_size, max_bitarray_size;		/* Bitarray memory in MB */
	int	numvals;			/* Number of gwnums we can allocate */
	struct pm1_stage2_cost_data cost_data;
	int	stop_reason;

/* Calculate (roughly) the memory required for the bit-array.  A typical stage 2 implementation is D=210, B2_start=B2/11, numrels=24, totrels=7*24. */
/* A typical optimal B2 is 40*B1.  This means bitarray memory = B2 * 10/11 / 210 * 24 * 7 / 8 bytes. */

	bitarray_size = divide_rounding_up ((unsigned int) ((pm1data->optimal_B2 ? 130 * pm1data->B : pm1data->C) * 1680 / 18480), 1 << 20);
	max_bitarray_size = IniGetInt (INI_FILE, "MaximumBitArraySize", 250);
	if (bitarray_size > max_bitarray_size) bitarray_size = max_bitarray_size;

/* Calculate the number of temporaries we can use in stage 2.  Base our decision on the current amount of memory available. */
/* If continuing from a stage 2 save file then we desire as many temporaries as the save file used.  Otherwise, assume 28 temporaries will */
/* provide us with a reasonable execution speed (D = 210).  We must have a minimum of 8 temporaries (D = 30). */

	min_memory = bitarray_size + cvt_gwnums_to_mem (&pm1data->gwdata, 8);
	desired_memory = bitarray_size + cvt_gwnums_to_mem (&pm1data->gwdata, pm1data->state >= PM1_STATE_STAGE2 ? pm1data->stage2_numvals : 28);
	stop_reason = avail_mem (pm1data->thread_num, min_memory, desired_memory, &memory);
	if (stop_reason) return (stop_reason);

/* Factor in the multiplier that we set to less than 1.0 when we get unexpected memory allocation errors. */
/* Make sure we can still allocate 8 temporaries. */

	memory = (unsigned int) (pm1data->pct_mem_to_use * (double) memory);
	if (memory < min_memory) return (avail_mem_not_sufficient (pm1data->thread_num, min_memory, desired_memory));
	if (memory < 8) memory = 8;

/* Output a message telling us how much memory is available */

	if (NUM_WORKER_THREADS > 1) {
		char	buf[100];
		sprintf (buf, "Available memory is %dMB.\n", memory);
		OutputStr (pm1data->thread_num, buf);
	}

/* Compute the number of gwnum temporaries we can allocate.  User nordi had over-allocating memory troubles on Linux testing M1277, presumably */
/* because our estimate genum size was too low.  As a work-around limit numvals to 100,000 by default. */

	numvals = cvt_mem_to_gwnums (&pm1data->gwdata, memory - bitarray_size);
	if (numvals < 8) numvals = 8;
	if (numvals > 100000) numvals = 100000;
	if (QA_TYPE) numvals = QA_TYPE;			/* Optionally override numvals for QA purposes */

/* Set C_done for future best_stage2_impl calls. */
/* Override B2 with optimal B2 based on amount of memory available. */

	if (pm1data->state == PM1_STATE_MIDSTAGE && pm1data->C_done == pm1data->B) {
		pm1data->C_done = pm1data->first_C_start;
		if (pm1data->optimal_B2) pm1_choose_B2 (pm1data, numvals);
	}

/* If are continuing from a save file that was in stage 2, check to see if we currently have enough memory to continue with the save file's */
/* stage 2 implementation.  Also check if we now have significantly more memory available and stage 2 is not near complete such that a new */
/* stage 2 implementation might give us a faster stage 2 completion. */

//GW: These are rather arbitrary heuristics
	if (pm1data->state >= PM1_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= pm1data->stage2_numvals &&				// We have enough memory and
	    (numvals < pm1data->stage2_numvals * 2 ||				// less than twice as much memory now available or
	     pm1data->Dsection >= pm1data->numDsections / 2))			// stage 2 more than half done
		return (0);							// Use old plan

/* Get the Brent-Suyama power.  We used to try various values of E (1, 2, 4, 6, 12). */
/* With Preda/Montgomery's optimization using gwnum temporaries to improve prime pairing we are better off only supporting E = 2. */

	cost_data.E = IniGetInt (INI_FILE, "BrentSuyama", 2);
	if (cost_data.E < 2) cost_data.E = 2;
	if (cost_data.E > 24) cost_data.E = 24;

/* Find the least costly stage 2 plan. */
/* Try various values of D until we find the best one.  E+1 gwnums are required for eQx calculations and one gwnum for gg. */

	best_stage2_impl (pm1data->C_done, pm1data->C, numvals - (cost_data.E + 2), &pm1_stage2_cost, &cost_data);

/* If are continuing from a save file that was in stage 2 and the new plan doesn't look significant better than the old plan, then */
/* we use the old plan and its partially completed bit array. */

	if (pm1data->state >= PM1_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= pm1data->stage2_numvals &&				// We have enough memory and
	    cost_data.stage2_numvals < pm1data->stage2_numvals * 2)		// new plan does not use significantly more memory
		return (0);							// Use old plan

/* If are continuing from a save file that was in stage 2, toss the save file's bit array. */

	if (pm1data->state >= PM1_STATE_STAGE2) {
		free (pm1data->bitarray);
		pm1data->bitarray = NULL;
	}

/* Set all the variables needed for this stage 2 plan */

	pm1data->interim_C = pm1data->C;
	pm1data->stage2_numvals = cost_data.stage2_numvals;
	pm1data->D = cost_data.D;
	pm1data->E = cost_data.E;
	pm1data->totrels = cost_data.totrels;
	pm1data->B2_start = cost_data.B2_start;
	pm1data->numDsections = cost_data.numDsections;
	pm1data->bitarraymaxDsections = cost_data.bitarraymaxDsections;
	pm1data->Dsection = 0;

	if (pm1data->state < PM1_STATE_STAGE2) {
		pm1data->first_relocatable = pm1data->first_C_start;
		pm1data->last_relocatable = pm1data->B2_start;
		pm1data->C_done = pm1data->B2_start;
	}

/* Create a bitmap maximizing the prime pairings */

	pm1data->bitarrayfirstDsection = pm1data->Dsection;
	stop_reason = fill_work_bitarray (pm1data->thread_num, &pm1data->sieve_info, pm1data->D, pm1data->totrels,
					  pm1data->first_relocatable, pm1data->last_relocatable, pm1data->C_done, pm1data->C,
					  pm1data->bitarraymaxDsections, &pm1data->bitarray);
	if (stop_reason) return (stop_reason);

	return (0);
}


/* Recursively compute exp used in initial 3^exp calculation of a P-1 or P+1 factoring run.  We include 2*n in exp, since P-1 factoring */
/* benefits from Mersenne number factors being 1 mod 2n.  Generalized Fermats also benefit in P-1 factoring. */

void calc_exp (
	void	*sieve_info,	/* Sieve to use calculating primes */
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c,		/* C in K*B^N+C */
	mpz_t	g,		/* Variable to accumulate multiplied small primes */
	uint64_t B1,		/* P-1 stage 1 bound */
	uint64_t *p,		/* Variable to fetch next small prime into */
	unsigned long lower,
	unsigned long upper)
{
	unsigned long len;

/* Compute the number of result bits we are to calculate */

	len = upper - lower;

/* Use recursion to compute the exponent.  This will perform better because mpz_mul will be handling arguments of equal size. */

	if (len >= 1024) {
		mpz_t	x;
		calc_exp (sieve_info, k, b, n, c, g, B1, p, lower, lower + (len >> 1));
		mpz_init (x);
		calc_exp (sieve_info, k, b, n, c, x, B1, p, lower + (len >> 1), upper);
		mpz_mul (g, x, g);
		mpz_clear (x);
		return;
	}

/* For Mersenne numbers, 2^n-1, make sure we include 2n in the calculated exponent (since factors */
/* are of the form 2kn+1).  For generalized Fermat numbers, b^n+1 (n is a power of 2), make sure n */
/* is included in the calculated exponent as factors are of the form kn+1 (actually forum posters */
/* have pointed out that Fermat numbers should include 4n and generalized Fermat should include 2n). */
/* Heck, maybe other forms may also need n included, so just always include 2n -- it is very cheap. */

	if (lower == 0) mpz_set_ui (g, 2*n);
	else mpz_set_ui (g, 1);

/* Find all the primes in the range and use as many powers as possible */

	for ( ; *p <= B1 && mpz_sizeinbase (g, 2) < len; *p = sieve (sieve_info)) {
		uint64_t val, max;
		val = *p;
		max = B1 / *p;
		while (val <= max) val *= *p;
		if (sizeof (unsigned int) == 4 && val > 0xFFFFFFFF) {
			mpz_t	mpz_val;
			mpz_init_set_d (mpz_val, (double) val);		/* Works for B1 up to 2^53 */
			mpz_mul (g, g, mpz_val);
			mpz_clear (mpz_val);
		} else
			mpz_mul_ui (g, g, (unsigned int) val);
	}
}

/*****************************************************************************/
/*                         Main P-1 routine				     */
/*****************************************************************************/

int pminus1 (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	pm1handle pm1data;
	giant	N;		/* Number being factored */
	giant	factor;		/* Factor found, if any */
	mpz_t	exp;
	int	exp_initialized;
	uint64_t stage0_limit;
	unsigned int memused;
	unsigned long SQRT_B;
	int	i, totrels;
	unsigned long len;
	unsigned long error_recovery_mode = 0;
	gwnum	tmp_gwnum;
	char	filename[32], buf[255], JSONbuf[4000], testnum[100];
	int	res, stop_reason, saving, near_fft_limit, echk;
	double	one_over_len, one_over_B, base_pct_complete, one_bit_pct;
	double	last_output, last_output_t, last_output_r;
	double	allowable_maxerr, output_frequency, output_title_frequency;
	int	first_iter_msg;
	int	msglen;
	char	*str, *msg;
	double	timers[2];

/* Output a blank line to separate multiple P-1 runs making the result more readable */

	OutputStr (thread_num, "\n");

/* Clear pointers to allocated memory (so common error exit code knows what to free) */

	N = NULL;
	factor = NULL;
	str = NULL;
	msg = NULL;
	exp_initialized = FALSE;

/* Begin initializing P-1 data structure */
/* Choose a default value for the second bound if none was specified */

	memset (&pm1data, 0, sizeof (pm1handle));
	pm1data.thread_num = thread_num;
	pm1data.w = w;
	pm1data.B = (uint64_t) w->B1;
	pm1data.C = (uint64_t) w->B2;
	if (pm1data.B < 30) {
		OutputStr (thread_num, "Using minimum bound #1 of 30\n");
		pm1data.B = 30;
	}
	if (pm1data.C == 0) pm1data.C = pm1data.B * 100;
	if (pm1data.C < pm1data.B) pm1data.C = pm1data.B;
	pm1data.first_C_start = (uint64_t) w->B2_start;
	if (pm1data.first_C_start < pm1data.B) pm1data.first_C_start = pm1data.B;
	pm1data.pct_mem_to_use = 1.0;				// Use as much memory as we can unless we get allocation errors

/* Decide if we will calculate an optimal B2 when stage 2 begins.  We do this by default for P-1 work where we know how much TF has been done. */

	pm1data.optimal_B2 = (!QA_IN_PROGRESS && w->work_type == WORK_PMINUS1 && pm1data.first_C_start == pm1data.B &&
			      w->sieve_depth > 50 && IniGetInt (INI_FILE, "Pminus1BestB2", 1));
	if (pm1data.optimal_B2) pm1data.C = 100 * pm1data.B;	// A guess to use for calling start_sieve_with_limit

/* Compute the number we are factoring */

	stop_reason = setN (thread_num, w, &N);
	if (stop_reason) goto exit;

/* Output startup message, but only if work type is P-1.  Pfactor work type has already output a startup message. */

	gw_as_string (testnum, w->k, w->b, w->n, w->c);
	sprintf (buf, "%s P-1", testnum);
	title (thread_num, buf);
	if (w->work_type == WORK_PMINUS1) {
		if (pm1data.C <= pm1data.B)
			sprintf (buf, "P-1 on %s with B1=%" PRIu64 "\n", testnum, pm1data.B);
		else if (pm1data.optimal_B2)
			sprintf (buf, "P-1 on %s with B1=%" PRIu64 ", B2=TBD\n", testnum, pm1data.B);
		else
			sprintf (buf, "P-1 on %s with B1=%" PRIu64 ", B2=%" PRIu64 "\n", testnum, pm1data.B, pm1data.C);
		OutputStr (thread_num, buf);
		if (w->sieve_depth > 0.0 && !pm1data.optimal_B2) {
			double prob = guess_pminus1_probability (w);
			sprintf (buf, "Chance of finding a factor is an estimated %.3g%%\n", prob * 100.0);
			OutputStr (thread_num, buf);
		}
	}

/* Clear all timers */

restart:
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init filename.  This is a little kludgy as we want to generate a P-1 save file that does not conflict with an LL or PRP save file name. */
/* Both save files can exist at the same time when stage 2 is delayed waiting for more memory. */

	tempFileName (w, filename);
	filename[0] = 'm';

/* Perform setup functions.  This includes decding how big an FFT to use, allocating memory, calling the FFT setup code, etc. */

/* Setup the assembly code */

	gwinit (&pm1data.gwdata);
	gwset_sum_inputs_checking (&pm1data.gwdata, SUM_INPUTS_ERRCHK);
	if (IniGetInt (LOCALINI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&pm1data.gwdata);
	if (IniGetInt (INI_FILE, "HyperthreadPrefetch", 0)) gwset_hyperthread_prefetch (&pm1data.gwdata);
	if (HYPERTHREAD_LL) {
		sp_info->normal_work_hyperthreads = IniGetInt (LOCALINI_FILE, "HyperthreadLLcount", CPU_HYPERTHREADS);
		gwset_will_hyperthread (&pm1data.gwdata, sp_info->normal_work_hyperthreads);
	}
	gwset_bench_cores (&pm1data.gwdata, NUM_CPUS);
	gwset_bench_workers (&pm1data.gwdata, NUM_WORKER_THREADS);
	if (ERRCHK) gwset_will_error_check (&pm1data.gwdata);
	else gwset_will_error_check_near_limit (&pm1data.gwdata);
	gwset_num_threads (&pm1data.gwdata, CORES_PER_TEST[thread_num] * sp_info->normal_work_hyperthreads);
	gwset_thread_callback (&pm1data.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&pm1data.gwdata, sp_info);
	gwset_safety_margin (&pm1data.gwdata, IniGetFloat (INI_FILE, "ExtraSafetyMargin", 0.0));
	gwset_minimum_fftlen (&pm1data.gwdata, w->minimum_fftlen);
	res = gwsetup (&pm1data.gwdata, w->k, w->b, w->n, w->c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}

/* A kludge so that the error checking code is not as strict. */

	pm1data.gwdata.MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);

/* More miscellaneous initializations */

	last_output = last_output_t = last_output_r = 0;
	gw_clear_fft_count (&pm1data.gwdata);
	first_iter_msg = TRUE;
	calc_output_frequencies (&pm1data.gwdata, &output_frequency, &output_title_frequency);

/* Output message about the FFT length chosen */

	{
		char	fft_desc[200];
		gwfft_description (&pm1data.gwdata, fft_desc);
		sprintf (buf, "Using %s\n", fft_desc);
		OutputStr (thread_num, buf);
	}

/* If we are near the maximum exponent this fft length can test, then we will roundoff check all multiplies */

	near_fft_limit = exponent_near_fft_limit (&pm1data.gwdata);
	gwsetnormroutine (&pm1data.gwdata, 0, ERRCHK || near_fft_limit, 0);

/* Figure out the maximum round-off error we will allow.  By default this is 27/64 when near the FFT limit and 26/64 otherwise. */
/* We've found that this default catches errors without raising too many spurious error messages.  We let the user override */
/* this default for user "Never Odd Or Even" who tests exponents well beyond an FFT's limit.  He does his error checking by */
/* running the first-test and double-check simultaneously. */

	allowable_maxerr = IniGetFloat (INI_FILE, "MaxRoundoffError", (float) (near_fft_limit ? 0.421875 : 0.40625));

/* Check for a save file and read the save file.  If there is an error */
/* reading the file then restart the P-1 factoring job from scratch. */
/* Limit number of backup files we try */
/* to read in case there is an error deleting bad save files. */

	readSaveFileStateInit (&pm1data.read_save_file_state, thread_num, filename);
	writeSaveFileStateInit (&pm1data.write_save_file_state, filename, 0);
	for ( ; ; ) {
		if (! saveFileExists (&pm1data.read_save_file_state)) {
			/* If there were save files, they are all bad.  Report a message */
			/* and temporarily abandon the work unit.  We do this in hopes that */
			/* we can successfully read one of the bad save files at a later time. */
			/* This sounds crazy, but has happened when OSes get in a funky state. */
			if (pm1data.read_save_file_state.a_non_bad_save_file_existed) {
				OutputBoth (thread_num, ALLSAVEBAD_MSG);
				return (0);
			}
			/* No save files existed, start from scratch. */
			break;
		}

		if (!pm1_restore (&pm1data)) {
			/* Close and rename the bad save file */
			saveFileBad (&pm1data.read_save_file_state);
			continue;
		}

/* Handle stage 0 save files.  If the B values do not match we use the bound from the save file for now.  Later we'll do more B if necessary. */

		if (pm1data.state == PM1_STATE_STAGE0) {
			goto restart0;
		}

/* To avoid an infinite loop of repeatable roundoff errors, we square the value read in from the P-1 save file. */
/* This won't affect our final results, but will change the FFT data. */

//GW: This wont work if stage 1 is changed to sliding window
		if (error_recovery_mode) {
			gwmul3_carefully (&pm1data.gwdata, pm1data.x, pm1data.x, pm1data.x, 0);
			pm1_save (&pm1data);
			error_recovery_mode = 0;
		}

/* Handle stage 1 save files.  If the save file that had a higher B1 target then we can reduce the target B1 to the desired B1. */

		if (pm1data.state == PM1_STATE_STAGE1) {
			if (pm1data.interim_B > pm1data.B) pm1data.interim_B = pm1data.B;
			goto restart1;
		}

/* Handle between stages save files */

		if (pm1data.state == PM1_STATE_MIDSTAGE) {
			if (pm1data.B > pm1data.B_done) {
				gwfree (&pm1data.gwdata, pm1data.gg), pm1data.gg = NULL;
				goto more_B;
			}
			goto restart3b;
		}

/* Handle stage 2 save files */

		if (pm1data.state == PM1_STATE_STAGE2) {

/* If B is larger than the one in the save file, then do more stage 1 processing.  Since this is very upsetting to */
/* an LL/PRP tester that has already begun stage 2 only do this for the non-LL/PRP tester. */

			if (pm1data.B > pm1data.B_done && w->work_type == WORK_PMINUS1) {
				gwfree (&pm1data.gwdata, pm1data.gg), pm1data.gg = NULL;
				free (pm1data.bitarray), pm1data.bitarray = NULL;
				goto more_B;
			}

/* If B is different than the one in the save file, then use the one in the save file rather than discarding all the work done thusfar in stage 2. */

			if (pm1data.B != pm1data.B_done) {
				pm1data.B = pm1data.B_done;
				sprintf (buf, "Ignoring suggested B1 value, using B1=%" PRIu64 " from the save file\n", pm1data.B);
				OutputStr (thread_num, buf);
			}

/* If LL testing and bound #2 has changed then use the original bound #2. */
/* If explicit P-1 testing and bound #2 is larger in the save file then use the original bound #2. */
/* The user doing explicit P-1 testing that wants to discard the stage 2 work he has done thusfar */
/* and reduce the stage 2 bound must manually delete the save file. */

			if ((w->work_type != WORK_PMINUS1 && pm1data.C != pm1data.interim_C) ||
			    (w->work_type == WORK_PMINUS1 && pm1data.C < pm1data.interim_C)) {
				pm1data.C = pm1data.interim_C;
				sprintf (buf, "Ignoring suggested B2 value, using B2=%" PRIu64 " from the save file\n", pm1data.C);
				OutputStr (thread_num, buf);
			}

/* Resume stage 2 */

			if (pm1data.optimal_B2) {
				pm1data.C = pm1data.interim_C;
				sprintf (buf, "Resuming P-1 in stage 2 with B2=%" PRIu64 "\n", pm1data.interim_C);
				OutputStr (thread_num, buf);
			}
			goto restart3b;
		}

/* Handle stage 2 GCD save files */

		if (pm1data.state == PM1_STATE_GCD) {
			if (pm1data.optimal_B2) pm1data.C = pm1data.C_done;
			goto restart4;
		}

/* Handle case where we have a completed save file (the PM1_STATE_DONE state) */
/* Note: if first_C_start != B then the user is using the undocumented feature of doing stage 2 in pieces.  Assume he knows what he is doing. */

		ASSERTG (pm1data.state == PM1_STATE_DONE);
		if (pm1data.B > pm1data.B_done) goto more_B;
		if (pm1data.C > pm1data.C_done) {
			pm1data.state = PM1_STATE_MIDSTAGE;
			if (pm1data.first_C_start == pm1data.B) pm1data.first_C_start = pm1data.C_done;
			goto restart3a;
		}

/* The save file indicates we've tested to these bounds already */

		sprintf (buf, "%s already tested to B1=%" PRIu64 " and B2=%" PRIu64 ".\n",
			 gwmodulo_as_string (&pm1data.gwdata), pm1data.B_done, pm1data.C_done);
		OutputBoth (thread_num, buf);
		goto done;
	}

/* Start this P-1 run from scratch starting with x = 3 */

	strcpy (w->stage, "S1");
	w->pct_complete = 0.0;
	pm1data.state = PM1_STATE_STAGE0;
	pm1data.interim_B = pm1data.B;
	pm1data.stage0_bitnum = 0;
	pm1data.x = gwalloc (&pm1data.gwdata);
	if (pm1data.x == NULL) goto oom;
	dbltogw (&pm1data.gwdata, 3.0, pm1data.x);

/* Stage 0 pre-calculates an exponent that is the product of small primes.  Our default uses only small */
/* primes below 40,000,000 (roughly 60 million bits).  This is configurable starting in version 29.8 build 8. */

	pm1data.max_stage0_prime = IniGetInt (INI_FILE, "MaxStage0Prime", 40000000);

/* First restart point.  Compute the big exponent (product of small primes).  Then compute 3^exponent. */
/* The exponent always contains 2*p.  We only use primes below (roughly) max_stage0_prime.  The rest of the */
/* exponentiation will be done one prime at a time in the second part of stage 1. */
/* This stage uses 2 transforms per exponent bit. */

restart0:
	pm1data.B_done = 0;
	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&pm1data.gwdata, 1));
	start_timer (timers, 0);
	start_timer (timers, 1);
	stop_reason = start_sieve_with_limit (thread_num, 2, (uint32_t) sqrt ((double) pm1data.C), &pm1data.sieve_info);
	if (stop_reason) goto exit;
	pm1data.stage1_prime = sieve (pm1data.sieve_info);
	stage0_limit = (pm1data.interim_B > pm1data.max_stage0_prime) ? pm1data.max_stage0_prime : pm1data.interim_B;
	mpz_init (exp);  exp_initialized = TRUE;
	calc_exp (pm1data.sieve_info, w->k, w->b, w->n, w->c, exp, pm1data.interim_B, &pm1data.stage1_prime, 0, (unsigned long) (stage0_limit * 1.5));

/* Find number of bits, ignoring the most significant bit.  Init variables used in calculating percent complete. */

	len = (unsigned long) mpz_sizeinbase (exp, 2) - 1;
	one_over_len = 1.0 / (double) len;
	if (pm1data.stage1_prime < pm1data.B) one_over_len *= (double) pm1data.stage1_prime / (double) pm1data.B;

/* Now take the exponent and raise x to that power */

	gwsetmulbyconst (&pm1data.gwdata, 3);
	while (pm1data.stage0_bitnum < len) {

/* To avoid an infinite loop of repeatable roundoff errors, carefully get us past the offending iteration. */

		if (error_recovery_mode && pm1data.stage0_bitnum == error_recovery_mode) {
			gwstartnextfft (&pm1data.gwdata, FALSE);
			gwsetnormroutine (&pm1data.gwdata, 0, 0, mpz_tstbit (exp, len - pm1data.stage0_bitnum - 1));
			gwsquare_carefully (&pm1data.gwdata, pm1data.x);
			error_recovery_mode = 0;
			saving = TRUE;
		}

/* Set various flags.  They control whether error-checking or the next FFT can be started. */

		else {
			stop_reason = stopCheck (thread_num);
			saving = testSaveFilesFlag (thread_num);
			echk = stop_reason || saving || ERRCHK || near_fft_limit || ((pm1data.stage0_bitnum & 127) == 64);

/* Either square x or square x and multiply it by three. */

#ifndef SERVER_TESTING
			gwstartnextfft (&pm1data.gwdata, !stop_reason && !saving && pm1data.stage0_bitnum+1 != error_recovery_mode && pm1data.stage0_bitnum+1 != len);
			gwsetnormroutine (&pm1data.gwdata, 0, echk, mpz_tstbit (exp, len - pm1data.stage0_bitnum - 1));
			if (pm1data.stage0_bitnum < 30) gwsquare_carefully (&pm1data.gwdata, pm1data.x);
			else gwsquare (&pm1data.gwdata, pm1data.x);
#endif
		}

/* Test for an error */

		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;

/* Calculate our stage 1 percentage complete */

		pm1data.stage0_bitnum++;
		w->pct_complete = (double) pm1data.stage0_bitnum * one_over_len;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P-1 stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata));
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pm1data.gwdata);
		}

/* Every N squarings, output a progress report */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			if (pm1data.stage0_bitnum > 1)
				OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pm1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Every N squarings, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pm1data.gwdata);
		}

/* Check for escape and/or if its time to write a save file */

		if (stop_reason || saving) {
			pm1_save (&pm1data);
			if (stop_reason) goto exit;
		}
	}

/* If roundoff error recovery returned to restart0, but the roundoff error */
/* occurs after the above loop, then square the value and create a save file. */
/* This won't affect our final results, but will change the FFT data. */

	if (error_recovery_mode) {
		gwmul3_carefully (&pm1data.gwdata, pm1data.x, pm1data.x, pm1data.x, 0);
		pm1_save (&pm1data);
		error_recovery_mode = 0;
	}

/* Do stage 0 cleanup */

	gwsetnormroutine (&pm1data.gwdata, 0, ERRCHK || near_fft_limit, 0);
	mpz_clear (exp), exp_initialized = FALSE;
	end_timer (timers, 0);
	end_timer (timers, 1);

/* Set up sieving tart point for next section */

	pm1data.stage1_prime--;				// Stage1_prime was not included by calc_exp.  Back up one for next sieve call.

/* Do the larger primes of stage 1.  This stage uses 2.5 transforms per exponent bit.  Proceed until interim_B is finished. */

//GW  more efficient would be to build another big exponent and do a sliding window exponentiate!  Far less than 2.5 transforms / bit.
// We could (possibly) eliminate interim_B and put it's primes in the big exponent.
// For really large B, repeat until done

restart1:
	strcpy (w->stage, "S1");
	one_over_B = 1.0 / (double) pm1data.B;
	w->pct_complete = (double) pm1data.stage1_prime * one_over_B;
	start_timer (timers, 0);
	start_timer (timers, 1);
	pm1data.state = PM1_STATE_STAGE1;
	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&pm1data.gwdata, 2));
	tmp_gwnum = gwalloc (&pm1data.gwdata);
	if (tmp_gwnum == NULL) goto oom;
	stop_reason = start_sieve_with_limit (thread_num, pm1data.stage1_prime + 1, (uint32_t) sqrt((double) pm1data.C), &pm1data.sieve_info);
	if (stop_reason) goto exit;
	SQRT_B = (unsigned long) sqrt ((double) pm1data.interim_B);
	for ( ; ; ) {
		pm1data.stage1_prime = sieve (pm1data.sieve_info);
		if (pm1data.stage1_prime > pm1data.interim_B) break;

/* Test for user interrupt, save files, and error checking */

		stop_reason = stopCheck (thread_num);
		saving = testSaveFilesFlag (thread_num);
		echk = stop_reason || saving || ERRCHK || near_fft_limit || ((pm1data.stage1_prime & 127) == 127);
		gwsetnormroutine (&pm1data.gwdata, 0, echk, 0);

/* Apply as many powers of prime as long as prime^n <= B */

		if (pm1data.stage1_prime > pm1data.B_done) {
			pm1_mul (&pm1data, pm1data.x, tmp_gwnum, pm1data.stage1_prime);
		}
		if (pm1data.stage1_prime <= SQRT_B) {
			uint64_t mult, max;
			mult = pm1data.stage1_prime;
			max = pm1data.interim_B / pm1data.stage1_prime;
			for ( ; ; ) {
				mult *= pm1data.stage1_prime;
				if (mult > pm1data.B_done) {
					pm1_mul (&pm1data, pm1data.x, tmp_gwnum, pm1data.stage1_prime);
				}
				if (mult > max) break;
			}
		}

/* Test for an error */

		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;

/* Calculate our stage 1 percentage complete */

		w->pct_complete = (double) pm1data.stage1_prime * one_over_B;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P-1 stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata));
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pm1data.gwdata);
		}

/* Every N primes, output a progress report */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pm1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Every N primes, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pm1data.gwdata);
		}

/* Check for escape and/or if its time to write a save file */

		if (stop_reason || saving) {
			pm1_save (&pm1data);
			if (stop_reason) goto exit;
		}
	}
	gwfree (&pm1data.gwdata, tmp_gwnum);
	pm1data.B_done = pm1data.interim_B;
	end_timer (timers, 0);
	end_timer (timers, 1);

/* Check for the rare case where we need to do even more stage 1.  This happens using a save file created with a smaller bound #1. */

	if (pm1data.B > pm1data.B_done) {
more_B:		pm1data.interim_B = pm1data.B;
		pm1data.stage1_prime = 2;
//GW - pct_complete resets  to 0 because of setting stage1_prime to 2
		goto restart1;
	}
	pm1data.C_done = pm1data.B;

/* Stage 1 complete, print a message */

	sprintf (buf, "%s stage 1 complete. %.0f transforms. Time: ", gwmodulo_as_string (&pm1data.gwdata), gw_get_fft_count (&pm1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));
	gw_clear_fft_count (&pm1data.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pm1data.gwdata);
	}

/* Check to see if we found a factor - do GCD (x-1, N) */

	strcpy (w->stage, "S1");
	w->pct_complete = 1.0;
	if (pm1data.C <= pm1data.B || (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "Stage1GCD", 1))) {
		if (w->work_type != WORK_PMINUS1) OutputStr (thread_num, "Starting stage 1 GCD - please be patient.\n");
		start_timer (timers, 0);
		gwaddsmall (&pm1data.gwdata, pm1data.x, -1);
		stop_reason = gcd (&pm1data.gwdata, thread_num, pm1data.x, N, &factor);
		gwaddsmall (&pm1data.gwdata, pm1data.x, 1);
		if (stop_reason) {
			pm1_save (&pm1data);
			goto exit;
		}
		end_timer (timers, 0);
		strcpy (buf, "Stage 1 GCD complete. Time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputStr (thread_num, buf);
		if (factor != NULL) goto bingo;
	}

/* Skip second stage if so requested */

	if (pm1data.C <= pm1data.B) goto msg_and_exit;

/* Since gwmul_carefully was called in stage 1, free the memory allocated for GW_RANDOM */

	gwfree_internal_memory (&pm1data.gwdata);

/*
   Stage 2:  Use ideas from Crandall, Zimmermann, Montgomery, and Preda on each prime below C.
   This code is more efficient the more memory you can give it.
   Inputs: x, the value at the end of stage 1
*/

/* Change state to between stage 1 and 2 */

	pm1data.state = PM1_STATE_MIDSTAGE;
	strcpy (w->stage, "S2");
	w->pct_complete = 0.0;

/* Initialize variables for second stage.  We set gg to x-1 in case the user opted to skip the GCD after stage 1. */

restart3a:
	pm1data.gg = gwalloc (&pm1data.gwdata);
	if (pm1data.gg == NULL) goto oom;
	gwcopy (&pm1data.gwdata, pm1data.x, pm1data.gg);
	gwaddsmall (&pm1data.gwdata, pm1data.gg, -1);

/* Restart here when in the middle of stage 2 */

restart3b:
more_C:
	start_timer (timers, 0);
	sprintf (buf, "%s P-1 stage 2 init", gwmodulo_as_string (&pm1data.gwdata));
	title (thread_num, buf);

/* Clear flag indicating we need to restart if the maximum amount of memory changes. */
/* Prior to this point we allow changing the optimal bounds of a Pfactor assignment. */

	clear_restart_if_max_memory_change (thread_num);

/* Test if we will ever have enough memory to do stage 2 based on the maximum available memory. */
/* Our minimum working set is one gwnum for gg, 4 for nQx, 3 for eQx. */

replan:	{
		unsigned long min_memory = cvt_gwnums_to_mem (&pm1data.gwdata, 8);
		if (max_mem (thread_num) < min_memory) {
			sprintf (buf, "Insufficient memory to ever run stage 2 -- %luMB needed.\n", min_memory);
			OutputStr (thread_num, buf);
			pm1data.C = pm1data.B_done;
			goto restart4;
		}
	}

/* Choose the best plan implementation given the currently available memory. */
/* This implementation could be "wait until we have more memory". */

	stop_reason = pm1_stage2_impl (&pm1data);
	if (stop_reason) {
		if (pm1data.state == PM1_STATE_MIDSTAGE) pm1_save (&pm1data);
		goto exit;
	}

/* Record the amount of memory this thread will be using in stage 2. */

	memused = cvt_gwnums_to_mem (&pm1data.gwdata, pm1data.stage2_numvals);
	memused += (_intmin (pm1data.numDsections - pm1data.bitarrayfirstDsection, pm1data.bitarraymaxDsections) * pm1data.totrels) >> 23;
	// To dodge possible infinite loop if pm1_stage2_impl allocates too much memory (it shouldn't), decrease the percentage of memory we are allowed to use
	// Beware that replaning may allocate a larger bitarray
	if (set_memory_usage (thread_num, MEM_VARIABLE_USAGE, memused)) {
		pm1data.pct_mem_to_use *= 0.99;
		free (pm1data.bitarray); pm1data.bitarray = NULL;
		goto replan;
	}

/* Output a useful message regarding memory usage */

	sprintf (buf, "Using %uMB of memory.\n", memused);
	OutputStr (thread_num, buf);

/* Initialize variables for second stage */

	// Calculate the percent completed by previous bit arrays
	base_pct_complete = (double) (pm1data.B2_start - pm1data.last_relocatable) / (double) (pm1data.C - pm1data.last_relocatable);
	// Calculate the percent completed by each bit in this bit array
	one_bit_pct = (1.0 - base_pct_complete) / (double) (pm1data.numDsections * pm1data.totrels);
	// Calculate the percent completed by previous bit arrays and the current bit array
	w->pct_complete = base_pct_complete + (double) (pm1data.Dsection * pm1data.totrels) * one_bit_pct;

/* Allocate arrays of pointers to stage 2 gwnums */

	pm1data.nQx = (gwnum *) malloc (pm1data.totrels * sizeof (gwnum));
	if (pm1data.nQx == NULL) goto lowmem;
	pm1data.eQx = (gwnum *) malloc ((pm1data.E + 1) * sizeof (gwnum));
	if (pm1data.eQx == NULL) goto lowmem;

/* Compute x^(first_relative_prime^e), x^(second_relative_prime^e), ..., x^(last_relative_prime^e) */

//gw is there a 2/4 increment optimization available to reduce nQx init costs by 33%
	stop_reason = fd_init (&pm1data, 1, 2, pm1data.x);			// This FFTs x as a side effect
	if (stop_reason) goto exit;
	for (totrels = 0, i = 1; ; i += 2) {		/* Compute x^(i^e) */
		if (relatively_prime (i, pm1data.D)) {
			pm1data.nQx[totrels] = gwalloc (&pm1data.gwdata);
			if (pm1data.nQx[totrels] == NULL) goto lowmem;
			gwfft (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[totrels]);
			totrels++;
			if (totrels == pm1data.totrels) break;
		}
		fd_next (&pm1data);
		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;
		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			fd_term (&pm1data);
			pm1_save (&pm1data);
			goto exit;
		}
	}
	fd_term (&pm1data);

/* Initialize for computing successive x^(m^e) where m corresponds to the first D section we are working on. */

//GW: Should we have nQx setup compute x^D (like ECM code does)?
	fd_init (&pm1data, pm1data.B2_start + pm1data.Dsection * pm1data.D + pm1data.D / 2, pm1data.D, pm1data.x);

/* Stage 2 init complete, change the title, output a message */

	sprintf (buf, "%.*f%% of %s P-1 stage 2 (using %uMB)",
		 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata), memused);
	title (thread_num, buf);

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 init complete. %.0f transforms. Time: ", gw_get_fft_count (&pm1data.gwdata));
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	gw_clear_fft_count (&pm1data.gwdata);

/* Since E >= 2, we do prime pairing with each loop iteration handling the range m-Q to m+Q where m is a multiple of D and Q is the */
/* Q-th relative prime to D.  Totrels is often much larger than the number of relative primes less than D.  This Preda/Montgomery */
/* optimization provides us with many more chances to find a prime pairing. */

	pm1data.state = PM1_STATE_STAGE2;
	start_timer (timers, 0);
//gw: ecm does not use timer1.  Change ECM to use timer 1 and P-1 to message "total time"
	start_timer (timers, 1);
	last_output = last_output_t = last_output_r = 0;
	for ( ; ; ) {
		uint64_t bitnum;

/* Test which relprimes in this D section need to be processed */

		bitnum = (pm1data.Dsection - pm1data.bitarrayfirstDsection) * pm1data.totrels;
		for (i = 0; i < pm1data.totrels; i++) {

/* Skip this relprime if the pairing routine did not set the corresponding bit in the bitarray. */

			if (! bittst (pm1data.bitarray, bitnum + i)) continue;

/* Multiply this eQx - nQx value into the gg accumulator */

			saving = testSaveFilesFlag (thread_num);
			stop_reason = stopCheck (thread_num);
			gwsubmul4 (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[i], pm1data.gg, pm1data.gg, (!stop_reason && !saving) ? GWMUL_STARTNEXTFFT : 0);

/* Clear pairing bit from the bitarray in case a save file is written.  Calculate stage 2 percentage. */

			bitclr (pm1data.bitarray, bitnum + i);
			w->pct_complete = base_pct_complete + (double) (pm1data.Dsection * pm1data.totrels + i) * one_bit_pct;

/* Test for errors */

			if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;

/* Output the title every so often */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
				sprintf (buf, "%.*f%% of %s P-1 stage 2 (using %uMB)",
					 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata), memused);
				title (thread_num, buf);
				last_output_t = gw_get_fft_count (&pm1data.gwdata);
			}

/* Write out a message every now and then */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
				sprintf (buf, "%s stage 2 is %.*f%% complete.",
					 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
				end_timer (timers, 0);
				if (first_iter_msg) {
					strcat (buf, "\n");
					clear_timer (timers, 0);
				} else {
					strcat (buf, " Time: ");
					print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
				}
				OutputStr (thread_num, buf);
				start_timer (timers, 0);
				last_output = gw_get_fft_count (&pm1data.gwdata);
				first_iter_msg = FALSE;
			}

/* Write out a message to the results file every now and then */

			if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
			    (NO_GUI && stop_reason)) {
				sprintf (buf, "%s stage 2 is %.*f%% complete.\n",
					 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
				writeResults (buf);
				last_output_r = gw_get_fft_count (&pm1data.gwdata);
			}

/* Periodicly write a save file.  If we escaped, free eQx memory so that pm1_save can reuse it to convert x and gg to binary. */

			if (stop_reason || saving) {
				if (stop_reason) fd_term (&pm1data);
				pm1_save (&pm1data);
				if (stop_reason) goto exit;
			}
		}

/* Move onto the next D section when we are done with all the relprimes */

		pm1data.Dsection++;
		pm1data.C_done = pm1data.B2_start + pm1data.Dsection * pm1data.D;
		if (pm1data.Dsection >= pm1data.numDsections) break;
		fd_next (&pm1data);
//GW: check for errors?

/* See if more bitarrays are required to get us to bound #2 */

		if (pm1data.Dsection < pm1data.bitarrayfirstDsection + pm1data.bitarraymaxDsections) continue;

		pm1data.first_relocatable = calc_new_first_relocatable (pm1data.D, pm1data.C_done);
		pm1data.bitarrayfirstDsection = pm1data.Dsection;
		stop_reason = fill_work_bitarray (pm1data.thread_num, &pm1data.sieve_info, pm1data.D, pm1data.totrels,
						  pm1data.first_relocatable, pm1data.last_relocatable, pm1data.C_done, pm1data.C,
						  pm1data.bitarraymaxDsections, &pm1data.bitarray);
		if (stop_reason) {
			pm1_save (&pm1data);
			goto exit;
		}
	}
	pm1data.C_done = pm1data.interim_C;
	fd_term (&pm1data);

/* Free up memory */

	for (i = 0; i < pm1data.totrels; i++) gwfree (&pm1data.gwdata, pm1data.nQx[i]);
	free (pm1data.nQx), pm1data.nQx = NULL;
	free (pm1data.eQx), pm1data.eQx = NULL;
	free (pm1data.bitarray), pm1data.bitarray = NULL;

/* Check for the rare cases where we need to do even more stage 2.  This happens when continuing a save file in the middle of stage 2 and */
/* the save file's target bound #2 was smaller than our target bound #2. */

	if (pm1data.C > pm1data.C_done) {
		pm1data.state = PM1_STATE_MIDSTAGE;
		pm1data.first_C_start = pm1data.C_done;
		goto more_C;
	}

/* Stage 2 is complete */

	end_timer (timers, 1);
	sprintf (buf, "%s stage 2 complete. %.0f transforms. Time: ", gwmodulo_as_string (&pm1data.gwdata), gw_get_fft_count (&pm1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pm1data.gwdata);
	}

/* See if we got lucky! */

restart4:
	pm1data.state = PM1_STATE_GCD;
	strcpy (w->stage, "S2");
	w->pct_complete = 1.0;
	if (w->work_type != WORK_PMINUS1) OutputStr (thread_num, "Starting stage 2 GCD - please be patient.\n");
	start_timer (timers, 0);
	stop_reason = gcd (&pm1data.gwdata, thread_num, pm1data.gg, N, &factor);
	if (stop_reason) {
		pm1_save (&pm1data);
		goto exit;
	}
	pm1data.state = PM1_STATE_DONE;
	end_timer (timers, 0);
	strcpy (buf, "Stage 2 GCD complete. Time: ");
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	if (factor != NULL) goto bingo;

/* Output line to results file indicating P-1 run */

msg_and_exit:
	sprintf (buf, "%s completed P-1, B1=%" PRIu64, gwmodulo_as_string (&pm1data.gwdata), pm1data.B);
	if (pm1data.C > pm1data.B) {
		if (pm1data.E <= 2)
			sprintf (buf+strlen(buf), ", B2=%" PRIu64, pm1data.C);
		else
			sprintf (buf+strlen(buf), ", B2=%" PRIu64 ", E=%d", pm1data.C, pm1data.E);
	}
	sprintf (buf+strlen(buf), ", Wi%d: %08lX\n", PORT, SEC5 (w->n, pm1data.B, pm1data.C));
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"NF", "exponent":45581713, "worktype":"P-1", "b1":50000, "b2":5000000, "brent-suyama":6, */
/* "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"NF\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P-1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64, pm1data.B);
	if (pm1data.C > pm1data.B) {
		sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%" PRIu64, pm1data.C);
		if (pm1data.E > 2) sprintf (JSONbuf+strlen(JSONbuf), ", \"brent-suyama\":%d", pm1data.E);
	}
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pm1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, pm1data.B, pm1data.C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send P-1 completed message to the server.  Although don't do it for puny B1 values as this is just the user tinkering with P-1 factoring. */

	if (!QA_IN_PROGRESS && (pm1data.B >= 50000 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type = PRIMENET_AR_P1_NOFACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.B1 = (double) pm1data.B;
		pkt.B2 = (double) pm1data.C;
		pkt.fftlen = gwfftlen (&pm1data.gwdata);
		pkt.done = (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR);
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* If this is pre-factoring for an LL or PRP test, then delete the large save file. */
/* Create save file so that we can expand bound 1 or bound 2 at a later date. */

	unlinkSaveFiles (&pm1data.write_save_file_state);
	if (!QA_IN_PROGRESS && w->work_type == WORK_PMINUS1 && IniGetInt (INI_FILE, "KeepPminus1SaveFiles", 1))
		pm1_save (&pm1data);

/* Return stop code indicating success or work unit complete */ 

done:	if (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR)
		stop_reason = STOP_WORK_UNIT_COMPLETE;
	else {
		w->pminus1ed = 1;		// Flag to indicate LL test has completed P-1
		w->tests_saved = 0.0;		// Variable to indicate PRP test has completed P-1
		stop_reason = updateWorkToDoLine (thread_num, w);
		if (stop_reason) goto exit;
	}

/* Free memory and return */

exit:	pm1_cleanup (&pm1data);
	free (N);
	free (factor);
	free (str);
	free (msg);
	if (exp_initialized) mpz_clear (exp);
	return (stop_reason);

/* Low on memory, reduce memory settings and try again */

lowmem:	fd_term (&pm1data);
	pm1_save (&pm1data);
	pm1_cleanup (&pm1data);
	OutputBoth (thread_num, "Memory allocation error.  Trying again using less memory.\n");
	pm1data.pct_mem_to_use *= 0.8;
	goto restart;

/* We've run out of memory.  Print error message and exit. */

oom:	stop_reason = OutOfMemory (thread_num);
	goto exit;

/* Print a message if we found a factor! */

bingo:	if (pm1data.state < PM1_STATE_MIDSTAGE)
		sprintf (buf, "P-1 found a factor in stage #1, B1=%" PRIu64 ".\n", pm1data.B);
	else if (pm1data.E <= 2)
		sprintf (buf, "P-1 found a factor in stage #2, B1=%" PRIu64 ", B2=%" PRIu64 ".\n", pm1data.B, pm1data.C);
	else
		sprintf (buf, "P-1 found a factor in stage #2, B1=%" PRIu64 ", B2=%" PRIu64 ", E=%d.\n", pm1data.B, pm1data.C, pm1data.E);
	OutputBoth (thread_num, buf);

/* Allocate memory for the string representation of the factor and for */
/* a message.  Convert the factor to a string.  Allocate lots of extra space */
/* as formatMsgForResultsFile can append a lot of text. */	

	msglen = factor->sign * 10 + 400;
	str = (char *) malloc (msglen);
	if (str == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	msg = (char *) malloc (msglen);
	if (msg == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	gtoc (factor, str, msglen);

/* Validate the factor we just found */

	if (!testFactor (&pm1data.gwdata, w, factor)) {
		sprintf (msg, "ERROR: Bad factor for %s found: %s\n", gwmodulo_as_string (&pm1data.gwdata), str);
		OutputBoth (thread_num, msg);
		unlinkSaveFiles (&pm1data.write_save_file_state);
		OutputStr (thread_num, "Restarting P-1 from scratch.\n");
		stop_reason = 0;
		goto error_restart;
	}

/* Output the validated factor */

	if (pm1data.state < PM1_STATE_MIDSTAGE)
		sprintf (msg, "%s has a factor: %s (P-1, B1=%" PRIu64 ")\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, pm1data.B);
	else if (pm1data.E <= 2)
		sprintf (msg, "%s has a factor: %s (P-1, B1=%" PRIu64 ", B2=%" PRIu64 ")\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, pm1data.B, pm1data.C);
	else
		sprintf (msg, "%s has a factor: %s (P-1, B1=%" PRIu64 ", B2=%" PRIu64 ", E=%d)\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, pm1data.B, pm1data.C, pm1data.E);
	OutputStr (thread_num, msg);
	formatMsgForResultsFile (msg, w);
	writeResults (msg);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"F", "exponent":45581713, "worktype":"P-1", "factors":["430639100587696027847"], */
/* "b1":50000, "b2":5000000, "brent-suyama":6, */
/* "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"F\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P-1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"factors\":[\"%s\"]", str);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64, pm1data.B);
	if (pm1data.state > PM1_STATE_MIDSTAGE) {
		sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%" PRIu64, pm1data.C);
		if (pm1data.E > 2) sprintf (JSONbuf+strlen(JSONbuf), ", \"brent-suyama\":%d", pm1data.E);
	}
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pm1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, pm1data.B, pm1data.C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send assignment result to the server.  To avoid flooding the server with small factors from users needlessly redoing */
/* factoring work, make sure the factor is more than 60 bits or so. */

	if (!QA_IN_PROGRESS && (strlen (str) >= 18 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, msg);
		pkt.result_type = PRIMENET_AR_P1_FACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		truncated_strcpy (pkt.factor, sizeof (pkt.factor), str);
		pkt.B1 = (double) pm1data.B;
		pkt.B2 = (double) (pm1data.state < PM1_STATE_MIDSTAGE ? 0 : pm1data.C);
		pkt.fftlen = gwfftlen (&pm1data.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Free save files.  If LL testing, free those save files too. */
/* Then create save file so that we can expand bound 1 or bound 2 at a later date. */

	unlinkSaveFiles (&pm1data.write_save_file_state);
	if (w->work_type != WORK_PMINUS1) {
		pm1data.write_save_file_state.base_filename[0] = 'p';
		unlinkSaveFiles (&pm1data.write_save_file_state);
	}
	if (!QA_IN_PROGRESS && w->work_type == WORK_PMINUS1 && IniGetInt (INI_FILE, "KeepPminus1SaveFiles", 1))
		pm1_save (&pm1data);

/* Since we found a factor, then we may have performed less work than */
/* expected.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	invalidateNextRollingAverageUpdate ();

/* Remove the exponent from the worktodo.ini file */

	stop_reason = STOP_WORK_UNIT_COMPLETE;
	goto exit;

/* Output an error message saying we are restarting. */
/* Sleep five minutes before restarting from last save file. */

err:	if (gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) {
		sprintf (buf, "Possible roundoff error (%.8g), backtracking to last save file.\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
	} else {
		OutputBoth (thread_num, "SUMOUT error occurred.\n");
		stop_reason = SleepFive (thread_num);
		if (stop_reason) goto exit;
	}
error_restart:
	pm1_cleanup (&pm1data);
	free (factor), factor = NULL;
	free (str), str = NULL;
	free (msg), msg = NULL;
	if (exp_initialized) mpz_clear (exp), exp_initialized = FALSE;
	goto restart;
}

/* Read a file of P-1 tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9992 to run the QA suite */

int pminus1_QA (
	int	thread_num,
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	FILE	*fd;

/* Set the title */

	title (thread_num, "QA");

/* Open QA file */

	fd = fopen ("qa_pm1", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa_pm1' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	QA_TYPE = 0;
	for ( ; ; ) {
		struct work_unit w;
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		int	stop_reason;

/* Read a line from the file */

		n = 0;
		(void) fscanf (fd, "%lf,%lu,%lu,%ld,%lu,%lu,%lu,%s\n", &k, &b, &n, &c, &B1, &B2_start, &B2_end, fac_str);
		if (n == 0) break;

/* If p is 1, set QA_TYPE */

		if (n == 1) {
			QA_TYPE = c;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		QA_FACTOR = allocgiant ((int) strlen (fac_str));
		ctog (fac_str, QA_FACTOR);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Do the P-1 */

		if (B2_start < B1) B2_start = B1;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_PMINUS1;
		w.k = k;
		w.b = b;
		w.n = n;
		w.c = c;
		w.B1 = B1;
		w.B2_start = B2_start;
		w.B2 = B2_end;
		QA_IN_PROGRESS = TRUE;
		stop_reason = pminus1 (0, sp_info, &w);
		QA_IN_PROGRESS = FALSE;
		free (QA_FACTOR);
		if (stop_reason != STOP_WORK_UNIT_COMPLETE) {
			fclose (fd);
			return (stop_reason);
		}
	}

/* Cleanup */

	fclose (fd);
	return (0);
}

/**************************************************************/
/* Routines to compute optimal and test to optimal P-1 bounds */
/**************************************************************/

struct global_pm1_cost_data {
	unsigned long n;				// Exponent being tested
	double	takeAwayBits;				// Bits we get for free in smoothness of P-1 factor
	double	how_far_factored;			// How far the number has been trial factored
	double	gcd_cost;				// Cost (in squarings) of running GCD
	double	ll_testing_cost;			// Cost (in squarings) of running LL/PRP tests should we fail to find a factor
	unsigned long vals;				// Number of temporaries we can allocate in pass 2
	int	E;					// Brent-Suyama
};

struct cost_pm1_data {
	unsigned long B1;
	unsigned long B2;
	double	prob;
	double	pass1_squarings;
	double	pass2_squarings;
	double	savings;
};

/* For a given B1,B2 calculate the costs and savings */

void cost_pm1 (
	struct global_pm1_cost_data *g,
	struct cost_pm1_data *c)
{
	struct pm1_stage2_cost_data cost_data;

/* Not sure this test is needed.  Handle no stage 2. */

	if (c->B2 <= c->B1) {
		c->pass2_squarings = 0.0;
	}

/* Compute how many squarings will be required in the best implementation of pass 2 given gg and (E+1) gwnum temporaries will be needed. */

	else {
		cost_data.E = g->E;
		c->pass2_squarings = best_stage2_impl (c->B1, c->B2, g->vals - (cost_data.E + 2), &pm1_stage2_cost, &cost_data);
	}

/* Calculate probability of finding a factor (courtesy of Mihai Preda) */

	c->prob = pm1prob (g->takeAwayBits, (unsigned int) g->how_far_factored, c->B1, c->B2);

/* Calculate our savings using this B1/B2.  Savings is success_probability * cost_of_LL_tests - cost_of_Pminus1. */

	c->savings = c->prob * g->ll_testing_cost - (c->pass1_squarings + c->pass2_squarings + g->gcd_cost);
}

/* For a given B1, find the best B2 */

void pminus1_choose_B2 (
	struct global_pm1_cost_data *g,
	struct cost_pm1_data *c)
{
	struct cost_pm1_data best[3];

/* Estimate how many squarings will be required in pass 1 */

	c->pass1_squarings = ceil (1.44 * c->B1);

/* Look for the best B2 somewhere between 1*B1 and 100*B1 */

	best[0] = *c;
	best[0].B2 = c->B1;
	cost_pm1 (g, &best[0]);
	best[1] = *c;
	best[1].B2 = 50*c->B1;
	cost_pm1 (g, &best[1]);
	best[2] = *c;
	best[2].B2 = 100*c->B1;
	cost_pm1 (g, &best[2]);

/* Handle case where midpoint has worse savings than the start point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[0].savings > best[1].savings) {
		best[2] = best[1];
		best[1].B2 = (best[0].B1 + best[2].B1) / 2;
		cost_pm1 (g, &best[1]);
	}

/* Handle case where midpoint has worse savings than the end point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[2].savings > best[1].savings) {
		best[0] = best[1];
		best[1] = best[2];
		best[2].B2 = best[1].B2 * 2;
		cost_pm1 (g, &best[2]);
	}

/* Find the best B2.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].B2 != best[2].B2) {
		struct cost_pm1_data midpoint;

		ASSERTG (best[1].savings >= best[0].savings);
		ASSERTG (best[1].savings >= best[2].savings);

/* Work on the bigger of the lower section and upper section */

		if (best[1].B2 - best[0].B2 > best[2].B2 - best[1].B2) {	// Work on lower section
			// If B2's are close together or the savings difference is real small, then we've searched this section enough
			if (best[1].B2 - best[0].B2 < 1000 || best[1].savings - best[0].savings < 100.0) {
				best[0] = best[1];
				continue;
			}
			midpoint = *c;
			midpoint.B2 = (best[0].B2 + best[1].B2) / 2;
			cost_pm1 (g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			// If B2's are close together or the savings difference is real small, then we've searched this section enough
			if (best[2].B2 - best[1].B2 < 1000 || best[1].savings - best[2].savings < 100.0) {
				best[2] = best[1];
				continue;
			}
			midpoint = *c;
			midpoint.B2 = (best[1].B2 + best[2].B2) / 2;
			cost_pm1 (g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Return the best B2 we found */

	c->B2 = best[1].B2;
	c->prob = best[1].prob;
	c->pass2_squarings = best[1].pass2_squarings;
	c->savings = best[1].savings;
}

/* Calculate the best B1 and B2 values to use in a P-1 factoring job. */
/* Return the B1 and B2 bounds, execution cost, and chance of success. */

void guess_pminus1_bounds (
	int	thread_num,
	double	k,			/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,		/* B in K*B^N+C. Must be two. */
	unsigned long n,		/* N in K*B^N+C. Exponent to test. */
	signed long c,			/* C in K*B^N+C. */
	double	how_far_factored,	/* Bit depth of trial factoring */
	double	tests_saved,		/* 1 if doublecheck, 2 if first test */
	unsigned long *bound1,
	unsigned long *bound2,
	unsigned long *squarings,
	double	*success_rate)
{
	struct global_pm1_cost_data g;
	struct cost_pm1_data best[3];

/* Copy exponent, how_far_factored to global costing data */

	g.n = n;
	g.takeAwayBits = isMersenne (k, b, n, c) ? log2 (n) + 1.0 : isGeneralizedFermat (k, b, n, c) ? log2 (n) : 0.0;
	g.how_far_factored = how_far_factored;

/* Guard against wild tests_saved values.  Huge values will cause this routine */
/* to run for a very long time.  This shouldn't happen as auxiliaryWorkUnitInit */
/* now has the exact same test. */

	if (tests_saved > 10) tests_saved = 10;

/* Balance P-1 against 1 or 2 LL/PRP tests (actually more since we get a */
/* corrupt result reported some of the time). */

	g.ll_testing_cost = (tests_saved + 2 * ERROR_RATE) * n;

/* The GCD cost comes from the timing code running ECM on M604 and a spreadsheet. */
/* Since GCDs are single-threaded we double the GCD cost for multi-threaded P-1 runs. */

	g.gcd_cost = 160.265 * log ((double) n) - 1651.0;
	if (g.gcd_cost < 50.0) g.gcd_cost = 50.0;
	if (CORES_PER_TEST[thread_num] > 1) g.gcd_cost *= 2.0;

/* Compute how many temporaries we can use given our memory constraints. */
/* Allow 1MB for code and data structures. */

	g.vals = cvt_mem_to_estimated_gwnums (max_mem (thread_num), k, b, n, c);
	if (g.vals < 1) g.vals = 1;

/* Let the user force a Brent-Suyama stage 2 */

	g.E = IniGetInt (INI_FILE, "BrentSuyama", 2);
	if (g.E < 2) g.E = 2;
	if (g.E > 24) g.E = 24;

/* Find the best B1 somewhere between n/3300 and 250*(n/3300). */

	best[0].B1 = n / 3300;
	pminus1_choose_B2 (&g, &best[0]);
	best[1].B1 = 125 * best[0].B1;
	pminus1_choose_B2 (&g, &best[1]);
	best[2].B1 = 250 * best[0].B1;
	pminus1_choose_B2 (&g, &best[2]);

/* Handle case where midpoint has worse savings than the start point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[0].savings > best[1].savings) {
		best[2] = best[1];
		best[1].B1 = (best[0].B1 + best[2].B1) / 2;
		pminus1_choose_B2 (&g, &best[1]);
	}

/* Handle case where midpoint has worse savings than the end point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (best[2].savings > best[1].savings) {
		best[0] = best[1];
		best[1] = best[2];
		best[2].B1 = best[1].B1 * 2;
		pminus1_choose_B2 (&g, &best[2]);
	}

/* Find the best B1.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].B1 != best[2].B1) {
		struct cost_pm1_data midpoint;

		ASSERTG (best[1].savings >= best[0].savings);
		ASSERTG (best[1].savings >= best[2].savings);

/* Work on the bigger of the lower section and upper section */

		if (best[1].B1 - best[0].B1 > best[2].B1 - best[1].B1) {	// Work on lower section
			// If B1's are close together or the savings difference is real small, then we've searched this section enough
			if (best[1].B1 - best[0].B1 < 1000 || best[1].savings - best[0].savings < 100.0) {
				best[0] = best[1];
				continue;
			}
			midpoint.B1 = (best[0].B1 + best[1].B1) / 2;
			pminus1_choose_B2 (&g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			// If B1's are close together or the savings difference is real small, then we've searched this section enough
			if (best[2].B1 - best[1].B1 < 1000 || best[1].savings - best[2].savings < 100.0) {
				best[2] = best[1];
				continue;
			}
			midpoint.B1 = (best[1].B1 + best[2].B1) / 2;
			pminus1_choose_B2 (&g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Round up B1 and B2 to nearest 1000 -- just to look pretty */

	best[1].B1 = round_up_to_multiple_of (best[1].B1, 1000);
	best[1].B2 = round_up_to_multiple_of (best[1].B2, 1000);
	cost_pm1 (&g, &best[1]);

/* Return the final best choice */

	if (best[1].savings > 0.0) {
		*bound1 = best[1].B1;
		*bound2 = best[1].B2;
		*squarings = (unsigned long) (best[1].pass1_squarings + best[1].pass2_squarings + g.gcd_cost);
		*success_rate = best[1].prob;
	} else {
		*bound1 = 0;
		*bound2 = 0;
		*squarings = 0;
		*success_rate = 0.0;
	}
}

/* Determine the probability of P-1 finding a factor.  Return the Mihai estimated P-1 success probability */

double guess_pminus1_probability (
	struct work_unit *w)
{
	double	takeAwayBits;
	takeAwayBits = isMersenne (w->k, w->b, w->n, w->c) ? log2 (w->n) + 1.0 :
		       isGeneralizedFermat (w->k, w->b, w->n, w->c) ? log2 (w->n) : 0.0;
	return (pm1prob (takeAwayBits, (unsigned int) w->sieve_depth, w->B1, w->B2));
}

/* Do the P-1 factoring step prior to a Lucas-Lehmer test */
/* Similar to the main P-1 entry point, except bounds are not known */

int pfactor (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	unsigned long bound1, bound2, squarings;
	double	prob;
	char	buf[120], testnum[120];
	int	stop_reason;

/* Choose the best FFT size */

	stop_reason = pick_fft_size (thread_num, w);
	if (stop_reason) return (stop_reason);

/* Set flag indicating we need to restart if the maximum amount of memory changes (as opposed to currently available memory!) */
/* If maximum memory changes we want to recompute the P-1 bounds. */

	set_restart_if_max_memory_change (thread_num);

/* Output a message that P-1 factoring is about to begin */

	gw_as_string (testnum, w->k, w->b, w->n, w->c);
	sprintf (buf, "Optimal P-1 factoring of %s using up to %luMB of memory.\n", testnum, max_mem (thread_num));
	OutputStr (thread_num, buf);
	sprintf (buf, "Assuming no factors below 2^%.2g and %.2g primality test%s saved if a factor is found.\n",
		 w->sieve_depth, w->tests_saved, w->tests_saved == 1.0 ? "" : "s");
	OutputStr (thread_num, buf);

/* Deduce the proper P-1 bounds */

	guess_pminus1_bounds (thread_num, w->k, w->b, w->n, w->c, w->sieve_depth, w->tests_saved, &bound1, &bound2, &squarings, &prob);
	if (bound1 == 0) {
		sprintf (buf, "%s does not need P-1 factoring.\n", testnum);
		OutputBoth (thread_num, buf);
		if (w->work_type == WORK_PFACTOR) {
			//bug - do we need to tell the server to cancel the
			//assignment?  In theory, server shouldn't ever send
			//this assignment out.
			return (STOP_WORK_UNIT_COMPLETE);
		} else {
			w->pminus1ed = 1;		// Flag to indicate LL test has completed P-1
			w->tests_saved = 0.0;		// Variable to indicate PRP test has completed P-1
			stop_reason = updateWorkToDoLine (thread_num, w);
			if (stop_reason) return (stop_reason);
			return (0);
		}
	}

/* Output a message that P-1 factoring is about to begin */

	sprintf (buf, "Optimal bounds are B1=%lu, B2=%lu\n", bound1, bound2);
	OutputStr (thread_num, buf);
	sprintf (buf, "Chance of finding a factor is an estimated %.3g%%\n", prob * 100.0);
	OutputStr (thread_num, buf);

/* Call the P-1 factoring code */

	w->B1 = bound1;
	w->B2_start = bound1;
	w->B2 = bound2;
	return (pminus1 (thread_num, sp_info, w));
}

/**************************************************************
 *	P+1 Functions
 **************************************************************/

/* Data maintained during P+1 process */

#define PP1_STATE_STAGE1	1	/* In stage 1 */
#define PP1_STATE_MIDSTAGE	2	/* Between stage 1 and stage 2 */
#define PP1_STATE_STAGE2	3	/* In middle of stage 2 (processing a bit array) */
#define PP1_STATE_GCD		4	/* Stage 2 GCD */
#define PP1_STATE_DONE		5	/* P+1 job complete */

typedef struct {
	gwhandle gwdata;	/* GWNUM handle */
	int	thread_num;	/* Worker thread number */
	struct work_unit *w;	/* Worktodo.txt entry */
	int	state;		/* One of the states listed above */
	double	takeAwayBits;	/* Bits we get for free in smoothness of P+1 factor */
	double	success_rate;	/* Percentage of smooth P+1 factors we expect to find.  Percentage goes down with each P+1 run. */
	uint32_t numerator;	/* Starting point numerator */
	uint32_t denominator;	/* Starting point denominator */
	uint64_t B;		/* Bound #1 (a.k.a. B1) */
	uint64_t C;		/* Bound #2 (a.k.a. B2 */
	uint64_t interim_B;	/* B1 we are currently calculating (equals B except when finishing stage 1 from a save file using a different B1). */
	uint64_t B_done;	/* We have completed calculating 3^e to this bound #1 */
	int	optimal_B2;	/* TRUE if we calculate optimal bound #2 given currently available memory.  FALSE for a fixed bound #2. */

	readSaveFileState read_save_file_state;	/* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	void	*sieve_info;	/* Prime number sieve */
	uint64_t stage1_prime;	/* Prime number being processed */

	int	D;		/* Stage 2 loop size */
	int	totrels;	/* Number relatively prime nQx values used */
	uint64_t B2_start;	/* Starting point of first D section to be processed in stage 2 (an odd multiple of D/2) */
	uint64_t numDsections;	/* Number of D sections to process in stage 2 */
	uint64_t Dsection;	/* Current D section being processed in stage 2 */
	uint64_t interim_C;	/* B2 we are currently calculating (equals C except when finishing stage 2 a save file using different B2) */

	char	*bitarray;	/* Bit array for prime pairings in each D section */
	uint64_t bitarraymaxDsections;	/* Maximum number of D sections per bit array */
	uint64_t bitarrayfirstDsection; /* First D section in the bit array */
	uint64_t first_relocatable; /* First relocatable prime (same as B1 unless bit arrays must be split or mem change caused a replan) */
	uint64_t last_relocatable; /* Last relocatable prime for filling bit arrays (unless mem change causes a replan) */
	uint64_t C_done;	/* Stage 2 completed thusfar (updated every D section that is completed) */

	int	stage2_numvals;	/* Number of gwnums used in stage 2 */
	gwnum	*nQx;		/* Array of relprime data used in stage 2 */

	double	pct_mem_to_use;	/* If we get memory allocation errors, we progressively try using less and less. */

	gwnum	V;		/* V_1 in a Lucas sequence */
	gwnum	Vn;		/* V_n in a Lucas sequence */
	gwnum	Vn1;		/* V_{n+1} in a Lucas sequence */
	gwnum	gg;		/* An accumulator in stage 2 */
} pp1handle;

/* Perform cleanup functions. */

void pp1_cleanup (
	pp1handle *pp1data)
{

/* Free memory */

	free (pp1data->nQx), pp1data->nQx = NULL;
	free (pp1data->bitarray), pp1data->bitarray = NULL;
	gwdone (&pp1data->gwdata);
	end_sieve (pp1data->sieve_info), pp1data->sieve_info = NULL;
}

/* Routines to create and read save files for a P+1 factoring job */

#define PP1_MAGICNUM	0x912a374a
#define PP1_VERSION	1				/* New in 30.6 */

void pp1_save (
	pp1handle *pp1data)
{
	int	fd;
	struct work_unit *w = pp1data->w;
	unsigned long sum = 0;

/* Create the intermediate file */

	fd = openWriteSaveFile (&pp1data->write_save_file_state);
	if (fd < 0) return;

/* Write the file header */

	if (!write_header (fd, PP1_MAGICNUM, PP1_VERSION, w)) goto writeerr;

/* Write the file data */

	if (! write_int (fd, pp1data->state, &sum)) goto writeerr;
	if (! write_int (fd, pp1data->numerator, &sum)) goto writeerr;
	if (! write_int (fd, pp1data->denominator, &sum)) goto writeerr;

	if (pp1data->state == PP1_STATE_STAGE1) {
		if (! write_longlong (fd, pp1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->interim_B, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->stage1_prime, &sum)) goto writeerr;
	}

	else if (pp1data->state == PP1_STATE_MIDSTAGE) {
		if (! write_longlong (fd, pp1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->C_done, &sum)) goto writeerr;
	}

	// Save everything necessary to restart stage 2 without calling pp1_stage2_impl again
	else if (pp1data->state == PP1_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! write_longlong (fd, pp1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->C_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->interim_C, &sum)) goto writeerr;
		if (! write_int (fd, pp1data->stage2_numvals, &sum)) goto writeerr;
		if (! write_int (fd, pp1data->totrels, &sum)) goto writeerr;
		if (! write_int (fd, pp1data->D, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->first_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->last_relocatable, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->B2_start, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->numDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->Dsection, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->bitarraymaxDsections, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->bitarrayfirstDsection, &sum)) goto writeerr;
		// Output the truncated bit array
		bitarray_numDsections = pp1data->numDsections - pp1data->bitarrayfirstDsection;
		if (bitarray_numDsections > pp1data->bitarraymaxDsections) bitarray_numDsections = pp1data->bitarraymaxDsections; 		
		bitarray_len = divide_rounding_up (bitarray_numDsections * pp1data->totrels, 8);
		bitarray_start = divide_rounding_down ((pp1data->Dsection - pp1data->bitarrayfirstDsection) * pp1data->totrels, 8);
		if (! write_array (fd, pp1data->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto writeerr;
	}

	else if (pp1data->state == PP1_STATE_GCD) {
		if (! write_longlong (fd, pp1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->C_done, &sum)) goto writeerr;
	}

	else if (pp1data->state == PP1_STATE_DONE) {
		if (! write_longlong (fd, pp1data->B_done, &sum)) goto writeerr;
		if (! write_longlong (fd, pp1data->C_done, &sum)) goto writeerr;
	}

/* Write the gwnum value used in stage 1 and stage 2.  There are occasions where x or gg may be in a partially FFTed state. */

	if (! write_gwnum (fd, &pp1data->gwdata, pp1data->V, &sum)) goto writeerr;
	if (pp1data->state >= PP1_STATE_MIDSTAGE && pp1data->state <= PP1_STATE_GCD) {
		if (! write_gwnum (fd, &pp1data->gwdata, pp1data->gg, &sum)) goto writeerr;
	}

/* Write the checksum, we're done */

	if (! write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (&pp1data->write_save_file_state, fd);
	return;

/* An error occurred.  Close and delete the current file. */

writeerr:
	deleteWriteSaveFile (&pp1data->write_save_file_state, fd);
}


/* Read a save file */

int pp1_restore (
	pp1handle *pp1data)
{
	int	fd, numerator, denominator;
	struct work_unit *w = pp1data->w;
	unsigned long version;
	unsigned long sum = 0, filesum;

/* Open the intermediate file */

	fd = _open (pp1data->read_save_file_state.current_filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto err;

/* Read the file header */

	if (! read_magicnum (fd, PP1_MAGICNUM)) goto readerr;
	if (! read_header (fd, &version, w, &filesum)) goto readerr;
	if (version < 1 || version > PP1_VERSION) goto readerr;

/* Read the first part of the save file */

	if (! read_int (fd, &pp1data->state, &sum)) goto readerr;
	if (! read_int (fd, &numerator, &sum)) goto readerr;
	if (! read_int (fd, &denominator, &sum)) goto readerr;
	if (numerator != pp1data->numerator || denominator != pp1data->denominator) {
		if (pp1data->w->nth_run <= 2) goto bad_nth_run;			// User wants to do 2/7 or 6/5 and save file does not match
		if (numerator == 2 && denominator == 7) goto bad_nth_run;	// User wants a random start, not 2/7
		if (numerator == 6 && denominator == 5) goto bad_nth_run;	// User wants a random start, not 6/5
		pp1data->numerator = numerator;					// Replace random start with random start from save file
		pp1data->denominator = denominator;
	}

/* Read state dependent data */

	if (pp1data->state == PP1_STATE_STAGE1) {
		if (! read_longlong (fd, &pp1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->interim_B, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->stage1_prime, &sum)) goto readerr;
	}

	else if (pp1data->state == PP1_STATE_MIDSTAGE) {
		if (! read_longlong (fd, &pp1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->C_done, &sum)) goto readerr;
	}

	else if (pp1data->state == PP1_STATE_STAGE2) {
		uint64_t bitarray_numDsections, bitarray_start, bitarray_len;
		if (! read_longlong (fd, &pp1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->C_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->interim_C, &sum)) goto readerr;
		if (! read_int (fd, &pp1data->stage2_numvals, &sum)) goto readerr;
		if (! read_int (fd, &pp1data->totrels, &sum)) goto readerr;
		if (! read_int (fd, &pp1data->D, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->first_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->last_relocatable, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->B2_start, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->numDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->Dsection, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->bitarraymaxDsections, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->bitarrayfirstDsection, &sum)) goto readerr;
		// Read the truncated bit array
		bitarray_numDsections = pp1data->numDsections - pp1data->bitarrayfirstDsection;
		if (bitarray_numDsections > pp1data->bitarraymaxDsections) bitarray_numDsections = pp1data->bitarraymaxDsections;
		bitarray_len = divide_rounding_up (bitarray_numDsections * pp1data->totrels, 8);
		bitarray_start = divide_rounding_down ((pp1data->Dsection - pp1data->bitarrayfirstDsection) * pp1data->totrels, 8);
		pp1data->bitarray = (char *) malloc ((size_t) bitarray_len);
		if (pp1data->bitarray == NULL) goto readerr;
		if (! read_array (fd, pp1data->bitarray + bitarray_start, (unsigned long) (bitarray_len - bitarray_start), &sum)) goto readerr;
	}

	else if (pp1data->state == PP1_STATE_GCD) {
		if (! read_longlong (fd, &pp1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->C_done, &sum)) goto readerr;
	}

	else if (pp1data->state == PP1_STATE_DONE) {
		if (! read_longlong (fd, &pp1data->B_done, &sum)) goto readerr;
		if (! read_longlong (fd, &pp1data->C_done, &sum)) goto readerr;
	}

/* Read the gwnum value used in stage 1 and stage 2 */

	if (! read_gwnum (fd, &pp1data->gwdata, pp1data->V, &sum)) goto readerr;

/* Read stage 2 accumulator gwnum */

	if (pp1data->state >= PP1_STATE_MIDSTAGE && pp1data->state <= PP1_STATE_GCD) {
		pp1data->gg = gwalloc (&pp1data->gwdata);
		if (pp1data->gg == NULL) goto readerr;
		if (! read_gwnum (fd, &pp1data->gwdata, pp1data->gg, &sum)) goto readerr;
	}

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);

/* All done */

	return (TRUE);

/* An error occurred.  Cleanup and return. */

bad_nth_run:
	OutputBoth (pp1data->thread_num, "P+1 starting point in save file does not match nth_run parameter from worktodo.txt\n");
readerr:
	_close (fd);
err:
	return (FALSE);
}


/* Compute the cost (in squarings) of a particular P+1 stage 2 implementation. */

struct pp1_stage2_cost_data {
				/* Data returned from cost function follows */
	int	D;		/* D value for big steps */
	int	totrels;	/* Total number of relative primes used for pairing */
	int	stage2_numvals;	/* Total number of gwnum temporaries to be used in stage 2 */
	uint64_t B2_start;	/* Stage 2 start */
	uint64_t numDsections;	/* Number of D sections to process */
	uint64_t bitarraymaxDsections; /* Maximum number of D sections per bit array (in case work bit array must be split) */
};

double pp1_stage2_cost (
	int	D,		/* Stage 2 big step */
	uint64_t B2_start,	/* Stage 2 start point */
	uint64_t numDsections,	/* Number of D sections to process */
	uint64_t bitarraymaxDsections, /* Maximum number of D sections per bit array (if work bit array must be split) */
	int	numrels,	/* Number of relative primes less than D */
	double	multiplier,	/* Multiplier * numrels relative primes will be used in stage 2 */
	double	numpairs,	/* Estimated number of paired primes in stage 2 */
	double	numsingles,	/* Estimated number of prime singles in stage 2 */
	void	*data)		/* P+1 specific costing data */
{
	int	totrels;	/* Total number of relative primes in use */
	struct pp1_stage2_cost_data *cost_data = (struct pp1_stage2_cost_data *) data;
	double	cost;

/* Calculate the number of temporaries used for relative primes.  This will leave some available for modinv pooling. */

	totrels = (int) (numrels * multiplier);

/* Compute the stage 2 setup costs.  They are tiny and if inaccurately estimated can be safely ignored. */

/* For the minimum number of relative primes for a D, there are D/6 luc_add calls at a cost of 1 multiply each call. */
/* If using more than the minimum number of relative primes, increase the cost proportionally. */
/* Also, each nQx value will be FFTed once but this happens as a by-product of the D / 4 luc_add calls. */

	cost = multiplier * ((double) D / 6) * (1.0) + totrels * 0.0 +

/* Compute the Vn, Vn1 setup costs in pp1_luc_prep_for_additions (2 squarings per bit) */

	       _log2 (B2_start / D) * 2.0;

/* Add in the cost of luc_add calls (#sections * 1 multiply) */

	cost += (double) (numDsections * 1);

/* Finally, each prime pair and prime single costs one multiply */

	cost += numpairs + numsingles;

/* Return data P+1 implementation will need */
/* Stage 2 memory is totrels gwnums for the nQx array, 3 gwnums to calculate multiples of D values, one gwnum for gg. */

	cost_data->stage2_numvals = totrels + 4;
	cost_data->D = D;
	cost_data->totrels = (int) (numrels * multiplier);
	cost_data->B2_start = B2_start;
	cost_data->numDsections = numDsections;
	cost_data->bitarraymaxDsections = bitarraymaxDsections;

/* Return the resulting cost */

	return (cost);
}

/* Choose the most effective B2 for a P+1 run with a fixed B1 given the number of gwnums we are allowed to allocate. */
/* That is, find the B2 such that investing a fixed cost in either a larger B1 or B2 results in the same increase in chance of finding a factor. */

void pp1_choose_B2 (
	pp1handle *pp1data,
	unsigned long numvals)
{
	struct pp1_stage2_cost_data cost_data;	/* Extra data passed to P+1 costing function */
	struct pp1_stage2_efficiency {
		int	i;
		double	B2_cost;		/* Cost of stage 2 in squarings */
		double	fac_pct;		/* Percentage chance of finding a factor */
	} best[3];
	int	sieve_depth;
	char	buf[160];

// P+1 probability of success
#define pp1prob(B1,B2)		(pm1prob (pp1data->takeAwayBits, sieve_depth, (double) (B1), (double) (B2)) * pp1data->success_rate)
// Cost out a B2 value
#define p1eval(x,B2mult)	x.i = B2mult; \
				x.B2_cost = best_stage2_impl (pp1data->B, x.i * pp1data->B, numvals - 4, &pp1_stage2_cost, &cost_data); \
				x.fac_pct = pp1prob (pp1data->B, x.i * pp1data->B);

// Return TRUE if x is better than y.  Determined by seeing if taking the increased cost of y's higher B2 and investing it in increasing x's bounds
// results in a higher chance of finding a factor.
#define B1increase(x,y)	(int) ((y.B2_cost - x.B2_cost) / (1.44*1.52))	// Each B1 increase 1.44x more bits to process at a cost of 1.52 squarings per bit
#define p1compare(x,y)	(pp1prob (pp1data->B + B1increase(x,y), x.i * pp1data->B + B1increase(x,y)) > y.fac_pct)

/* Look for the best B2 which is likely between 10*B1 and 80*B1.  If optimal is not between these bounds, don't worry we'll locate the optimal spot anyway. */

	sieve_depth = (int) pp1data->w->sieve_depth;
	p1eval (best[0], 10);
	p1eval (best[1], 40);
	p1eval (best[2], 80);

/* Handle case where midpoint is worse than the start point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (p1compare (best[0], best[1])) {
		best[2] = best[1];
		p1eval (best[1], (best[0].i + best[2].i) / 2);
	}

/* Handle case where midpoint is worse than the end point */
/* The search code requires best[1] is better than best[0] and best[2] */

	while (!p1compare (best[1], best[2])) {
		best[0] = best[1];
		best[1] = best[2];
		p1eval (best[2], (int) (best[1].i + 20));
	}

/* Find the best B2.  We use a binary-like search to speed things up. */

	while (best[0].i + 2 != best[2].i) {
		struct pp1_stage2_efficiency midpoint;

		// Work on the bigger of the lower section and upper section
		if (best[1].i - best[0].i > best[2].i - best[1].i) {		// Work on lower section
			p1eval (midpoint, (best[0].i + best[1].i) / 2);
			if (p1compare (midpoint, best[1])) {			// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			p1eval (midpoint, (best[1].i + best[2].i) / 2);
			if (!p1compare (best[1], midpoint)) {			// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Return the best B2 */

	pp1data->C = best[1].i * pp1data->B;
	sprintf (buf, "With trial factoring done to 2^%d, optimal B2 is %d*B1 = %" PRIu64 ".\n", sieve_depth, best[1].i, pp1data->C);
	OutputStr (pp1data->thread_num, buf);
	sprintf (buf, "Chance of a new factor assuming no ECM has been done is %.3g%%\n", best[1].fac_pct * 100.0);
	OutputStr (pp1data->thread_num, buf);
}
#undef p1eval
#undef B1increase
#undef p1compare

/* Choose the best implementation of P+1 stage 2 given the current memory settings.  We may decide there will never be enough memory. */
/* We may decide to wait for more memory to be available. */
/* We choose the best value for D that reduces the number of multiplications with the current memory constraints. */

int pp1_stage2_impl (
	pp1handle *pp1data)
{
	unsigned int memory, min_memory, desired_memory;	/* Memory is in MB */
	unsigned int bitarray_size, max_bitarray_size;		/* Bitarray memory in MB */
	int	numvals;			/* Number of gwnums we can allocate */
	struct pp1_stage2_cost_data cost_data;
	int	stop_reason;

/* Calculate (roughly) the memory required for the bit-array.  A typical stage 2 implementation is D=210, B2_start=B2/11, numrels=24, totrels=7*24. */
/* A typical optimal B2 is 40*B1.  This means bitarray memory = B2 * 10/11 / 210 * 24 * 7 / 8 bytes. */

	bitarray_size = divide_rounding_up ((unsigned int) ((pp1data->optimal_B2 ? 130 * pp1data->B : pp1data->C) * 1680 / 18480), 1 << 20);
	max_bitarray_size = IniGetInt (INI_FILE, "MaximumBitArraySize", 250);
	if (bitarray_size > max_bitarray_size) bitarray_size = max_bitarray_size;

/* Calculate the number of temporaries we can use in stage 2.  Base our decision on the current amount of memory available. */
/* If continuing from a stage 2 save file then we desire as many temporaries as the save file used.  Otherwise, assume 28 temporaries will */
/* provide us with a reasonable execution speed (D = 210).  We must have a minimum of 8 temporaries (D = 30). */

	min_memory = bitarray_size + cvt_gwnums_to_mem (&pp1data->gwdata, 8);
	desired_memory = bitarray_size + cvt_gwnums_to_mem (&pp1data->gwdata, pp1data->state >= PP1_STATE_STAGE2 ? pp1data->stage2_numvals : 28);
	stop_reason = avail_mem (pp1data->thread_num, min_memory, desired_memory, &memory);
	if (stop_reason) return (stop_reason);

/* Factor in the multiplier that we set to less than 1.0 when we get unexpected memory allocation errors. */
/* Make sure we can still allocate 8 temporaries. */

	memory = (unsigned int) (pp1data->pct_mem_to_use * (double) memory);
	if (memory < min_memory) return (avail_mem_not_sufficient (pp1data->thread_num, min_memory, desired_memory));
	if (memory < 8) memory = 8;

/* Output a message telling us how much memory is available */

	if (NUM_WORKER_THREADS > 1) {
		char	buf[100];
		sprintf (buf, "Available memory is %dMB.\n", memory);
		OutputStr (pp1data->thread_num, buf);
	}

/* Compute the number of gwnum temporaries we can allocate.  User nordi had over-allocating memory troubles on Linux testing M1277, presumably */
/* because our estimate genum size was too low.  As a work-around limit numvals to 100,000 by default. */

	numvals = cvt_mem_to_gwnums (&pp1data->gwdata, memory - bitarray_size);
	if (numvals < 8) numvals = 8;
	if (numvals > 100000) numvals = 100000;
	if (QA_TYPE) numvals = QA_TYPE;			/* Optionally override numvals for QA purposes */

/* Set C_done for future best_stage2_impl calls. */
/* Override B2 with optimal B2 based on amount of memory available. */

	if (pp1data->state == PP1_STATE_MIDSTAGE) {
		pp1data->C_done = pp1data->B;
		if (pp1data->optimal_B2) pp1_choose_B2 (pp1data, numvals);
	}

/* If are continuing from a save file that was in stage 2, check to see if we currently have enough memory to continue with the save file's */
/* stage 2 implementation.  Also check if we now have significantly more memory available and stage 2 is not near complete such that a new */
/* stage 2 implementation might give us a faster stage 2 completion. */

//GW: These are rather arbitrary heuristics
	if (pp1data->state >= PP1_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= pp1data->stage2_numvals &&				// We have enough memory and
	    (numvals < pp1data->stage2_numvals * 2 ||				// less than twice as much memory now available or
	     pp1data->Dsection >= pp1data->numDsections / 2))			// stage 2 more than half done
		return (0);							// Use old plan

/* Find the least costly stage 2 plan. */
/* Try various values of D until we find the best one.  3 gwnums are required for multiples of V_D calculations and one gwnum for gg. */

	best_stage2_impl (pp1data->C_done, pp1data->C, numvals - 4, &pp1_stage2_cost, &cost_data);

/* If are continuing from a save file that was in stage 2 and the new plan doesn't look significant better than the old plan, then */
/* we use the old plan and its partially completed bit array. */

	if (pp1data->state >= PP1_STATE_STAGE2 &&				// Continuing a stage 2 save file and
	    numvals >= pp1data->stage2_numvals &&				// We have enough memory and
	    cost_data.stage2_numvals < pp1data->stage2_numvals * 2)		// new plan does not use significantly more memory
		return (0);							// Use old plan

/* If are continuing from a save file that was in stage 2, toss the save file's bit array. */

	if (pp1data->state >= PP1_STATE_STAGE2) {
		free (pp1data->bitarray);
		pp1data->bitarray = NULL;
	}

/* Set all the variables needed for this stage 2 plan */

	pp1data->interim_C = pp1data->C;
	pp1data->stage2_numvals = cost_data.stage2_numvals;
	pp1data->D = cost_data.D;
	pp1data->totrels = cost_data.totrels;
	pp1data->B2_start = cost_data.B2_start;
	pp1data->numDsections = cost_data.numDsections;
	pp1data->bitarraymaxDsections = cost_data.bitarraymaxDsections;
	pp1data->Dsection = 0;

	if (pp1data->state < PP1_STATE_STAGE2) {
		pp1data->first_relocatable = pp1data->B;
		pp1data->last_relocatable = pp1data->B2_start;
		pp1data->C_done = pp1data->B2_start;
	}

/* Create a bitmap maximizing the prime pairings */

	pp1data->bitarrayfirstDsection = pp1data->Dsection;
	stop_reason = fill_work_bitarray (pp1data->thread_num, &pp1data->sieve_info, pp1data->D, pp1data->totrels,
					  pp1data->first_relocatable, pp1data->last_relocatable, pp1data->C_done, pp1data->C,
					  pp1data->bitarraymaxDsections, &pp1data->bitarray);
	if (stop_reason) return (stop_reason);

	return (0);
}

/* Calculate the P+1 start point.  Peter Montgomery suggests 2/7 or 6/5. */

void pp1_calc_start (
	gwhandle *gwdata,
	int	numerator,
	int	denominator,
	giant	N,		/* Number we are factoring */
	gwnum	V)		/* Returned starting point */
{
	mpz_t	__frac, __inv, __N;
	giant	g_frac;

/* Convert number we are factoring to mpz_t */

	mpz_init (__N);
	gtompz (N, __N);

/* Compute inverse */

	mpz_init_set_ui (__inv, denominator);
	mpz_invert (__inv, __inv, __N);

/* Compute fraction and convert to gwnum */

	mpz_init (__frac);
	mpz_mul_ui (__frac, __inv, numerator);
	mpz_mod (__frac, __frac, __N);
	g_frac = allocgiant ((int) divide_rounding_up (mpz_sizeinbase (__frac, 2), 32));
	mpztog (__frac, g_frac);
	gianttogw (gwdata, g_frac, V);
	free (g_frac);

/* Cleanup and return */

	mpz_clear (__N);
	mpz_clear (__inv);
	mpz_clear (__frac);
}

/* Handy macros for Lucas doubling and adding */

#define luc_dbl(h,s,d)			gwsquare2 (h, s, d, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT)
#define luc_add(h,s1,s2,diff,d)		gwmulsub4 (h, s1, s2, diff, d, GWMUL_FFT_S1 | GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT)
#define luc_add_last(h,s1,s2,diff,d)	gwmulsub4 (h, s1, s2, diff, d, GWMUL_FFT_S1 | GWMUL_FFT_S2)

/* Helper routines for implementing a PRAC-like lucas multiply */

/* The costing function assigns a doubling operation a cost of 2 FFTs and an addition operation a cost of 2 FFTs. */
/* The addition should cost a little more as it does more memory accesses. */

#define swap(a,b)	{t=a;a=b;b=t;}
#define PP1_ADD_COST	2
#define PP1_DBL_COST	2

int pp1_lucas_cost (
	uint64_t n,
	uint64_t d)
{
	uint64_t e, t;
	unsigned long c;

	if (d >= n || d <= n/2) return (999999999);		/* Catch invalid costings */

	c = 0;
	e = n - d;
	d = d - e;

	c += PP1_ADD_COST;

	while (d != e) {
		if (d < e) {
			swap (d,e);
		}
		if (100 * d <= 296 * e) {
			d = d-e;
			c += PP1_ADD_COST;
		} else if ((d&1) == (e&1)) {
			d = (d-e) >> 1;
			c += PP1_DBL_COST + PP1_ADD_COST;
		} else if ((d&1) == 0) {
			d = d >> 1;
			c += PP1_DBL_COST + PP1_ADD_COST;
//		} else if (d > 4 * e && d%3 == 0) {
//			d = d/3-e;
//			c += PP1_DBL_COST + 3*PP1_ADD_COST;
//		} else if (d > 4 * e && d%3 == 3 - e%3) {
//			d = (d-e-e)/3;
//			c += PP1_DBL_COST + 3*PP1_ADD_COST;
		} else {
			d = d-e;
			c += PP1_ADD_COST;
//			e = e >> 1;
//			c += PP1_DBL_COST + PP1_ADD_COST;
		}
	}

	return (c);
}

__inline int pp1_lucas_cost_several (uint64_t n, uint64_t *d) {
	int	i, c, min;
	uint64_t testd;
	for (i = 0, testd = *d - PRAC_SEARCH / 2; i < PRAC_SEARCH; i++, testd++) {
		c = pp1_lucas_cost (n, testd);
		if (i == 0 || c < min) min = c, *d = testd;
	}
	return (min);
}

void pp1_lucas_mul (
	pp1handle *pp1data,
	uint64_t n,
	uint64_t d,
	int	last_mul)	// TRUE if the last multiply should not use the GWMUL_STARTNEXTFFT option
{
	uint64_t e, t;
	gwnum	A, B, C;
	A = pp1data->V;
	B = pp1data->Vn;
	C = pp1data->Vn1;

	luc_dbl (&pp1data->gwdata, A, B);				/* B = 2*A */
									/* C = A (but we delay setting that up) */

	e = n - d;
	d = d - e;

	// To save a gwcopy setting C=A, we handle the most common case for the first iteration of the following loop.
	// I've only seen three cases that end up doing the gwcopy, n=3, n=11 and n=17.  With change to 2.96 there are a couple more.

	if (e > d && (100 * e <= 296 * d)) {
		swap (d, e);
		gwswap (A, B);							/* swap A & B, thus diff C = B */
		luc_add (&pp1data->gwdata, A, B, B, C);				/* B = A+B */
		gwswap (B, C);							/* C = B */
		d = d-e;
	}
	else if (d > e && (100 * d <= 296 * e)) {
		luc_add (&pp1data->gwdata, A, B, A, C);				/* B = A+B */
		gwswap (B, C);							/* C = B */
		d = d-e;
	} else {
		gwcopy (&pp1data->gwdata, A, C);				/* C = A */
	}

	while (d != e) {
		if (d < e) {
			swap (d, e);
			gwswap (A, B);
		}
		if (100 * d <= 296 * e) {					/* d <= 2.96 * e (Montgomery used 4.00) */
			luc_add (&pp1data->gwdata, A, B, C, C);			/* B = A+B */
			gwswap (B, C);						/* C = B */
			d = d-e;
		} else if ((d&1) == (e&1)) {
			luc_add (&pp1data->gwdata, A, B, C, B);			/* B = A+B */
			luc_dbl (&pp1data->gwdata, A, A);			/* A = 2*A */
			d = (d-e) >> 1;
		} else if ((d&1) == 0) {
			luc_add (&pp1data->gwdata, A, C, B, C);			/* C = A+C */
			luc_dbl (&pp1data->gwdata, A, A);			/* A = 2*A */
			d = d >> 1;
//		} else if (d > 4 * e && d%3 == 0) {		These provided little benefit, more research needed (adjust d>4*e)? Try remaining PRAC rules?
//			gwnum S = gwalloc (&pp1data->gwdata);
//			gwnum T = gwalloc (&pp1data->gwdata);
//			luc_dbl (&pp1data->gwdata, A, S);			/* S = 2*A */
//			luc_add (&pp1data->gwdata, A, B, C, T);			/* T = A+B */
//			luc_add (&pp1data->gwdata, S, T, C, C);			/* B = S+T */
//			luc_add (&pp1data->gwdata, S, A, A, A);			/* A = S+A */
//			gwswap (B, C);						/* C = B */
//			d = d/3-e;
//			gwfree (&pp1data->gwdata, S);
//			gwfree (&pp1data->gwdata, T);
//		} else if (d > 4 * e && d%3 == 3 - e%3) {
//			gwnum S = gwalloc (&pp1data->gwdata);
//			luc_add (&pp1data->gwdata, A, B, C, S);			/* S = A+B */
//			luc_add (&pp1data->gwdata, A, S, B, B);			/* B = A+S */
//			luc_dbl (&pp1data->gwdata, A, S);			/* S = 2*A */
//			luc_add (&pp1data->gwdata, S, A, A, A);			/* A = S+A */
//			d = (d-e-e)/3;
//			gwfree (&pp1data->gwdata, S);
		} else {
			luc_add (&pp1data->gwdata, A, B, C, C);			/* B = A+B */
			gwswap (B, C);						/* C = B */
			d = d-e;
//			luc_add (&pp1data->gwdata, B, C, A, C);			/* C = C-B */		Montgomery's default is worse than the above
//			luc_dbl (&pp1data->gwdata, B, B);			/* B = 2*B */
//			e = e >> 1;
		}
	}

	if (!last_mul)
		luc_add (&pp1data->gwdata, A, B, C, A);				/* A = A+B */
	else
		luc_add_last (&pp1data->gwdata, A, B, C, A);			/* A = A+B */

	pp1data->V = A;
	pp1data->Vn = B;
	pp1data->Vn1 = C;
}
#undef swap
#undef PP1_ADD_COST
#undef PP1_DBL_COST

/* Compute V_i from V_1 and i.  Routine limited to a multiplier of 2 or an odd number above 2. */

void pp1_mul (
	pp1handle *pp1data,
	uint64_t multiplier,
	int	last_mul)		// TRUE if the last multiply should not use the GWMUL_STARTNEXTFFT option
{
	ASSERTG (multiplier == 2 || (multiplier > 2 && (multiplier & 1)));

/* Perform a P+1 multiply using a modified PRAC algorithm developed by Peter Montgomery.  Basically, we try to find a near optimal Lucas */
/* chain of additions that generates the number we are multiplying by. */

	if (multiplier > 12) {
		int	c, min;
		uint64_t n, d, mind;

/* Cost a series of Lucas chains to find the cheapest.  First try v = (1+sqrt(5))/2, then (2+v)/(1+v), then (3+2*v)/(2+v), then (5+3*v)/(3+2*v), etc. */

		n = multiplier;
		mind = (uint64_t) ceil((double) 0.6180339887498948 * n);		/*v=(1+sqrt(5))/2*/
		min = pp1_lucas_cost_several (n, &mind);

		d = (uint64_t) ceil ((double) 0.7236067977499790 * n);			/*(2+v)/(1+v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.5801787282954641 * n);			/*(3+2*v)/(2+v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6328398060887063 * n);			/*(5+3*v)/(3+2*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6124299495094950 * n);			/*(8+5*v)/(5+3*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6201819808074158 * n);			/*(13+8*v)/(8+5*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6172146165344039 * n);			/*(21+13*v)/(13+8*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6183471196562281 * n);			/*(34+21*v)/(21+13*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6179144065288179 * n);			/*(55+34*v)/(34+21*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

		d = (uint64_t) ceil ((double) 0.6180796684698958 * n);			/*(89+55*v)/(55+34*v)*/
		c = pp1_lucas_cost_several (n, &d);
		if (c < min) min = c, mind = d;

/* Execute the cheapest Lucas chain */

		pp1_lucas_mul (pp1data, n, mind, last_mul);
	}

/* Otherwise, use a standard 2 operations (4 FFTs) per bit algorithm */

	else {
		gwnum	V, Vn, Vn1;
		uint64_t mask;
		int	last_mul_options = last_mul ? GWMUL_STARTNEXTFFT : 0;

		V = pp1data->V;
		Vn = pp1data->Vn;
		Vn1 = pp1data->Vn1;

/* Handle multiplier of 2 */

		if (multiplier == 2) {
			gwsquare2 (&pp1data->gwdata, V, V, GWMUL_ADDINCONST | last_mul_options);
			return;
		}
		if (multiplier == 3) {
			gwsquare2 (&pp1data->gwdata, V, Vn, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);	// V2
			gwmulsub4 (&pp1data->gwdata, V, Vn, V, V, last_mul_options);					// V3 = V1 + V2, diff 1
			return;
		}

/* Find top bit in multiplier */

		for (mask = 1ULL << 63; (multiplier & mask) == 0; mask >>= 1);

/* Init V values processing second mask bit */

		mask >>= 1;
		if (multiplier & mask) {
			gwsquare2 (&pp1data->gwdata, V, Vn1, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);	// V2
			gwmulsub4 (&pp1data->gwdata, V, Vn1, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);			// next n   = V3 = V1 + V2, diff 1
			gwsquare2 (&pp1data->gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);			// next n+1 = V4 = V2 + V2
		}
		else {
			gwsquare2 (&pp1data->gwdata, V, Vn, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);	// next n   = V2
			gwmulsub4 (&pp1data->gwdata, V, Vn, V, Vn1, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);			// next n+1 = V3 = V1 + V2, diff 1
		}

/* Loop processing rest of the mask one multiplier bit at a time.  Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1} */

		gwfft_for_fma (&pp1data->gwdata, V, V);	//GW: option that does this as part of mulsub4 would be better: GWMUL_FFT_FOR_FMA_S3!  Use it during init above.
		for (mask >>= 1; ; mask >>= 1) {
			if (mask == 1) {
				gwmulsub4 (&pp1data->gwdata, Vn, Vn1, V, V, last_mul_options);				// final n   = n + (n+1), diff 1
				break;
			}
			if (multiplier & mask) {
				gwmulsub4 (&pp1data->gwdata, Vn, Vn1, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// next n   = n + (n+1), diff 1
				gwsquare2 (&pp1data->gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);		// next n+1 = (n+1) + (n+1)
			}
			else {
				gwmulsub4 (&pp1data->gwdata, Vn, Vn1, V, Vn1, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	// next n+1 = n + (n+1), diff 1
				gwsquare2 (&pp1data->gwdata, Vn, Vn, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);		// next n   = n + n
			}
		}
	}
}

/* Compute V_{D*i} and V_{D*(i+1)} from V_D and i.  Use this to prepare for a series of lucas additions. */

void pp1_luc_prep_for_additions (
	pp1handle *pp1data,
	uint64_t multiplier)
{
	gwnum	V, Vn, Vn1;
	uint64_t mask;

/* Find top bit in multiplier */

	for (mask = 1ULL << 63; (multiplier & mask) == 0; mask >>= 1);

/* Init V values processing second mask bit */

	V = pp1data->V;
	Vn = pp1data->Vn;
	Vn1 = pp1data->Vn1;
	mask >>= 1;
	if (multiplier & mask) {
		gwsquare2 (&pp1data->gwdata, V, Vn1, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);	// V2
		gwmulsub4 (&pp1data->gwdata, V, Vn1, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);			// next n   = V3 = V1 + V2, diff 1
		gwsquare2 (&pp1data->gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);			// next n+1 = V4 = V2 + V2
	}
	else {
		gwsquare2 (&pp1data->gwdata, V, Vn, GWMUL_FFT_S1 | GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);	// next n   = V2
		gwmulsub4 (&pp1data->gwdata, V, Vn, V, Vn1, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);			// next n+1 = V3 = V1 + V2, diff 1
	}

/* Loop processing one multiplier bit at a time.  Each iteration computes V_n = V_{2n+bit}, V_{n+1} = V_{2n+bit+1} */

	for (mask >>= 1; mask; mask >>= 1) {
	        if (multiplier & mask) {
			gwmulsub4 (&pp1data->gwdata, Vn, Vn1, V, Vn, GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);	// next n   = n + (n+1), diff 1
			gwsquare2 (&pp1data->gwdata, Vn1, Vn1, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);		// next n+1 = (n+1) + (n+1)
		}
		else {
			gwmulsub4 (&pp1data->gwdata, Vn, Vn1, V, Vn1, GWMUL_FFT_S1 | GWMUL_STARTNEXTFFT);	// next n+1 = n + (n+1), diff 1
			gwsquare2 (&pp1data->gwdata, Vn, Vn, GWMUL_ADDINCONST | GWMUL_STARTNEXTFFT);		// next n   = n + n
		}
	}
}

/*****************************************************************************/
/*                         Main P+1 routine				     */
/*****************************************************************************/

int pplus1 (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	pp1handle pp1data;
	giant	N;		/* Number being factored */
	giant	factor;		/* Factor found, if any */
	unsigned int memused;
	int	i, res, stop_reason, first_iter_msg, saving, near_fft_limit, echk;
	char	filename[32], buf[255], JSONbuf[4000], testnum[100];
	uint64_t sieve_start, next_prime;
	double	one_over_B, base_pct_complete, one_bit_pct;
	double	last_output, last_output_t, last_output_r;
	double	allowable_maxerr, output_frequency, output_title_frequency;
	char	*str, *msg;
	int	msglen;
	double	timers[2];

/* Output a blank line to separate multiple P+1 runs making the result more readable */

	OutputStr (thread_num, "\n");

/* Clear pointers to allocated memory (so common error exit code knows what to free) */

	N = NULL;
	factor = NULL;
	str = NULL;
	msg = NULL;

/* Begin initializing P+1 data structure */
/* Choose a default value for the second bound if none was specified */

	memset (&pp1data, 0, sizeof (pp1handle));
	pp1data.thread_num = thread_num;
	pp1data.w = w;
	pp1data.B = (uint64_t) w->B1;
	pp1data.C = (uint64_t) w->B2;
	if (pp1data.B < 30) {
		OutputStr (thread_num, "Using minimum bound #1 of 30\n");
		pp1data.B = 30;
	}
	if (pp1data.C == 0) pp1data.C = pp1data.B * 100;
	if (pp1data.C < pp1data.B) pp1data.C = pp1data.B;
	pp1data.pct_mem_to_use = 1.0;				// Use as much memory as we can unless we get allocation errors

/* Convert run number into starting point */

	if (w->nth_run <= 1) {				// Per Peter Montgomery, make 2/7 the first starting point to try
		pp1data.numerator = 2;
		pp1data.denominator = 7;
		pp1data.takeAwayBits = log2 (6.0);	// 2/7 finds all factors that are 5 mod 6.  That is, p+1 is divisible by 6.
		pp1data.success_rate = 1.0;
	} else if (w->nth_run == 2) {			// Per Peter Montgomery, make 6/5 the second starting point to try
		pp1data.numerator = 6;
		pp1data.denominator = 5;
		pp1data.takeAwayBits = log2 (4.0);	// 6/5 finds all factors that are 3 mod 4.  That is, p+1 is divisible by 4.
		pp1data.success_rate = 0.5;		// Assume 2/7 has already been run.  It will have found half our 3 mod 4 factors.
	} else {
		unsigned char small_primes[16] = {11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
		srand ((unsigned) time (NULL));
		pp1data.numerator = 72 + (rand () & 0x7F);
		pp1data.denominator = small_primes[rand() & 0xF];
		pp1data.takeAwayBits = log2 (2.0);	// Random start point has no special characteristics.  That is, p+1 is divisible by 2.
		pp1data.success_rate = pow (0.5, (double) w->nth_run); // This run and each previous P+1 run find only half the smooth factors.
	}

/* Decide if we will calculate an optimal B2 when stage 2 begins.  We do this by default for P+1 work where we know how much TF has been done. */

	pp1data.optimal_B2 = (!QA_IN_PROGRESS && w->sieve_depth > 50 && IniGetInt (INI_FILE, "Pplus1BestB2", 1));
	if (pp1data.optimal_B2) pp1data.C = 100 * pp1data.B;	// A guess to use for calling start_sieve_with_limit

/* Compute the number we are factoring */

	stop_reason = setN (thread_num, w, &N);
	if (stop_reason) goto exit;

/* Other initialization */

	PRAC_SEARCH = IniGetInt (INI_FILE, "PracSearch", 7);
	if (PRAC_SEARCH < 1) PRAC_SEARCH = 1;
	if (PRAC_SEARCH > 50) PRAC_SEARCH = 50;

/* Output startup message */

	gw_as_string (testnum, w->k, w->b, w->n, w->c);
	sprintf (buf, "%s P+1", testnum);
	title (thread_num, buf);
	if (pp1data.C <= pp1data.B)
		sprintf (buf, "P+1 on %s, start=%" PRIu32 "/%" PRIu32 ", B1=%" PRIu64 "\n", testnum, pp1data.numerator, pp1data.denominator, pp1data.B);
	else if (pp1data.optimal_B2)
		sprintf (buf, "P+1 on %s, start=%" PRIu32 "/%" PRIu32 ", B1=%" PRIu64 ", B2=TBD\n", testnum, pp1data.numerator, pp1data.denominator, pp1data.B);
	else
		sprintf (buf, "P+1 on %s, start=%" PRIu32 "/%" PRIu32 ", B1=%" PRIu64 ", B2=%" PRIu64 "\n", testnum, pp1data.numerator, pp1data.denominator, pp1data.B, pp1data.C);
	OutputStr (thread_num, buf);
	if (w->sieve_depth > 0.0 && !pp1data.optimal_B2) {
		double prob = pm1prob (pp1data.takeAwayBits, (unsigned int) w->sieve_depth, (double) pp1data.B, (double) pp1data.C) * pp1data.success_rate;
		sprintf (buf, "Chance of finding a factor assuming no ECM has been done is an estimated %.3g%%\n", prob * 100.0);
		OutputStr (thread_num, buf);
	}

/* Clear all timers */

restart:
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init filename */

	tempFileName (w, filename);
	filename[0] = 'n';

/* Perform setup functions.  This includes decding how big an FFT to use, allocating memory, calling the FFT setup code, etc. */

/* Setup the assembly code */

	gwinit (&pp1data.gwdata);
	gwset_sum_inputs_checking (&pp1data.gwdata, SUM_INPUTS_ERRCHK);
	if (IniGetInt (LOCALINI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&pp1data.gwdata);
	if (IniGetInt (INI_FILE, "HyperthreadPrefetch", 0)) gwset_hyperthread_prefetch (&pp1data.gwdata);
	if (HYPERTHREAD_LL) {
		sp_info->normal_work_hyperthreads = IniGetInt (LOCALINI_FILE, "HyperthreadLLcount", CPU_HYPERTHREADS);
		gwset_will_hyperthread (&pp1data.gwdata, sp_info->normal_work_hyperthreads);
	}
	gwset_bench_cores (&pp1data.gwdata, NUM_CPUS);
	gwset_bench_workers (&pp1data.gwdata, NUM_WORKER_THREADS);
	if (ERRCHK) gwset_will_error_check (&pp1data.gwdata);
	else gwset_will_error_check_near_limit (&pp1data.gwdata);
	gwset_num_threads (&pp1data.gwdata, CORES_PER_TEST[thread_num] * sp_info->normal_work_hyperthreads);
	gwset_thread_callback (&pp1data.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&pp1data.gwdata, sp_info);
	gwset_safety_margin (&pp1data.gwdata, IniGetFloat (INI_FILE, "ExtraSafetyMargin", 0.0));
	gwset_minimum_fftlen (&pp1data.gwdata, w->minimum_fftlen);
	res = gwsetup (&pp1data.gwdata, w->k, w->b, w->n, w->c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}
	gwsetaddin (&pp1data.gwdata, -2);

/* A kludge so that the error checking code is not as strict. */

	pp1data.gwdata.MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);

/* Allocate gwnums for computing Lucas sequences */

	pp1data.V = gwalloc (&pp1data.gwdata);
	if (pp1data.V == NULL) goto oom;
	pp1data.Vn = gwalloc (&pp1data.gwdata);
	if (pp1data.Vn == NULL) goto oom;
	pp1data.Vn1 = gwalloc (&pp1data.gwdata);
	if (pp1data.Vn1 == NULL) goto oom;

/* More miscellaneous initializations */

	last_output = last_output_t = last_output_r = 0;
	gw_clear_fft_count (&pp1data.gwdata);
	first_iter_msg = TRUE;
	calc_output_frequencies (&pp1data.gwdata, &output_frequency, &output_title_frequency);

/* Output message about the FFT length chosen */

	{
		char	fft_desc[200];
		gwfft_description (&pp1data.gwdata, fft_desc);
		sprintf (buf, "Using %s\n", fft_desc);
		OutputStr (thread_num, buf);
	}

/* If we are near the maximum exponent this fft length can test, then we will roundoff check all multiplies */

	near_fft_limit = exponent_near_fft_limit (&pp1data.gwdata);
	gwerror_checking (&pp1data.gwdata, ERRCHK || near_fft_limit);

/* Figure out the maximum round-off error we will allow.  By default this is 27/64 when near the FFT limit and 26/64 otherwise. */
/* We've found that this default catches errors without raising too many spurious error messages.  We let the user override */
/* this default for user "Never Odd Or Even" who tests exponents well beyond an FFT's limit.  He does his error checking by */
/* running the first-test and double-check simultaneously. */

	allowable_maxerr = IniGetFloat (INI_FILE, "MaxRoundoffError", (float) (near_fft_limit ? 0.421875 : 0.40625));

/* Check for a save file and read the save file.  If there is an error reading the file then restart the P+1 factoring job from scratch. */
/* Limit number of backup files we try to read in case there is an error deleting bad save files. */

	readSaveFileStateInit (&pp1data.read_save_file_state, thread_num, filename);
	writeSaveFileStateInit (&pp1data.write_save_file_state, filename, 0);
	for ( ; ; ) {
		if (! saveFileExists (&pp1data.read_save_file_state)) {
			/* If there were save files, they are all bad.  Report a message */
			/* and temporarily abandon the work unit.  We do this in hopes that */
			/* we can successfully read one of the bad save files at a later time. */
			/* This sounds crazy, but has happened when OSes get in a funky state. */
			if (pp1data.read_save_file_state.a_non_bad_save_file_existed) {
				OutputBoth (thread_num, ALLSAVEBAD_MSG);
				return (0);
			}
			/* No save files existed, start from scratch. */
			break;
		}

		if (!pp1_restore (&pp1data)) {
			/* Close and rename the bad save file */
			saveFileBad (&pp1data.read_save_file_state);
			continue;
		}

/* Handle stage 1 save files.  If the save file that had a higher B1 target then we can reduce the target B1 to the desired B1. */

		if (pp1data.state == PP1_STATE_STAGE1) {
			sieve_start = pp1data.stage1_prime + 1;
			if (pp1data.interim_B > pp1data.B) pp1data.interim_B = pp1data.B;
			goto restart1;
		}

/* Handle between stages save files */

		if (pp1data.state == PP1_STATE_MIDSTAGE) {
			if (pp1data.B > pp1data.B_done) {
				gwfree (&pp1data.gwdata, pp1data.gg), pp1data.gg = NULL;
				goto more_B;
			}
			goto restart3b;
		}

/* Handle stage 2 save files */

		if (pp1data.state == PP1_STATE_STAGE2) {

/* If B is larger than the one in the save file, then do more stage 1 processing. */

			if (pp1data.B > pp1data.B_done) {
				gwfree (&pp1data.gwdata, pp1data.gg), pp1data.gg = NULL;
				free (pp1data.bitarray), pp1data.bitarray = NULL;
				goto more_B;
			}

/* If B is different than the one in the save file, then use the one in the save file rather than discarding all the work done thusfar in stage 2. */

			if (pp1data.B != pp1data.B_done) {
				pp1data.B = pp1data.B_done;
				sprintf (buf, "Ignoring suggested B1 value, using B1=%" PRIu64 " from the save file\n", pp1data.B);
				OutputStr (thread_num, buf);
			}

/* If bound #2 is larger in the save file then use the original bound #2.  The user that wants to discard the stage 2 work he */
/* has done thusfar and reduce the stage 2 bound must manually delete the save file. */

			if (pp1data.C < pp1data.interim_C) {
				pp1data.C = pp1data.interim_C;
				sprintf (buf, "Ignoring suggested B2 value, using B2=%" PRIu64 " from the save file\n", pp1data.C);
				OutputStr (thread_num, buf);
			}

/* Resume stage 2 */

			if (pp1data.optimal_B2) {
				pp1data.C = pp1data.interim_C;
				sprintf (buf, "Resuming P+1 in stage 2 with B2=%" PRIu64 "\n", pp1data.interim_C);
				OutputStr (thread_num, buf);
			}
			goto restart3b;
		}

/* Handle stage 2 GCD save files */

		if (pp1data.state == PP1_STATE_GCD) {
			if (pp1data.optimal_B2) pp1data.C = pp1data.C_done;
			goto restart4;
		}

/* Handle case where we have a completed save file (the PP1_STATE_DONE state) */

		ASSERTG (pp1data.state == PP1_STATE_DONE);
		if (pp1data.B > pp1data.B_done) goto more_B;
		if (pp1data.C > pp1data.C_done) {
			pp1data.state = PP1_STATE_MIDSTAGE;
			goto restart3a;
		}

/* The save file indicates we've tested to these bounds already */

		sprintf (buf, "%s already tested to B1=%" PRIu64 " and B2=%" PRIu64 ".\n",
			 gwmodulo_as_string (&pp1data.gwdata), pp1data.B_done, pp1data.C_done);
		OutputBoth (thread_num, buf);
		goto done;
	}

/* Start this P+1 run from scratch starting with V = 2/7 or 6/5 or random */

	pp1_calc_start (&pp1data.gwdata, pp1data.numerator, pp1data.denominator, N, pp1data.V);	/* V = 2/7 or 6/5 or random mod N */
	pp1data.B_done = 0;
	pp1data.interim_B = pp1data.B;
	sieve_start = 2;

/* The stage 1 restart point */

restart1:
	strcpy (w->stage, "S1");
	pp1data.state = PP1_STATE_STAGE1;
	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&pp1data.gwdata, 3));
	start_timer (timers, 0);
	start_timer (timers, 1);
	stop_reason = start_sieve_with_limit (thread_num, sieve_start, (uint32_t) sqrt ((double) pp1data.C), &pp1data.sieve_info);
	if (stop_reason) goto exit;
	pp1data.stage1_prime = sieve (pp1data.sieve_info);
	one_over_B = 1.0 / (double) pp1data.B;
	w->pct_complete = (double) pp1data.stage1_prime * one_over_B;
	for ( ; pp1data.stage1_prime <= pp1data.interim_B; pp1data.stage1_prime = next_prime) {
		uint64_t mult, max;
		int	count;

		next_prime = sieve (pp1data.sieve_info);

/* Test for user interrupt, save files, and error checking */

		stop_reason = stopCheck (thread_num);
		saving = testSaveFilesFlag (thread_num);
		echk = stop_reason || saving || ERRCHK || near_fft_limit || ((pp1data.stage1_prime & 255) == 255);
		gwerror_checking (&pp1data.gwdata, echk);

/* Count the number of prime powers where B_done < prime^n <= interim_B */

		count = 0;
		max = pp1data.interim_B / pp1data.stage1_prime;
		for (mult = pp1data.stage1_prime; ; mult *= pp1data.stage1_prime) {
			if (mult > pp1data.B_done) count++;
			if (mult > max) break;
		}

/* Adjust the count of prime powers.  All p+1 factors are a multiple of 2.  For primes 2 and 3, we increase the count by one since a */
/* starting fraction of 2/7 has p+1 factors that are a multiple of 6 and fraction of 6/5 has p+1 factors that are a multiple of 4. */

		if (pp1data.stage1_prime <= 3 && pp1data.B_done == 0) count += (pp1data.stage1_prime == 2) ? 2 : 1;

/* For prime 2, we use square_carefully as the fraction used in calculating the initial V value can result in pathological data. */

		if (pp1data.stage1_prime == 2 && pp1data.B_done == 0) {
			if (count < (int) log2 (w->n) + 3) count = (int) log2 (w->n) + 3;
			for ( ; count; count--) {
				gwmul3_carefully (&pp1data.gwdata, pp1data.V, pp1data.V, pp1data.V, GWMUL_ADDINCONST);
			}
			/* Free the memory allocated for GW_RANDOM */
			gwfree_internal_memory (&pp1data.gwdata);
			/* Include the Mersenne (or generalized Fermat) exponent.  I was against this (encourages using P+1 to find P-1 factors */
			/* rather than the more efficient P-1 factoring algorithm).  But just in case a P-1 run was done and a hardware error */
			/* occurred, give P+1 a chance to find a missed factor. */
			if (w->n > pp1data.B && (isMersenne (w->k, w->b, w->n, w->c) || isGeneralizedFermat (w->k, w->b, w->n, w->c)))
				pp1_mul (&pp1data, w->n, FALSE);
		}

/* Apply the proper number of primes */

		else {
			for ( ; count; count--) {
				int	last_mul = (count == 1 && (saving || stop_reason || next_prime > pp1data.B));
				pp1_mul (&pp1data, pp1data.stage1_prime, last_mul);
			}
		}

/* Test for an error */

		if (gw_test_for_error (&pp1data.gwdata) || gw_get_maxerr (&pp1data.gwdata) > allowable_maxerr) goto err;

/* Calculate our stage 1 percentage complete */

		w->pct_complete = (double) pp1data.stage1_prime * one_over_B;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P+1 stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pp1data.gwdata));
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pp1data.gwdata);
		}

/* Every N squarings, output a progress report */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.",
				 gwmodulo_as_string (&pp1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			if (pp1data.stage1_prime != 2 || pp1data.B_done != 0) OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pp1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Every N squarings, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pp1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pp1data.gwdata);
		}

/* Check for escape and/or if its time to write a save file */

		if (stop_reason || saving) {
			pp1_save (&pp1data);
			if (stop_reason) goto exit;
		}
	}

	pp1data.B_done = pp1data.interim_B;
	end_timer (timers, 0);
	end_timer (timers, 1);

/* Check for the rare case where we need to do even more stage 1.  This happens using a save file created with a smaller bound #1. */

	if (pp1data.B > pp1data.B_done) {
more_B:		pp1data.interim_B = pp1data.B;
		sieve_start = 2;
//GW - pct_complete resets  to 0 because of setting stage1_prime to 2
		goto restart1;
	}
	pp1data.C_done = pp1data.B;

/* Do stage 1 cleanup, resume "standard" error checking */

	gwerror_checking (&pp1data.gwdata, ERRCHK || near_fft_limit);

/* Stage 1 complete, print a message */

	sprintf (buf, "%s stage 1 complete. %.0f transforms. Time: ", gwmodulo_as_string (&pp1data.gwdata), gw_get_fft_count (&pp1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));
	gw_clear_fft_count (&pp1data.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pp1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pp1data.gwdata);
	}

/* Check to see if we found a factor - do GCD (V-2, N) */

	strcpy (w->stage, "S1");
	w->pct_complete = 1.0;
	if (pp1data.C <= pp1data.B || (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "Stage1GCD", 1))) {
		start_timer (timers, 0);
		gwaddsmall (&pp1data.gwdata, pp1data.V, -2);
		stop_reason = gcd (&pp1data.gwdata, thread_num, pp1data.V, N, &factor);
		gwaddsmall (&pp1data.gwdata, pp1data.V, 2);
		if (stop_reason) {
			pp1_save (&pp1data);
			goto exit;
		}
		end_timer (timers, 0);
		strcpy (buf, "Stage 1 GCD complete. Time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputStr (thread_num, buf);
		if (factor != NULL) goto bingo;
	}

/* Skip second stage if so requested */

	if (pp1data.C <= pp1data.B) goto msg_and_exit;

/*
   Stage 2:  Use ideas from Crandall, Zimmermann, Montgomery, and Preda on each prime below C.
   This code is more efficient the more memory you can give it.
   Inputs: x, the value at the end of stage 1
*/

/* Change state to between stage 1 and 2 */

	pp1data.state = PP1_STATE_MIDSTAGE;
	strcpy (w->stage, "S2");
	w->pct_complete = 0.0;

/* Initialize variables for second stage.  We set gg to Vn-2 in case the user opted to skip the GCD after stage 1. */

restart3a:
	pp1data.gg = gwalloc (&pp1data.gwdata);
	if (pp1data.gg == NULL) goto oom;
//GW: We could start gg with VD-2 later on -- reduces peak by one gwnum (ECM too?)  BUT VD is FFTed possibly FFTED_FOR_FMA
	gwcopy (&pp1data.gwdata, pp1data.V, pp1data.gg);
	gwaddsmall (&pp1data.gwdata, pp1data.gg, -2);

/* Restart here when in the middle of stage 2 */

restart3b:
more_C:
	start_timer (timers, 0);
	sprintf (buf, "%s P+1 stage 2 init", gwmodulo_as_string (&pp1data.gwdata));
	title (thread_num, buf);

/* Clear flag indicating we need to restart if the maximum amount of memory changes. */
/* Prior to this point we allow changing the optimal bounds of a Pfactor assignment. */

	clear_restart_if_max_memory_change (thread_num);

/* Test if we will ever have enough memory to do stage 2 based on the maximum available memory. */
/* Our minimum working set is one gwnum for gg, 4 for nQx, 3 for VD, Vn, Vn1. */

replan:	{
		unsigned long min_memory = cvt_gwnums_to_mem (&pp1data.gwdata, 8);
		if (max_mem (thread_num) < min_memory) {
			sprintf (buf, "Insufficient memory to ever run stage 2 -- %luMB needed.\n", min_memory);
			OutputStr (thread_num, buf);
			pp1data.C = pp1data.B_done;
			goto restart4;
		}
	}

/* Choose the best plan implementation given the currently available memory. */
/* This implementation could be "wait until we have more memory". */

	stop_reason = pp1_stage2_impl (&pp1data);
	if (stop_reason) {
		if (pp1data.state == PP1_STATE_MIDSTAGE) pp1_save (&pp1data);
		goto exit;
	}

/* Record the amount of memory this thread will be using in stage 2. */

	memused = cvt_gwnums_to_mem (&pp1data.gwdata, pp1data.stage2_numvals);
	memused += (_intmin (pp1data.numDsections - pp1data.bitarrayfirstDsection, pp1data.bitarraymaxDsections) * pp1data.totrels) >> 23;
	// To dodge possible infinite loop if pp1_stage2_impl allocates too much memory (it shouldn't), decrease the percentage of memory we are allowed to use
	// Beware that replaning may allocate a larger bitarray
	if (set_memory_usage (thread_num, MEM_VARIABLE_USAGE, memused)) {
		pp1data.pct_mem_to_use *= 0.99;
		free (pp1data.bitarray); pp1data.bitarray = NULL;
		goto replan;
	}

/* Output a useful message regarding memory usage */

	sprintf (buf, "Using %uMB of memory.\n", memused);
	OutputStr (thread_num, buf);

/* Initialize variables for second stage */

	// Calculate the percent completed by previous bit arrays
	base_pct_complete = (double) (pp1data.B2_start - pp1data.last_relocatable) / (double) (pp1data.C - pp1data.last_relocatable);
	// Calculate the percent completed by each bit in this bit array
	one_bit_pct = (1.0 - base_pct_complete) / (double) (pp1data.numDsections * pp1data.totrels);
	// Calculate the percent completed by previous bit arrays and the current bit array
	w->pct_complete = base_pct_complete + (double) (pp1data.Dsection * pp1data.totrels) * one_bit_pct;

/* Allocate arrays of pointers to stage 2 gwnums */

	pp1data.nQx = (gwnum *) malloc (pp1data.totrels * sizeof (gwnum));
	if (pp1data.nQx == NULL) goto lowmem;

/* A fast 1,5 mod 6 nQx initialization */

	{
		gwnum	t1, t2, t3, t4, t5;
		gwnum	*V1, *V2, *V3, *V4, *V5, *V6, *V7, *V11, *V1mod6, *V1mod6minus6, *V5mod6, *V5mod6minus6;
		int	totrels, have_VD;

/* Allocate memory and init values for computing nQx.  We need V_1, V_5, V_6, V_7, V_11. */
/* We also need V_4 to compute V_D in the middle of the nQx init loop. */
/* NOTE:  At end of the loop V_D is stored in pp1data.V */

		t1 = gwalloc (&pp1data.gwdata);
		if (t1 == NULL) goto lowmem;
		t2 = gwalloc (&pp1data.gwdata);
		if (t2 == NULL) goto lowmem;
		t3 = gwalloc (&pp1data.gwdata);
		if (t3 == NULL) goto lowmem;
		t4 = pp1data.Vn;
		t5 = pp1data.Vn1;
//GW  we could save more temps (perhaps) by doing all the 1 mod 6, then all the 5 mod 6

/* Compute V_2 and place in pp1data.V.  This is the value we will save if init is aborted and we create a save file. */
/* V_2 will be replaced by V_D later on in this loop.  At all times pp1data.V will be a save-able value! */

		luc_dbl (&pp1data.gwdata, pp1data.V, t1);			// V2 = 2 * V1
		gwswap (t1, pp1data.V);
		V1 = &t1;
		V2 = &pp1data.V;
		pp1data.nQx[0] = *V1; totrels = 1;

/* Compute V_5, V_6, V_7, V_11 */

		V3 = &t2;
		luc_add (&pp1data.gwdata, *V2, *V1, *V1, *V3);			// V3 = V2 + V1 (diff V1)

		V5 = &t3;
		luc_add (&pp1data.gwdata, *V3, *V2, *V1, *V5);			// V5 = V3 + V2 (diff V1)
		if (relatively_prime (5, pp1data.D)) pp1data.nQx[totrels++] = *V5;

		V6 = V3;
		luc_dbl (&pp1data.gwdata, *V3, *V6);				// V6 = 2 * V3, V3 no longer needed

		V7 = &t4;
		luc_add (&pp1data.gwdata, *V6, *V1, *V5, *V7);			// V7 = V6 + V1 (diff V5)
		if (relatively_prime (7, pp1data.D)) pp1data.nQx[totrels++] = *V7;

		V11 = &t5;
		luc_add (&pp1data.gwdata, *V6, *V5, *V1, *V11);			// V11 = V6 + V5 (diff V1)
		if (relatively_prime (11, pp1data.D)) pp1data.nQx[totrels++] = *V11;

/* Compute the rest of the nQx values (V_i for i >= 7) */

		V1mod6 = V7;
		V1mod6minus6 = V1;
		V5mod6 = V11;
		V5mod6minus6 = V5;
		have_VD = FALSE;
		for (i = 7; ; i += 6) {

/* Compute V_D which we will need in pp1_luc_prep_for_additions.  Do this with a single luc_add call when we reach two values that */
/* are 4 apart that add to D.  This only works for D = 2 mod 4. */

			if (i + i+4 == pp1data.D) {
				ASSERTG (V2 == &pp1data.V);
				V4 = V2;
				luc_dbl (&pp1data.gwdata, *V2, *V4);					// V4 = 2 * V2, V2 no longer needed
				luc_add (&pp1data.gwdata, *V1mod6, *V5mod6, *V4, pp1data.V);		// VD = Vi + V{i+4} (diff V4), V4 no longer needed
				have_VD = TRUE;
			}

/* Break out of loop when we have all our nQx values */

			if (have_VD && totrels >= pp1data.totrels) break;

/* Get next 1 mod 6 value.  Don't destroy V1mod6minus6 if it is an nQx value. */

			if (relatively_prime (i-6, pp1data.D)) {
				gwnum	next1mod6 = gwalloc (&pp1data.gwdata);
				if (next1mod6 == NULL) goto lowmem;
				luc_add (&pp1data.gwdata, *V1mod6, *V6, *V1mod6minus6, next1mod6);	// Next 1 mod 6 value
				*V1mod6minus6 = *V1mod6;
				*V1mod6 = next1mod6;
			} else {
				luc_add (&pp1data.gwdata, *V1mod6, *V6, *V1mod6minus6, *V1mod6minus6);	// Next 1 mod 6 value
				gwswap (*V1mod6, *V1mod6minus6);
			}
			if (relatively_prime (i+6, pp1data.D) && totrels < pp1data.totrels) pp1data.nQx[totrels++] = *V1mod6;

/* Compute V_D which we will need in pp1_luc_prep_for_additions.  Do this with a single luc_add call when we reach the two values that */
/* are 2 apart that add to D.  This only works for D = 0 mod 4. */

			if (i+6 + i+4 == pp1data.D) {
				ASSERTG (V2 == &pp1data.V);
				luc_add (&pp1data.gwdata, *V1mod6, *V5mod6, *V2, pp1data.V);		// VD = V{i+6} + V{i+4} (diff V2), V2 no longer needed
				have_VD = TRUE;
			}

/* Break out of loop when we have VD and all our nQx values */

			if (have_VD && totrels >= pp1data.totrels) break;

/* Get next 5 mod 6 value.  Don't destroy V5mod6minus6 if it is an nQx value */

			if (relatively_prime (i+4-6, pp1data.D)) {
				gwnum	next5mod6 = gwalloc (&pp1data.gwdata);
				if (next5mod6 == NULL) goto lowmem;
				luc_add (&pp1data.gwdata, *V5mod6, *V6, *V5mod6minus6, next5mod6);	// Next 5 mod 6 value
				*V5mod6minus6 = *V5mod6;
				*V5mod6 = next5mod6;
			} else {
				luc_add (&pp1data.gwdata, *V5mod6, *V6, *V5mod6minus6, *V5mod6minus6);	// Next 5 mod 6 value
				gwswap (*V5mod6, *V5mod6minus6);
			}
			if (relatively_prime (i+4+6, pp1data.D) && totrels < pp1data.totrels) pp1data.nQx[totrels++] = *V5mod6;

/* Check for errors and save file creation */

			if (gw_test_for_error (&pp1data.gwdata)) goto err;

			stop_reason = stopCheck (thread_num);
			if (stop_reason) goto possible_lowmem;
		}

/* Free memory used in computing nQx values, complicated by the fact some of these could be among the final four nQx values. */

		for (i = 1; i <= 4; i++) {
			if (t1 == pp1data.nQx[totrels-i]) t1 = NULL;
			if (t2 == pp1data.nQx[totrels-i]) t2 = NULL;
			if (t3 == pp1data.nQx[totrels-i]) t3 = NULL;
			if (t4 == pp1data.nQx[totrels-i]) t4 = NULL;
			if (t5 == pp1data.nQx[totrels-i]) t5 = NULL;
		}
		gwfree (&pp1data.gwdata, t1);
		gwfree (&pp1data.gwdata, t2);
		gwfree (&pp1data.gwdata, t3);
		gwfree (&pp1data.gwdata, t4);
		gwfree (&pp1data.gwdata, t5);

/* Replace the Vn and Vn1 gwnums that we used as temporaries. */

		pp1data.Vn = gwalloc (&pp1data.gwdata);
		pp1data.Vn1 = gwalloc (&pp1data.gwdata);
	}

/* Initialize for computing successive multiples of V_D */

	uint64_t D_multiplier = (pp1data.B2_start + pp1data.D / 2) / pp1data.D + pp1data.Dsection;
	pp1_luc_prep_for_additions (&pp1data, D_multiplier - 1);	// Vn = V_{Dsection-1}, Vn1 = V_{Dsection}

/* Stage 2 init complete, change the title, output a message */

	sprintf (buf, "%.*f%% of %s P+1 stage 2 (using %uMB)",
		 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pp1data.gwdata), memused);
	title (thread_num, buf);

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 init complete. %.0f transforms. Time: ", gw_get_fft_count (&pp1data.gwdata));
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	gw_clear_fft_count (&pp1data.gwdata);

/* We do prime pairing with each loop iteration handling the range m-Q to m+Q where m is a multiple of D and Q is the */
/* Q-th relative prime to D.  Totrels is often much larger than the number of relative primes less than D.  This Preda/Montgomery */
/* optimization provides us with many more chances to find a prime pairing. */

	pp1data.state = PP1_STATE_STAGE2;
	start_timer (timers, 0);
	start_timer (timers, 1);
	last_output = last_output_t = last_output_r = 0;
	for ( ; ; ) {
		uint64_t bitnum;

/* Test which relprimes in this D section need to be processed */

		bitnum = (pp1data.Dsection - pp1data.bitarrayfirstDsection) * pp1data.totrels;
		for (i = 0; i < pp1data.totrels; i++) {

/* Skip this relprime if the pairing routine did not set the corresponding bit in the bitarray. */

			if (! bittst (pp1data.bitarray, bitnum + i)) continue;

/* Multiply this Vn1 - nQx value into the gg accumulator */

			saving = testSaveFilesFlag (thread_num);
			stop_reason = stopCheck (thread_num);
			gwsubmul4 (&pp1data.gwdata, pp1data.Vn1, pp1data.nQx[i], pp1data.gg, pp1data.gg, (!stop_reason && !saving) ? GWMUL_STARTNEXTFFT : 0);

/* Clear pairing bit from the bitarray in case a save file is written.  Calculate stage 2 percentage. */

			bitclr (pp1data.bitarray, bitnum + i);
			w->pct_complete = base_pct_complete + (double) (pp1data.Dsection * pp1data.totrels + i) * one_bit_pct;

/* Test for errors */

			if (gw_test_for_error (&pp1data.gwdata) || gw_get_maxerr (&pp1data.gwdata) > allowable_maxerr) goto err;

/* Output the title every so often */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
				sprintf (buf, "%.*f%% of %s P+1 stage 2 (using %uMB)",
					 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pp1data.gwdata), memused);
				title (thread_num, buf);
				last_output_t = gw_get_fft_count (&pp1data.gwdata);
			}

/* Write out a message every now and then */

			if (first_iter_msg ||
			    (ITER_OUTPUT != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
				sprintf (buf, "%s stage 2 is %.*f%% complete.",
					 gwmodulo_as_string (&pp1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
				end_timer (timers, 0);
				if (first_iter_msg) {
					strcat (buf, "\n");
					clear_timer (timers, 0);
				} else {
					strcat (buf, " Time: ");
					print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
				}
				OutputStr (thread_num, buf);
				start_timer (timers, 0);
				last_output = gw_get_fft_count (&pp1data.gwdata);
				first_iter_msg = FALSE;
			}

/* Write out a message to the results file every now and then */

			if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pp1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
			    (NO_GUI && stop_reason)) {
				sprintf (buf, "%s stage 2 is %.*f%% complete.\n",
					 gwmodulo_as_string (&pp1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
				writeResults (buf);
				last_output_r = gw_get_fft_count (&pp1data.gwdata);
			}

/* Periodicly write a save file */

			if (stop_reason || saving) {
				pp1_save (&pp1data);
				if (stop_reason) goto exit;
			}
		}

/* Move onto the next D section when we are done with all the relprimes */

		pp1data.Dsection++;
		pp1data.C_done = pp1data.B2_start + pp1data.Dsection * pp1data.D;
		if (pp1data.Dsection >= pp1data.numDsections) break;

		luc_add (&pp1data.gwdata, pp1data.Vn1, pp1data.V, pp1data.Vn, pp1data.Vn);	// Vn2 = Vn1 + V, diff = Vn
		gwswap (pp1data.Vn, pp1data.Vn1);						// Vn = Vn1, Vn1 = Vn2
//GW: check for errors?

/* See if more bitarrays are required to get us to bound #2 */

		if (pp1data.Dsection < pp1data.bitarrayfirstDsection + pp1data.bitarraymaxDsections) continue;

		pp1data.first_relocatable = calc_new_first_relocatable (pp1data.D, pp1data.C_done);
		pp1data.bitarrayfirstDsection = pp1data.Dsection;
		stop_reason = fill_work_bitarray (pp1data.thread_num, &pp1data.sieve_info, pp1data.D, pp1data.totrels,
						  pp1data.first_relocatable, pp1data.last_relocatable, pp1data.C_done, pp1data.C,
						  pp1data.bitarraymaxDsections, &pp1data.bitarray);
		if (stop_reason) {
			pp1_save (&pp1data);
			goto exit;
		}
	}
	pp1data.C_done = pp1data.interim_C;

/* Free up memory */

	for (i = 0; i < pp1data.totrels; i++) gwfree (&pp1data.gwdata, pp1data.nQx[i]);
	free (pp1data.nQx), pp1data.nQx = NULL;
	free (pp1data.bitarray), pp1data.bitarray = NULL;

/* Check for the rare cases where we need to do even more stage 2.  This happens when continuing a save file in the middle of stage 2 and */
/* the save file's target bound #2 is smaller than our target bound #2. */

	if (pp1data.C > pp1data.C_done) {
		pp1data.state = PP1_STATE_MIDSTAGE;
		goto more_C;
	}

/* Stage 2 is complete */

	end_timer (timers, 1);
	sprintf (buf, "%s stage 2 complete. %.0f transforms. Time: ", gwmodulo_as_string (&pp1data.gwdata), gw_get_fft_count (&pp1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pp1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pp1data.gwdata);
	}

/* See if we got lucky! */

restart4:
	pp1data.state = PP1_STATE_GCD;
	strcpy (w->stage, "S2");
	w->pct_complete = 1.0;
	start_timer (timers, 0);
	stop_reason = gcd (&pp1data.gwdata, thread_num, pp1data.gg, N, &factor);
	if (stop_reason) {
		pp1_save (&pp1data);
		goto exit;
	}
	pp1data.state = PP1_STATE_DONE;
	end_timer (timers, 0);
	strcpy (buf, "Stage 2 GCD complete. Time: ");
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	if (factor != NULL) goto bingo;

/* Output line to results file indicating P+1 run */

msg_and_exit:
	sprintf (buf, "%s completed P+1, B1=%" PRIu64, gwmodulo_as_string (&pp1data.gwdata), pp1data.B);
	if (pp1data.C > pp1data.B) sprintf (buf+strlen(buf), ", B2=%" PRIu64, pp1data.C);
	sprintf (buf+strlen(buf), ", Wi%d: %08lX\n", PORT, SEC5 (w->n, pp1data.B, pp1data.C));
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"NF", "exponent":45581713, "worktype":"P+1", "b1":50000, "b2":5000000, "start":"2/7", */
/* "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"30.6", "build":"1"}, "timestamp":"2021-04-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"NF\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P+1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64, pp1data.B);
	if (pp1data.C > pp1data.B) sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%" PRIu64, pp1data.C);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"start\":\"%" PRIu32 "/%" PRIu32 "\"", pp1data.numerator, pp1data.denominator);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pp1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, pp1data.B, pp1data.C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send P+1 completed message to the server.  Although don't do it for puny B1 values as this is just the user tinkering with P+1 factoring. */

	if (!QA_IN_PROGRESS && (pp1data.B >= 50000 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type = PRIMENET_AR_PP1_NOFACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.pp1_numerator = pp1data.numerator;
		pkt.pp1_denominator = pp1data.denominator;
		pkt.B1 = (double) pp1data.B;
		pkt.B2 = (double) pp1data.C;
		pkt.fftlen = gwfftlen (&pp1data.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Free save files (especially big ones created during stage 2) */
/* Create save file so that we can expand bound 1 or bound 2 at a later date.  Only allow this for 2/7 and 6/5 start. */
/* If we don't delete random start save files, then a new random start attempt would use the existing random start save file. */ 

	unlinkSaveFiles (&pp1data.write_save_file_state);
	if (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "KeepPplus1SaveFiles", 0) && w->nth_run <= 2)
		pp1_save (&pp1data);

/* Return stop code indicating success or work unit complete */ 

done:	stop_reason = STOP_WORK_UNIT_COMPLETE;

/* Free memory and return */

exit:	pp1_cleanup (&pp1data);
	free (N);
	free (factor);
	free (str);
	free (msg);
	return (stop_reason);

/* Low or possibly low on memory in stage 2 init, create save file, reduce memory settings, and try again */

lowmem:	stop_reason = OutOfMemory (thread_num);
possible_lowmem:
	if (pp1data.state == PP1_STATE_MIDSTAGE) pp1_save (&pp1data);
	if (stop_reason != STOP_OUT_OF_MEM) goto exit;
	pp1_save (&pp1data);
	pp1_cleanup (&pp1data);
	OutputBoth (thread_num, "Memory allocation error.  Trying again using less memory.\n");
	pp1data.pct_mem_to_use *= 0.8;
	goto restart;

/* We've run out of memory.  Print error message and exit. */

oom:	stop_reason = OutOfMemory (thread_num);
	goto exit;

/* Print a message if we found a factor! */

bingo:	if (pp1data.state < PP1_STATE_MIDSTAGE)
		sprintf (buf, "P+1 found a factor in stage #1, B1=%" PRIu64 ".\n", pp1data.B);
	else
		sprintf (buf, "P+1 found a factor in stage #2, B1=%" PRIu64 ", B2=%" PRIu64 ".\n", pp1data.B, pp1data.C);
	OutputBoth (thread_num, buf);

/* Allocate memory for the string representation of the factor and for */
/* a message.  Convert the factor to a string.  Allocate lots of extra space */
/* as formatMsgForResultsFile can append a lot of text. */	

	msglen = factor->sign * 10 + 400;
	str = (char *) malloc (msglen);
	if (str == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	msg = (char *) malloc (msglen);
	if (msg == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	gtoc (factor, str, msglen);

/* Validate the factor we just found */

	if (!testFactor (&pp1data.gwdata, w, factor)) {
		sprintf (msg, "ERROR: Bad factor for %s found: %s\n", gwmodulo_as_string (&pp1data.gwdata), str);
		OutputBoth (thread_num, msg);
		unlinkSaveFiles (&pp1data.write_save_file_state);
		OutputStr (thread_num, "Restarting P+1 from scratch.\n");
		stop_reason = 0;
		goto error_restart;
	}

/* Output the validated factor */

	if (pp1data.state < PP1_STATE_MIDSTAGE)
		sprintf (msg, "%s has a factor: %s (P+1, B1=%" PRIu64 ")\n",
			 gwmodulo_as_string (&pp1data.gwdata), str, pp1data.B);
	else
		sprintf (msg, "%s has a factor: %s (P+1, B1=%" PRIu64 ", B2=%" PRIu64 ")\n",
			 gwmodulo_as_string (&pp1data.gwdata), str, pp1data.B, pp1data.C);
	OutputStr (thread_num, msg);
	formatMsgForResultsFile (msg, w);
	writeResults (msg);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"F", "exponent":45581713, "worktype":"P+1", "factors":["430639100587696027847"], */
/* "b1":50000, "b2":5000000, "start":"2/7", "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"30.6", "build":"1"}, "timestamp":"2021-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"F\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P+1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"factors\":[\"%s\"]", str);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%" PRIu64, pp1data.B);
	if (pp1data.state > PP1_STATE_MIDSTAGE) sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%" PRIu64, pp1data.C);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"start\":\"%" PRIu32 "/%" PRIu32 "\"", pp1data.numerator, pp1data.denominator);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pp1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, pp1data.B, pp1data.C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send assignment result to the server.  To avoid flooding the server with small factors from users needlessly redoing */
/* factoring work, make sure the factor is more than 60 bits or so. */

	if (!QA_IN_PROGRESS && (strlen (str) >= 18 || IniGetInt (INI_FILE, "SendAllFactorData", 0))) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, msg);
		pkt.result_type = PRIMENET_AR_PP1_FACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.pp1_numerator = pp1data.numerator;
		pkt.pp1_denominator = pp1data.denominator;
		truncated_strcpy (pkt.factor, sizeof (pkt.factor), str);
		pkt.B1 = (double) pp1data.B;
		pkt.B2 = (double) (pp1data.state < PP1_STATE_MIDSTAGE ? 0 : pp1data.C);
		pkt.fftlen = gwfftlen (&pp1data.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Free save files (especially big ones created during stage 2) */
/* Then create save file so that we can expand bound 1 or bound 2 at a later date. */

	unlinkSaveFiles (&pp1data.write_save_file_state);
	if (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "KeepPplus1SaveFiles", 0) && w->nth_run <= 2)
		pp1_save (&pp1data);

/* Since we found a factor, then we may have performed less work than */
/* expected.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	invalidateNextRollingAverageUpdate ();

/* Remove the exponent from the worktodo.ini file */

	stop_reason = STOP_WORK_UNIT_COMPLETE;
	goto exit;

/* Output an error message saying we are restarting. */
/* Sleep five minutes before restarting from last save file. */

err:	if (gw_get_maxerr (&pp1data.gwdata) > allowable_maxerr) {
		sprintf (buf, "Possible roundoff error (%.8g), backtracking to last save file.\n", gw_get_maxerr (&pp1data.gwdata));
		OutputStr (thread_num, buf);
	} else {
		OutputBoth (thread_num, "SUMOUT error occurred.\n");
		stop_reason = SleepFive (thread_num);
		if (stop_reason) goto exit;
	}
error_restart:
	pp1_cleanup (&pp1data);
	free (factor), factor = NULL;
	free (str), str = NULL;
	free (msg), msg = NULL;
	goto restart;
}

/* Read a file of P+1 tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9993 to run the QA suite */

int pplus1_QA (
	int	thread_num,
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	FILE	*fd;

/* Set the title */

	title (thread_num, "QA");

/* Open QA file */

	fd = fopen ("qa_pp1", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa_pp1' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	QA_TYPE = 0;
	for ( ; ; ) {
		struct work_unit w;
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		int	stop_reason;

/* Read a line from the file */

		n = 0;
		(void) fscanf (fd, "%lf,%lu,%lu,%ld,%lu,%lu,%lu,%s\n", &k, &b, &n, &c, &B1, &B2_start, &B2_end, fac_str);
		if (n == 0) break;

/* If p is 1, set QA_TYPE */

		if (n == 1) {
			QA_TYPE = c;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		QA_FACTOR = allocgiant ((int) strlen (fac_str));
		ctog (fac_str, QA_FACTOR);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Do the P+1 */

		if (B2_start < B1) B2_start = B1;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_PPLUS1;
		w.k = k;
		w.b = b;
		w.n = n;
		w.c = c;
		w.B1 = B1;
		w.B2_start = B2_start;
		w.B2 = B2_end;
		QA_IN_PROGRESS = TRUE;
		stop_reason = pplus1 (0, sp_info, &w);
		QA_IN_PROGRESS = FALSE;
		free (QA_FACTOR);
		if (stop_reason != STOP_WORK_UNIT_COMPLETE) {
			fclose (fd);
			return (stop_reason);
		}
	}

/* Cleanup */

	fclose (fd);
	return (0);
}
