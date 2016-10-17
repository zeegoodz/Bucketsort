
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
/*
 * Sequential bucketsort for randomly generated integers.
 *
 */

const int DEBUG_LEVEL = DEBUG;
int count_array_grows = 0;

double getMilliSeconds();
int compareTo(const void *, const void *);

const int TRUE = 1;
const int FALSE = 0;

int *A;
int **bucket; // bucket[i] = array holding the ith bucket
int *capacity; // capacity[i] = capacity of ith bucket
int *size; // size[i] = next free location in ith bucket




void checkIfSorted(int *array, int n)
{
	int i;
	int sorted;

	sorted = TRUE;
	for (i=0; i<n-1; i++)
		if (array[i] > array[i+1]) {
				sorted = FALSE;
				break;
		}

	if (sorted) {
		if (DEBUG_LEVEL >= 1)
			fprintf(stderr, "array is sorted\n");
	} else {
			fprintf(stderr, "Error: array is not sorted!\n");
	}
}



/*
 * Generate n numbers using the given seed
 */
void generateInput(int *A, int n, long int seed)
{
	int i;

	srandom(seed);
	for (i=0; i<n; i++) {
		A[i] = random();
	}
}


/*
 * Print the array, one element per line
 */
void printArray(int *A, int n)
{
	int i;

	for (i=0; i<n; i++)
		printf(" %16d \n", A[i]);
}


/*
 *
 * Insert a given value into the specified bucket
 */

void insertInBucket(int value, int bucketIndex)
{
	int *tmp;

	if (size[bucketIndex] == capacity[bucketIndex])
	{
		//grow the bucket array
		tmp = (int *) malloc(sizeof(int)*(2*capacity[bucketIndex]));
		memcpy(tmp, bucket[bucketIndex], capacity[bucketIndex]*sizeof(int));
		free(bucket[bucketIndex]);
		bucket[bucketIndex] = tmp;
		capacity[bucketIndex] = 2 * capacity[bucketIndex];
		count_array_grows++;
		if (DEBUG_LEVEL >= 1) {
			fprintf(stderr, "Growing bucket %d from %d to %d elements\n",
					bucketIndex, capacity[bucketIndex]/2, capacity[bucketIndex]);
		}
	}

	bucket[bucketIndex][size[bucketIndex]] = value;
	size[bucketIndex]++;
}


/*
 * compareTo function for using qsort
 * returns  -ve if *x < *y, 0 if *x == *y, +ve if *x > *y
 */
int compareTo(const void *x, const void *y)
{
		return ((*(int *)x) - (*(int *)y));
}


/*
 * Sort indiviual bucket using quick sort from std C library
 */

void sortEachBucket(int numBuckets)
{
	int i;

	for (i=0; i<numBuckets; i++)
	{
		qsort(bucket[i], size[i], sizeof(int), compareTo);
	}
	if (DEBUG_LEVEL >= 2)
		for (i=0; i<numBuckets; i++) {
			fprintf(stderr, "bucket %d has %d elements\n", i, size[i]);
		}
}


/*
 * Combine all buckets back into the original array to finish the sorting
 *
 */
void combineBuckets(int *A, int n, int numBuckets)
{
	int i;

	int start = 0;
	for (i=0; i<numBuckets; i++) {
		memcpy(A+start, bucket[i], sizeof(int)*size[i]);
		start = start + size[i];
		free(bucket[i]);
	}
	free(bucket);
}

/*
 * Use bucketsort to sort n uniformly distributed numbers in the range [0..2^31-1].
 * Input: int *A: array of ints A[0..n-1]
 *        int n: number of elements in the input array
 *        int numBuckets: number of buckets to use
 *
 */

void sequentialBucketsort(int *A, int n, int numBuckets)
{
	int share;
	int i;
	int bucketRange;
	int bucketIndex;

	share = n / numBuckets;
	share = share + (share * 11)/100; // 11% extra for overflow

	capacity = (int *) malloc(sizeof(int)*numBuckets);
	size = (int *) malloc(sizeof(int)*numBuckets);
	bucket = (int **) malloc(sizeof(int *)* numBuckets);

	for (i=0; i<numBuckets; i++) {
		bucket[i] = (int *) malloc(sizeof(int)*share);
		capacity[i] = share;
		size[i] = 0;
	}

	bucketRange = RAND_MAX/numBuckets;
	for (i=0; i<n; i++)
	{
		bucketIndex = A[i]/bucketRange;
		if (bucketIndex > numBuckets - 1)
				bucketIndex = numBuckets - 1;
		insertInBucket(A[i], bucketIndex);
	}

	sortEachBucket(numBuckets);
	combineBuckets(A, n, numBuckets);
	free(capacity);
	free(size);
}


void print_usage(char * program)
{
		fprintf(stderr, "Usage %s <n, must be > 1> <#buckets, must between 1 and n> <random seed>\n", program);
}



int main(int argc, char **argv)
{
	int n, myid, size, elementsPerProc;
	int numBuckets;
	unsigned int seed;
	double startTime;
	double totalTime;

	if (argc != 4) {
		print_usage(argv[0]);
		exit(1);
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	n = atoi(argv[1]);
	numBuckets = atoi(argv[2]);
	seed = atoi(argv[3]);
	elementsPerProc = (n/size);

	if ((numBuckets < 1) || (n < 1) || (n < numBuckets)) {
		print_usage(argv[0]);
		exit(1);
	}

	A = (int *) malloc(sizeof(int) * n);
	int subA = (int *) malloc(sizeof(int) * (n/size));

	if (myid == 0){
		generateInput(A, n, seed);
 	if (DEBUG_LEVEL >= 3)
		printArray(A,n);
	}

	startTime = getMilliSeconds();


	MPI_Scatter(A, elementsPerProc, MPI_INT, subA, elementsPerProc, MPI_INT, 0, MPI_COMM_WORLD);
	sequentialBucketsort(subA, elementsPerProc, numBuckets);

	int *sortedA = (int *) malloc(sizeof(int)*size);

	MPI_Allgather(&subA, 1, MPI_INT, sortedA, 1, MPI_INT, MPI_COMM_WORLD);
	sequentialBucketsort(sortedA, n, numBuckets);


	totalTime = getMilliSeconds() - startTime;

	checkIfSorted(A, n);

	if (DEBUG_LEVEL >= 1)
		printf("Number of array grows is %d\n", count_array_grows);
 	if (DEBUG_LEVEL >= 3)
 		printArray(A,n);
	printf("bucketsort: n = %d  m = %d buckets seed = %d time = %lf seconds\n",
		   n, numBuckets, seed,	totalTime/1000.0);

	free(A);
	exit(0);
}

/* vim: set ts=4: */
