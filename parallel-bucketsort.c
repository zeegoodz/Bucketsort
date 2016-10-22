// Afton and Zach

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <prand.h>
#include <mpi.h>

/*
 * Sequential bucketsort for randomly generated integers.
 *
 */
void parallelBucketsort(int *, int, int, int, int);

const int DEBUG_LEVEL = DEBUG;
int count_array_grows = 0;

double getMilliSeconds();
int compareTo(const void *, const void *);

const int TRUE = 1;
const int FALSE = 0;

int *A;
int *newA;
int *recvBuff;
int **parBucket; // bucket[i] = array holding the ith bucket
int *parCapacity; // capacity[i] = capacity of ith bucket
int *parSize; // size[i] = next free location in ith bucket
int **bucket;
int *capacity;
int *size;
int myId;
int count;
int *allRecvBuff;


void checkIfSorted(int *array, int n) {
	int i;
	int sorted;

	sorted = TRUE;
	for (i=0; i<n-1; i++)
		if (array[i] > array[i+1]) {
				printf("array[i] %d\narray[i+1] %d\n",array[i],array[i+1]);
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
void generateInput(int *A, int fillN, long int seed, int myId, long long int offset) {
	int i;
    long long int startingPosition = 0;

	srandom(seed);
    startingPosition = myId * offset;
    unrankRand(startingPosition);

	for (i=0; i<fillN; i++) {
		A[i] = random();
	}
}


/*
 * Print the array, one element per line
 */
void printArray(int *A, int n) {
	int i;

	for (i=0; i<n; i++)
		printf(" %16d \n", A[i]);
}


/*
 *
 * Insert a given value into the specified bucket
 */

void insertInBucket(int value, int bucketIndex) {
	int *tmp;

	if (size[bucketIndex] == capacity[bucketIndex]) {
		//grow the bucket array
		tmp = (int *) malloc(sizeof(int)*(2*capacity[bucketIndex]));
		memcpy(tmp, bucket[bucketIndex], capacity[bucketIndex]*sizeof(int));
		free(bucket[bucketIndex]);
		bucket[bucketIndex] = tmp;
		capacity[bucketIndex] = 2 * capacity[bucketIndex];
		count_array_grows++;
		if (DEBUG_LEVEL >= 1) {
	//		fprintf(stderr, "Growing bucket %d from %d to %d elements\n",
	//				bucketIndex, capacity[bucketIndex]/2, capacity[bucketIndex]);
		}
	}

	bucket[bucketIndex][size[bucketIndex]] = value;
	size[bucketIndex]++;
}

void parInsertInBucket(int value, int bucketIndex) {
	int *tmp;

	if (parSize[bucketIndex] == parCapacity[bucketIndex]) {
		//grow the bucket array
		tmp = (int *) malloc(sizeof(int)*(2*parCapacity[bucketIndex]));
		memcpy(tmp, parBucket[bucketIndex], parCapacity[bucketIndex]*sizeof(int));
		free(parBucket[bucketIndex]);
		parBucket[bucketIndex] = tmp;
		parCapacity[bucketIndex] = 2 * parCapacity[bucketIndex];
		count_array_grows++;
		if (DEBUG_LEVEL >= 1) {
//			fprintf(stderr, "Growing bucket %d from %d to %d elements\n",
//					bucketIndex, parCapacity[bucketIndex]/2, parCapacity[bucketIndex]);
		}
	}

	parBucket[bucketIndex][parSize[bucketIndex]] = value;
	parSize[bucketIndex]++;

}
/*
 * compareTo function for using qsort
 * returns  -ve if *x < *y, 0 if *x == *y, +ve if *x > *y
 */
int compareTo(const void *x, const void *y) {
		return ((*(int *)x) - (*(int *)y));
}


/*
 * Sort indiviual bucket using quick sort from std C library
 */

void sortEachBucket(int numBuckets) {
	int i;

	for (i=0; i<numBuckets; i++) {
		qsort(bucket[i], size[i], sizeof(int), compareTo);
	}
	if (DEBUG_LEVEL >= 2)
		for (i=0; i<numBuckets; i++) {
	//		fprintf(stderr, "bucket %d has %d elements\n", i, size[i]);
		}
}

/*
 * Combine all buckets back into the original array to finish the sorting
 *
 */
void combineBuckets(int *A, int n, int numBuckets) {
	int i;

	int start = 0;
	for (i=0; i<numBuckets; i++) {
		memcpy(A+start, bucket[i], sizeof(int)*size[i]);
		start = start + size[i];
	//	free(bucket[i]);
	}
//	free(bucket);
}

void parCombineBuckets(int *A, int n, int numBuckets) {
	int i;

	int start = 0;
	for (i=0; i<numBuckets; i++) {
		memcpy(A+start, parBucket[i], sizeof(int)*parSize[i]);
		start = start + parSize[i];
	//	free(parBucket[i]);
	}
//	free(parBucket);
}

/*
 * Use bucketsort to sort n uniformly distributed numbers in the range [0..2^31-1].
 * Input: int *A: array of ints A[0..n-1]
 *        int n: number of elements in the input array
 *        int numBuckets: number of buckets to use
 *
 */

void sequentialBucketsort(int *A, int n, int numBuckets, int myid) {
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
	for (i=0; i<n; i++) {
		bucketIndex = recvBuff[i]/bucketRange;
		if (bucketIndex > numBuckets - 1)
				bucketIndex = numBuckets - 1;
		insertInBucket(recvBuff[i], bucketIndex);
	}

	sortEachBucket(numBuckets);
	combineBuckets(A, n, numBuckets);
//	free(capacity);
//	free(size);
}

void parallelBucketsort(int *A, int myN, int numBuckets, int numProcs, int myId) {
	int share;
	int i;
	int bucketRange;
	int bucketIndex;
    int *sdispls;
    int *rdispls;

	share = myN / numBuckets;
	share = share + (share * 11)/100; // 11% extra for overflow

    sdispls = (int *) malloc(sizeof(int)*numBuckets);
    rdispls = (int *) malloc(sizeof(int)*numBuckets);
    allRecvBuff = (int *) malloc(sizeof(int)*numBuckets);
	parCapacity = (int *) malloc(sizeof(int)*numBuckets);
	parSize = (int *) malloc(sizeof(int)*numBuckets);
	parBucket = (int **) malloc(sizeof(int *)* numBuckets);

	for (i=0; i<numBuckets; i++) {
		parBucket[i] = (int *) malloc(sizeof(int)*share);
		parCapacity[i] = share;
		parSize[i] = 0;
		allRecvBuff[i] = 0;
	}

	bucketRange = RAND_MAX/numBuckets;
	for (i=0; i<myN; i++) {
		bucketIndex = A[i]/bucketRange;
		if (bucketIndex > numBuckets - 1) {
				bucketIndex = numBuckets - 1;
        }
		parInsertInBucket(A[i], bucketIndex);
	}

	parCombineBuckets(A, myN, numBuckets);
    sdispls[0] = 0;

    for (i=0;i<numBuckets;i++) {
        sdispls[i+1] = sdispls[i] + parSize[i];
    }

	int count = numBuckets/numProcs;
	MPI_Alltoall(parSize, count, MPI_INT, allRecvBuff, count, MPI_INT, MPI_COMM_WORLD);


	for(int i = 0; i<numBuckets; i++){
		rdispls[i]=0;
	}
	for (int i=0;i<(numBuckets);i++) {
        rdispls[i+1] = rdispls[i] + allRecvBuff[i];
    }

	if(myId == 1){

	}

	MPI_Alltoallv(A, parSize, sdispls, MPI_INT, recvBuff, allRecvBuff, rdispls, MPI_INT, MPI_COMM_WORLD);

	for(int i = 0; i<numProcs; i++){
		count += allRecvBuff[i];
	}
}

void print_usage(char * program) {
	fprintf(stderr, "Usage %s <n, must be > 1> <#buckets, must between 1 and n> <random seed>\n",
           program);
}

void checkIfSortedParallel(int myId, int numProcs, int count, int *recvBuff){

		int sendVal, recvVal;

		sendVal = recvBuff[count-1];

		if(myId == 0){
			MPI_Send(&sendVal, 1, MPI_INT, 1, myId, MPI_COMM_WORLD);
		}else if (myId == numProcs-1){
			MPI_Recv(&recvVal, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}else{
			MPI_Send(&sendVal, 1, MPI_INT, (myId+1), myId, MPI_COMM_WORLD);
			MPI_Recv(&recvVal, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}

		if(myId != 0){
			if(recvVal<recvBuff[count-1]){
				printf("All sorted values in process %d are leff than all sorted values of process %d.\n",myId-1, myId);
			}else{
				fprintf(stderr, "ERROR: Arrays are not sorted properly\n");
			}

		}
}

int main(int argc, char **argv) {
	int n, numProcs, myN;
	int numBuckets;
    int myNumBuckets;
    int fillN;
    long long int offset;
	unsigned int seed;
	double startTime=0;
	double totalTime=0;

	if (argc != 4) {
		print_usage(argv[0]);
		exit(1);
	}

	numBuckets = atoi(argv[2]);
	seed = atoi(argv[3]);
	n = atoi(argv[1]);

	if ((numBuckets < 1) || (n < 1) || (n < numBuckets)) {
		print_usage(argv[0]);
		exit(1);
	}

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    myN = n / numProcs;
   if(myId == numProcs-1){
		myN = myN + (n % numProcs);
   }
    fillN = n / numProcs;
    offset = (long long int)fillN;
    myNumBuckets = numProcs;


	newA = (int *) malloc(numProcs*sizeof(int) * n);
	A = (int *) malloc(sizeof(int) * myN);
    recvBuff = (int *) malloc(sizeof(int) * myN);

	if(myId==0){
	startTime = MPI_Wtime();
	}
    if (myId == numProcs - 1) {
	    generateInput(A, myN, seed, myId, offset);
    } else {
        generateInput(A, fillN,seed, myId, offset);
    }

 	if (DEBUG_LEVEL >= 3)
		printArray(A, myN);

	parallelBucketsort(A, myN, myNumBuckets, numProcs, myId);
	int newCount=0;

	for(int i=0;i<numProcs; i++){
		newCount += allRecvBuff[i];
	}
	sequentialBucketsort(recvBuff, newCount, numBuckets, myId);
	checkIfSorted(recvBuff, newCount);
	checkIfSortedParallel(myId, numProcs, newCount, recvBuff);
	checkIfSorted(newA, n);

	if (DEBUG_LEVEL >= 1)
		printf("Number of array grows is %d\n", count_array_grows);
 	if (DEBUG_LEVEL >= 3)
 		printArray(A,n);

	MPI_Barrier(MPI_COMM_WORLD);
	if(myId==0){

	totalTime = MPI_Wtime() - startTime;
	printf("bucketsort: n = %d  m = %d buckets seed = %d time = =%.2f seconds\n",
		   n, numBuckets, seed,	totalTime);
	}
//	free(A);
 //   free(recvBuff);
    MPI_Finalize();
	exit(0);
}

/* vim: set ts=4: */
