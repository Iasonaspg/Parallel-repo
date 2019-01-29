#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cuda.h>
#include <curand.h>
#include <math.h>


#define CHECK(call) \
{                    \
    const cudaError_t error = call; \
    if (error != cudaSuccess){       \
        printf("Error: %s:%d, ", __FILE__, __LINE__);  \
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
        exit(1); \
    } \
} \

// Kernel to change 1 with 0 because we want values inside [0,1) and currand produces (0,1]
__global__ void corr(float* ptr, size_t N){

    int index = threadIdx.x + blockDim.x*blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i=index; i<N; i+=stride){
        if (ptr[i] == 1){
            ptr[i] = 0;
        }
    }
}

// Kernel to print points from a big array
__global__ void check(float* ptr, size_t N, int len){

    int index = threadIdx.x + blockDim.x*blockIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i=index; i<len; i+=stride){
        printf("Point[%d]: (%f,%f,%f)\n",i,ptr[i],ptr[i+N],ptr[i+2*N]);
    }

}

// Kernel to initialize a big array with a specific number. I eventually didn't use it
__global__ void init(float* ptr, int len, int num){

    int index = threadIdx.x + blockDim.x*blockIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i=index; i<len; i+=stride){
        ptr[i] = num;
    }
}

// Kernel to find the nearest neighbor of a point
__global__ void findNN(float* queries, int len, float* boxesC, int len1, int elements, float side, int* aux, int d, float* distances){

    int index = threadIdx.x + blockDim.x*blockIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i=index; i<len; i+=stride){
        int coX,coY,coZ;
        int coXtmp,coYtmp,coZtmp;
        int start,boxMin;
        float distance, sum;
        coX = queries[i]/side;
        coY = queries[i+len]/side;
        coZ = queries[i+2*len]/side;
        int box = coX + coY*d + coZ*d*d;    // calculate the box in which the point belongs to
        // printf("Box inside kernel: %d\n",box);
        // printf("Queries inside kernel: (%f,%f,%f)\n",queries[i],queries[i+len],queries[i+2*len]);
        float min = 10000000;
        int pos = -1;
        boxMin = -1;
        start = box*elements;   // Starting point inside the big 1-d array for the desired box
        for (int j=0; j<aux[box];j++){
            sum = 0;
            for (int k=0; k<3; k++){
                sum += powf(queries[i+k*len]-boxesC[start+j+k*len1],2);
            }
            distance = powf(sum,0.5);
            //printf("Distance: %f\n", distance);
            if (distance < min){
                min = distance;     // Calculate the min distance inside the box
                pos = j;
                boxMin = box;
            }
        }

        // Check if I have to search adjacent boxes
        if ((fabsf(queries[i+2*len]-(coZ+1)*side)<min) || (fabsf(queries[i+2*len]-coZ*side)<min) || (fabsf(queries[i+len]-(coY+1)*side)<min) || (fabsf(queries[i+len]-coY*side)<min) || (fabsf(queries[i]-(coX+1)*side)<min) || (fabsf(queries[i]-coX*side)<min)){
            for (int dx=-1;dx<=1;dx++){
                coXtmp = coX + dx;
                if (coXtmp>=d || coXtmp<0){
                    continue;
                }
                for (int dy=-1;dy<=1;dy++){
                    coYtmp = coY + dy;
                    if (coYtmp>=d || coYtmp<0){
                        continue;
                    }
                    for (int dz=-1;dz<=1;dz++){
                        coZtmp = coZ + dz;
                        if (coZtmp>=d || coZtmp<0){
                            continue;
                        }
                        if (dz==0 && dx==0 && dy==0){
                            continue;
                        }
                        box = coXtmp + coYtmp*d + coZtmp*d*d;
                        // printf("Box number: %d\n",box);
                        start = box*elements;
                        for (int j=0; j<aux[box];j++){
                            sum = 0;
                            for (int k=0; k<3; k++){
                                sum += powf(queries[i+k*len]-boxesC[start+j+k*len1],2);
                            }
                            distance = powf(sum,0.5);
                            //printf("Distance: %f\n", distance);
                            if (distance < min){
                                min = distance;
                                pos = j;
                                boxMin = box;
                            }
                        }
                    }
                }
            }
        } 
        // printf("Min distance from (%f,%f,%f): %f at box: %d and pos: %d\n",queries[i],queries[i+len],queries[i+2*len],min,boxMin,pos);
        int start1 = boxMin*elements;
        // printf("Nearest neighbor of (%f,%f,%f): (%f,%f,%f)\n",queries[i],queries[i+len],queries[i+2*len],boxesC[start1+pos],boxesC[start1+pos+len1],boxesC[start1+pos+2*len1]);
        distances[i] = min;  // Store the final min distance
    }
}

// Kernel that calculates minimum distance in a brute-force way
__global__ void validate(float* queries, float* corpus, int len, int* flag, float* distances){

    int index = threadIdx.x + blockDim.x*blockIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int j=index;j<len;j+=stride){
        float min = 10000;
        for (int i=0;i<len;i++){
            float sum = 0;
            for (int k=0; k<3; k++){
                sum += powf(queries[j+k*len]-corpus[i+k*len],2);
            }
            float distance = powf(sum,0.5);
            if (distance < min){
                min = distance;
            }
        }
        // printf("Minimum distance for query[%d]: %f\n",j,min);
        if (min != distances[j]) {
            *flag = 1;
        }
    }

}

// Function that returns current time in seconds
double cpuSecond() {
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}



int main(int argc, char** argv){

    size_t Nq = (1<<atoi(argv[1]));     // Nq=Nc
    int d = (1<<atoi(argv[2]));
    float side = 1/((float) d);
    printf("Side length of each box: %f\n",side);
    
    
    curandGenerator_t gen;
    float *queries, *corpus;

    // Generate Uniform Points
    CHECK(cudaMallocManaged(&queries, 3*Nq*sizeof(float)));
    CHECK(cudaMallocManaged(&corpus, 3*Nq*sizeof(float)));
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen,time(NULL));
    
    double start = cpuSecond();
    curandGenerateUniform(gen, queries, 3*Nq);
    curandGenerateUniform(gen, corpus, 3*Nq);
    corr<<<((Nq+1023)/1024),1024>>>(queries,3*Nq);
    corr<<<((Nq+1023)/1024),1024>>>(corpus,3*Nq);
    CHECK(cudaDeviceSynchronize());
    double elaps = cpuSecond() - start;
    printf("Time elapsed for generating points: %f seconds\n",elaps);
    //check<<<1,32>>>(queries,Nq,Nq);
    //check<<<1,32>>>(corpus,Nq,Nq);


    int newNq = 8*Nq;
    float* boxesC;
    CHECK(cudaMallocManaged(&boxesC, 3*newNq*sizeof(float)));   // Array that is going to store the points
    
    int* auxC;                                                  // Array that will store the number of points inside each box
    CHECK(cudaMallocManaged(&auxC, d*d*d*sizeof(int)));
    CHECK(cudaMemset(auxC,0,d*d*d*sizeof(int)));
    

    int elements = newNq/(d*d*d);       // Number of points that a box can contain
    printf("Elements: %d\n",elements);

    start = cpuSecond();
    for (int i=0; i<Nq; i++){
        int coX_C,coY_C,coZ_C;
        coX_C = corpus[i]/side;
        coY_C = corpus[i+Nq]/side;
        coZ_C = corpus[i+2*Nq]/side;
        int boxC = coX_C + coY_C*d + coZ_C*d*d;     // Find the corresponding box of the point
        if ((boxC >= d*d*d) || (boxC < 0)){
            printf("Error, non existing box\n");
            exit(1);
        }
        boxesC[boxC*elements + auxC[boxC]] = corpus[i];     // --> Store the point inside the box 
        boxesC[boxC*elements + auxC[boxC] + newNq] = corpus[i + Nq];
        boxesC[boxC*elements + auxC[boxC] + 2*newNq] = corpus[i + 2*Nq];  // <--
        auxC[boxC]++;
        if (auxC[boxC]>=elements){
            printf("We have problem\n");
            exit(1);
        }    
    }
    elaps = cpuSecond() - start;
    printf("Time elapsed for hashing the total number of points on host: %f seconds\n",elaps);

    // for (int i=0; i<newNq; i++){
    //     printf("BoxesC[%d]: (%f,%f,%f)\n",i, boxesC[i], boxesC[i+newNq], boxesC[i+2*newNq]);
    // }

    // for (int i=0; i<d*d*d; i++){
    //     printf("AuxC: %d\n",auxC[i]);
    // }

    float* distances;
    CHECK(cudaMallocManaged(&distances,Nq*sizeof(float)));

    start = cpuSecond();
    findNN<<<((Nq+1023)/1024),1024>>>(queries,Nq,boxesC,newNq,elements,side,auxC,d,distances);
    // findNN<<<1,1>>>(queries,Nq,boxesC,newNq,elements,side,auxC,d);
    CHECK(cudaDeviceSynchronize());
    elaps = cpuSecond() - start;
    printf("Time elapsed for knn on device: %f seconds\n",elaps);

    int* flag;
    CHECK(cudaMallocManaged(&flag,sizeof(int)));

    start = cpuSecond();
    validate<<<((Nq+1023)/1024),1024>>>(queries,corpus,Nq,flag,distances);
    CHECK(cudaDeviceSynchronize());
    elaps = cpuSecond() - start;
    printf("Time elapsed for brute force knn on device: %f seconds\n",elaps);

    if (*flag != 1){
        printf("Validation PASSED\n");
    }
    else{
        printf("Validation FAILED\n");
    }

    if (Nq < (1<<16)){
        start = cpuSecond();
        for (int j=0;j<Nq;j++){
            float min = 10000;
            for (int i=0;i<Nq;i++){
                float sum = 0;
                for (int k=0; k<3; k++){
                    sum += powf(queries[j+k*Nq]-corpus[i+k*Nq],2);
                }
                float distance = powf(sum,0.5);
                if (distance < min){
                    min = distance;
                }
            }
            // printf("Minimum distance for query[%d]: %f\n",j,min);
            if (fabsf(min-distances[j]) > 0.000001){
                printf("Validation failed\n");
                exit(1);
            }
        }
        elaps = cpuSecond() - start;
        printf("Time elapsed for brute force knn on host: %f seconds\n",elaps);
        printf("Validation PASSED\n");
    }
   


    cudaFree(boxesC);
    cudaFree(auxC);
    cudaFree(queries);
    cudaFree(corpus);
    cudaFree(distances);
    curandDestroyGenerator(gen);
    return 0;
}
