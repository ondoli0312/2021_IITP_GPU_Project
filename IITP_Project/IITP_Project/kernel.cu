#include "type.cuh"
#include "SHA512.cuh"

//int main()
//{
//	cudaEvent_t start, stop;
//	uint8_t pt[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
//	uint8_t salt[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
//	uint64_t DK[8];
//	uint8_t* cu_salt = NULL;
//	uint64_t* cu_DK = NULL;
//	float elapsed_time_ms = 0.0f;
//
//	cudaMalloc((void**)&cu_salt, 8);
//	cudaMemcpy(cu_salt, salt, 8 * sizeof(uint64_t), cudaMemcpyHostToDevice);
//	cudaMalloc((void**)&cu_DK, 8 * 8);
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);
//	cudaEventRecord(start, 0);
//	for (int i = 0; i < 1; i++) {
//		KMU_PBKDF2v2 << <BLOCKFULL, THREADFULL >> > (129973, cu_salt, 8, cu_DK);
//	}
//	cudaEventRecord(stop, 0);
//	cudaDeviceSynchronize();
//	cudaEventSynchronize(start);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
//	printf("Performance : %4.2f GB/s\n", elapsed_time_ms);
//}