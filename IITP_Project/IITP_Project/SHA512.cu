#include "type.cuh"
#include "SHA512.cuh"

//CUDA C SHA512 version
__device__ bool KMU_SHA256(uint_8 plaintext[], uint_64 ptLen, uint_8 output[64])
{
	if (ptLen > 56)
		return FAILURE;
	uint_8 buffer[64] = { 0, };
	for (int i = 0; i < ptLen; i++) {
		buffer[i] = plaintext[i];
	}
	buffer[ptLen] = 0x80;
	buffer[62] = (ptLen * 8) / 256;
	buffer[63] = (ptLen * 8) % 256
}

__device__ bool KMU_SHA256_Process(uint_32 word[16],)
{

}
