#include "type.cuh"
#include "SHA512.cuh"
//CPU part

__constant__ uint64_t cont_512[80] =
{
	0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
	0x3956c25bf348b538, 0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
	0xd807aa98a3030242, 0x12835b0145706fbe, 0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
	0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235, 0xc19bf174cf692694,
	0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
	0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
	0x983e5152ee66dfab, 0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4,
	0xc6e00bf33da88fc2, 0xd5a79147930aa725, 0x06ca6351e003826f, 0x142929670a0e6e70,
	0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
	0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
	0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30,
	0xd192e819d6ef5218, 0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8,
	0x19a4c116b8d2d0c8, 0x1e376c085141ab53, 0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
	0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
	0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
	0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915,0xc67178f2e372532b,
	0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
	0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
	0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
	0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817
};

__device__ void KMU_PBKDF_SHA_BLOCK(const uint64_t* PT, uint64_t* output)
{
	volatile uint64_t w0_t = (PT[0]);
	volatile uint64_t w1_t = (PT[1]);
	volatile uint64_t w2_t = (PT[2]);
	volatile uint64_t w3_t = (PT[3]);
	volatile uint64_t w4_t = (PT[4]);
	volatile uint64_t w5_t = (PT[5]);
	volatile uint64_t w6_t = (PT[6]);
	volatile uint64_t w7_t = (PT[7]);
	volatile uint64_t w8_t = (PT[8]);
	volatile uint64_t w9_t = (PT[9]);
	volatile uint64_t wa_t = (PT[10]);
	volatile uint64_t wb_t = (PT[11]);
	volatile uint64_t wc_t = (PT[12]);
	volatile uint64_t wd_t = (PT[13]);
	volatile uint64_t we_t = (PT[14]);
	volatile uint64_t wf_t = (PT[15]);

	uint64_t a, b, c, d, e, f, g, h = 0;
	a = 0x6a09e667f3bcc908;
	b = 0xbb67ae8584caa73b;
	c = 0x3c6ef372fe94f82b;
	d = 0xa54ff53a5f1d36f1;
	e = 0x510e527fade682d1;
	f = 0x9b05688c2b3e6c1f;
	g = 0x1f83d9abfb41bd6b;
	h = 0x5be0cd19137e2179;

	#define ROUND_EXPAND()								\
	{													\
		w0_t = SHA512_EXPAND (we_t, w9_t, w1_t, w0_t);  \
		w1_t = SHA512_EXPAND (wf_t, wa_t, w2_t, w1_t);  \
		w2_t = SHA512_EXPAND (w0_t, wb_t, w3_t, w2_t);  \
		w3_t = SHA512_EXPAND (w1_t, wc_t, w4_t, w3_t);  \
		w4_t = SHA512_EXPAND (w2_t, wd_t, w5_t, w4_t);  \
		w5_t = SHA512_EXPAND (w3_t, we_t, w6_t, w5_t);  \
		w6_t = SHA512_EXPAND (w4_t, wf_t, w7_t, w6_t);  \
		w7_t = SHA512_EXPAND (w5_t, w0_t, w8_t, w7_t);  \
		w8_t = SHA512_EXPAND (w6_t, w1_t, w9_t, w8_t);  \
		w9_t = SHA512_EXPAND (w7_t, w2_t, wa_t, w9_t);  \
		wa_t = SHA512_EXPAND (w8_t, w3_t, wb_t, wa_t);  \
		wb_t = SHA512_EXPAND (w9_t, w4_t, wc_t, wb_t);  \
		wc_t = SHA512_EXPAND (wa_t, w5_t, wd_t, wc_t);  \
		wd_t = SHA512_EXPAND (wb_t, w6_t, we_t, wd_t);  \
		we_t = SHA512_EXPAND (wc_t, w7_t, wf_t, we_t);  \
		wf_t = SHA512_EXPAND (wd_t, w8_t, w0_t, wf_t);  \
	}

	#define ROUND_STEP(i)																	\
	{																						\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512[i +  0]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512[i +  1]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512[i +  2]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512[i +  3]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512[i +  4]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512[i +  5]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512[i +  6]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512[i +  7]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512[i +  8]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512[i +  9]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512[i + 10]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512[i + 11]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512[i + 12]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512[i + 13]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512[i + 14]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512[i + 15]);	\
	}
	ROUND_STEP(0);
	for (int i = 16; i < 80; i += 16) {
		ROUND_EXPAND(); 
		ROUND_STEP(i);
	}
	output[0] = a + 0x6a09e667f3bcc908;
	output[1] = b + 0xbb67ae8584caa73b;
	output[2] = c + 0x3c6ef372fe94f82b;
	output[3] = d + 0xa54ff53a5f1d36f1;
	output[4] = e + 0x510e527fade682d1;
	output[5] = f + 0x9b05688c2b3e6c1f;
	output[6] = g + 0x1f83d9abfb41bd6b;
	output[7] = h + 0x5be0cd19137e2179;
}

__device__ void KMU_PBKDF2_BLOCK_CORE(uint64_t* PRE_IN, uint64_t* hash, uint64_t ptLen, uint64_t* OUT)
{
	uint64_t a, b, c, d, e, f, g, h = 0;
	a = PRE_IN[0];
	b = PRE_IN[1];
	c = PRE_IN[2];
	d = PRE_IN[3];
	e = PRE_IN[4];
	f = PRE_IN[5];
	g = PRE_IN[6];
	h = PRE_IN[7];

	//Padding
	volatile uint64_t w0_t = hash[0];
	volatile uint64_t w1_t = hash[1];
	volatile uint64_t w2_t = hash[2];
	volatile uint64_t w3_t = hash[3];
	volatile uint64_t w4_t = hash[4];
	volatile uint64_t w5_t = hash[5];
	volatile uint64_t w6_t = hash[6];
	volatile uint64_t w7_t = hash[7];
	volatile uint64_t w8_t = 0x8000000000000000;
	volatile uint64_t w9_t = 0;
	volatile uint64_t wa_t = 0;
	volatile uint64_t wb_t = 0;
	volatile uint64_t wc_t = 0;
	volatile uint64_t wd_t = 0;
	volatile uint64_t we_t = 0;
	volatile uint64_t wf_t = ptLen  << 3;
#define ROUND_EXPAND()									\
	{													\
		w0_t = SHA512_EXPAND (we_t, w9_t, w1_t, w0_t);  \
		w1_t = SHA512_EXPAND (wf_t, wa_t, w2_t, w1_t);  \
		w2_t = SHA512_EXPAND (w0_t, wb_t, w3_t, w2_t);  \
		w3_t = SHA512_EXPAND (w1_t, wc_t, w4_t, w3_t);  \
		w4_t = SHA512_EXPAND (w2_t, wd_t, w5_t, w4_t);  \
		w5_t = SHA512_EXPAND (w3_t, we_t, w6_t, w5_t);  \
		w6_t = SHA512_EXPAND (w4_t, wf_t, w7_t, w6_t);  \
		w7_t = SHA512_EXPAND (w5_t, w0_t, w8_t, w7_t);  \
		w8_t = SHA512_EXPAND (w6_t, w1_t, w9_t, w8_t);  \
		w9_t = SHA512_EXPAND (w7_t, w2_t, wa_t, w9_t);  \
		wa_t = SHA512_EXPAND (w8_t, w3_t, wb_t, wa_t);  \
		wb_t = SHA512_EXPAND (w9_t, w4_t, wc_t, wb_t);  \
		wc_t = SHA512_EXPAND (wa_t, w5_t, wd_t, wc_t);  \
		wd_t = SHA512_EXPAND (wb_t, w6_t, we_t, wd_t);  \
		we_t = SHA512_EXPAND (wc_t, w7_t, wf_t, we_t);  \
		wf_t = SHA512_EXPAND (wd_t, w8_t, w0_t, wf_t);  \
	}

#define ROUND_STEP(i)																		\
	{																						\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512[i +  0]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512[i +  1]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512[i +  2]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512[i +  3]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512[i +  4]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512[i +  5]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512[i +  6]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512[i +  7]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512[i +  8]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512[i +  9]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512[i + 10]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512[i + 11]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512[i + 12]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512[i + 13]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512[i + 14]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512[i + 15]);	\
	}
	ROUND_STEP(0);
	for (int i = 16; i < 80; i += 16) {
		ROUND_EXPAND();
		ROUND_STEP(i);
	}

	OUT[0] = a + PRE_IN[0];
	OUT[1] = b + PRE_IN[1];
	OUT[2] = c + PRE_IN[2];
	OUT[3] = d + PRE_IN[3];
	OUT[4] = e + PRE_IN[4];
	OUT[5] = f + PRE_IN[5];
	OUT[6] = g + PRE_IN[6];
	OUT[7] = h + PRE_IN[7];
}

__device__ void KMU_PRE_HMAC_SHA512(uint8_t* pt, uint64_t ptLen, uint64_t* OPAD_out, uint64_t* IPAD_out)
{
	uint8_t IPAD[128] = { 0x36, };
	uint8_t OPAD[128] = { 0x5c, };
	for (int i = 0; i < 128; i++) {
		IPAD[i] = 0x36;
		OPAD[i] = 0x5c;
	}
	uint64_t buffer[16];
	for (int i = 0; i < ptLen; i++) {
		IPAD[i] = IPAD[i] ^ pt[i];
		OPAD[i] = OPAD[i] ^ pt[i];
	}
	for (int i = 0; i < 16; i++) {
		*(uint64_t*)(IPAD + (i << 3)) = ENDIAN_CHANGE(*(uint64_t*)(IPAD + (i << 3)));
		*(uint64_t*)(OPAD + (i << 3)) = ENDIAN_CHANGE(*(uint64_t*)(IPAD + (i << 3)));
	}
	KMU_PBKDF_SHA_BLOCK((uint64_t*)IPAD, IPAD_out);
	KMU_PBKDF_SHA_BLOCK((uint64_t*)OPAD, OPAD_out);


}

__device__ void KMU_PBKDF2_CORE(uint8_t* pt, uint64_t ptLen, uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* DK) {
	uint64_t dk[8] = { 0, };
	uint64_t IPAD_out[8] = { 0, };
	uint64_t OPAD_out[8] = { 0, };
	uint64_t out[8] = { 0, };
	uint64_t Buffer[8] = { 0, };

	KMU_PRE_HMAC_SHA512(pt, ptLen, IPAD_out, OPAD_out);

	//U0 operation
	uint8_t buffer[64] = { 0, };
	for (int i = 0; i < saLen; i++) {
		buffer[i] = salt[i];
	}
	
	buffer[saLen + 0] = (1 << 24) & 0xff;
	buffer[saLen + 1] = (1 << 16) & 0xff;
	buffer[saLen + 2] = (1 << 8) & 0xff;
	buffer[saLen + 3] = (1) & 0xff;

	for (int i = 0; i < 8; i++) {
		*(uint64_t*)(buffer + (i << 3)) = ENDIAN_CHANGE(*(uint64_t*)(buffer + (i << 3)));
	}

	KMU_PBKDF2_BLOCK_CORE(IPAD_out, (uint64_t*)buffer, 64 + saLen + 4, dk);
	KMU_PBKDF2_BLOCK_CORE(OPAD_out, dk, 64, Buffer);

	for (int i = 1; i < iteration; i++) {
		KMU_PBKDF2_BLOCK_CORE(IPAD_out, Buffer, 128, dk);
		KMU_PBKDF2_BLOCK_CORE(OPAD_out, dk, 128, Buffer);
		for (int i = 0; i < 8; i++)
			out[i] ^= Buffer[i];
	}
	
	for (int i = 0; i < 8; i++) {
		DK[i] = out[i];
	}
}

__global__ void KMU_PBKDF2(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* DK)
{
	uint8_t pt[4] = { 0, };
	pt[0] = blockIdx.x;
	pt[1] = threadIdx.x;
	uint64_t out[8];
	for (int i = 0; i < 256 * 256; i++) {
		pt[2] = i / 256;
		pt[3] = i % 256;
	}
	KMU_PBKDF2_CORE(pt, 4, 129973, salt, saLen, out);

	for (int i = 0; i < 8; i++) {
		DK[i] = out[i];
	}
	if (blockIdx.x == 0 && threadIdx.x == 0) {
		for (int i = 0; i < 8; i++) {
			printf("%016llx ", DK[i]);
		}
		printf("\n");
	}
}

int main()
{
	cudaEvent_t start, stop;
	uint8_t pt[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	uint8_t salt[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	uint64_t DK[8];
	uint8_t* cu_salt = NULL;
	uint64_t* cu_DK = NULL;
	float elapsed_time_ms = 0.0f;

	cudaMalloc((void**)&cu_salt, 8);
	cudaMemcpy(cu_salt, salt, 8, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&cu_DK, 8 * 8);
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	for (int i = 0; i < 1; i++) {
		KMU_PBKDF2 << <256, 256 >> > (129973, cu_salt, 8, cu_DK);
	}
	cudaEventRecord(stop, 0);
	cudaDeviceSynchronize();
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	printf("Performance : %4.2f GB/s\n", elapsed_time_ms);
}

//int main()
//{
//	uint64_t A[16] = {0, };
//	A[0] = 0x6162636465666768;
//	A[1] = 0x6263646566676869;
//	A[2] = 0x636465666768696a;
//	A[3] = 0x6465666768696a6b;
//	A[4] = 0x65666768696a6b6c;
//	A[5] = 0x666768696a6b6c6d;
//	A[6] = 0x6768696a6b6c6d6e;
//	A[7] = 0x68696a6b6c6d6e6f;
//	A[8] = 0x696a6b6c6d6e6f70;
//	A[9] = 0x6a6b6c6d6e6f7071;
//	A[10] = 0x6b6c6d6e6f707172;
//	A[11] = 0x6c6d6e6f70717273;
//	A[12] = 0x6d6e6f7071727374;
//	A[13] = 0x6e6f707172737475;
//	A[14] = 0x8000000000000000;
//	A[15] = 0x0000000000000000;
//
//	uint64_t B[8] = { 0, };
//	uint64_t* cu_A = NULL;
//	uint64_t* cu_out = NULL;
//	uint64_t* cu_B = NULL;
//	uint64_t* cu_pre_IN = NULL;
//	uint64_t Pre[8];
//	Pre[0] = 0x4319017a2b706e69;
//	Pre[1] = 0xcd4b05938bae5e89;
//	Pre[2] = 0x0186bf199f30aa95;
//	Pre[3] = 0x6ef8b71d2f810585;
//	Pre[4] = 0xd787d6764b20bda2;
//	Pre[5] = 0xa260144709736920;
//	Pre[6] = 0x00ec057f37d14b8e;
//	Pre[7] = 0x06add5b50e671c72;
//	cudaMalloc((void**)&cu_A, 128);
//	cudaMalloc((void**)&cu_B, 64);
//	cudaMalloc((void**)&cu_out, 64);
//	cudaMalloc((void**)&cu_pre_IN, 64);
//	cudaMemcpy(cu_A, A, 128, cudaMemcpyHostToDevice);
//	cudaMemcpy(cu_B, B, 64, cudaMemcpyHostToDevice);
//	cudaMemcpy(cu_pre_IN, Pre, 64, cudaMemcpyHostToDevice);
//	//KMU_PBKDF_SHA_BLOCK << <1, 1 >> > (cu_A, cu_out);
//	//KMU_PBKDF2_BLOCK_CORE << <1, 1 >> > (cu_pre_IN, cu_B, 112, cu_out);
//	KMU_PRE_HMAC_SHA512 << <1, 1 >> > (cu_out);
//	cudaMemcpy(B, cu_out, 64, cudaMemcpyDeviceToHost);
//	for (int i = 0; i < 8; i++)
//		printf("%016llx ", B[i]);
//
//}

//
//__device__ void KMU_PBKDF2_BLOCK_CORE(uint64_t* PRE_IN, uint64_t* hash, uint64_t ptLen, uint64_t* OUT)
//{
//	uint64_t a, b, c, d, e, f, g, h = 0;
//	a = PRE_IN[0];
//	b = PRE_IN[1];
//	c = PRE_IN[2];
//	d = PRE_IN[3];
//	e = PRE_IN[4];
//	f = PRE_IN[5];
//	g = PRE_IN[6];
//	h = PRE_IN[7];
//
//	//Padding
//	uint64_t w0_t = hash[0];
//	uint64_t w1_t = hash[1];
//	uint64_t w2_t = hash[2];
//	uint64_t w3_t = hash[3];
//	uint64_t w4_t = hash[4];
//	uint64_t w5_t = hash[5];
//	uint64_t w6_t = hash[6];
//	uint64_t w7_t = hash[7];
//	uint64_t w8_t = 0x8000000000000000;
//	uint64_t w9_t = 0;
//	uint64_t wa_t = 0;
//	uint64_t wb_t = 0;
//	uint64_t wc_t = 0;
//	uint64_t wd_t = 0;
//	uint64_t we_t = 0;
//	uint64_t wf_t = ptLen  << 3;
//
//#define ROUND_EXPAND()									\
//	{													\
//		w0_t = SHA512_EXPAND (we_t, w9_t, w1_t, w0_t);  \
//		w1_t = SHA512_EXPAND (wf_t, wa_t, w2_t, w1_t);  \
//		w2_t = SHA512_EXPAND (w0_t, wb_t, w3_t, w2_t);  \
//		w3_t = SHA512_EXPAND (w1_t, wc_t, w4_t, w3_t);  \
//		w4_t = SHA512_EXPAND (w2_t, wd_t, w5_t, w4_t);  \
//		w5_t = SHA512_EXPAND (w3_t, we_t, w6_t, w5_t);  \
//		w6_t = SHA512_EXPAND (w4_t, wf_t, w7_t, w6_t);  \
//		w7_t = SHA512_EXPAND (w5_t, w0_t, w8_t, w7_t);  \
//		w8_t = SHA512_EXPAND (w6_t, w1_t, w9_t, w8_t);  \
//		w9_t = SHA512_EXPAND (w7_t, w2_t, wa_t, w9_t);  \
//		wa_t = SHA512_EXPAND (w8_t, w3_t, wb_t, wa_t);  \
//		wb_t = SHA512_EXPAND (w9_t, w4_t, wc_t, wb_t);  \
//		wc_t = SHA512_EXPAND (wa_t, w5_t, wd_t, wc_t);  \
//		wd_t = SHA512_EXPAND (wb_t, w6_t, we_t, wd_t);  \
//		we_t = SHA512_EXPAND (wc_t, w7_t, wf_t, we_t);  \
//		wf_t = SHA512_EXPAND (wd_t, w8_t, w0_t, wf_t);  \
//	}
//
//#define ROUND_STEP(i)																		\
//	{																						\
//		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512[i +  0]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512[i +  1]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512[i +  2]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512[i +  3]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512[i +  4]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512[i +  5]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512[i +  6]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512[i +  7]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512[i +  8]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512[i +  9]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512[i + 10]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512[i + 11]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512[i + 12]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512[i + 13]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512[i + 14]);	\
//		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512[i + 15]);	\
//	}
//	ROUND_STEP(0);
//	for (int i = 16; i < 80; i += 16) {
//		ROUND_EXPAND();
//		ROUND_STEP(i);
//	}
//
//	OUT[0] = a + 0x6a09e667f3bcc908;
//	OUT[1] = b + 0xbb67ae8584caa73b;
//	OUT[2] = c + 0x3c6ef372fe94f82b;
//	OUT[3] = d + 0xa54ff53a5f1d36f1;
//	OUT[4] = e + 0x510e527fade682d1;
//	OUT[5] = f + 0x9b05688c2b3e6c1f;
//	OUT[6] = g + 0x1f83d9abfb41bd6b;
//	OUT[7] = h + 0x5be0cd19137e2179;
//}
//
//__device__ void KMU_PRE_HMAC_SHA512(uint8_t* pt, uint64_t ptLen, uint64_t* IPAD_out, uint64_t* OPAD_out)
//{
//	uint8_t IPAD[128] = { 0x36, };
//	uint8_t OPAD[128] = { 0x5c, };
//	uint64_t buffer[16];
//	for (int i = 0; i < ptLen; i++) {
//		IPAD[i] = IPAD[i] ^ pt[i];
//		OPAD[i] = OPAD[i] ^ pt[i];
//	}
//	KMU_PBKDF_SHA_BLOCK((uint64_t*)IPAD, IPAD_out);
//	KMU_PBKDF_SHA_BLOCK((uint64_t*)OPAD, OPAD_out);
//
//}
//
////Fixed dkLen = 64byte
//__device__ void KMU_PBKDF2_CORE(uint8_t* pt, uint64_t ptLen, uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* DK) {
//	uint64_t dk[8] = { 0, };
//	uint64_t IPAD_out[8] = { 0, };
//	uint64_t OPAD_out[8] = { 0, };
//	uint64_t out[8] = { 0, };
//	KMU_PRE_HMAC_SHA512(pt, ptLen, IPAD_out, OPAD_out);
//
//	//U0 operation
//	uint8_t buffer[64] = { 0, };
//	for (int i = 0; i < saLen; i++) {
//		buffer[i] = salt[i];
//	}
//	buffer[saLen + 0] = (1 >> 24) & 0xff;
//	buffer[saLen + 1] = (1 >> 16) & 0xff;
//	buffer[saLen + 2] = (1 >> 8) & 0xff;
//	buffer[saLen + 3] = (1) & 0xff;
//	KMU_PBKDF2_BLOCK_CORE(IPAD_out, (uint64_t*)buffer, 64 + saLen + 4, dk);
//	KMU_PBKDF2_BLOCK_CORE(OPAD_out, dk, 64, dk);
//	for (int i = 0; i < 8; i++)
//		out[i] = dk[i];
//	for (int i = 1; i < iteration; i++) {
//		KMU_PBKDF2_BLOCK_CORE(IPAD_out, dk, 128, dk);
//		KMU_PBKDF2_BLOCK_CORE(OPAD_out, dk, 128, dk);
//		for (int i = 0; i < 8; i++)
//			out[i] ^= dk[i];
//	}
//	
//	for (int i = 0; i < 8; i++) {
//		DK[i] = out[i];
//	}
//}
//
//
////passwordLen = 4
//__global__ void KMU_PBKDF2(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* DK)
//{
//	uint8_t pt[4] = { 0, };
//	pt[0] = blockIdx.x;
//	pt[1] = threadIdx.x;
//	uint64_t out[8];
//	for (int i = 0; i < 256 * 256; i++) {
//		pt[2] = i / 256;
//		pt[3] = i % 256;
//	}
//	KMU_PBKDF2_CORE(pt, 4, 129973, salt, saLen, out);
//
//	for (int i = 0; i < 8; i++) {
//		DK[i] = out[i];
//	}
//	if (blockIdx.x == 0 && threadIdx.x == 0) {
//		for (int i = 0; i < 8; i++) {
//			printf("%016llu ", DK[i]);
//		}
//		printf("\n");
//	}
//
//}
//
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
//	cudaMemcpy(cu_salt, salt, 8, cudaMemcpyHostToDevice);
//	cudaMalloc((void**)&cu_DK, 8 * 8);
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);
//	cudaEventRecord(start, 0);
//	for (int i = 0; i < 1; i++) {
//		KMU_PBKDF2 << <256, 256 >> > (129973, cu_salt, 8, cu_DK);
//	}
//	cudaEventRecord(stop, 0);
//	cudaDeviceSynchronize();
//	cudaEventSynchronize(start);
//	cudaEventSynchronize(stop);
//	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
//	printf("Performance : %4.2f GB/s\n", elapsed_time_ms);
//}