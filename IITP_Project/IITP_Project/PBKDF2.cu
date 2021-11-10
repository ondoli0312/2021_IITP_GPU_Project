#define CUDA_API_PER_THREAD_DEFAULT_STEAM
#include "type.cuh"
#include "SHA512.cuh"

//Clear
__constant__ uint64_t cont_512v2[80] =
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
	0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b,
	0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
	0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
	0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
	0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817
};

__constant__ uint64_t TT[1] = { 0x2527B51365102B38 };

//Clear
__device__ void KMU_PBKDF_SHA_BLOCKv2(uint64_t* PT, uint64_t* output)
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

#define ROUND_STEP(i)																			\
	{																							\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512v2[i +  0]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512v2[i +  1]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512v2[i +  2]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512v2[i +  3]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512v2[i +  4]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512v2[i +  5]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512v2[i +  6]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512v2[i +  7]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512v2[i +  8]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512v2[i +  9]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512v2[i + 10]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512v2[i + 11]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512v2[i + 12]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512v2[i + 13]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512v2[i + 14]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512v2[i + 15]);	\
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

__device__ void KMU_PBKDF_SHA_BLOCK_COREv2(uint64_t* PREIN, uint64_t* hash, uint64_t ptLen, uint64_t* OUT) {
	
	uint64_t a = PREIN[0];
	uint64_t b = PREIN[1];
	uint64_t c = PREIN[2];
	uint64_t d = PREIN[3];
	uint64_t e = PREIN[4];
	uint64_t f = PREIN[5];
	uint64_t g = PREIN[6];
	uint64_t h = PREIN[7];

	//word setting
	volatile uint64_t w0_t = (hash[0]);
	volatile uint64_t w1_t = (hash[1]);
	volatile uint64_t w2_t = (hash[2]);
	volatile uint64_t w3_t = (hash[3]);
	volatile uint64_t w4_t = (hash[4]);
	volatile uint64_t w5_t = (hash[5]);
	volatile uint64_t w6_t = (hash[6]);
	volatile uint64_t w7_t = (hash[7]);
	volatile uint64_t w8_t = (0x8000000000000000);
	volatile uint64_t w9_t = 0;
	volatile uint64_t wa_t = 0;
	volatile uint64_t wb_t = 0;
	volatile uint64_t wc_t = 0;
	volatile uint64_t wd_t = 0;
	volatile uint64_t we_t = 0;
	volatile uint64_t wf_t = (ptLen << 3);


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

#define ROUND_STEP(i)																			\
	{																							\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512v2[i +  0]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512v2[i +  1]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512v2[i +  2]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512v2[i +  3]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512v2[i +  4]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512v2[i +  5]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512v2[i +  6]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512v2[i +  7]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512v2[i +  8]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512v2[i +  9]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512v2[i + 10]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512v2[i + 11]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512v2[i + 12]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512v2[i + 13]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512v2[i + 14]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512v2[i + 15]);	\
	}
	ROUND_STEP(0);
	for (int i = 16; i < 80; i += 16) {
		ROUND_EXPAND();
		ROUND_STEP(i);
	}

	OUT[0] = a + PREIN[0];
	OUT[1] = b + PREIN[1];
	OUT[2] = c + PREIN[2];
	OUT[3] = d + PREIN[3];
	OUT[4] = e + PREIN[4];
	OUT[5] = f + PREIN[5];
	OUT[6] = g + PREIN[6];
	OUT[7] = h + PREIN[7];
}

__device__ void KMU_PBKDF_SALTv1(uint64_t* PREIN, uint64_t* salt, uint64_t* OUT) {

	uint64_t a = PREIN[0];
	uint64_t b = PREIN[1];
	uint64_t c = PREIN[2];
	uint64_t d = PREIN[3];
	uint64_t e = PREIN[4];
	uint64_t f = PREIN[5];
	uint64_t g = PREIN[6];
	uint64_t h = PREIN[7];

	//word setting
	volatile uint64_t w0_t = ENDIAN_CHANGE(salt[0]);
	volatile uint64_t w1_t = ENDIAN_CHANGE(salt[1]);
	volatile uint64_t w2_t = ENDIAN_CHANGE(salt[2]);
	volatile uint64_t w3_t = ENDIAN_CHANGE(salt[3]);
	volatile uint64_t w4_t = ENDIAN_CHANGE(salt[4]);
	volatile uint64_t w5_t = ENDIAN_CHANGE(salt[5]);
	volatile uint64_t w6_t = ENDIAN_CHANGE(salt[6]);
	volatile uint64_t w7_t = ENDIAN_CHANGE(salt[7]);
	volatile uint64_t w8_t = ENDIAN_CHANGE(salt[8]);
	volatile uint64_t w9_t = ENDIAN_CHANGE(salt[9]);
	volatile uint64_t wa_t = ENDIAN_CHANGE(salt[10]);
	volatile uint64_t wb_t = ENDIAN_CHANGE(salt[11]);
	volatile uint64_t wc_t = ENDIAN_CHANGE(salt[12]);
	volatile uint64_t wd_t = ENDIAN_CHANGE(salt[13]);
	volatile uint64_t we_t = ENDIAN_CHANGE(salt[14]);
	volatile uint64_t wf_t = (salt[15]);

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

#define ROUND_STEP(i)																			\
	{																							\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w0_t, cont_512v2[i +  0]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w1_t, cont_512v2[i +  1]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, w2_t, cont_512v2[i +  2]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, w3_t, cont_512v2[i +  3]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, w4_t, cont_512v2[i +  4]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, w5_t, cont_512v2[i +  5]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, w6_t, cont_512v2[i +  6]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, w7_t, cont_512v2[i +  7]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, a, b, c, d, e, f, g, h, w8_t, cont_512v2[i +  8]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, h, a, b, c, d, e, f, g, w9_t, cont_512v2[i +  9]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, g, h, a, b, c, d, e, f, wa_t, cont_512v2[i + 10]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, f, g, h, a ,b, c, d, e, wb_t, cont_512v2[i + 11]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, e, f, g, h, a ,b, c, d, wc_t, cont_512v2[i + 12]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, d, e, f, g, h, a ,b, c, wd_t, cont_512v2[i + 13]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, c, d, e, f, g, h, a, b, we_t, cont_512v2[i + 14]);	\
		SHA512_STEP(SHA512_F0, SHA512_F1, b, c, d, e, f, g, h, a, wf_t, cont_512v2[i + 15]);	\
	}
	ROUND_STEP(0);
	for (int i = 16; i < 80; i += 16) {
		ROUND_EXPAND();
		ROUND_STEP(i);
	}

	OUT[0] = a + PREIN[0];
	OUT[1] = b + PREIN[1];
	OUT[2] = c + PREIN[2];
	OUT[3] = d + PREIN[3];
	OUT[4] = e + PREIN[4];
	OUT[5] = f + PREIN[5];
	OUT[6] = g + PREIN[6];
	OUT[7] = h + PREIN[7];
}

__device__ void KMU_PRE_HMAC_SHA512v2(uint64_t pt, uint64_t ptLen, uint64_t* OPAD_out, uint64_t* IPAD_out)
{
	
	uint64_t IPAD[16];
	uint64_t OPAD[16];
	for (int i = 0; i < 16; i++) {
		IPAD[i] = 0x3636363636363636;
		OPAD[i] = 0x5c5c5c5c5c5c5c5c;
	}
	
	//여기 password word가 여러개 오는 경우 수정하기
	for (int i = 0; i < 1; i++) {
		IPAD[i] = IPAD[i] ^ pt;
		OPAD[i] = OPAD[i] ^ pt;
	}
	
	KMU_PBKDF_SHA_BLOCKv2(IPAD, IPAD_out);
	KMU_PBKDF_SHA_BLOCKv2(OPAD, OPAD_out);
}

__device__ void KMU_PKBDF2_Core_v2(uint64_t pt, uint64_t ptLen, uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* DK, uint8_t* cracking) {
	uint64_t OPAD_out[16];
	uint64_t IPAD_out[16];
	uint64_t dk1[8];
	uint64_t dk2[8];
	uint64_t out[8];
	uint8_t buffer[128] = { 0, };

	//pre_computation
	uint64_t crackpt = pt;
	KMU_PRE_HMAC_SHA512v2(crackpt, ptLen, OPAD_out, IPAD_out);

	//salt value gen
	for (int i = 0; i < saLen; i++)
		buffer[i] = salt[i];
	buffer[saLen + 3] = (1) & 0xff;
	buffer[saLen + 4] = 0x80;

	//U0 Gen
	//8 -> passwordLen
	((uint64_t*)(buffer))[15] = (128 + 8) << 3;
	KMU_PBKDF_SALTv1(IPAD_out, (uint64_t*)buffer, dk1);
	KMU_PBKDF_SHA_BLOCK_COREv2(OPAD_out, dk1, 192, dk2);

	for (int j = 0; j < 8; j++)
		out[j] = dk2[j];
	for (int i = 1; i < iteration; i++) {
		KMU_PBKDF_SHA_BLOCK_COREv2(IPAD_out, dk2, 192, dk1);
		KMU_PBKDF_SHA_BLOCK_COREv2(OPAD_out, dk1, 192, dk2);
		
		for (int j = 0; j < 8; j++) {
			out[j] ^= dk2[j];
		}
	}
	if (out[0] == DK[0]) {
		cracking[0] = 1;
	}
}

//6byte cracking Sample
//__global__ void KMU_PBKDF2v2(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* Target)
//{
//	//Password Word Setting
//	uint64_t pt = 0;
//	uint64_t pt1 = 0;
//	uint32_t flag = blockDim.x * blockIdx.x + threadIdx.x;
//	uint8_t cracking[1];
//	cracking[0] = 1;
//	pt = (((flag) % 128)) | (((flag / 128) % 128) << 8) | (((flag / 16384) % 128) << 16);
//	pt1 = pt;
//	flag = (flag / 2097152) % 128; //(0 <= flag < 16)
//	for (uint64_t i = 0; i < 16384; i++) {
//		pt = pt + ((i % 128) << 24);
//		pt = pt + ((i / 128) << 32);
//		KMU_PKBDF2_Core_v2(pt, 1, iteration, salt, saLen, Target, cracking);
//		pt = pt1;
//	}
//	
//}

//2byte cracking Sample
__global__ void KMU_PBKDF2v3(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* Target)
{
	//Password Word Setting
	uint64_t pt = 0;
	uint64_t pt1 = 0;
	//uint64_t tt = Target[0];
	uint8_t cracking[1];
	cracking[0] = 0;
	pt = blockIdx.x;
	pt = pt << 56;
	pt1 = threadIdx.x;
	pt = pt + (pt1 << 48);

	KMU_PKBDF2_Core_v2(pt, 1, iteration, salt, saLen, Target, cracking);
	if (cracking[0] == 1) {
		printf("%016llx ", pt);
	}

}

__global__ void KMU_PBKDF2v4(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* Target)
{
	//Password Word Setting
	uint64_t pt = 0;
	uint64_t pt1 = 0;
	//uint64_t tt = Target[0];
	uint8_t cracking[1];
	cracking[0] = 0;
	uint64_t flag = (blockDim.x * blockIdx.x) + threadIdx.x;

	pt = (flag % 128) << 56;
	pt = pt | (((flag >> 7) & 0x7F) << 48);
	pt = pt | (((flag / 16384) & 0x7F) << 40);
	for (int j = 0; j < 2; j++) {
		pt = pt | ((j & 0x7F) << 32);
		KMU_PKBDF2_Core_v2(pt, 1, iteration, salt, saLen, Target, cracking);
	}
	if (cracking[0] == 1) {
		printf("%016llx ", pt);
	}

}


__global__ void KMU_PBKDF2v5(uint64_t iteration, uint8_t* salt, uint64_t saLen, uint64_t* Target, int stream)
{
	//Password Word Setting
	uint64_t pt = 0;
	uint64_t pt1 = 0;
	//uint64_t tt = Target[0];
	uint8_t cracking[1];
	cracking[0] = 0;
	uint64_t flag = (blockDim.x * blockIdx.x) + threadIdx.x;
	pt = (flag % 128) << 56;
	pt = pt | (((flag >> 7) & 0x7F) << 48);
	pt = pt | (((flag / 16384) & 0x7F) << 40);
	pt = pt | ((stream & 0x7F) << 32);
	KMU_PKBDF2_Core_v2(pt, 1, iteration, salt, saLen, Target, cracking);
	if (cracking[0] == 1) {
		printf("%016llx ", pt);
	}
}

__global__ void Core()
{

}

int main()
{
	uint8_t salt[4];
	uint64_t Target[1];
	cudaStream_t stream[2];

	cudaEvent_t start, stop;
	cudaError_t err;
	float elapsed_time_ms = 0.0f;
	
	salt[0] = 0x73;
	salt[1] = 0x61;
	salt[2] = 0x6c;
	salt[3] = 0x74;
	Target[0] = 0x9D536145420D4242;
	TT[0] = Target[0];
	uint8_t* cuda_salt = NULL;
	uint64_t* cuda_Target = NULL;

	cudaMalloc((void**)&cuda_salt, 4);
	cudaMalloc((void**)&cuda_Target, 8);

	cudaMemcpy(cuda_salt, salt, 4, cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_Target, Target, 8, cudaMemcpyHostToDevice);
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	for (int i = 0; i < 1; i++) {
		//KMU_PBKDF2v4 << <4096, 512 >> > (129937, cuda_salt, 4, cuda_Target);
	}cudaEventRecord(stop, 0);
	cudaDeviceSynchronize();
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	printf("Performance : %4.2f/s\n", elapsed_time_ms/1000);

	uint8_t* Salt = NULL;
	uint64_t* TTarget = NULL;
	uint8_t* cuda_salt_s0 = NULL;
	uint8_t* cuda_salt_s1 = NULL;
	uint64_t* cuda_Target0 = NULL;
	uint64_t* cuda_Target1 = NULL;

	cudaMalloc((void**)&cuda_salt_s0, 4);
	cudaMalloc((void**)&cuda_salt_s1, 4);
	cudaMalloc((void**)&cuda_Target0, 8);
	cudaMalloc((void**)&cuda_Target1, 8);
	cudaMallocHost((void**)&Salt, 4);
	cudaMallocHost((void**)&TTarget, 8);

	cudaStreamCreate(&stream[0]);
	cudaStreamCreate(&stream[1]);

	Salt[0] = 0x73;
	Salt[1] = 0x61;
	Salt[2] = 0x6c;
	Salt[3] = 0x74;
	TTarget[0] = 0x9D536145420D4242;


	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	cudaMemcpyAsync(cuda_salt_s0, Salt, 4, cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(cuda_Target0, TTarget, 8, cudaMemcpyHostToDevice, stream[0]);
	KMU_PBKDF2v5 << <4096, 512, 0, stream[0] >> > (129937, cuda_salt, 4, cuda_Target, 0);
	cudaMemcpyAsync(cuda_salt_s1, Salt, 4, cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(cuda_Target1, TTarget, 8, cudaMemcpyHostToDevice, stream[1]);
	KMU_PBKDF2v5 << <4096, 512, 0, stream[1] >> > (129937, cuda_salt, 4, cuda_Target, 1);
	cudaEventRecord(stop, 0);
	cudaDeviceSynchronize();
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time_ms, start, stop);
	printf("Performance : %4.2f/s\n", elapsed_time_ms / 1000);

}