#include "type.cuh"

typedef struct {
	uint64_t hash[8];
	uint64_t ptLen;
}SHA512_INFO;

#define ROTL64(x, n)			(((x) << (n)) | ((x) >> (64 - (n))))
#define ROTR64(x, n)			(((x) >> (n)) | ((x) << (64 - (n))))
#define	SF(x, n)				(x >> (n))

//THETA
#define WE0_512(x)				(ROTR64(x,  1) ^ ROTR64(x, 8) ^ SF(x, 7))
#define WE1_512(x)				(ROTR64(x,  19) ^ ROTR64(x, 61) ^ SF(x, 6))

//SIGMA
#define BS0_512(x)				((ROTR64(x,  28)) ^ ROTR64(x, 34) ^ ROTR64(x,  39))
#define BS1_512(x)				(ROTR64(x,  14) ^ ROTR64(x, 18) ^ ROTR64(x,  41))

//OPERATOR
#define CH(x, y, z)			((x & y) ^ (~(x) & (z)))
#define MAJ(x, y, z)		(((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define SHA512_F0(x,y,z) ((z) ^ ((x) & ((y) ^ (z))))
#define SHA512_F1(x,y,z) (((x) & (y)) | ((z) & ((x) ^ (y))))
//ENDIAN
#define ENDIAN_CHANGE(val)	(\
(((val) >> 56) & 0x00000000000000FF) | (((val) >> 40) & 0x000000000000FF00) | \
(((val) >> 24) & 0x0000000000FF0000) | (((val) >>  8) & 0x00000000FF000000) | \
(((val) <<  8) & 0x000000FF00000000) | (((val) << 24) & 0x0000FF0000000000) | \
(((val) << 40) & 0x00FF000000000000) | (((val) << 56) & 0xFF00000000000000))

//CORE OPERATION
#define SHA512_STEP(F0, F1, a, b, c ,d ,e ,f ,g ,h, x, K)	\
{															\
	h += K;													\
	h += x;													\
	h += BS1_512(e);										\
	h += F0(e, f, g);										\
	d += h;													\
	h += BS0_512(a);										\
	h += F1(a, b, c);										\
}

//Word Expansion
#define SHA512_EXPAND(x, y, z ,w) (WE1_512(x) + y + WE0_512(z) + w)
