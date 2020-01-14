/*
 ============================================================================
 Name        : RSA_Montgomery.c
 Author      : Tim Salomonsson
 Version     :
 Copyright   :
 Description : RSA using Montgomery Reduction
 ============================================================================
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "bn.h"

#define LUT_SIZE 2049

static struct bn lut[LUT_SIZE];

/*
 * Performs bitwise Montgomery modular multiplication ( X*Y*R^(-1) mod M)
 *
 * Parameters:
 * 		x,y,m - bignums
 * 		mBits - # of bits in m
 * 		out	  - bignum result
 */

void montMult(struct bn*  x, struct bn*  y, struct bn*  m, int mBits, struct bn*  out){

	struct bn t;
	bignum_init(&t);

	int i;
	for(i = mBits; i > 0 ; i--){					//efficient loop exit

		int t0Bit = bignum_getbit(&t,0);
		int xiBit = bignum_getbit(x, mBits - i);	//loop exit requires subtraction here
		int y0Bit = bignum_getbit(y,0);
		int op = t0Bit + (xiBit * y0Bit);

		if(xiBit == 1){
			bignum_add(&t, y, &t);
		}

		if(op == 1){
			bignum_add(&t, m, &t);
		}

		bignum_rshift(&t,&t, 1);
	}

	if(bignum_cmp(&t, m) >= 0){
		bignum_sub(&t,m,&t);
	}

	bignum_assign(out,&t);
}



/*mod exp, no LUT */

void modExp(struct bn*  x, struct bn*   e, int eBits, struct bn*  m, int mBits, struct bn*  r2m,  struct bn*   out){

	struct bn z,one;
	struct bn parr[3];
	struct bn zarr[3];

	//reduce z?
	bignum_from_int(&z, 1);
	montMult(&z,r2m,m, mBits, &zarr[1]);

	//reduce x, assign to p
	montMult(x,r2m,m, mBits,&parr[1]);

	struct bn tm;

	int i = 0;
	for(; i < eBits; i++){

		bignum_assign(&tm, &parr[1]);
		montMult(&tm,&parr[1],m, mBits, &parr[2]);

		if(bignum_getbit(e, i) == 1){
			montMult(&zarr[1],&parr[1],m,mBits,&zarr[2]);
		}else{
			bignum_assign(&zarr[2],&zarr[1]);
		}

		//printf("num bits p: %d, num bits z: %d\n", bignum_numbits(&parr[1]), bignum_numbits(&zarr[1]));
		bignum_assign(&parr[1], &parr[2]);
		bignum_assign(&zarr[1], &zarr[2]);
	}

	bignum_from_int(&one, 1);
	montMult(&zarr[1], &one, m, mBits, out);
}


/* Mod Exp using precomputed LUT */

void modExpLUT(struct bn*  x, struct bn*  e, int eBits, struct bn*  m, int mBits, struct bn*  r2m, struct bn*   out){

	struct bn z,one;
	struct bn parr[3];
	struct bn zarr[3];

	//reduce z?
	bignum_from_int(&z, 1);
	montMult(&z,r2m,m, mBits, &zarr[1]);

	bignum_assign(&parr[1],&lut[0]);

	int b = 1;
	int i = 0;
	for(; i < eBits; i++){
		bignum_assign(&parr[2],&lut[i+1]);

		if(bignum_getbit(e, i) == 1){
			montMult(&zarr[1],&lut[i],m,mBits,&zarr[2]);
		}else{
			bignum_assign(&zarr[2],&zarr[1]);
		}

		//printf("num bits p: %d, num bits z: %d\n", bignum_numbits(&parr[1]), bignum_numbits(&zarr[1]));
		bignum_assign(&parr[1], &parr[2]);
		bignum_assign(&zarr[1], &zarr[2]);
		b++;
	}

	bignum_from_int(&one, 1);
	montMult(&zarr[1], &one, m, mBits, out);
}

/*
void genLUT(struct bn* valToDec, struct bn* m, int mBits, struct bn* r2m) {

	struct bn two;
	bignum_from_int(&two, 2);

	bignum_assign(&lut[0],valToDec);

	int i;
	for (i = 1; i < (LUT_SIZE + 1); i++) {
		struct bn temp;
		bignum_assign(&temp,&lut[i-1]);
		struct bn tmp, tmp1, tmp2;
		bignum_assign(&tmp2, &two);
		bignum_pow(&temp, &tmp2, &tmp);
		bignum_mod(&tmp, m, &tmp1);
		bignum_assign(&lut[i],&tmp1);
	}

	int a = 0;
	FILE *f;
	f = fopen("LUT.txt", "a");

	for(; a < (mBits + 1); a++){
		montMult(&lut[a],r2m,m, mBits, &lut[a]);
		int size = 8192;
		char str[size];
		bignum_to_string(&lut[a],str, size);
		fprintf(f,"%s\n", str);
	}

	fclose(f);
}
*/

int parseLUT(int start, int num){

	FILE *f;
	f = fopen("LUT.txt", "r");
	char str[num+1];
	static char* zpad[8] = {"", "0", "00", "000", "0000", "00000", "000000", "0000000"};
	//printf("Parse LUT: %d to %d\n", start, num + start);

	//puts("");
	int i = 0;
	int a = 0;
	for(; i < (num + start); i++){
		fscanf(f, "%s\n", str);

		if(i >= start){

			int len = strlen(str);
			char* dup;
			if((len & 1) == 1){
				sprintf(str, "%s%s", zpad[1], (dup = strdup(str)));
				len++;
			}

			int lenMod8 = len -((len >> 3) << 3);
			if(lenMod8 != 0){
				sprintf(str, "%s%s", zpad[lenMod8], (dup = strdup(str)));
				len += lenMod8;
			}
			bignum_from_string(&lut[i-start],str,len);
		}
		a++;
	}
	return i;
}

int main(void) {

	/* ----------- 12-bit Test -------------- */

	struct bn n,e,d,r2m;
	bignum_from_int(&n, 3233); 			//modulus
	bignum_from_int(&e, 17);			//public
	bignum_from_int(&d, 2753);			//private
	bignum_from_int(&r2m, 1179);		//R^2m mod M

	struct bn valToDec;					//value to decrypt/encrypt
	bignum_from_int(&valToDec, 855);

	int nBits = bignum_numbits(&n);
	int dBits = bignum_numbits(&d);

	struct bn result;
	bignum_init(&result);

	clock_t before = clock();
	modExp(&valToDec, &d, dBits, &n, nBits, &r2m, &result);
	clock_t after = clock();

	double msec = (double)(after - before) / CLOCKS_PER_SEC;

	//print result and timing
	printf("-------Test 1--------\n");
	printf(" RSA Keysize: %3d [bits]\n",dBits);
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("Time(no LUT):  %.5f [sec]\n", msec);

	bignum_init(&result);

	//genLUT(&valToDec, &n, nBits, &r2m);
	int lutSeek = parseLUT(0,dBits+1);
	//printf("Lutseek: %d\n", lutSeek);

	before = clock();
	modExpLUT(&valToDec, &d, dBits, &n, nBits, &r2m, &result);
	after = clock();

	msec = (double)(after - before) / CLOCKS_PER_SEC;

	//print result and timing
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("   Time(LUT):  %.5f [sec]\n", msec);

	/* ----------- 512-bit Modulus Test -------------- */


	struct bn n512, pub, priv, v2Dec;
	int e1 = 65537;
	bignum_from_int(&pub,e1);
	char str1[] = "758463d46999c11496449db8dddd1e407de2e9a8f33612f454866acddd759da8173d4e3fe8c4eaf121f86f87ac8e1d58f54e2c6a80bcf8c404884795252224ad";
	bignum_from_string(&n512,str1, 128);
	char str2[] = "68827b718d1452d4e72a5085f6b14dd516df34e3ae9fb94d96da0fa3d33e651cc244b0275a24ab0753b5c01eac2f8f0d700c587bbd6d8aeb6a4e99e1a9372655";
	bignum_from_string(&priv,str2, 128);
	char str3[] = "45462476f31c3dfde5ac5fde4862d33d917f52255d80555b543584a32b71762a1fc719a341c0e925e9fff02a657764ae78b143d324cfc8892695c55801237885";
	bignum_from_string(&v2Dec,str3,128);
	char r2ms[] = "47395beb0ae85106f9f8548040a9b165d9a37499d0d98a14a5bcd0b943d0549be18b2ced65bfc42db40331f3ec67faf9cccf19e51d3ef7a09e03ebb1855d5e5e";
	bignum_from_string(&r2m,r2ms,128);

	nBits = bignum_numbits(&n512);
	dBits = bignum_numbits(&priv);

	bignum_init(&result);

	clock_t before1 = clock();
	modExp(&v2Dec, &priv,dBits, &n512, nBits, &r2m, &result);
	clock_t after1 = clock();

	double msec1 = (float)(after1 - before1) / CLOCKS_PER_SEC;

	//print result and timing
	printf("-------Test 2--------\n");
	printf(" RSA Keysize: %4d [bits]\n",dBits);
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("Time(no LUT):  %.5f [sec]\n", msec1);

	bignum_init(&result);

	//genLUT(&v2Dec, &n512, nBits, &r2m);
	lutSeek = parseLUT(lutSeek,dBits+1);


	before1 = clock();
	modExpLUT(&v2Dec, &priv,dBits, &n512, nBits, &r2m, &result);
	after1 = clock();

	msec1 = (float)(after1 - before1) / CLOCKS_PER_SEC;

	//print result and timing
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("   Time(LUT):  %.5f [sec]\n", msec1);


	/* ----------- 1024-bit Modulus Test -------------- */

	struct bn n1024, pub1, priv1, v2Dec1;

	bignum_from_int(&pub1,e1);
	char str4[] = "79eec1e33a41bf4592557bb1991b1830d4b445f55e3c9e683afc7a7f4abf05549a5e7ea811f8c3faf58450c2eafce1a25c5eb49821d0f930247ef2c6a6e426f01f91a6090292a433d84b93a1e6c5ba933c48f48923aa727f3de18c5fa4f1c0f7cce43cf407f94ee1d316d572b4428c7399158b76fa15f8b3dfbb36bd5f4bc5d1";
	bignum_from_string(&n1024,str4, 256);
	char str5[] = "233c05371e4c85731b382c88438ffacb918b8e73bb099554d546c43728684ea805fbac69f0d78bfa671c17225c393b1269d2cc28f20cab1568566edd4cb8bd2f59e4b25f4b3787af54e002216bc42a34a2bdbd7bfe4ddab35dde5256fc7bfbc1b39f641c86e99950768214e69b18f806b0d200908484eb7cf6e817ab57400861";
	bignum_from_string(&priv1,str5, 256);
	char str6[] = "4e29e645da6efddda068a8dcfceea970a5e86f7b518655cd3fba103d6899618a6b7caa86df16f28f7bdadbe2ad250794c9f20c9c42338624ab077f9f9ae3733a5c3bf8b4686b56cfe635be0010bf734fdc2a4f2ce5cf920fd4e79c6b7330a8fc2025e61d33dd8b3056390a2226d9d9eaec37f7aea1682f25120c260ecb165823";
	bignum_from_string(&v2Dec1,str6,256);
	char r2ms1[] = "1a32ca1d9343f9ac08567501d91b0b29540e5e6914aaf46c460b92007b6264ca7a4be15e5346933dd2865022a2535729ea817c215f80714384b8235705b88bc3a295fe00ae789bd241d5816e5d617c362a2ed1bdd8b45ca26f558a987de829afe0253c33b6a7bab59c35429c29c4ab63a0ab16c7f8c4b9319f6f1947266522a5";
	bignum_from_string(&r2m,r2ms1,256);

	dBits = bignum_numbits(&priv1);
	nBits = bignum_numbits(&n1024);

	bignum_init(&result);

	clock_t before2 = clock();
	modExp(&v2Dec1, &priv1,dBits, &n1024, nBits, &r2m, &result);
	clock_t after2 = clock();

	msec1 = (float)(after2 - before2) / CLOCKS_PER_SEC;

	//print result and timing
	printf("-------Test 3--------\n");
	printf(" RSA Keysize: %5d [bits]\n",dBits);
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("Time(no LUT):  %.5f [sec]\n", msec1);
	bignum_init(&result);

	//genLUT(&v2Dec1, &n1024, nBits, &r2m);
	lutSeek = parseLUT(lutSeek,dBits+2);
	before2 = clock();
	modExpLUT(&v2Dec1, &priv1,dBits, &n1024, nBits,&r2m, &result);
	after2 = clock();

	msec1 = (float)(after2 - before2) / CLOCKS_PER_SEC;

	//print result and timing
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("   Time(LUT):  %.5f [sec]\n", msec1);


	/* ----------- 2048-bit Modulus Test -------------- */

	struct bn n2048, pub2, priv2, v2Dec2;

	bignum_from_int(&pub2,e1);
	char str7[] = "bc07d529450214ef63a8d61966987e8ca0594d9a7ec4f1881117b4f8ecbdc74b8769f6c98bfe931c9474116be8bd36527acfd95f6633d12cc8a960ab3d3e7a0b4b3e4990b594ee61af3b56315337501225525fb997b65c38118d614601dcb8bd631673a510498f2c3dab44d723d8b6daa697d0108e7fcb4d27525f386e7fcd9ce29c4ab12c4258aa77872259a25804791a1eaef54b65226ec84765442ac839db30467d86910e700d802807de1f4fef5235738d66359cb0a2707cb9cd90e90bb1f2d0d807aafbd048b1ddbb156d4984cfbbaa9a435b9230d213140dd5be64b5e594945d474665eaf5267fc598a5f75b99f83b029971b80c4149891d43abe62b95";
	bignum_from_string(&n2048,str7, 512);
	char str8[] = "9b78ea83264133684182400d5eaca6aec68330cc97176712f7f71f3758210f61df44f9beead78372753987922f2e0c75a480aa1edc95e9d65ad0da529ce044ef83b6ac03507125ae75c2dd61098ac9d54730d65fd21702278633dd8392549c18548f22ee100a92aca50d316da68131a897691dac22f77df57c96fa8ee1a7212db313396410a5c9c8a31f6f940724cff2b2db5eb078eedad92b6ff29a8636fcd370e99773e96168f34839693f84b7a083597bfbe0f674c79b2348b038ca730cada30bcf2dd9cde27dd555891d3cc10b7831b23e7cda163570635727f11d569492a201f55c56d9a92d46f71b6ecea30f28f8c040f834a2da43f72a1ec927df9441";
	bignum_from_string(&priv2,str8, 512);
	char str9[] = "4bd1a139b2ae5bcae58410ae32ce65e41ef226d30bb2d020e1cb02387f13985d4f18d154444954a4831a26870d1671f54f3ff87efe3adc66cac098160a524674740948d6ba8466054cce1d27018bbd5dd9c4b58def7d62cf8c0d6621ad846324b72a0414c56843075c4199c5963f55977ea6437a501afe3eebf150b1e2cbbbcde4be89e1ce8e72f8a334297224418e29ad44882a99ba59c5eb481b5faeeebd423b5bdab6c7edb288e8ab42f01a18cef3521c3cbcad8fde05e2f189070725c13716112b7497bb27250fe4141b41de67e0b0fe1763e0831ecf692ad1ff18e5f0186a6e7729ec84e7b9e2483838be73cbe4fea67fa186b329bfd1434dad528524a4";
	bignum_from_string(&v2Dec2,str9,512);
	char r2ms2[]= "10e2f70a5f5ea34371cd7f6d36ce95604746f2aa503bc45369201212a4006df4433827b085890ed3a614058df7af4caa9a988bb5cbc49179e0a4e76b046926b3f700532e1ed1d191985176c2cd7f9600f45eb96323d975060c44f06ef3bdaa220957e7905c5641276e7752e7e503f930cb49a4abe90cae46270a41e17964206bd6edaa7943b32237d2bfa4063060b388424944ec21c7c2f3bc29554214dee86c848116fc1fd28b60b0b438aa8bc8303c0788ea216bace026f78c09aa10b139a5ee415aa73888ac15157ab9a355eda90b7838e8cddb44a626d1c17a203eaf3c64be524f077df6892984a7198b9c3ba31228bb49259162572747ec51e5b49849ae";
	bignum_from_string(&r2m,r2ms2,512);

	dBits = bignum_numbits(&priv2);
	nBits = bignum_numbits(&n2048);

	bignum_init(&result);

	clock_t before3 = clock();
	modExp(&v2Dec2, &priv2,dBits, &n2048, nBits, &r2m, &result);
	clock_t after3 = clock();

	msec1 = (float)(after3 - before3) / CLOCKS_PER_SEC;

	//print result and timing
	printf("-------Test 4--------\n");
	printf(" RSA Keysize: %6d [bits]\n",dBits);
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("Time(no LUT):  %.5f [sec]\n", msec1);
	bignum_init(&result);

	//genLUT(&v2Dec2, &n2048, nBits, &r2m);
	parseLUT(lutSeek, dBits+1);
	before3 = clock();
	modExpLUT(&v2Dec2, &priv2,dBits, &n2048, nBits, &r2m, &result);
	after3 = clock();

	msec1 = (float)(after3 - before3) / CLOCKS_PER_SEC;

	//print result and timing
	printf("  RSA Result: ");
	bignum_print(&result);
	printf("   Time(LUT):  %.5f [sec]\n", msec1);



	return 0;
}

