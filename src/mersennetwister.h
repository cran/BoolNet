#ifndef MERSENNETWISTER_H_
#define MERSENNETWISTER_H_

extern void init_genrand(unsigned long s);

extern void init_by_array(unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
extern unsigned long genrand_int32(void);

/* generates a random number on [0,0x7fffffff]-interval */
extern long genrand_int31(void);

/* generates a random number on [0,1]-real-interval */
extern double genrand_real1(void);

/* generates a random number on [0,1)-real-interval */
extern double genrand_real2(void);

/* generates a random number on (0,1)-real-interval */
extern double genrand_real3(void);

/* generates a random number on [0,1) with 53-bit resolution*/
extern double genrand_res53(void);

#endif /*MERSENNETWISTER_H_*/
