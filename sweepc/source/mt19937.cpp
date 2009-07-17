
#include "mt19937.h"

#include <cassert>

/* initializes mt[N] with a seed */
void Sweep::Mt19937::init_genrand(unsigned long s)
{
    //std::cout << "Seeding with " << s;
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<MTN; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    //std::cout << " to give " << &mti << ' ' << mti << ' ' << mt[0] << ' ' << mt[1] << std::endl;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void Sweep::Mt19937::init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (MTN>key_length ? MTN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=MTN) { mt[0] = mt[MTN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MTN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=MTN) { mt[0] = mt[MTN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long Sweep::Mt19937::genrand_int32(void)
{
    //std::cout << "genrand_int32() " << &mti << ' ' << mti << ' ' << mt[0] << ' ' << mt[1] << '\n';;
    unsigned long y;
    const static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MTN) { /* generate N words at one time */
        //std::cout << "updating rng state\n";
        int kk;

        // Removed code that performs default initialisation as there
        // is now a default constructor, but check the condition hat
        // used to be here is still met.
        assert(mti <= MTN);

        //std::cout << "from " << mti << ' ' << mt[0] << ' ' << mt[1] << std::endl;
        for (kk=0;kk<MTN-MTM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MTM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<MTN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MTM-MTN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MTN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MTN-1] = mt[MTM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    //std::cout << "random int is " << y << std::endl;
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long Sweep::Mt19937::genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double Sweep::Mt19937::genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double Sweep::Mt19937::genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double Sweep::Mt19937::genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double Sweep::Mt19937::genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
}


// Create a static object to provide previous functionality via global methods
static Sweep::Mt19937 generator;

/**
 * Global method for backward compatibility.  Operates on global instance of
 * Mersenne Twister class defined above.
 *
 * \param[in]   s       New seed for generator
 */
void Sweep::init_genrand(unsigned long s) {
    generator.init_genrand(s);
}

/**
 * Global method for backward compatibility.  Uses on global instance of
 * Mersenne Twister class defined above.
 *
 *\return       U[0,1] random sample
 */
double Sweep::genrand_real1(void) {
    return generator.genrand_real1();
}
