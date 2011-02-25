#ifndef RNGLIB_H

extern "C" void __stdcall MERSENNETWISTER_mp_SGRND(int* seed);
extern "C" double __stdcall MERSENNETWISTER_mp_GRND(void);
extern "C" int __stdcall MERSENNETWISTER_mp_IGRND(void);
extern "C" int __stdcall MERSENNETWISTER_mp_IGNPOI(float *mu);

namespace RNGLIB 
{
inline void sgrnd(int *seed)
{
    MERSENNETWISTER_mp_SGRND(seed);
}

inline double grnd(void)
{
    return MERSENNETWISTER_mp_GRND();
}

inline int igrnd(void)
{
    return MERSENNETWISTER_mp_IGRND();
}

inline int ignpoi(float *mu)
{
    return MERSENNETWISTER_mp_IGNPOI(mu);
}
}

#endif
