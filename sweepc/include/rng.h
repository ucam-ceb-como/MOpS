/*
  Author(s):      Matthew Celnik (msc37)
  Project:        flopsy (population balance solver)

  File purpose:
    Random number generators.
*/

#ifndef SWEEP_RNG_H
#define SWEEP_RNG_H

#include "swp_params.h"
#include "mt19937ar.h"

namespace Sweep
{
    inline void srnd(unsigned long seed)
    {
        init_genrand(seed);
    };

    /* Returns a uniform random number [0,1]. */
    inline real rnd()
    {
        //return (real)RNGLIB::grnd();
        return genrand_real1();
        //return (real)(rand()+0.5) / (real)RAND_MAX;
    };

    /* Returns a uniform random integer in the given range. */
    inline int irnd(int min, int max)
    {
        //return (RNGLIB::igrnd() & (max-min)) + min;
        return (int)(rnd() * (real)(max-min)) + min;
    };

    /* From Numerical Recipes 6.1. */
    inline real gammln(real xx)
    {
        real x, y, tmp, ser;
        static real cof[6] = { 76.18009172947146,
                              -86.50532032941677,
                               24.01409824083091,
                              -1.231739572450155,
                               0.1208650973866179e-2,
                              -0.5395239384953e-5};
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x+0.5)*log(tmp);
        ser = 1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }

    /* From Numerical Recipes 7.3. */
    inline int ignpoi(real mu)
    {
        static real oldmu=(-1.0), g;
        real em=0.0, t=0.0, y=0.0, lnmu=0.0, sq=0.0;

        if (mu<12.0) {
            // Direct method.
            if (mu != oldmu) {
                oldmu = mu;
                g = exp(-mu);
            }
            em = -1;
            t = 1.0;
            do {
                ++em;
                t *= rnd();
            } while (t > g);
        } else {
            // Rejection method.
            if (mu!=oldmu) {
                oldmu = mu;
                sq = sqrt(2.0*mu);
                lnmu = log(mu);
                g = (mu * lnmu) - gammln(mu+1.0);
            }
            do {
                do {
                    // y = deviate from a Lorentzian comparison
                    // function.
                    y = tan(PI*rnd());
                    em = (sq*y) + mu;
                } while (em < 0.0);
                em=floor(em);
                t = 0.9 * (1.0+(y*y)) * exp(em*lnmu-gammln(em+1.0)-g);
            } while (rnd() > t);
        }
        return (int)em;
    }

/*
    inline int ignpoi2(real mu)
    {
        int ans;
        real b1,b2,del,difmuk,e,fk,fx,fy,g,px,py,t,u,v,x,xx,pp[35];
        static real c,c0,c1,c2,c3,s,d,omega,p,p0,q;
        int j,k,kflag;
        static int ll,l,m;
        static real a0=-0.5,a1=0.3333333,a2=-0.2500068,a3=0.2000118,
                    a4=-0.1661269,a5=0.1421878,a6=-0.1384794,a7=0.1250060;
        static real muprev=-1.0e37, muold=-1.0e37;
        static real fact[10] = {1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0};

        if (mu==muprev) {
            // continue.
        } else {
            if (mu<10.0) {
                //goto 420;
            } else {
                // mu has changed so recalculate s, d and ll.
                muprev = mu;
                s = sqrt(mu);
                d = 6.0 * mu * mu; 
                ll = (int)floor(mu-1.1484);
            }
        }

        // 310.

        // Normal sample - snorm() for standard normal deviate.
        g = mu + (s*snorm());
        if (g<0.0) {
            //goto 320;
        } else {
            // continue.
        }
        ans = (int)floor(g);

        // Return immediately if large enough.
        if (ans > ll) return ans;

        // Squeeze acceptance.
        fk = (real)ans;
        difmuk = mu - fk;
        u = rnd();
        if ((d*u) > (difmuk*difmuk*difmuk)) return ans;

        // 320.


//C     STEP P. PREPARATIONS FOR STEPS Q AND H.
//C             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
//C             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
//C             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
//C             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
//C             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

        if (mu==muold) {
            // goto 330.
        }
        
        muold = mu;
        omega = 0.3989423/s;
        b1    = 0.4166667e-1/mu;
        b2    = 0.3*b1*b1;
        c3    = 0.1428571*b1*b2;
        c2    = b2 - 15.0*c3;
        c1    = b1 - 6.0*b2 + 45.0*c3;
        c0    = 1.0 - b1 + 3.0*b2 - 15.0*c3;
        c     = 0.1069/mu;

        // 330.

        if (g < 0.0) {
            //goto 350.
        }

        kflag = 0;
        //goto 370.

        // 340.

        // Quotient acceptance.
        if ((fy - (u*fy)) < (py * exp(px-fx))) return ans;

        // 350.

        // Exponential sample - sexpo() for standard exponential deviatee and sample
        // t from the laplace "hat".
        e = sexpo();
        u = rnd();
        u = u + u - 1.0;
        t = 1.8 + (u>=0.0 ? e : -e);
        if (t < -0.6744) {
            //goto 350.
        }
        fk = floor(mu + (s*t));
        ans = (int)fk;
        difmuk = mu - fk;

        kflag = 1;
        //goto 370.

        // 360.

        // HAT acceptance (e is repeated on rejection).
        if ((c * abs(u)) > (py*exp(px+e) - fy*exp(fx+e))) {
            //goto 350.
        }

        return ans;

        // 370.

        // If ans is less than 10 then use factorials from array fact.
        if (ans > 10) {
            // goto 380.
        }

        px = -mu;
        py = pow(mu, ans/fact[ans]);
        //goto 410.

        // 380.

        // If ans is greater than 10 then use polynomial approximation
        // a0-a7 for accuracy when advisable.
        del = 0.8333333e-1/fk;
        del = del - (4.8*del*del*del);
        v = difmuk/fk;
        if (abs(v) < 0.25) {
            // goto 390.
        }
        px = fk*log(1.0+v) - difmuk - del;
        //goto 400.

        // 390.

        px = (fk * v * v * (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)) - del;

        // 400.

        py = 0.3989423 / sqrt(fk);

        // 410.

        x = (0.5 - difmuk) / s;
        xx = x * x;
        fx = -0.5*xx;
        fy = omega * (((c3*xx+c2)*xx+c1)*xx+c0);
        if (kflag<=0) {
            // goto 340.
        } else {
            // goto 360.
        }

        // 420.
        
        // Case B:  Start new table and calculate p0 if necessary.

        muprev = -1.0e37;
        if (mu!=muold) {
            if (mu>0.0) {
                muold = mu;
                l = 0;
                p = exp(-mu);
                q = p;
                p0 = p;
            } else {
                // Return an error.
                return -1;
            }
        }

        // 430.
        do {
            // Uniform sample for inversion method.
            u = rnd();
            if (u<p0) return 0;

            // Table comparison until the end pp(l) of the pp-table
            // of cumulative poisson probabilities (pp(8)==0.458 for mu==10).

            if (l!=0) {
                j = 1;

                if (u>0.458) j = min<int>(l,m);

                for (k=j-1; k<l; k++) {
                    if (u<pp(k)) {
                        return k;
                    }
                }
                if (l==34) {
                    continue;
                }
            }

            // Creation of new poisson probabilities p and their
            // cumulatives q=pp(k).
            for (k=l; k<35; k++) {
                p = p * mu / (real)k;
                q = q + p;
                pp(k) = q;
                if (u<q) {
                    l = k;
                    return k;
                }
            }

            l = 34;
        }
    }
    */
};

#endif
