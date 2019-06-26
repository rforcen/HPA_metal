//
//  hpa.hpp
//  high precission aritmethics
//
//  Created by asd on 24/06/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//  based in http://www.nongnu.org/hpalib/

#ifndef hpa_hpp
#define hpa_hpp

#include <float.h>
#include <stdio.h>
#include <string>

using std::string;

#define XLITTLE_ENDIAN 1
#define XULONG_BITSIZE 64

static const int XDIM = 7; // 7 default value -> long double 

typedef unsigned short ushort;
typedef unsigned int uint;

// consts for XDIM==7
typedef struct _xpr { // hpa binary representation
    ushort nmm[XDIM + 1];
} xpr;

static const xpr
xZero = {{0x0, 0x0}}, xOne = {{0x3fff, 0x8000}}, xTwo = {{0x4000, 0x8000}},
xTen = {{0x4002, 0xa000}},xPinf = {{0x7fff, 0x0}},xMinf = {{0xffff, 0x0}},
xVSV = {{0x3ff2, 0x8000}},xVGV = {{0x4013, 0x8000}},
xEmax = {{0x400c, 0xb16c}}, xEmin = {{0xc00c, 0xb16c}},
xE2min = {{0xc00c, 0xfffb}}, /* -16382.75 */ xE2max = {{0x400c, 0xfffb}}, /* +16382.75 */

xPi4 = {{0x3ffe, 0xc90f, 0xdaa2, 0x2168, 0xc234, 0xc4c6, 0x628b, 0x80dc}},
xPi2 = {{0x3fff, 0xc90f, 0xdaa2, 0x2168, 0xc234, 0xc4c6, 0x628b, 0x80dc}},
xPi = {{0x4000, 0xc90f, 0xdaa2, 0x2168, 0xc234, 0xc4c6, 0x628b, 0x80dc}},
xEe = {{0x4000, 0xadf8, 0x5458, 0xa2bb, 0x4a9a, 0xafdc, 0x5620, 0x273d}},
xLn2 = {{0x3ffe, 0xb172, 0x17f7, 0xd1cf, 0x79ab, 0xc9e3, 0xb398, 0x03f3}},
xLn10 = {{0x4000, 0x935d, 0x8ddd, 0xaaa8, 0xac16, 0xea56, 0xd62b, 0x82d3}},
xSqrt2 = {{0x3fff, 0xb504, 0xf333, 0xf9de, 0x6484, 0x597d, 0x89b3, 0x754b}},
xLog2_e = {{0x3fff, 0xb8aa, 0x3b29, 0x5c17, 0xf0bb, 0xbe87, 0xfed0, 0x691d}},
xLog2_10 = {{0x4000, 0xd49a, 0x784b, 0xcd1b, 0x8afe, 0x492b, 0xf6ff, 0x4db0}},
xLog10_e = {{0x3ffd, 0xde5b, 0xd8a9, 0x3728, 0x7195, 0x355b, 0xaaaf, 0xad34}},
xRndcorr = {{0x3ffe, 0x8000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0068}},
xFixcorr = {{0x3f97, 0xc000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000}},
xNaN = {{0x0000, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff}};

static const int XERR_DFL = 1,XENONE = 0,XEDIV = 1,XEDOM = 2,XEBADEXP = 3,XFPOFLOW = 4; // Floating point overflow
static const int XNERR = 4, XEINV = 5, XMAX_10EX = 4931;
static const ushort xM_sgn = 0x8000, xM_exp = 0x7fff;
static const short xBias = 16383,xD_bias = 15360 ,xD_max = 2047,xD_lex = 12,xF_bias = 16256,xF_max = 255,xF_lex = 9;
static const short xMax_p = 16 * XDIM,xK_lin = -8 * XDIM;
static const int xItt_div = 2, xK_tanh = 5, xMS_exp = 21, xMS_hyp = 25,xMS_trg = 31;
static const double __Log10_2__ = .3010299956639812;



class HPA {
    
public:
    // constructors
    HPA(const xpr *px = &xZero) : br(*px) {}
    HPA(xpr x) : br(x) {}
    HPA(double x) { br = dbltox(x); }
    HPA(float x) { br = flttox(x); }
    HPA(int n) { br = inttox(n); }
    HPA(long n) { br = inttox(n); }
    HPA(uint u) { br = uinttox(u); }
    HPA(unsigned long u) { br = uinttox(u); }
    HPA(const char *str, char **endptr = 0) { br = strtox(str, endptr); }
    HPA(string str) { br = atox(str.c_str()); }
    HPA(const HPA &x) { br = x.br; }
    
    // assign
    HPA &operator=(const double d) {
        br = dbltox(d);
        return *this;
    }
    HPA &operator=(const float d) {
        br = flttox(d);
        return *this;
    }
    HPA &operator=(const int d) {
        br = inttox(d);
        return *this;
    }
    HPA &operator=(const long d) {
        br = inttox(d);
        return *this;
    }
    HPA &operator=(char *d) {
        br = strtox(d, 0);
        return *this;
    }
    HPA &operator=(string &d) {
        br = strtox(d.c_str(), 0);
        return *this;
    }
    
    HPA &operator=(HPA x) {
        br = x.br;
        return *this;
    }
    
    // aritmetic
    HPA operator+(const HPA &x) { return HPA(xsum(br, x.br)); }
    HPA operator-(const HPA &x) { return HPA(xsub(br, x.br)); }
    HPA operator*(const HPA &x) { return HPA(xmul(br, x.br)); }
    HPA operator/(const HPA &x) { return HPA(xdiv(br, x.br)); }
    HPA operator-() { return xneg(br); }
    
    HPA &operator+=(const HPA &x) {    br = xsum(br, x.br);    return *this;    }
    HPA &operator-=(const HPA &x) {    br = xsub(br, x.br);    return *this;    }
    HPA &operator*=(const HPA &x) {    br = xmul(br, x.br);    return *this;    }
    HPA &operator/=(const HPA &x) {    br = xdiv(br, x.br);    return *this;    }
    
    // compare
    bool operator==(const HPA &x) { return xeq(br, x.br); }
    bool operator!=(const HPA &x) { return xneq(br, x.br); }
    bool operator>(const HPA &x) { return xgt(br, x.br); }
    bool operator<(const HPA &x) { return xlt(br, x.br); }
    bool operator>=(const HPA &x) { return xge(br, x.br); }
    bool operator<=(const HPA &x) { return xle(br, x.br); }
    
    // tras. func's
    HPA exp() { return HPA(xexp(br)); }
    HPA exp10() { return HPA(xexp10(br)); }
    HPA exp2() { return HPA(xexp2(br)); }
    HPA log() { return HPA(xlog(br)); }
    HPA log10() { return HPA(xlog10(br)); }
    HPA log2() { return HPA(xlog2(br)); }
    
    HPA pow(HPA &y) { return HPA(xpow(br, y.br)); }
    HPA pwr(int n) { return HPA(xpwr(br, n)); }
    
    HPA fmod(HPA &y) {
        xpr q;
        return HPA(xfmod(br, y.br, &q));
    }
    // trigs.
    HPA sin() { return xsin(br); }
    HPA cos() { return xcos(br); }
    HPA tan() { return xtan(br); }
    
    // conversion
    int sign() { return xsgn(&br); }
    float tofloat() { return xtoflt(br); }
    double todouble() { return xtodbl(br); }
    char *tochar(int lim = 9, int sc_not=1, int sgn=0) { return xpr_asprint(br, sc_not, sgn, lim); }
    bool isNaN() { return xisNaN(&br); }
    bool isInf() { return xisPinf(&br); }
    
    // misc
    HPA abs() { return xabs(br); }
    
private:
    xpr xsum(xpr a, xpr b) { return xadd(a, b, 0); }
    xpr xsub(xpr a, xpr b) { return xadd(a, b, 1); }
    
    xpr xadd(xpr s, xpr t, int f) {
        ushort pe[XDIM + 1], h, u;
        ushort *pa, *pb, *pc, *pf = pe;
        uint n = 0;
        short e;
        short k;
        
        pa = (ushort *)&s;
        pb = (ushort *)&t;
        e = *pa & xM_exp;
        k = *pb & xM_exp;
        if (f)
            *pb ^= xM_sgn;
        u = (*pb ^ *pa) & xM_sgn;
        f = 0;
        if (e > k) {
            if ((k = e - k) >= xMax_p)
                return s;
            xrshift(k, pb + 1, XDIM);
        } else if (e < k) {
            if ((e = k - e) >= xMax_p)
                return t;
            xrshift(e, pa + 1, XDIM);
            e = k;
            pc = pa;
            pa = pb;
            pb = pc;
        } else if (u) {
            for (pc = pa, pf = pb; *(++pc) == *(++pf) && f < XDIM; ++f)
                ;
            if (f >= XDIM)
                return xZero;
            if (*pc < *pf) {
                pc = pa;
                pa = pb;
                pb = pc;
            }
            pf = pe + f;
        }
        h = *pa & xM_sgn;
        if (u) {
            for (pc = pb + XDIM; pc > pb; --pc)
                *pc = ~(*pc);
            n = 1L;
        }
        for (pc = pe + XDIM, pa += XDIM, pb += XDIM; pc > pf;) {
            n += *pa;
            pa--;
            n += *pb;
            pb--;
            *pc = n;
            pc--;
            n >>= 16;
        }
        if (u) {
            for (; *(++pc) == 0; ++f)
                ;
            for (k = 0; !((*pc << k) & xM_sgn); ++k)
                ;
            if ((k += 16 * f)) {
                if ((e -= k) <= 0)
                    return xZero;
                xlshift(k, pe + 1, XDIM);
            }
        } else {
            if (n) {
                ++e;
                if ((xsigerr(e == (short)xM_exp, XFPOFLOW, NULL)))
                    return (!h ? xPinf : xMinf);
                ++pf;
                xrshift(1, pf, XDIM);
                *pf |= xM_sgn;
            }
        }
        *pe = e;
        *pe |= h;
        return *(xpr *)pe;
    }
    
    xpr xmul(xpr s, xpr t) {
        ushort pe[XDIM + 2], *q0, *q1, h;
        ushort *pa, *pb, *pc;
        uint m, n, p;
        short e;
        short k;
        
        q0 = (ushort *)&s;
        q1 = (ushort *)&t;
        e = (*q0 & xM_exp) - xBias;
        k = (*q1 & xM_exp) + 1;
        if ((xsigerr(e > (short)xM_exp - k, XFPOFLOW, NULL)))
            return (((s.nmm[0] & xM_sgn) ^ (t.nmm[0] & xM_sgn)) ? xMinf : xPinf);
        if ((e += k) <= 0)
            return xZero;
        h = (*q0 ^ *q1) & xM_sgn;
        for (++q1, k = XDIM, p = n = 0L, pc = pe + XDIM + 1; k > 0; --k) {
            for (pa = q0 + k, pb = q1; pa > q0;) {
                m = *pa--;
                m *= *pb++;
                n += (m & 0xffffL);
                p += (m >> 16);
            }
            *pc-- = n;
            n = p + (n >> 16);
            p = 0L;
        }
        *pc = n;
        if (!(*pc & xM_sgn)) {
            --e;
            if (e <= 0)
                return xZero;
            xlshift(1, pc, XDIM + 1);
        }
        if ((xsigerr(e == (short)xM_exp, XFPOFLOW, NULL)))
            return (!h ? xPinf : xMinf);
        *pe = e;
        *pe |= h;
        return *(xpr *)pe;
    }
    
    xpr xdiv(xpr s, xpr t) {
        xpr a;
        ushort *pc, e, i;
        
        pc = (ushort *)&t;
        e = *pc;
        *pc = xBias;
        if ((xsigerr(xprcmp(&t, &xZero) == 0, XEDIV, "xdiv()")))
            return xZero;
        else {
            a = dbltox(1 / xtodbl(t));
            *pc = e;
            pc = (ushort *)&a;
            *pc += xBias - (e & xM_exp);
            *pc |= e & xM_sgn;
            for (i = 0; i < xItt_div; ++i)
                a = xmul(a, xadd(xTwo, xmul(a, t), 1));
            return xmul(s, a);
        }
    }
    
    int xprcmp(const xpr *pa, const xpr *pb) {
        ushort e, k, *p, *q, p0, q0;
        int m;
        
        p = (ushort *)pa;
        q = (ushort *)pb;
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        if (m > XDIM && (*p & xM_exp) < xM_exp)
        /* *pa is actually zero */
            p0 = 0;
        else
            p0 = *p;
        for (m = 1; m <= XDIM && q[m] == 0; m++)
            ;
        if (m > XDIM && (*q & xM_exp) < xM_exp)
        /* *pb is actually zero */
            q0 = 0;
        else
            q0 = *q;
        e = p0 & xM_sgn;
        k = q0 & xM_sgn;
        if (e && !k)
            return -1;
        else if (!e && k)
            return 1;
        else /* *pa and *pb have the same sign */
        {
            m = (e) ? -1 : 1;
            e = p0 & xM_exp;
            k = q0 & xM_exp;
            if (e > k)
                return m;
            else if (e < k)
                return -m;
            else {
                for (e = 0; *(++p) == *(++q) && e < XDIM; ++e)
                    ;
                if (e < XDIM)
                    return (*p > *q ? m : -m);
                else
                    return 0;
            }
        }
    }
    
    int xeq(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) == 0); }
    
    int xneq(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) != 0); }
    
    int xgt(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) > 0); }
    
    int xge(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) >= 0); }
    
    int xlt(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) < 0); }
    
    int xle(xpr x1, xpr x2) { return (xprcmp(&x1, &x2) <= 0); }
    
    /*
     xisNaN (&x) returns 1 if and only if x is not a valid number
     */
    
    int xisNaN(const xpr *u) {
        ushort *p = (ushort *)u;
        
        if ((*p))
            return 0;
        else {
            int i;
            
            for (i = 1; i <= XDIM && p[i] == 0x0; i++)
                ;
            return (i <= XDIM ? 1 : 0);
        }
    }
    
    int xis0(const xpr *u) {
        ushort *p = (ushort *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        return (m > XDIM && (*p & xM_exp) < xM_exp ? 1 : 0);
    }
    
    int xnot0(const xpr *u) {
        ushort *p = (ushort *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        return (m > XDIM && (*p & xM_exp) < xM_exp ? 0 : 1);
    }
    
    int xsgn(xpr *u) {
        ushort *p = (ushort *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        if ((m > XDIM && (*p & xM_exp) < xM_exp) || !*p)
            return 0;
        else
            return ((*p & xM_sgn) ? -1 : 1);
    }
    
    int xisPinf(const xpr *u) { return (*u->nmm == xM_exp ? 1 : 0); }
    
    int xisMinf(const xpr *u) { return (*u->nmm == (xM_exp | xM_sgn) ? 1 : 0); }
    
    int xisordnumb(const xpr *u) {
        int isNaN, isfinite;
        ushort *p = (ushort *)u;
        
        if ((*p))
            isNaN = 0;
        else {
            int i;
            
            for (i = 1; i <= XDIM && p[i] == 0x0; i++)
                ;
            isNaN = i <= XDIM;
        }
        isfinite = (*p & xM_exp) < xM_exp;
        return (!isNaN && (isfinite) ? 1 : 0);
    }
    
    int xsigerr(int errcond, int errcode, const char *where) {
        if (!errcond)
            errcode = 0;
        if (errcode < 0 || errcode > XNERR)
            errcode = XEINV;
        return errcode;
    }
    void xlshift(int n, ushort *pm, int m) {
        ushort *pa, *pc;
        
        pc = pm + m - 1;
        if (n < 16 * m) {
            pa = pm + n / 16;
            m = n % 16;
            n = 16 - m;
            while (pa < pc) {
                *pm = (*pa++) << m;
                *pm++ |= *pa >> n;
            }
            *pm++ = *pa << m;
        }
        while (pm <= pc)
            *pm++ = 0;
    }
    
    void xrshift(int n, ushort *pm, int m) {
        ushort *pa, *pc;
        
        pc = pm + m - 1;
        if (n < 16 * m) {
            pa = pc - n / 16;
            m = n % 16;
            n = 16 - m;
            while (pa > pm) {
                *pc = (*pa--) >> m;
                *pc-- |= *pa << n;
            }
            *pc-- = *pa >> m;
        }
        while (pc >= pm)
            *pc-- = 0;
    }
    
    xpr dbltox(double y) {
        ushort pe[XDIM + 1], *pc, u;
        short i, e;
        
        if (y < DBL_MIN && y > -DBL_MIN)
            return xZero;
        pc = (ushort *)&y;
        /* Change made by Ivano Primi - 11/19/2004 */
#ifdef XLITTLE_ENDIAN
        pc += 3;
#endif
        u = *pc & xM_sgn;
        e = xD_bias + ((*pc & xM_exp) >> (16 - xD_lex));
        /* Change made by Ivano Primi - 11/19/2004 */
#ifdef XLITTLE_ENDIAN
        for (i = 1; i < 5; pe[i] = *pc--, i++)
            ;
#else
        for (i = 1; i < 5; pe[i] = *pc++, i++)
            ;
#endif
        while (i <= XDIM)
            pe[i++] = 0;
        pc = pe + 1;
        xlshift(xD_lex - 1, pc, 4);
        *pc |= xM_sgn;
        *pe = e;
        *pe |= u;
        return *(xpr *)pe;
    }
    
    double xtodbl(xpr s) {
        ushort pe[4], *pc, u;
        short i, e;
        
        pc = (ushort *)&s;
        u = *pc & xM_sgn;
        e = (*pc & xM_exp) - xD_bias;
        if (e >= xD_max)
            return (!u ? DBL_MAX : -DBL_MAX);
        if (e <= 0)
            return 0.;
        for (i = 0; i < 4; pe[i] = *++pc, i++)
            ;
        pe[0] &= xM_exp;
        xrshift(xD_lex - 1, pe, 4);
        pe[0] |= (e << (16 - xD_lex));
        pe[0] |= u;
        /* Change made by Ivano Primi - 11/19/2004 */
#ifdef XLITTLE_ENDIAN
        u = pe[3];
        pe[3] = pe[0];
        pe[0] = u;
        u = pe[2];
        pe[2] = pe[1];
        pe[1] = u;
#endif
        return *(double *)pe;
    }
    
    inline unsigned long _abs(long n) { return ((n) >= 0 ? (n) : -(n)); }
    
    xpr inttox(long n) {
        ushort pe[XDIM + 1], *pc;
        short e;
        unsigned long k, h;
        
        k = _abs(n);
        pc = (ushort *)&k;
        for (e = 0; e <= XDIM; pe[e++] = 0)
            ;
        if (n == 0)
            return *(xpr *)pe;
        
#if (XULONG_BITSIZE == 64)
        
#ifdef XLITTLE_ENDIAN
        pe[1] = *(pc + 3);
        pe[2] = *(pc + 2);
        pe[3] = *(pc + 1);
        pe[4] = *pc;
#else
        pe[1] = *pc;
        pe[2] = *(pc + 1);
        pe[3] = *(pc + 2);
        pe[4] = *(pc + 3);
#endif
        
#else /* (XULONG_BITSIZE == 32) */
        
#ifdef XLITTLE_ENDIAN
        pe[1] = *(pc + 1);
        pe[2] = *pc;
#else
        pe[1] = *pc;
        pe[2] = *(pc + 1);
#endif
        
#endif
        
        for (e = 0, h = 1; h <= k && e < (XULONG_BITSIZE - 1); h <<= 1, ++e)
            ;
        if (h <= k)
            e += 1;
        *pe = xBias + e - 1;
        if (n < 0)
            *pe |= xM_sgn;
        xlshift(XULONG_BITSIZE - e, pe + 1, XDIM);
        return *(xpr *)pe;
    }
    
    xpr uinttox(unsigned long n) {
        ushort pe[XDIM + 1], *pc;
        short e;
        unsigned long h;
        
        pc = (ushort *)&n;
        for (e = 0; e <= XDIM; pe[e++] = 0)
            ;
        if (n == 0)
            return *(xpr *)pe;
        
#if (XULONG_BITSIZE == 64)
        
#ifdef XLITTLE_ENDIAN
        pe[1] = *(pc + 3);
        pe[2] = *(pc + 2);
        pe[3] = *(pc + 1);
        pe[4] = *pc;
#else
        pe[1] = *pc;
        pe[2] = *(pc + 1);
        pe[3] = *(pc + 2);
        pe[4] = *(pc + 3);
#endif
        
#else /* (XULONG_BITSIZE == 32) */
        
#ifdef XLITTLE_ENDIAN
        pe[1] = *(pc + 1);
        pe[2] = *pc;
#else
        pe[1] = *pc;
        pe[2] = *(pc + 1);
#endif
        
#endif
        
        for (e = 0, h = 1; h <= n && e < (XULONG_BITSIZE - 1); h <<= 1, ++e)
            ;
        if (h <= n)
            e += 1;
        *pe = xBias + e - 1;
        xlshift(XULONG_BITSIZE - e, pe + 1, XDIM);
        return *(xpr *)pe;
    }
    
    float xtoflt(xpr s) {
        ushort pe[2], *pc, u;
        short i, e;
        
        pc = (ushort *)&s;
        u = *pc & xM_sgn;
        e = (*pc & xM_exp) - xF_bias;
        /*
         u is the sign of the number s.
         e == (exponent of s) + 127
         */
        if (e >= xF_max)
            return (!u ? FLT_MAX : -FLT_MAX);
        if (e <= 0)
            return 0.;
        for (i = 0; i < 2; pe[i] = *++pc, i++)
            ;
        
        /* In the IEEE 754 Standard the leading 1 */
        /* is not represented.                    */
        pe[0] &= xM_exp;
        /* Now in pe[0],pe[1] we have 31 bits of mantissa. */
        /* But only the first 23 ones must be put in the   */
        /* final float number.                             */
        xrshift(xF_lex - 1, pe, 2);
        /* We have just loaded the mantissa and now we */
        /* are going to load exponent and sign.        */
        pe[0] |= (e << (16 - xF_lex));
        pe[0] |= u;
#ifdef XLITTLE_ENDIAN
        u = pe[0];
        pe[0] = pe[1];
        pe[1] = u;
#endif
        return *(float *)pe;
    }
    
    xpr flttox(float y) {
        ushort pe[XDIM + 1], *pc, u;
        short i, e;
        
        if (y < FLT_MIN && y > -FLT_MIN)
            return xZero;
        pc = (ushort *)&y;
#ifdef XLITTLE_ENDIAN
        pc += 1;
#endif
        u = *pc & xM_sgn;
        e = xF_bias + ((*pc & xM_exp) >> (16 - xF_lex));
        /*
         Now u is the sign of y and e is the
         biased exponent (exponent + bias).
         */
#ifdef XLITTLE_ENDIAN
        for (i = 1; i < 3; pe[i] = *pc--, i++)
            ;
#else
        for (i = 1; i < 3; pe[i] = *pc++, i++)
            ;
#endif
        while (i <= XDIM)
            pe[i++] = 0;
        pc = pe + 1;
        xlshift(xF_lex - 1, pc, 2);
        *pc |= xM_sgn;
        /* We have just put in pe[1],pe[2] the whole */
        /* mantissa of y with a leading 1.           */
        /* Now we have only to put exponent and sign */
        /* in pe[0].                                 */
        *pe = e;
        *pe |= u;
        return *(xpr *)pe;
    }
    
#define UPPER_BOUND (XMAX_10EX + 100)
    
    xpr strtox(const char *q, char **endptr) {
        xpr s, f;
        ushort pc[XDIM + 1];
        ushort *pn, *pf, *pa, *pb;
        ushort sfg, ibex, fbex;
        unsigned long n;
        long j;
        int idex = 0, fdex = 0, c, m;
        short noip, nofp;
        const char *ptr;
        
        pn = (ushort *)&s;
        pf = (ushort *)&f;
        for (j = 0; j <= XDIM; pn[j] = pf[j] = pc[j] = 0, ++j)
            ;
        sfg = 0;
        m = XDIM + 1;
        if ((endptr))
            *endptr = (char *)q;
        /* Skip the leading spaces if there are some */
        while ((isspace(*q)))
            q++;
        /* Sign */
        if (*q == '+')
            ++q;
        else if (*q == '-') {
            sfg = 0x8000;
            ++q;
        }
        /* Integer part */
        for (ptr = q; (c = *q - '0') >= 0 && c <= 9; ++q) {
            if (pn[0])
                ++idex;
            else {
                xlshift(1, pn, m);
                for (j = 0; j < m; ++j)
                    pc[j] = pn[j];
                xlshift(2, pn, m);
                for (n = (uint)c, pa = pn + XDIM, pb = pc + XDIM; pa >= pn;
                     pa--, pb--) {
                    n += *pa + *pb;
                    *pa = n;
                    n >>= 16;
                }
            }
        }
        for (j = 0; j < m && pn[j] == 0; ++j)
            ;
        if (j == m)
            ibex = 0;
        else {
            ibex = xBias + xMax_p - 1;
            if (j) {
                j <<= 4;
                ibex -= j;
                xlshift((int)j, pn, m);
            }
            while (pn[0]) {
                xrshift(1, pn, m);
                ++ibex;
            }
            pn[0] = ibex | sfg;
        }
        noip = ptr == q;
        /* End Integer part */
        if (*q == '.') {
            /* Fractionary part */
            for (j = 0; j <= XDIM; ++j)
                pc[j] = 0;
            for (ptr = ++q; (c = *q - '0') >= 0 && c <= 9 && pf[0] == 0;
                 --fdex, ++q) {
                xlshift(1, pf, m);
                for (j = 0; j < m; ++j)
                    pc[j] = pf[j];
                xlshift(2, pf, m);
                for (n = (uint)c, pa = pf + XDIM, pb = pc + XDIM; pa >= pf;
                     pa--, pb--) {
                    n += *pa + *pb;
                    *pa = n;
                    n >>= 16;
                }
            }
            for (j = 0; j < m && pf[j] == 0; ++j)
                ;
            if (j == m)
                fbex = 0;
            else {
                fbex = xBias + xMax_p - 1;
                if (j) {
                    j <<= 4;
                    fbex -= j;
                    xlshift((int)j, pf, m);
                }
                while (pf[0]) {
                    xrshift(1, pf, m);
                    ++fbex;
                }
                pf[0] = fbex | sfg;
            }
            nofp = ptr == q;
        } /* end if (*q == '.') */
        else
            nofp = 1;
        if ((noip) && (nofp))
        /* Error ! */
            return xNaN;
        else {
            /*
             Added on August 4th 2007:
             If the decimal digits after the radix character ('.')
             are too much with respect to the given precision,
             all the meaningless ones must be neglected.
             The absence of this loop produced a VERY hazardous bug
             (thanks to Shina Tan <tansn ~at~ phys ~dot~ washington ~dot~ edu>
             for the bug notification)
             */
            for (; *q >= '0' && *q <= '9'; q++)
                ;
            /*
             End addition
             */
            if ((endptr))
                *endptr = (char *)q;
        }
        /* Exponent */
        if (*q == 'e' || *q == 'E') {
            ++q;
            sfg = 0;
            if (*q == '+')
                ++q;
            else if (*q == '-') {
                sfg = 1;
                ++q;
            }
            for (ptr = q, j = 0; (c = *q - '0') >= 0 && c <= 9 && j <= UPPER_BOUND;
                 ++q) {
                j <<= 1;
                m = (int)j;
                j <<= 2;
                j += c + m;
            }
            if (ptr != q && (endptr))
                *endptr = (char *)q;
            if (sfg)
                j = -j;
            idex += j;
            fdex += j;
        }
        /*
         Remark: s and f have the same sign (see above).
         */
        if (idex > XMAX_10EX || fdex > XMAX_10EX)
            return ((s.nmm[0] & xM_sgn) ? xMinf : xPinf);
        else {
            if (idex)
                s = xmul(s, xpwr(xTen, idex));
            if (fdex)
                f = xmul(f, xpwr(xTen, fdex));
            return xadd(s, f, 0);
        }
    }
    
    xpr atox(const char *q) { return strtox(q, NULL); }
    
    xpr xpwr(xpr s, int n) {
        xpr t;
        unsigned k, m;
        
        t = xOne;
        if (n < 0) {
            m = -n;
            if ((xsigerr(xprcmp(&s, &xZero) == 0, XEBADEXP, "xpwr()")))
                return xZero;
            s = xdiv(xOne, s);
        } else
            m = n;
        if ((m)) {
            k = 1;
            while (1) {
                if ((k & m))
                    t = xmul(s, t);
                if ((k <<= 1) <= m)
                    s = xmul(s, s);
                else
                    break;
            }
        } else
            xsigerr(xprcmp(&s, &xZero) == 0, XEBADEXP, "xpwr()");
        return t;
    }
    
    xpr xpr2(xpr s, int m) {
        ushort *p = (ushort *)&s;
        long e;
        
        for (e = 1; e <= XDIM && p[e] == 0; e++)
            ;
        if (e <= XDIM) {
            e = *p & xM_exp; /* biased exponent */
            if (e + m < 0)
                return xZero;
            else if ((xsigerr(e + m >= xM_exp, XFPOFLOW, NULL)))
                return ((s.nmm[0] & xM_sgn) ? xMinf : xPinf);
            else {
                *p += m;
                return s;
            }
        } else /* s is zero or +-Inf */
            return s;
    }
    
    xpr xpow(xpr x, xpr y) {
        if (xsigerr((x_neg(&x)) || x_exp(&x) == -xBias, XEDOM, "xpow()"))
            return xZero;
        else
            return xexp2(xmul(xlog2(x), y));
    }
    
    xpr xneg(xpr s) {
        ushort *p = (ushort *)&s;
        
        *p ^= xM_sgn;
        return s;
    }
    
    xpr xabs(xpr s) {
        ushort *p = (ushort *)&s;
        
        *p &= xM_exp;
        return s;
    }
    
    int x_exp(const xpr *ps) {
        ushort *q = (ushort *)ps;
        
        return (*q & xM_exp) - xBias;
    }
    
    int x_neg(const xpr *ps) {
        ushort *q = (ushort *)ps;
        
        return (*q & xM_sgn);
    }
    
    xpr xlog(xpr z) {
        xpr f, h;
        int k, m;
        
        if ((xsigerr((x_neg(&z)) || x_exp(&z) == -xBias, XEDOM, "xlog()")))
            return xMinf;
        else if (xprcmp(&z, &xOne) == 0)
            return xZero;
        else {
            z = xfrexp(z, &m);
            z = xmul(z, xSqrt2);
            z = xdiv(xadd(z, xOne, 1), xadd(z, xOne, 0));
            h = xpr2(z, 1);
            z = xmul(z, z);
            for (f = h, k = 1; x_exp(&h) > -xMax_p;) {
                h = xmul(h, z);
                f = xadd(f, xdiv(h, inttox(k += 2)), 0);
            }
            return xadd(f, xmul(xLn2, dbltox(m - .5)), 0);
        }
    }
    
    xpr xlog2(xpr z) {
        xpr f, h;
        int k, m;
        
        if ((xsigerr((x_neg(&z)) || x_exp(&z) == -xBias, XEDOM, "xlog2()")))
            return xMinf;
        else if (xprcmp(&z, &xOne) == 0)
            return xZero;
        else {
            z = xfrexp(z, &m);
            z = xmul(z, xSqrt2);
            z = xdiv(xadd(z, xOne, 1), xadd(z, xOne, 0));
            h = xpr2(z, 1);
            z = xmul(z, z);
            for (f = h, k = 1; x_exp(&h) > -xMax_p;) {
                h = xmul(h, z);
                f = xadd(f, xdiv(h, inttox(k += 2)), 0);
            }
            return xadd(xmul(f, xLog2_e), dbltox(m - .5), 0);
        }
    }
    
    xpr xlog10(xpr z) {
        xpr w = xlog(z);
        
        if (xprcmp(&w, &xMinf) <= 0)
            return xMinf;
        else
            return xmul(w, xLog10_e);
    }
    
    xpr xsfmod(xpr s, int *p) {
        ushort *pa, *pb;
        short e, k;
        
        pa = (ushort *)&s;
        pb = pa + 1;
        e = (*pa & xM_exp) - xBias;
        if ((xsigerr(e >= 15, XFPOFLOW, NULL))) {
            *p = -1;
            return s;
        } else if (e < 0) {
            *p = 0;
            return s;
        }
        *p = *pb >> (15 - e);
        xlshift(++e, pb, XDIM);
        *pa -= e;
        for (e = 0; *pb == 0 && e < xMax_p; ++pb, e += 16)
            ;
        if (e == xMax_p)
            return xZero;
        for (k = 0; !((*pb << k) & xM_sgn); ++k)
            ;
        if ((k += e)) {
            xlshift(k, pa + 1, XDIM);
            *pa -= k;
        }
        return s;
    }
    
    xpr xexp2(xpr x) {
        xpr s, d, f;
        ushort *pf = (ushort *)&x;
        int m, k;
        
        if (xprcmp(&x, &xE2min) < 0)
            return xZero;
        else if ((xsigerr(xprcmp(&x, &xE2max) > 0, XFPOFLOW, NULL)))
            return xPinf;
        else {
            m = (*pf & xM_sgn) ? 1 : 0;
            x = xsfmod(x, &k);
            if ((m))
                k *= -1;
            /* -xBias <= k <= +xBias */
            x = xmul(x, xLn2);
            if (x_exp(&x) > -xBias) {
                x = xpr2(x, -1);
                s = xmul(x, x);
                f = xZero;
                for (d = inttox(m = xMS_exp); m > 1; m -= 2, d = inttox(m))
                    f = xdiv(s, xadd(d, f, 0));
                f = xdiv(x, xadd(d, f, 0));
                f = xdiv(xadd(d, f, 0), xadd(d, f, 1));
            } else
                f = xOne;
            pf = (ushort *)&f;
            if (-k > *pf)
                return xZero;
            else {
                *pf += k;
                if ((xsigerr(*pf >= xM_exp, XFPOFLOW, NULL)))
                    return xPinf;
                else
                    return f;
            }
        }
    }
    
    xpr xexp(xpr z) { return xexp2(xmul(z, xLog2_e)); }
    
    xpr xexp10(xpr z) { return xexp2(xmul(z, xLog2_10)); }
    
    xpr xfmod(xpr s, xpr t, xpr *q) {
        if ((xsigerr(xprcmp(&t, &xZero) == 0, XEDIV, "xfmod()")))
            return xZero;
        else {
            ushort *p, mask = 0xffff;
            short e, i;
            int u;
            
            *q = xdiv(s, t);
            p = (ushort *)q;
            u = (*p & xM_sgn) ? 0 : 1;
            e = (*p &= xM_exp); /* biased exponent of *q */
            e = e < xBias ? 0 : e - xBias + 1;
            for (i = 1; e / 16 > 0; i++, e -= 16)
                ;
            if (i <= XDIM) {
                /* e = 0, ..., 15 */
                mask <<= 16 - e;
                p[i] &= mask;
                for (i++; i <= XDIM; p[i] = 0, i++)
                    ;
            }
            /* Now *q == abs(quotient of (s/t)) */
            return xadd(s, xmul(t, *q), u);
        }
    }
    
    xpr xfrexp(xpr s, int *p) {
        ushort *ps = (ushort *)&s, u;
        
        *p = (*ps & xM_exp) - xBias + 1;
        u = *ps & xM_sgn;
        *ps = xBias - 1;
        *ps |= u;
        return s;
    }
    
    /*
     Remark: 's' must have a size >= 5 !
     */
    
    int copied_special_value(char *s, xpr u, int sign) {
        if ((xisPinf(&u))) {
            if ((sign)) {
                *s = '+';
                s++;
            }
            strcpy(s, "Inf");
            return 1;
        } else if ((xisMinf(&u))) {
            strcpy(s, "-Inf");
            return 1;
        } else if ((xisNaN(&u))) {
            if ((sign)) {
                *s = '\?';
                s++;
            }
            strcpy(s, "NaN");
            return 1;
        } else
            return 0;
    }
    
    void xstrputc(char c, char *buffer) {
        char *ptr;
        
        for (ptr = buffer; *ptr != '\0'; ptr++)
            ;
        *ptr = c;
    }
    
    void xsprintfmt(char *buffer, const char *fmt, ...) {
        char ibuff[1024];
        va_list ap;
        
        va_start(ap, fmt);
        vsprintf(ibuff, fmt, ap);
        va_end(ap);
        strcat(buffer, ibuff);
    }
    
    char *xpr_asprint(xpr u, int sc_not, int sign, int lim) {
        char q[5 * XDIM + 4], *buffer, *ptr;
        char *p = q;
        int k, m;
        int dig;
        ushort *pa = (ushort *)&u;
        const int BUFF_SIZE = 5120; /* 5 Kb */
        
        if (lim < 0)
            lim = 0;
        if (lim > 5 * XDIM + 2)
            lim = 5 * XDIM + 2;
        if (!(buffer = (char *)calloc(BUFF_SIZE, sizeof(char))))
            return NULL;
        else if ((copied_special_value(buffer, u, sign))) {
            for (k = 0; buffer[k] != '\0'; k++)
                ;
            /* Now k is the length of the buffer. */
            /* We shrink the buffer so that it has the exact */
            /* size to contain all its non null chars.       */
            ptr = (char *)realloc(buffer, k + 1);
            return (ptr != NULL) ? ptr : buffer;
        } else {
            if ((*pa & xM_sgn)) {
                *pa ^= xM_sgn;
                xstrputc('-', buffer);
            } else {
                if ((sign))
                    xstrputc('+', buffer);
            }
            if ((xis0(&u))) {
                xsprintfmt(buffer, "0.");
                for (k = 0; k < lim; ++k)
                    xstrputc('0', buffer);
                if ((sc_not))
                    xsprintfmt(buffer, "e+0");
            } else {
                m = ((*pa & xM_exp) - xBias);
                m = (int)((double)(m + 1) * __Log10_2__);
                if ((m))
                    u = xmul(u, xpwr(xTen, -m));
                while ((*pa & xM_exp) < xBias) {
                    --m;
                    u = xmul(u, xTen);
                }
                for (*p = 0, k = 0; k <= lim; ++k) {
                    u = xsfmod(u, &dig);
                    ++p;
                    *p = (char)dig;
                    if (*pa == 0)
                        break;
                    u = xmul(xTen, u);
                }
                for (; k <= lim; ++k)
                    *++p = 0;
                if ((*pa)) {
                    u = xsfmod(u, &dig);
                    if (dig >= 5)
                        ++(*p);
                    while (*p == 10) {
                        *p = 0;
                        ++(*--p);
                    }
                }
                p = q;
                if (*p == 0)
                    ++p;
                else
                    ++m;
                /* Now has come the moment to print */
                if (m > XMAX_10EX)
                    xsprintfmt(buffer, "Inf");
                else if ((sc_not)) {
                    xsprintfmt(buffer, "%c.", '0' + *p++);
                    for (k = 0; k < lim; ++k)
                        xstrputc('0' + *p++, buffer);
                    if (m >= 0)
                        xsprintfmt(buffer, "e+%d", m);
                    else
                        xsprintfmt(buffer, "e%d", m);
                } else {
                    if (m >= 0) {
                        for (k = 0; k <= m; k++) {
                            if (k <= lim)
                                xstrputc('0' + p[k], buffer);
                            else
                                xstrputc('0', buffer);
                        }
                        if (k <= lim) {
                            xstrputc('.', buffer);
                            for (; k <= lim; k++)
                                xstrputc('0' + p[k], buffer);
                        }
                    } else {
                        xsprintfmt(buffer, "0.");
                        for (k = 1; k < -m; k++)
                            xstrputc('0', buffer);
                        for (k = 0; k <= lim; ++k)
                            xstrputc('0' + *p++, buffer);
                    }
                }
            } /* End of *pa != 0 */
            for (k = 0; buffer[k] != '\0'; k++)
                ;
            /* Now k is the length of the buffer. */
            /* We shrink the buffer so that it has the exact */
            /* size to contain all its non null chars.       */
            ptr = (char *)realloc(buffer, k + 1);
            return (ptr != NULL) ? ptr : buffer;
        } /* End of buffer != 0 */
    }
    
    // trigs
    int xodd(xpr x) {
        ushort *p = (ushort *)&x;
        short e, i;
        
        e = (*p & xM_exp) - xBias; /* exponent of x */
        if (e < 0)
            return 0;
        else {
            for (i = 1; e / 16 > 0; i++, e -= 16)
                ;
            /* Now e = 0, ..., 15 */
            return (i <= XDIM ? p[i] & 0x8000 >> e : 0);
        }
    }
    
    xpr xtan(xpr z) {
        int k, m;
        
        z = rred(z, 't', &k);
        if ((xsigerr(xprcmp(&z, &xPi2) >= 0, XEDOM, "xtan()")))
            return (!k ? xPinf : xMinf);
        else {
            if (xprcmp(&z, &xPi4) == 1) {
                m = 1;
                z = xadd(xPi2, z, 1);
            } else
                m = 0;
            if ((k))
                z = xneg(c_tan(z));
            else
                z = c_tan(z);
            if (m)
                return xdiv(xOne, z);
            else
                return z;
        }
    }
    
    xpr xcos(xpr z) {
        int k;
        
        z = rred(z, 'c', &k);
        if (x_exp(&z) < xK_lin) {
            if ((k))
                return xneg(xOne);
            else
                return xOne;
        }
        z = c_tan(xpr2(z, -1));
        z = xmul(z, z);
        z = xdiv(xadd(xOne, z, 1), xadd(xOne, z, 0));
        if ((k))
            return xneg(z);
        else
            return z;
    }
    
    xpr xsin(xpr z) {
        int k;
        
        z = rred(z, 's', &k);
        if (x_exp(&z) >= xK_lin) {
            z = c_tan(xpr2(z, -1));
            z = xdiv(xpr2(z, 1), xadd(xOne, xmul(z, z), 0));
        }
        if ((k))
            return xneg(z);
        else
            return z;
    }
    
    xpr c_tan(xpr z) {
        xpr s, f, d;
        int m;
        ushort k;
        
        if (x_exp(&z) < xK_lin)
            return z;
        s = xneg(xmul(z, z));
        for (k = 1; k <= XDIM && s.nmm[k] == 0; k++)
            ;
        if ((xsigerr(s.nmm[0] == 0xffff && k > XDIM, XFPOFLOW, NULL)))
            return xZero;
        else {
            f = xZero;
            for (d = inttox(m = xMS_trg); m > 1;) {
                f = xdiv(s, xadd(d, f, 0));
                d = inttox(m -= 2);
            }
            return xdiv(z, xadd(d, f, 0));
        }
    }
    
    xpr rred(xpr z, int kf, int *ps) {
        xpr is, q;
        
        if (x_neg(&z)) {
            z = xneg(z);
            is = xOne;
        } else
            is = xZero;
        z = xfmod(z, xPi, &q);
        if (kf == 't')
            q = is;
        else if (kf == 's')
            q = xadd(q, is, 0);
        if (xprcmp(&z, &xPi2) == 1) {
            z = xadd(xPi, z, 1);
            if (kf == 'c' || kf == 't')
                q = xadd(q, xOne, 0);
        }
        *ps = (xodd(q)) ? 1 : 0;
        return z;
    }
    
    xpr br; // binary representation
};

// func's in .cpp
HPA sin(HPA x), cos(HPA x), tan(HPA x);
HPA exp(HPA x), exp10(HPA x), exp2(HPA x);
HPA log(HPA x), log10(HPA x), log2(HPA x);
HPA pow(HPA x, HPA y), pwr(HPA x, int n);
HPA abs(HPA x);

#endif /* hpa_hpp */
