//
//  hpa_metal.h
//  hpa
//
//  Created by asd on 25/06/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//  based in http://www.nongnu.org/hpalib/
//  with XDIM==7 -> long double precission

#ifndef hpa_metal_h
#define hpa_metal_h


#define XLITTLE_ENDIAN 1
#define XULONG_BITSIZE 64
#define UPPER_BOUND (XMAX_10EX + 100)
#define NULL 0

constant int XDIM = 7;
typedef struct _xpr { // hpa binary representation
    unsigned short nmm[XDIM + 1];
} xpr;

// constants for XDIM==7
constant xpr xZero = {{0x0, 0x0}}, xOne = {{0x3fff, 0x8000}},xTwo = {{0x4000, 0x8000}},xTen = {{0x4002, 0xa000}},xPinf = {{0x7fff, 0x0}},xMinf = {{0xffff, 0x0}},xVSV = {{0x3ff2, 0x8000}},xVGV = {{0x4013, 0x8000}},xEmax = {{0x400c, 0xb16c}},xEmin = {{0xc00c, 0xb16c}},xE2min = {{0xc00c, 0xfffb}}, /* -16382.75 */ xE2max = {{0x400c, 0xfffb}}; /* +16382.75 */
constant xpr xPi4 = {{0x3ffe, 0xc90f, 0xdaa2, 0x2168, 0xc234, 0xc4c6, 0x628b, 0x80dc}},
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

constant int XERR_DFL = 1,XENONE = 0,XEDIV = 1,XEDOM = 2,XEBADEXP = 3,XFPOFLOW = 4; // Floating point overflow
constant const int XNERR = 4,XEINV = 5, XMAX_10EX = 4931;
constant unsigned short xM_sgn = 0x8000, xM_exp = 0x7fff;
constant short xBias = 16383,xD_bias = 15360 ,xD_max = 2047,xD_lex = 12,xF_bias = 16256,xF_max = 255,xF_lex = 9;
constant short xMax_p = 16 * XDIM,xK_lin = -8 * XDIM;
constant int xItt_div = 2, xK_tanh = 5, xMS_exp = 21, xMS_hyp = 25,xMS_trg = 31;

constant thread char _sInf[4]="Inf", _smInf[5]="-Inf", _sNaN[4]="NaN";

class HPA {
public:
    xpr br; // binary representation
    
    // constructors
    HPA() : br(xZero) {}
    HPA(xpr x) : br(x) {}
    HPA(float x) { br = flttox(x); }
    HPA(int n) { br = inttox(n); }
    HPA(uint n) { br = inttox(n); }
    HPA(thread const char *str, thread char **endptr = 0) { br = strtox(str, endptr); }
    HPA(const thread HPA &x) { br = x.br; }
    HPA(const device HPA &x) { br = x.br; }
    
    // assign
    HPA operator=(const float d) {
        br = flttox(d);
        return *this;
    }
    HPA operator=(const int d) {
        br = inttox(d);
        return *this;
    }
    HPA operator=(const uint d) {
        br = inttox(d);
        return *this;
    }
    HPA operator=(thread char *d) {
        br = strtox(d, 0);
        return *this;
    }
    
    HPA operator=(thread HPA &x) {   br = x.br;   return *this;    }
    HPA operator=(device HPA &x) {   br = x.br;   return *this;    }
    HPA operator=(HPA x) {   br = x.br;   return *this;    }
    
    void setDevice(device HPA &x) { br=x.br; }
    void xxx() {}
    
    // aritmetic
    HPA operator+(const HPA x) { return HPA(xsum(br, x.br)); }
    HPA operator-(const HPA x) { return HPA(xsub(br, x.br)); }
    HPA operator*(const HPA x) { return HPA(xmul(br, x.br)); }
    HPA operator/(const HPA x) { return HPA(xdiv(br, x.br)); }
    HPA operator-() { return xneg(br); }
    
    HPA thread &operator+=(const HPA x) {
        br = xsum(br, x.br);
        return *this;
    }
    HPA thread &operator-=(const HPA x) {
        br = xsub(br, x.br);
        return *this;
    }
    HPA thread &operator*=(const HPA x) {
        br = xmul(br, x.br);
        return *this;
    }
    HPA thread &operator/=(const HPA x) {
        br = xdiv(br, x.br);
        return *this;
    }
    
    // compare
    bool operator==(const HPA x) { return xeq(br, x.br); }
    bool operator!=(const HPA x) { return xneq(br, x.br); }
    bool operator>(const HPA x) { return xgt(br, x.br); }
    bool operator<(const HPA x) { return xlt(br, x.br); }
    bool operator>=(const HPA x) { return xge(br, x.br); }
    bool operator<=(const HPA x) { return xle(br, x.br); }
    
    // tras. func's
    HPA exp() { return HPA(xexp(br)); }
    HPA exp10() { return HPA(xexp10(br)); }
    HPA exp2() { return HPA(xexp2(br)); }
    HPA log() { return HPA(xlog(br)); }
    HPA log10() { return HPA(xlog10(br)); }
    HPA log2() { return HPA(xlog2(br)); }
    
    HPA pow(thread HPA &y) { return HPA(xpow(br, y.br)); }
    HPA pwr(int n) { return HPA(xpwr(br, n)); }
    
    HPA fmod(thread HPA &y) {
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
    bool isNaN() { return xisNaN(&br); }
    bool isInf() { return xisPinf(&br); }
    
    // misc
    HPA abs() { return xabs(br); }
    
    //private:
    
    xpr xneg(xpr s) {
        thread unsigned short *p = (thread unsigned short *)&s;
        
        *p ^= xM_sgn;
        return s;
    }
    xpr xsum(xpr a, xpr b) { return xadd(a, b, 0); }
    xpr xsub(xpr a, xpr b) { return xadd(a, b, 1); }
    
    xpr xadd(xpr s, xpr t, int f) {
        unsigned short pe[XDIM + 1], h, u;
        thread unsigned short *pa, *pb, *pc, *pf = pe;
        unsigned int n = 0;
        short e, k;
        
        pa = (thread unsigned short *)&s;
        pb = (thread unsigned short *)&t;
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
        return *(thread xpr *)pe;
    }
    
    xpr xmul(xpr s, xpr t) {
        thread unsigned short pe[XDIM + 2], *q0, *q1, h;
        thread unsigned short *pa, *pb, *pc;
        unsigned int m, n, p;
        short e;
        short k;
        
        q0 = (thread unsigned short *)&s;
        q1 = (thread unsigned short *)&t;
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
        return *(thread xpr *)pe;
    }
    
    xpr xdiv(xpr s, xpr t) {
        xpr a;
        thread unsigned short *pc, e, i;
        auto zero=xZero;
        
        pc = (thread unsigned short *)&t;
        e = *pc;
        *pc = xBias;
        if ((xsigerr(xprcmp(&t, &zero) == 0, XEDIV, "xdiv()")))
            return xZero;
        else {
            a = flttox(1 / xtoflt(t));
            *pc = e;
            pc = (thread unsigned short *)&a;
            *pc += xBias - (e & xM_exp);
            *pc |= e & xM_sgn;
            for (i = 0; i < xItt_div; ++i)
                a = xmul(a, xadd(xTwo, xmul(a, t), 1));
            return xmul(s, a);
        }
    }
    
    float xtoflt(xpr s) {
        thread unsigned short pe[2], *pc, u;
        short i, e;
        
        pc = (thread unsigned short *)&s;
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
        return *(thread float *)pe;
    }
    
    xpr flttox(float y) {
        unsigned short pe[XDIM + 1], u;
        unsigned short thread *pc;
        short i, e;
        
        if (y < FLT_MIN && y > -FLT_MIN)
            return xZero;
        pc = (unsigned short thread *)&y;
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
        return *(xpr thread*)pe;
    }
    
    xpr inttox(int n) {
        thread unsigned short pe[XDIM + 1],  *pc;
        short e;
        unsigned int k, h;
        
        k = _abs(n);
        pc = (thread unsigned short *)&k;
        for (e = 0; e <= XDIM; pe[e++] = 0)
            ;
        if (n == 0)
            return *(thread xpr *)pe;
        
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
        return *(thread xpr *)pe;
    }
    
    void xlshift(int n, thread unsigned short *pm, int m) {
        thread unsigned short *pa, *pc;
        
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
    
    void xrshift(int n, thread unsigned short *pm, int m) {
        thread unsigned short *pa, *pc;
        
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
    int xprcmp(const thread xpr *pa, const thread xpr *pb) {
        thread unsigned short e, k, *p, *q, p0, q0;
        int m;
        
        p = (thread unsigned short *)pa;
        q = (thread unsigned short *)pb;
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
    
    int xisNaN(thread const xpr *u) {
        thread unsigned short *p = (thread unsigned short *)u;
        
        if ((*p))
            return 0;
        else {
            int i;
            
            for (i = 1; i <= XDIM && p[i] == 0x0; i++)
                ;
            return (i <= XDIM ? 1 : 0);
        }
    }
    
    int xis0(thread const xpr *u) {
        thread unsigned short *p = (thread unsigned short *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        return (m > XDIM && (*p & xM_exp) < xM_exp ? 1 : 0);
    }
    
    int xnot0(thread const xpr *u) {
        thread unsigned short *p = (thread unsigned short *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        return (m > XDIM && (*p & xM_exp) < xM_exp ? 0 : 1);
    }
    
    int xsgn(thread xpr *u) {
        thread unsigned short *p = (thread unsigned short *)u;
        int m;
        
        for (m = 1; m <= XDIM && p[m] == 0; m++)
            ;
        if ((m > XDIM && (*p & xM_exp) < xM_exp) || !*p)
            return 0;
        else
            return ((*p & xM_sgn) ? -1 : 1);
    }
    
    int xisPinf(thread const xpr *u) { return (*u->nmm == xM_exp ? 1 : 0); }
    
    int xisMinf(thread const xpr *u) { return (*u->nmm == (xM_exp | xM_sgn) ? 1 : 0); }
    
    int xisordnumb(thread const xpr *u) {
        int isNaN, isfinite;
        thread unsigned short *p = (thread unsigned short *)u;
        
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
    
    unsigned int _abs(int n) { return ((n) >= 0 ? (n) : -(n)); }
    
    bool isspace(char ch) { return ch<=' '; }
    
    xpr strtox(thread const char *q, thread char **endptr) {
        xpr s, f;
        unsigned short pc[XDIM + 1];
        thread unsigned short *pn, *pf, *pa, *pb;
        unsigned short sfg, ibex, fbex;
        unsigned int n;
        int j;
        int idex = 0, fdex = 0, c, m;
        short noip, nofp;
        thread const char *ptr;
        
        pn = (thread unsigned short *)&s;
        pf = (thread unsigned short *)&f;
        for (j = 0; j <= XDIM; pn[j] = pf[j] = pc[j] = 0, ++j)
            ;
        sfg = 0;
        m = XDIM + 1;
        if ((endptr))
            *endptr = (thread char *)q;
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
                for (n = (unsigned int)c, pa = pn + XDIM, pb = pc + XDIM; pa >= pn;
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
                for (n = (unsigned int)c, pa = pf + XDIM, pb = pc + XDIM; pa >= pf;
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
                *endptr = (thread char *)q;
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
                *endptr = (thread char *)q;
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
    
    xpr xpwr(xpr s, int n) {
        xpr t;
        unsigned k, m;
        auto zero = xZero;
        
        t = xOne;
        if (n < 0) {
            m = -n;
            
            if ((xsigerr(xprcmp(&s, &zero) == 0, XEBADEXP, "xpwr()")))
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
            xsigerr(xprcmp(&s, &zero) == 0, XEBADEXP, "xpwr()");
        return t;
    }
    
    
    int xsigerr(int errcond, int errcode, const constant char *where) {
        if (!errcond)  errcode = 0;
        if (errcode < 0 || errcode > XNERR) errcode = XEINV;
        return errcode;
    }
    
    
    // trigs
    int xodd(xpr x) {
        thread unsigned short *p = (thread unsigned short *)&x;
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
        xpr _xPi2=xPi2, _xPi4=xPi4;
        
        z = rred(z, 't', &k);
        if ((xsigerr(xprcmp(&z, &_xPi2) >= 0, XEDOM, "xtan()")))
            return (!k ? xPinf : xMinf);
        else {
            if (xprcmp(&z, &_xPi4) == 1) {
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
        if (k)
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
        if (k)
            return xneg(z);
        else
            return z;
    }
    
    xpr c_tan(xpr z) {
        xpr s, f, d;
        int m;
        unsigned short k;
        
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
    
    xpr rred(xpr z, int kf, thread int *ps) {
        xpr is, q;
        auto _xPi2=xPi2;
        
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
        if (xprcmp(&z, &_xPi2) == 1) {
            z = xadd(xPi, z, 1);
            if (kf == 'c' || kf == 't')
                q = xadd(q, xOne, 0);
        }
        *ps = (xodd(q)) ? 1 : 0;
        return z;
    }
    
    xpr xpr2(xpr s, int m) {
        thread unsigned short *p = (thread unsigned short *)&s;
        int e;
        
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
    
    xpr xabs(xpr s) {
        thread unsigned short *p = (thread unsigned short *)&s;
        
        *p &= xM_exp;
        return s;
    }
    
    int x_exp(thread const xpr *ps) {
        thread unsigned short *q = (thread unsigned short *)ps;
        
        return (*q & xM_exp) - xBias;
    }
    
    int x_neg(thread const xpr *ps) {
        thread  unsigned short *q = (thread unsigned short *)ps;
        
        return (*q & xM_sgn);
    }
    
    xpr xlog(xpr z) {
        xpr f, h;
        int k, m;
        auto _xOne=xOne;
        
        if ((xsigerr((x_neg(&z)) || x_exp(&z) == -xBias, XEDOM, "xlog()")))
            return xMinf;
        else if (xprcmp(&z, &_xOne) == 0)
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
            return xadd(f, xmul(xLn2, flttox(m - .5)), 0);
        }
    }
    
    xpr xlog2(xpr z) {
        xpr f, h;
        int k, m;
        auto _xOne=xOne;
        
        if ((xsigerr((x_neg(&z)) || x_exp(&z) == -xBias, XEDOM, "xlog2()")))
            return xMinf;
        else if (xprcmp(&z, &_xOne) == 0)
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
            return xadd(xmul(f, xLog2_e), flttox(m - .5), 0);
        }
    }
    
    xpr xlog10(xpr z) {
        auto _xMinf=xMinf;
        xpr w = xlog(z);
        
        if (xprcmp(&w, &_xMinf) <= 0)
            return xMinf;
        else
            return xmul(w, xLog10_e);
    }
    
    xpr xsfmod(xpr s, thread int *p) {
        thread unsigned short *pa, *pb;
        short e, k;
        
        pa = (thread unsigned short *)&s;
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
        thread unsigned short *pf = (thread unsigned short *)&x;
        int m, k;
        auto _xE2min=xE2min, _xE2max=xE2max;
        
        if (xprcmp(&x, &_xE2min) < 0)
            return xZero;
        else if ((xsigerr(xprcmp(&x, &_xE2max) > 0, XFPOFLOW, NULL)))
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
            pf = (thread unsigned short *)&f;
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
    
    xpr xfmod(xpr s, xpr t, thread xpr *q) {
        auto _xZero=xZero;
        if ((xsigerr(xprcmp(&t, &_xZero) == 0, XEDIV, "xfmod()")))
            return xZero;
        else {
            thread unsigned short *p, mask = 0xffff;
            short e, i;
            int u;
            
            *q = xdiv(s, t);
            p = (thread unsigned short *)q;
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
    
    xpr xfrexp(xpr s, thread int *p) {
        thread unsigned short *ps = (thread unsigned short *)&s, u;
        
        *p = (*ps & xM_exp) - xBias + 1;
        u = *ps & xM_sgn;
        *ps = xBias - 1;
        *ps |= u;
        return s;
    }
    
    /*
     Remark: 's' must have a size >= 5 !
     */
    void strcpy(thread char*dst, constant thread char*org) {  while (*dst) *(dst++)=*(org++);   }
    
    int copied_special_value(thread char *s, xpr u, int sign) {
        if ((xisPinf(&u))) {
            if ((sign)) {
                *s = '+';
                s++;
            }
            strcpy(s, _sInf);
            return 1;
        } else if ((xisMinf(&u))) {
            strcpy(s, _smInf);
            return 1;
        } else if ((xisNaN(&u))) {
            if ((sign)) {
                *s = '\?';
                s++;
            }
            strcpy(s, _sNaN);
            return 1;
        } else
            return 0;
    }
    
};

static void assign(device HPA&hd, HPA h) {  hd.br=h.br; }
static void assign(device HPA&hd, float f) {  hd.br=HPA(f).br; }
static void assign(device HPA&hd, int i) {  hd.br=HPA(i).br; }

HPA sin(HPA x) { return HPA(x).sin(); }
HPA cos(HPA x) { return HPA(x).cos(); }
HPA tan(HPA x) { return HPA(x).tan(); }
HPA exp(HPA x) { return HPA(x).exp(); }
HPA exp10(HPA x) { return HPA(x).exp10(); }
HPA exp2(HPA x)  { return HPA(x).exp2(); }
HPA log(HPA x)   { return HPA(x).log(); }
HPA log10(HPA x) { return HPA(x).log10(); }
HPA log2(HPA x)  { return HPA(x).log2(); }
HPA pow(HPA x, HPA y) { return HPA(x).pow(y); }
HPA pwr(HPA x, int n) { return HPA(x).pwr(n); }
HPA fmod(HPA x, HPA y) { return HPA(x).fmod(y); }
bool isNaN(HPA x) { return HPA(x).isNaN(); }
bool isInf(HPA x) { return HPA(x).isInf();}
HPA abs(HPA x) { return HPA(x).abs();}

#endif /* hpa_metal_h */
