########################################################################
# Author & Copyright : Peter Luschny
# Created: 2014-12-12
# License: LGPL version 2.1 or (at your option)
# Creative Commons Attribution-ShareAlike 3.0
########################################################################

def OEIS(T, num, info=false):  # needs internet

    def is_dead(id):
        for kw in oeis(id).keywords():
            if kw == 'dead':
                return true
        return false

    if len(T) < 3: return []

    # disregard the first entry and trailing zeros
    i, l = 1, len(T)
    while i < l and T[i] == 0: i += 1
    S = T[i:min(l, 20)] if (min(l, 20) - i) > 6 else T

    R = oeis(S, max_results=num)

    if len(R) == 0:
        # what if S is too long?
        Sx = S[:]
        while len(", ".join([str(a) for a in Sx])) > 270:
            del Sx[-1]
        else:
            del Sx[-1]
        R = oeis(Sx, max_results=num)

    if len(R) == 0:
        U = [z for z in S if z != 0]
        del U[0]
        R = oeis(U, max_results=num)

        if len(R) == 0:
            # what if U is too long?
            Ux = U[:]
            while len(", ".join([str(a) for a in Ux])) > 270:
                del Ux[-1]
            else:
                del Ux[-1]
            R = oeis(Ux, max_results=num)

    L = []
    for r in R:
        s = str(r)[0:7]
        i = int(s[1:7])

        if is_dead(i):
            continue

        L.append(s)
        if info:
            K = oeis(i).first_terms(10)
            print(s, K)

    return L


########################################################################
def padded_dict(T, n):
    """

    List[Pair] --> List

    Assume T a list of pairs [v, m] with values v and strictly increasing
    in the second component with nonnegative integers m.

    Return L a list of values v in the order given, padded with 0 for
    those values of m which are missing in T, for 0 <= m < n.

    In use T is often a dictionary of coefficients.

    EXAMPLES::
        in:  [[a, 1], [b, 4]], 6
        out: [0, a, 0, 0, b, 0]

        in : [[-2, 1], [-6, 2]], 2
        out: [0, -2]
    """
    L = [0] * (n + 1)

    for t in T:
        if t[1] > n: break
        L[t[1]] = t[0]
    return L


###########################################################
def is_intseq(seq):
    for s in seq:
        try:
            Integer(s)
        except:
            return false
    return true


###########################################################
def plot_gfun(B, genfun):
    GF = ['ogf', 'egf', 'lgf', 'bgf']
    OP = [' ', 'rev', 'inv', 'exp', 'der', 'log', 'logder', 'int']
    rgb = (171 / 256, 141 / 256, 63 / 256)
    L = [[0.8 * cos(pi * i / 4), 0.8 * sin(pi * i / 4)] for i in range(8)]
    a = polygon2d(L, rgbcolor=rgb)
    b = circle((0, 0), 0.8)
    c = circle((0, 0), 2)
    m = 0
    for n in range(8):
        t = n * pi / 4
        z = (0.8 * cos(t), 0.8 * sin(t))
        u = (0.65 * cos(t), 0.65 * sin(t))
        c += text(OP[n], u)
        for k in range(4):
            t = (2 * m - pi) * pi / 32
            y = (2 * cos(t), 2 * sin(t))
            if B[m] == 1:
                c += line([z, y], thickness=3, color='green')
            elif B[m] == -1:
                c += line([z, y], thickness=3, color='red')
            else:
                c += line([z, y])
            v = (2.1 * cos(t), 2.1 * sin(t))
            c += text(GF[k], v)
            m += 1
    c += text(genfun, (0, 2.25), fontsize=11, color='black')
    c += text('GFun', (0, 0), fontsize=22, color='black')
    (b + c + a).show(figsize=8, axes=false)


########################################################################
########################################################################
class gfun(object):
    """
    SR --> PSR

    Convert a symbolic expression in 0 or more variables into a
    Power Series Ring in x over the Symbolic Ring by applying a
    formal Taylor expansion in 0 to g up to degree n.
    """

    ###############################################
    @staticmethod
    def reversion_ext(g):
        """
        Reversion extended to the non-revertable case.
        """
        r = g / g[0]
        r = r.shift(1)
        r = r.reverse()
        r = r.shift(-1)
        return r * g[0]

    ###############################################
    @staticmethod
    def transform(P, typ, strict):
        if typ == 'rev':
            if P[0] != 0 and not strict:
                return gfun.reversion_ext(P)
            if P.valuation() != 1:
                raise ValueError('Series must have valuation one for reversion.')
            return P.reverse()
        if typ == 'inv':
            if P.valuation() != 0:
                raise ValueError('Series must have valuation zero for inversion.')
            return P.__invert__()
        if typ == 'int':
            return P.integral()
        if typ == 'der':
            return P.derivative()
        if typ == 'logder':
            if P.valuation() != 0:
                raise ValueError('Series must have valuation zero for inversion.')
            return P.derivative() / P
        if typ == 'log':
            if P[0] != 1:
                raise ValueError('Series must have constant term one for logarithm.')
            return P.log()
        if typ == 'exp':
            return P.exp()

        return P

    ###############################################
    @staticmethod
    def taylor_series(gf, n, typ, strict=false):
        g = SR(gf)  # make sure that we get an symbolic expression
        R.<x> = PowerSeriesRing(SR, default_prec=n)
        if g.variables() == (): return R(g)
        T = taylor(g, g.variables()[0], 0, n).coefficients()
        P = R(padded_dict(T, n))
        return gfun.transform(P, typ, strict)

    ###############################################
    def __init__(self, g, op=None):
        self.gfun = g
        self.prec = 20
        self.typ = op
        self.strict_rev = false
        self.PS = gfun.taylor_series(g, self.prec, op)

    ###############################################
    def strict_reversion(self, onoff):
        self.strict_rev = onoff
        gfun.taylor_series(self.gfun, self.prec, self.typ, onoff)

    ###############################################
    def check_prec(self, n):
        if n > self.prec:
            self.PS = gfun.taylor_series(self.gfun, n, self.typ, self.strict_rev)
            self.prec = n

    ###############################################
    def dump(self):
        print('Type ord/rev/inv etc.: ' + str(self.typ))
        print('Strict reversion on: ' + str(self.strict_rev))
        print('Formal power series: ' + str(self.gfun))
        print('Coefficients: ' + str(self.PS.padded_list(6)))

    ###############################################
    @staticmethod
    def lgf_to_ogf(s):
        """
        Returns the ordinary generating function power series,
        assuming self is a logarithmic generating function power series.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t - 3/2*t^2 + 4/3*t^3
            sage: f.lgf_to_ogf()
            t + 3*t^2 + 4*t^3
        """
        return s.parent()([s[i] * (-1)^(i + 1) * i for i in range(s.degree() + 1)])

    ###############################################
    @staticmethod
    def bgf_to_ogf(s):
        """
        Returns the ordinary generating function power series,
        assuming self is a binary generating function power series.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(QQ)
            sage: f = t/2 - 3/4*t^2 + 7/8*t^3
            sage: f.lgf_to_ogf()
            t - 3*t^2 + 7*t^3
        """
        return s.parent()([s[i] * 2 ^ i for i in range(s.degree() + 1)])

    ###############################################
    def as_series(self, prec=12):
        self.check_prec(prec)
        return self.PS.add_bigoh(prec)

    ###############################################
    ## Symbolic Ring -> List
    def as_seq(self, n, search, typ='ogf', coeffs=false):
        self.check_prec(n)

        if typ == 'egf':
            P = self.PS.egf_to_ogf()
        elif typ == 'lgf':
            P = gfun.lgf_to_ogf(self.PS)
        elif typ == 'bgf':
            P = gfun.bgf_to_ogf(self.PS)
        else:
            P = self.PS

        if coeffs:
            R.<x> = PolynomialRing(QQ)
            # padded_list() does not work here!
            return [[R(P[n])[i] for i in (0..n)] for n in (0..n)]

        L = P.padded_list(n)
        if search: print(OEIS(L, 4, info=true))
        return L

    ###############################################
    ## Symbolic Ring -> List
    def as_ogf(self, n=12, search=false):
        return self.as_seq(n, search, 'ogf')

    ## Symbolic Ring -> List
    def as_egf(self, n=12, search=false):
        return self.as_seq(n, search, 'egf')

    ## Symbolic Ring -> List
    def as_lgf(self, n=12, search=false):
        return self.as_seq(n, search, 'lgf')

    ## Symbolic Ring -> List
    def as_bgf(self, n=12, search=false):
        return self.as_seq(n, search, 'bgf')

    ## Symbolic Ring -> Triangular List of Lists
    def coeffs(self, typ='ogf', rows=8):
        return self.as_seq(rows, false, typ, true)

    ## Symbolic Ring -> Rectangular List of Lists
    def values(self, typ='ogf', rows=8, cols=8):
        L = self.as_seq(rows, false, typ)
        return [[L[n](x=k) for n in range(rows)] for k in range(cols)]

    ###############################################
    @staticmethod
    def explore_genfun(fps):
        print('Formal power series: ', fps)
        print()

        B = [0 for _ in range(32)]
        c = 0

        for op in [None, 'rev', 'inv', 'exp', 'der', 'log', 'logder', 'int']:
            for funcname in ['as_ogf', 'as_egf', 'as_lgf', 'as_bgf']:

                try:
                    G = gfun(fps, op)
                    func = getattr(G, funcname, None)
                except:
                    pass
                else:
                    L = func()
                    if is_intseq(L):
                        B[c] = true
                        print(op, funcname, L[:9])
                        O = OEIS(L, 4, info=true)
                        if O != []:
                            print(O)
                            B[c] = 1
                        else:
                            B[c] = -1
                        print()
                c += 1

        plot_gfun(B, fps)
        print('Done!')

    ###############################################
    @staticmethod
    def explore_values(gen, typ='ogf', rows=6, cols=6):
        L = gen.values(typ, 9, cols)
        for n in range(rows):
            print("row", [n], L[n])
            print(OEIS(L[n], 4, info=true))
            print()
        for n in range(cols):
            gf = gen.colgenerators(n, typ)
            M = gfun(gf).as_ogf(9)
            print("col", [n], M)
            print(OEIS(M, 4, info=true))
            print()
        print('Done!')

    ###############################################
    def colgenerators(self, n, typ='ogf'):
        from sage.calculus.calculus import symbolic_sum

        x, j = SR.var('x, j')
        assume(abs(x) < 1)
        P = self.as_seq(n + 1, false, typ)
        p = symbolic_sum(x^j * P[n].substitute(x=j), j, 0, oo)
        return p.partial_fraction()

    ###############################################
    # Thanks to Ralf Stephan for help with the implementation.
    def succ(self, typ='ogf', rows=8):
        from sage.calculus.calculus import symbolic_sum
        x, j = SR.var('x, j')
        assume(abs(x) < 1)
        L = []
        P = self.as_seq(rows, false, typ)
        for k, m in enumerate(P):
            p = symbolic_sum(x^j * m.substitute(x=j), j, 0, oo)
            q = p.numerator().polynomial(ZZ) / p.denominator().polynomial(ZZ)
            f = q.partial_fraction_decomposition()
            C = [0 for _ in range(k + 1)]
            for h in f[1]:
                i = h.denominator().degree() - 1
                C[i] = Integer((-1)^(k + 1) * h.numerator() / factorial(i))
            L.append(C)
        return L

    ###############################################
    def pred(self, typ = 'ogf', rows=8):
        x = SR.var('x')
        PR = PolynomialRing(QQ, 'x')
        P = self.coeffs(typ, rows)
        A = []
        for n, p in enumerate(P):
            W = sum(factorial(j) * p[j] / (x - 1)^(j + 1) for j in (0..n))
            T = SR(W).taylor(x, 0, n + 2).coefficients()
            C = padded_dict(T, n)
            xy = [(k, C[k]) for k in (0..n)]
            L = PR(PR.lagrange_polynomial(xy))
            A.append([(-1)^(n + 1) * L[i] for i in (0..n)]) 
        return A

###############################################
### Transformations
###############################################
def Succ(T, rows):
    from sage.calculus.calculus import symbolic_sum
    x, v = SR.var('x, v')
    PR = PolynomialRing(QQ, 'x')
    assume(abs(x) < 1)
    L = []
    for k in range(rows):
        h = PR(sum(T(k,j) * x^j for j in (0..k)))
        p = symbolic_sum(x^v * h.substitute(x=v), v, 0, oo)
        q = p.numerator().polynomial(ZZ) / p.denominator().polynomial(ZZ)
        f = q.partial_fraction_decomposition()
        C = [0 for _ in range(k + 1)]
        for h in f[1]:
            i = h.denominator().degree() - 1
            C[i] = Integer((-1)^(k + 1) * h.numerator() / factorial(i))
        L.append(C)
    return L

###############################################
def Pred(T, rows):
    PR = PolynomialRing(QQ, 'x')
    A = []
    for n in range(rows):
        W = sum(factorial(j) * T(n,j) / (x - 1)^(j + 1) for j in (0..n))
        S = SR(W).taylor(x, 0, n + 2).coefficients()
        C = [0] * (n + 1)
        for s in S:  # padding
            if s[1] > n: break
            C[s[1]] = s[0]
        xy = [(k, C[k]) for k in (0..n)]
        L = PR(PR.lagrange_polynomial(xy))
        A.append([(-1)^(n + 1) * L[i] for i in (0..n)]) 
    return A
