import itertools
import random
import matplotlib.pyplot as plt
import numpy as np

def T_Length(v,t):
    
    # Computes the t-length factorization of
    # v := (m_1,m_2,...,m_k), where m_i are exponents of atoms a_i
    
    if(t == oo):
        return max(v)
    if(t == 0):
        return sum((m != 0) for m in v)
    s = sum(m ** t for m in v)
    if(t <= 1):
        return s
    return s ** (t ** -1)
def T_Length_Integer(v,t: int):
    return sum([m**t for m in v])
def T_LengthSetIntegerT(S,x,t):
    return sorted([T_Length_Integer(v,t) for v in S.Factorizations(x)])

def T_LengthSet(S,x,t):
    return sorted(list(Set([T_Length(f,t) for f in S.Factorizations(x)])))
def T_ElementsUpToN(S,n):
    return list(i for i in [1..n] if S.Contains(i))
def T_MaxLength(S,x,t):
    return T_LengthSet(S,x,t)[-1]
def T_MinLength(S,x,t):
    return T_LengthSet(S,x,t)[0]
def T_Elasticity(S,x,t):
    L = T_LengthSet(S,x,t)
    return L[-1]/L[0]
def T_FirstNTerms(S,n):
    n_terms = []
    i = 1
    while(len(n_terms)!=n):
        if(S.Contains(i)):
            n_terms.append(i)
        i+=1
    return n_terms
def FirstTermsUpToN(S,n):
    return list(i for i in [0..n] if S.Contains(i))
def T_MaxPowersOfElement(S,x,t,n):
    F = S 
    F.FactorizationsUpToElement(n*x)
    PowerMap = [T_MaxLength(F,i*x,t) if F.Contains(i*x) else 0 for i in [0..n]]
    del F
    return PowerMap
def T_MinPowersOfElement(S,x,t,n):
    F = S 
    F.FactorizationsUpToElement(n*x)
    PowerMap = [T_MinLength(F,i*x,t) if F.Contains(i*x) else 0 for i in [0..n]]
    del F
    return PowerMap
def T_FullMaxPowers(S,m,t,n):
    SemigroupTerms = T_FirstNTerms(S,m)
    F = S
    F.FactorizationsUpToElement(n*SemigroupTerms[-1])
    FullPowers = [[T_MaxLength(F,i*x,t) for i in [0..n]] for x in SemigroupTerms]
    del F
    return FullPowers
def T_FullMinPowers(S,m,t,n):
    SemigroupTerms = T_FirstNTerms(S,m)
    F = S
    F.FactorizationsUpToElement(n*SemigroupTerms[-1])
    FullPowers = [[T_MinLength(F,i*x,t) for i in [0..n]] for x in SemigroupTerms]
    del F
    return FullPowers
def T_PandasFullMaxPowers(S,m,t,n):
    df = pd.DataFrame(T_FullMaxPowers(S,m,t,n))
    df = df[df.columns[1:]].T
    df.columns = ['L_t('+str(i)+'n)' for i in T_FirstNTerms(S,m)]
    df.index.name = 'n'
    return df
def T_PandasFullMinPowers(S,m,t,n):
    df = pd.DataFrame(T_FullMinPowers(S,m,t,n))
    df = df[df.columns[1:]].T
    df.columns = ['l_'+str(t)+'('+str(i)+'n)' for i in T_FirstNTerms(S,m)]
    df.index.name = 'n'
    return df
def T_CSVFullMaxPowers(S,m,t,n):
    df = T_PandasFullMaxPowers(S,m,t,n)
    try:
        os.makedirs('csv_exports',exist_ok=True)
        s = 'csv_exports/'+'_'.join(str(x) for x in S.gens)+'T_'+str(t)+'.csv'
        df.to_csv(s)
        print('exported to '+s)
    except:
        print('error')
def T_CoordinateMaxPowers(S,x,t,n):
    return list(zip([0..n],T_MaxPowersOfElement(S,x,t,n)))
def T_CoordinateMinPowers(S,x,t,n):
    return list(zip([0..n],T_MinPowersOfElement(S,x,t,n)))
def t0(v):
    return sum((a != 0) for a in v)
# Thank you Eric for these! VVV
def SupportUpToX(S,x):
    G = S.gens
    Basis = [2^i for i in range(len(G))]
    SupportDict = {0:set([0])}
    for m in range(1,x+1):
        mfact = set()
        for i in range(len(G)):
            k = m-G[i]
            if k>=0:
                mfact = mfact.union(set([Basis[i] | f for f in SupportDict[k]]))
            #print(m,mfact)
        SupportDict[m] = mfact
    return SupportDict
def LZUpToX(S,x):
    LDict = []
    U = SupportUpToX(S,x)
    for i in range(x+1):
        LDict.append( set([t0([1 if ch == "1" else 0 for ch in bin(s)[2:]]) for s in U[i]]) )
    return LDict
def T_ZeroPowerLengths(S,x,n):
    FullLengths = LZUpToX(S,x*n)
    PowerLengths = [FullLengths[i*x] for i in range(1,n+1)]
    return PowerLengths


####################################################################################################

import itertools
import random

def norm(t,vec): #from James
    v = list(vec)
    if t!= 0:
        powers = [a^t for a in v]
    else:
        powers = [0 if a==0 else 1 for a in v]
    if t >= 1:
        return((sum(powers))^(1/t))
    else:
        return(sum(powers))

def D_delta(vec):
    svec = sorted(vec)
    return set([svec[i+1]-svec[i] for i in range(len(svec)-1)])
    #for i in range(len(svec)-1):
     #   dset.add(svec[i+1]-svec[i])
    #return dset

def D_DeltaCheck(S,r=1000,t=0):
    S.FactorizationsUpToElement(r)
    F = [S.Factorizations(i) for i in range(r+1)]
    L = [e if e==[] else set([norm(t,vec) for vec in e]) for e in F]
    D = [D_delta(vec) for vec in L]
    SDelta = set().union(*D)
    pts = []
    for i in range(r+1):
        for d in D[i]:
            pts.append([i,d])
    return (SDelta,pts,D,L)

def D_ImprovedDeltaCheck(S,r=1000):
    L = LZUpToX(S,r)
    D = [D_delta(vec) for vec in L]
    SDelta = set().union(*D)
    pts = []
    for i in range(r+1):
        for d in D[i]:
            pts.append([i,d])
    return (SDelta,pts,D,L)

def D_Nablas(S):
    B = S.BettiElements()
    X = SupportUpToX(S, max(B))
    N = {}
    for b in B:
        G = Graph()
        G.add_vertices(X[b])
        for e in itertools.combinations(X[b],2):
            if (e[0] & e[1] > 0):
                G.add_edge(e[0],e[1])
        N[b] = G
    return B,N

def D_NBI(bitwiseInteger):
    if bitwiseInteger == 0:
        return 0
    else:
        return norm(0,[1 if ch == "1" else 0 for ch in bin(bitwiseInteger)[2:]])
def D_M(S):
    N = D_Nablas(S)
    B = N[0]
    N = N[1]
    return max(max([min([NBI(v) for v in C]) for C in N[b].connected_components()]) for b in B)
# Start A-Team Section
# --------------------
def tLengthSet(S: NumericalSemigroup, n: int, t: float) -> list[float]:
    """Return a list of factorization t-lengths of n in S.
    
    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    t : float
        t-length
    
    Returns
    -------
    list[float]
        The list of factorization t-lengths of n in S.
    """
    if t:
        if t > 1:
            return list(sum(num**t for num in factorization)**(1/t) for factorization in S.Factorizations(n))
        return list(sum(num**t for num in factorization) for factorization in S.Factorizations(n))
    return list(len(factorization) - factorization.count(0) for factorization in S.Factorizations(n))

def tMaxLength(S: NumericalSemigroup, n: int, t: float) -> float:
    """Return the maximum t-length factorization of n in S.

    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    t : float
        t-length
    
    Returns
    -------
    float
        The maximum t-length factorization of n in S.
    """
    return max(tLengthSet(S,n,t))

def tMinLength(S: NumericalSemigroup, n: int, t: float) -> float:
    """Return the minimum t-length factorization of n in S.
    
    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    t : float
        t-length
    
    Returns
    -------
    float
        The minimum t-length factorization of n in S.
    """
    return min(tLengthSet(S,n,t))

def tElasticity(S: NumericalSemigroup, n: int, t: float) -> float:
    """Return the t-elasticity of n in S.

    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    t : float
        t-length
    
    Returns
    -------
    float
        The t-elasticity of n in S.
    """
    lengths = tLengthSet(S,n,t)
    return max(lengths)/min(lengths)
    
def tFactorizationLengths(factorization: list[int], t_vals: list[float]) -> list[float]:
    """Return a list of the t-length of factorization for t values in t_vals.
    
    Parameters
    ----------
    factorization : list[int]
        A factorization
    t_vals : list[float]
        A list of t-values

    Returns
    -------
    list[float]
        A list of factorization t-lengths for each t-value in the same order as t_vals.
    """
    lengths = []
    for t in t_vals:
        if t:
            if t > 1:
                lengths.append(sum(num**t for num in factorization)**(1/t))
                continue
            lengths.append(sum(num**t for num in factorization))
            continue
        lengths.append(len(factorization) - factorization.count(0))
    return lengths

def tLengthPlot(S: NumericalSemigroup, n: int, start: float=0, stop: float=3, samples: int=200, vrange: list[float,float]=None) -> None:
    """Plot the factorization t-lengths of n as t varies from start to stop.
    
    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    start : float
        The minimum t-value
    stop : float
        The maximum t-value
    samples : int
        The number of t-values between start and stop (inclusive) sampled.
    vrange : list[float,float]
        An interval of y-values to show on the plot
    
    Returns
    -------
    None
    """
    t_vals = np.linspace(start, stop, num=samples)
    fig, ax = plt.subplots()
    for factorization in S.Factorizations(n):
        ax.plot(t_vals, tFactorizationLengths(factorization,t_vals), label=str(factorization))
    if vrange is not None:
        ax.set_ylim(vrange[0], vrange[1])
    ax.legend()
    ax.set_title(f'Factorization $t$-length for {n} in {S.gens}')
    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$t$-length')
    plt.show()

def tElasticityPlot(S: NumericalSemigroup, n: int, start: float=0, stop: float=3, sample: int=200, vrange: list[float,float]=None) -> None:
    """Plot the t-elasticity of n as t varies from start to stop.
    
    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    start : float
        The minimum t-value
    stop : float
        The maximum t-value
    samples : int
        The number of t-values between start and stop (inclusive) sampled.
    vrange : list[float,float]
        An interval of y-values to show on the plot

    Returns
    -------
        None
    """
    t_vals = np.linspace(start, stop, num=sample)
    t_elasticities = []
    for t in t_vals:
        lengths = tLengthSet(S, n, t)
        t_elasticities.append(max(lengths)/min(lengths))
    fig, ax = plt.subplots()
    ax.plot(t_vals, t_elasticities)
    if vrange is not None:
        ax.set_ylim(vrange[0], vrange[1])
    ax.set_title(f'$t$-Elastisity of {n} in {S.gens}')
    ax.set_xlabel('$t$')
    ax.set_ylabel(f'$\\rho_t({n})$')
    plt.show()

def tFactorizationRatioPlot(factorization_1: list[int], factorization_2: list[int], start: float=0, stop: float=3, samples: int=200, vrange: list[float,float]=None) -> None:
    """Plot the ratios of factorization_1 and factorization_2 t-lengths as t varies from start to stop.
    
    Parameters
    ----------
    factorization_1 : list[int]
    factorization_2 : list[int]
    start : float
        The minimum t-value
    stop : float
        The maximum t-value
    samples : int
        The number of t-values between start and stop (inclusive) sampled.
    vrange : list[float,float]
        An interval of y-values to show on the plot
    
    Returns
    -------
        None
    """
    t_vals = np.linspace(start,stop, num=samples)
    lengths_1 = tFactorizationLengths(factorization_1, t_vals)
    lengths_2 = tFactorizationLengths(factorization_2, t_vals)
    ratio_1 = [pair[0]/pair[1] for pair in zip(lengths_1,lengths_2)]
    ratio_2 = [pair[1]/pair[0] for pair in zip(lengths_1,lengths_2)]
    fig, ax = plt.subplots()
    ax.plot(t_vals, ratio_1, label=f'$t$-length of {factorization_1} / $t$-length of {factorization_2}')
    ax.plot(t_vals, ratio_2, label=f'$t$-length of {factorization_2} / $t$-length of {factorization_1}')
    if vrange is not None:
        ax.set_ylim(vrange[0], vrange[1])
    ax.set_title(f'Ratio of $t$-length of {factorization_1} and {factorization_2}')
    ax.set_xlabel('$t$')
    ax.set_ylabel('Ratio')
    ax.legend()
    plt.show()

def tAllFactorizationsRatioPlot(S: NumericalSemigroup, n: int, start: float=0, stop: float=3, samples: int=200, vrange: list[float,float]=None) -> None:
    """Plot the maximum ratio of all pairs of t-length factorizations of n as t varies from start to stop.
    
    Parameters
    ----------
    S : NumericalSemigroup
    n : int
        An element of S
    start : float
        The minimum t-value
    stop : float
        The maximum t-value
    samples : int
        The number of t-values between start and stop (inclusive) sampled.
    vrange : list[float,float]
        An interval of y-values to show on the plot
    
    Returns
    -------
        None
    """
    t_vals = np.linspace(start, stop, num=samples)
    fig, ax = plt.subplots()
    for factorization_1, factorization_2 in itertools.combinations(S.Factorizations(n), 2):
        lengths_1 = tFactorizationLengths(factorization_1, t_vals)
        lengths_2 = tFactorizationLengths(factorization_2, t_vals)
        ax.plot(t_vals, [max(pair)/min(pair) for pair in zip(lengths_1,lengths_2)], label=f'{factorization_1} and {factorization_2}')
    ax.plot(t_vals, [1 for i in range(samples)], label='factorization and itself')
    if vrange is not None:
        ax.set_ylim(vrange[0], vrange[1])
    ax.set_title(f'Maximum Factorization Ratios for {n} in {S.gens}')
    ax.set_xlabel('$t$')
    ax.set_ylabel('Maximum Factorization Ratio')
    ax.legend(bbox_to_anchor=(1,1))
    plt.show()
# End A-Team Section
# ------------------
