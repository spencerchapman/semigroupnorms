import itertools
import random

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
