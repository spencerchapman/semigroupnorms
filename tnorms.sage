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
