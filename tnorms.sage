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
    return sorted(list(Set([T_Length(f,t) for f in S.Factorizations(n)])))
def T_MaxLength(S,x,t):
    return T_LengthSet(S,x,t)[-1]
def T_MinLength(S,x,t):
    return T_LengthSet(S,x,t)[0]
def T_Elasticity(S,x,t):
    return T_MaxLength(S,x,t)/T_MinLength(S,x,t)
