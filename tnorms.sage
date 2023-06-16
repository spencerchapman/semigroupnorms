def t_length(v,t):
    
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
def T_Length_Set(x,t):
    F = S.Factorizations(x)
    Lengths = []
    for v in F:
        Lengths.append(t_length(v,t))
    return sorted(Lengths)
def T_Max_Length(x,t):
    return T_Length_Set(x,t)[-1]
def T_Min_Length(x,t):
    return T_Length_Set(x,t)[0]
def T_Elasticity(x,t):
    return T_Max_Length(x,t)/T_Min_Length(x,t)
