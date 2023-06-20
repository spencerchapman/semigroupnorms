# semigroupnorms

## Usage

### T_Length(v,t)
Takes a vector v and a value t, returns the t-length factorization of v. For infinity let t = oo

### T_LengthSet(S,x,t)
Takes a semigroup S, an element x and a value t, generates the length set wrt t

### T_MaxLength(S,x,t)
Returns the max of T_LengthSet(S,x,t)

### T_MinLength(S,x,t)
Analogous to above

### T_Elasticity(S,x,t)
Returns the elasticity of x wrt t

### T_FirstNTerms(S,n)
Returns the first n terms in S in increasing order

### T_MaxPowersOfElement(S,x,t,n)
Returns the list [L_t(x^0), L_t(x^1), L_t(x^2), ... L_t(x^n)]

### T_MinPowersOfElement(S,x,t,n)
Analogous to above

### T_FullMaxPowers(S,m,t,n)
Returns a matrix of the powers of the first m elements of S, from 0 up to n
Each column is T_MaxPowersOfElement(S,x,t,n), for the first m terms x in S
