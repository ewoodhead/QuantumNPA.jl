using Test



## Tests of basic operator properties

PA = projector(1, 1:2, 1:2)
PB = projector(2, 1:2, 1:2)



@test PA[1,1]*PB[1,1] == PB[1,1]*PA[1,1]

@test iszero(PA[1,1]*PA[2,1])

@test PA[1,1]*PA[1,1] == PA[1,1]



UA = unitary(1, 1:3)
PA11 = PA[1,1]
(UA1, UA2, UA3) = UA
(UAc1, UAc2, UAc3) = conj.(UA)

V = UA[1] * conj(UA[2]) * UA[3]

@test conj(V) == UAc3 * UA2 * UAc1



P = PA[1,1] + V

result = PA11 + UAc3 * UA2 * UAc1
@test conj(P) == result

result = (PA11 + PA11 * UA1 * UAc2 * UA3
          + UA1 * UAc2 * UA3 * PA11
          + UA1 * UAc2 * UA3 * UA1 * UAc2 * UA3)
@test P*P == result

result = (Id + PA11 + PA11 * UA1 * UAc2 * UA3
          + UAc3 * UA2 * UAc1 * PA11)
@test conj(P)*P == result

Q = Id + V*PA[1,1]

result = (Id + 2 * UA1 * UAc2 * UA3 * PA11
          + UA1 * UAc2 * UA3 * PA11 * UA1 * UAc2 * UA3 * PA11)
@test Q*Q == result

result = (Id + PA11 + PA11 * UAc3 * UA2 * UAc1
          + UA1 * UAc2 * UA3 * PA11)
@test conj(Q)*Q == result



ZE = zbff(5, 1:2)

R = PA[1,1] * PB[2,2] * ZE[1] * ZE[2]

result = PA[1,1] * PB[2,2] * conj(ZE[2]) * conj(ZE[1])
@test conj(R) == result



FA = fourier(1, 9, 1, 5)

@test isidentity(FA^0)
@test conj(FA^3) == FA^2
@test isidentity(FA^5)
@test FA^6 == FA^1
@test isidentity(conj(FA^3)*FA^3)
@test FA*FA == FA^2
@test conj((FA*FA)^4) == FA^2



A1, A2 = dichotomic(1, 1:2)
B1, B2 = dichotomic(2, 1:2)
S = A1*(B1 + B2) + A2*(B1 - B2)

@test S^2 == 4*Id - comm(A1, A2) * comm(B1, B2)



## Tests with dichotomic & multiparty operators

A1, A2 = dichotomic(1, 1:2)
B1, B2 = dichotomic(2, 1:2)
A_C1 = dichotomic([1,3], 1)

@test B1 * A1 == A1 * B1
@test A1 * A_C1 != A_C1 * A1
@test A1 * B1 * A_C1 == A1 * A_C1 * B1
@test A1 * B1 * A_C1 * A2 * B2 == A1 * B1 * B2 * A_C1 * A2
@test A1 * B1 * A_C1 * A1 * B1 == A1 * A_C1 * A1

X = A1 * B1
Y = A_C1 * X

@test Y == B1 * A_C1 * A1
@test A_C1 * Y == A1 * B1



## Some simple type checks

P = projector(1, 1, 1)

@test P isa Monomial
@test !(P isa Polynomial)

Q = 1*P

@test Q isa Polynomial
@test !(Q isa Monomial)

S = A1*B1 + A1*B2 + A2*B1 - A2*B2

@test S isa Polynomial

x = Polynomial(1)
y = Polynomial(Id)
z = Polynomial(S)

@test x == Id
@test y == Id
@test z === S
@test typeof.([x, y, z]) == [Polynomial, Polynomial, Polynomial]




## Tests of npa_max()

# Test CHSH

@dichotomic A1 A2 B1 B2

S = A1*(B1 + B2) + A2*(B1 - B2)

@test npa_max(S, 2) ≈ sqrt(8) atol=1e-3



# Test Svetlichny

@dichotomic A[1:2] B[1:2] C[1:2]

E(x,y,z) = A[x]*B[y]*C[z]

S = (-E(1,1,1) + E(1,1,2) + E(1,2,1) + E(1,2,2)
     + E(2,1,1) + E(2,1,2) + E(2,2,1) - E(2,2,2))

@test npa_max(S, "1 + A B + A C + B C") ≈ sqrt(32) atol=1e-3



# Test modified CHSH

b = 0.3
a = 0.6

S = b * A1 + a * A1*(B1 + B2) + A2*(B1 - B2)
f = ((b, a) -> 2*sqrt((1 + a^2)*(1 + 0.25*b^2)))

@test npa_max(S, "1 + A B + A^2 B") ≈ f(b, a) atol=1e-3



# Test max of <A1> with constraint that <A1 (B1 + B2)> = <A2 (B1 - B2)>

S1 = A1*(B1 + B2)
S2 = A2*(B1 - B2)
eq_constraints = [S1 - 1.4*Id, S2 - 1.4*Id]

f = ((s1, s2) -> sqrt(1 - s2^2/(4 - s1^2)))

@test npa_max(A1, 2, eq=eq_constraints) ≈ f(1.4, 1.4) atol=1e-3
 


# Test CH74 form of CHSH

PA11, PA12 = projector(1,1,1:2)
PB11, PB12 = projector(2,1,1:2)

S = -PA11 - PB11 + PA11*(PB11 + PB12) + PA12*(PB11 - PB12)

@test npa_max(S, 1) ≈ 0.5*(sqrt(2.0) - 1) atol=1e-3



# Test CGLMP with built-in cglmp() function

@test npa_max(cglmp(3), "1 + A B") ≈ 1 + sqrt(11/3) atol=1e-3



# Guessing probability test

PA = projector(1, 1:2, 1:2, full=true)
PB = projector(2, 1:2, 1:2, full=true)
PE = projector(5, 1:4, 1, full=true)

# CHSH = 2*sqrt(2) * p
p = 0.9

# Expectation value of G is the probability that Eve correctly guesses
# Alice's and Bob's joint outcome.
G = sum(PA[a,1] * PB[b,1] * PE[2*(a-1) + b]
        for a in 1:2 for b in 1:2)

# Ideal CHSH-violating correlations mixed with noise. N.B., the actual
# constraints imposed are on the expectation values of the operators
# in the array.
constraints = [PA[1,1] - 0.5*Id,
               PA[1,2] - 0.5*Id,
               PB[1,1] - 0.5*Id,
               PB[1,2] - 0.5*Id,
               PA[1,1]*PB[1,1] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,1]*PB[1,2] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,2]*PB[1,1] - 0.25*(1 + p/sqrt(2))*Id,
               PA[1,2]*PB[1,2] - 0.25*(1 - p/sqrt(2))*Id]

# This returns about 0.7467 for p = 0.9 at level 2 using the default SCS
# solver. Solution for comparison was computed separately using my older
# npa-hierarchy library in Lisp and the arbitrary precision SDP solver
# SDPA-GMP

pguess_sol = 0.7461756908058874
@test npa_max(G, 2, eq=constraints) ≈ pguess_sol atol=1e-3
