# SO3_singlets
Code to compute SO(3) singlets inside a Fock space of SO(3) oscillators

Consider a quantum mechanical theory of N real scalar fields (or flavors) on the sphere S^2. The full Hilbert space
of the theory is a Fock space of creation operators a^\dagger_{l,m}, b^\dagger_{l,m}, ... where the oscillators
{a, b, ...} correspond to the N different fields in the theory.

The Hilbert space falls into irreducible representations of SO(3): scalars (total spin L=0), vectors (L=1) and so forth.
In many computations, it is advantageous to organize states in this manner. However, Fock space states live are tensor products
of oscillators: a generic Fock space state

  a^\dagger_{l1 m1} b^\dagger_{l2 m2} ... | 0 >
  
does not transform in any spin-L representation.

It is easy to see SO(3) only mixes states with the same spins (l1, l2, ...). Given a tuple (l1, ...), 
this subspace has dimension \prod_i (2 l_i + 1). The algorithm presented here performs two steps:

1) given a tuple of positive integers (l1, ...), compute the list of all allowed states that have L_3 = 0
   (meaning that m1 + m2 + ... = 0)
2) inside this subspace, construct the set of all actual SO(3) singlets, having L = 0.

The second step can be done using various brute-force methods (e.g. finding states that obey L_{-} = 0)
but such algorithms tend to be very time-consuming.
