# Galois-Field-Project
<< trying to do working arithmetic within finite fields with polynomials included >>

1. Finite Field Class file :
  contains the class of integers modulo p, for any prime p supposedly
  the class GaloisFq in the polynomial_class.py file inherits this class as it needs a bit more polynomial technology to exist

2. Irreducible file :
  raw algorithms file making up the idea from which the project emerged; the need to compute irreducible polynomials
  in the file 2 methods of sieving are used to compute irreducible polynomials in $Z/pZ[X]$

3. Polynomial Class file :
  main file, to some degree
  the polynomial class started from the heuristic decision of representing polynomials as dictionaries, and not lists, to circumvent the unnecessary amount of zeroes in memory
  this was done after trying to compute the euclidean division of $X^{p^n}-X$ by some polynomial $P$ in $Z/pZ[X]$
  -> if this results in a clean division (remainder zero), a theorem states that $P$ is irreducible, the problem is that the first polynomial has a very large degree and contains useless zeroes mostly
  2 algorithms are added at the end of the polynomial class, which are FFT and multiplication via FFT
  FFT is the Fast Fourier Transform, a method which is supposed to have complexity $O(n\cdot log(n))$ compared to the usual $O(n^2)$ complexity of normal polynomial multiplication
  (polynomial multiplication is just the polynomial resulting in the convolution of the coefficients of its products)
  -> instead of using complex numbers we use integers extended by $\omega = \exp(2i\pi/n)$ for some n (a power of 2 in the case of FFT), so the Z-algebra $Z[\omega]$ of dimension $\deg(\Theta_n(X))$, so that we dont have to accumulate floating points errors due to roots of unity often having irrational, if not transcendental, cartesian coordinates in the complex plane
  -> sadly this makes us forced to compute $A(1)\times B(1)$, $A(\omega)\times B(\omega)$, etc, which are products taking place in the algebra $Z[w]$, and need distributive multiplication, coefficient by coefficient, which is the exact thing we are trying to avoid (remember polynomial coefficients)
  -> /!\ (?) this could be avoided if before doing these multiplications we converted each element of the algebra $Z[w]$ back to its complex form, with enough acuracy as to avoid errors when doing the inverse FFT

4. Root of Unity Class file :
   add-on to the couple of FFT algorithms side project
   $\omega$ is a root of unity $\omega = \exp(2i\pi/n)$, n often being a power of 2
   this represents a number in $Z[w]$ (could be $Q[w]$) as a list of coordinates which follow expected arithmetic rules
   for example, if $\omega^4 = 1$ (this is $\omega = i$), then $(1+\omega)\times (1-\omega) = 1-\omega + \omega-\omega^2 = 1-\omega^2 = 1-(-1) = 2$, which this class works out if you say, compute Cyclo([1,1], 4) * Cyclo([1,-1], 4) (this should return [2,0], of the class Cyclo of order 4 and dimension 2)
   /!\ still in construction because working on a way to compute cyclotomic polynomials so that dimension is accurate even when n is not a power of 2
