
# Galois-Field-Project
> *trying to do working arithmetic within finite fields with polynomials included*

## 1. Finite Field Class file :
contains the class of integers modulo p, for any prime p supposedly  
the class GaloisFq in the `polynomial_class.py` file inherits this class as it needs a bit more polynomial technology to exist

## 2. Irreducible file :
raw algorithms file making up the idea from which the project emerged; the need to compute irreducible polynomials  
in the file 2 methods of sieving are used to compute irreducible polynomials in $\displaystyle \mathbf{Z}/p\mathbf{Z}[X]$

## 3. Polynomial Class file :
main file, to some degree  
the polynomial class started from the heuristic decision of representing polynomials as dictionaries, and not lists, to circumvent the unnecessary amount of zeroes in memory  
this was done after trying to compute the euclidean division of $\displaystyle X^{p^n}-X$ by some polynomial $P$ in $\displaystyle \mathbf{Z}/p\mathbf{Z}[X]$  
-> if this results in a clean division (remainder zero), a theorem states that $P$ is irreducible, the problem is that the first polynomial has a very large degree and contains useless zeroes mostly  
2 algorithms are added at the end of the polynomial class, which are FFT and multiplication via FFT  
FFT is the Fast Fourier Transform, a method which is supposed to have complexity $\displaystyle O(n\cdot log(n))$ compared to the usual \displaystyle $O(n^2)$ complexity of normal polynomial multiplication  
(polynomial multiplication is just the polynomial resulting in the convolution of the coefficients of its products)  
-> instead of using complex numbers we use integers extended by $\omega = \exp(2i\pi/n)$ for some n (a power of 2 in the case of FFT), so the $\mathbf{Z}$-algebra $\displaystyle \mathbf{Z}[\omega]$ of dimension $\displaystyle \deg(\Phi_n(X))$, so that we dont have to accumulate floating points errors due to roots of unity often having irrational, if not transcendental, cartesian coordinates in the complex plane  
-> sadly this makes us forced to compute $A(1)\times B(1)$, $A(\omega)\times B(\omega)$, etc, which are products taking place in the algebra $\mathbf{Z}[\omega]$, and need distributive multiplication, coefficient by coefficient, which is the exact thing we are trying to avoid (remember polynomial coefficients)  
The other implementation of FFT aka the `ComplexFourierMult` function begins computations in the `Cyclo` class, then converts to complex numbers before doing the inverse FFT of each $C(\omega^k) = A(\omega^k)\times B(\omega^k)$ and rounds each coefficient to the closest integer in complex algebraic form, when the imaginary part is rounded to 0 (happens when $\Im(z) \in (-0.5,0.5)$ we consider that the coefficient is an integer, as expected.  
By doing that, we multiply floating point numbers (negligeable complexity compared to polynomial multiplication taking place in `Cyclo` class, used in `FourierMult` function) still taking advantage of the simplicity of basic arithmetic operations taking place in `Cyclo` class in the first call of the `FFT` function (which calculates evaluations of $A$ and $B$ at roots of unity), then knowing the result of the inverse FFT of the floating points complex coefficients is going to be very close to the expected integer coefficients of $C$.  
We know that from **(W. M. Gentleman & G. Sande, 1966, page 570)** derived from a result in **(Wilkinson, 1963)** which gives an upper bound on the floating point error on the multiplication of two matrices given their matrix norm, a b-bit mantissa (or significand) and dimensions.
```math
\displaystyle |A(\omega^k)| = \left|\sum_{0\leq k\leq n-1} a_k\omega^k\right| \leq \sum_{0\leq k\leq n-1}|a_k| =: |A|_1
```
thus,
```math
|\hat{C}|_2 = \left( \sum_{0\leq k\leq n-1} |A(\omega^k)\times B(\omega^k)|^2 \right)^{1/2} \leq \sqrt{n} |A|_ 1\times |B|_1
```
So, if $M$ represents the matrix multiplication which takes the Fourier Transform of $C$, namely $\hat{C}$, to $C$, then we can write :  
```math
|\text{error}(C)|_2 = |\text{error}(M\times \hat{C})|_2 \leq 1.06\cdot n\cdot 2^{-b}\cdot |M|_2\cdot |\hat{C}|_2 \leq  1.06\cdot n^2\cdot 2^{-b}\cdot |A|_ 1\times |B|_1
```
Since we know coefficients of $A$ and $B$, and the power of 2 that comes after the degree of $C$ (namely, $n$), we can mitigate this error precisely since it depends only on the variable $b$.

## 4. Root of Unity Class file :
add-on to the couple of FFT algorithms side project  
$\omega$ is a root of unity $\omega = \exp(2i\pi/n)$, n often being a power of 2  
this represents a number in $\mathbf{Z}[\omega]$ (could be $\mathbf{Q}[\omega]$) as a list of coordinates which follow expected arithmetic rules  
for example, if $\omega^4 = 1$ (this is $\omega = i$), then $(1+\omega)\times (1-\omega) = 1-\omega + \omega-\omega^2 = 1-\omega^2 = 1-(-1) = 2$, which this class works out if you say, compute `Cyclo([1,1], 4) * Cyclo([1,-1], 4)` (this should return `[2,0]`, of the class `Cyclo` of order 4 and dimension 2)  
I make use of the article **(Bosma, 1990)** to define another class `Zeta` anterior to `Cyclo` which represents successive powers of a root of unity $\zeta_n=\exp(2i\pi/n)$ of degree $n$ as `Zeta(k,n)`, in this file we compare the efficiency of the canonical basis constructed by W. Bosma at computing the remaining powers of $\zeta_n$ outside of this very basis, compared to, say, another easy way of computing these; set the basis as the powers of $\zeta_n$ strictly under the degree $d$ of the $n$-th cyclotomic polynomial, so $\zeta_n^0, \zeta_n^1, \cdots ,\zeta_n^{d-1}$, then compute the remaining powers $\zeta_n^k$ for $d\leq k\leq n-1$ using polynomial divisions and evaluations at $\zeta_n$.  
We observe for example that for $n = 4005$ the first method takes roughly 37 ms and the latter takes 2800 ms, at least for this rough implementation of both techniques, the computation of the $n$-th cyclotomic polynomial may not be the best and is derived from this formula :
```math
\Phi_n(x) = \prod_{d|n} \left(x^{n/d} - 1\right)^{\mu(d)}
```
Where $\mu$ is the Möbius arithmetic function.
