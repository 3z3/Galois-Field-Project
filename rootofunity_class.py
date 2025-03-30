import cmath
from itertools import product
from itertools import combinations
import time

class Cyclo:
    #implementation of the A-algebra A[w] (w being a root of unity) needed to compute the FFT of a poly P in A[X]

    def __init__(self, number, order):
        #can expect a w string (omega)
        if number == 'w':
            self.value = [0 for k in range(order)]
            self.value[1] = 1
        else:
            self.value = number #this is a list

        #whenever the order of w is divisible by 2, w^(o/2) = -1 so the dimension of the A-vs A[w] is not the order
        if order%2 == 0:
            self.dimension = int(order/2)
            self.even = True
        else:
            self.dimension = order
            self.even = False
            #TODO implement a program that calculates cyclotomic polynomials and deduces what Z[w] looks like algebraically
            #dimension, expressions of powers of w above the cyclo poly's degree, etc

        self.order = order
        if len(number) == self.dimension:
            pass
        else:
            raise ValueError('dimension doesn\'t match the length of the given list')
        
    def __repr__(self):
        #how the number appears when called in command line (as seen inside a list or whatnot)
        return str(self.value)

    def __mul__(self,other):
        if isinstance(other, Cyclo) and self.even:
            try:
                res = [0]*self.dimension
                for i in range(self.dimension):
                    for j in range(self.dimension):
                        res[(i+j)%self.dimension] += self.minus(i+j)*self.value[i]*other.value[j]    #loops back around
                return Cyclo(res,self.order)
            except:
                raise IndexError('dimensions on both side of the multiplication don\'t match')
        elif isinstance(other, Cyclo) and (not self.even):
            try:
                res = [0]*self.dimension
                for i in range(self.dimension):
                    for j in range(self.dimension):
                        res[(i+j)%self.dimension] += self.value[i]*other.value[j]    #loops back around
                return Cyclo(res,self.order)
            except:
                raise IndexError('dimensions on both side of the multiplication don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            res = [x*other for x in self.value]
            return Cyclo(res,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for *")
        
    def __rmul__(self,other):
        if isinstance(other, Cyclo):
            try:
                res = [0]*self.dimension
                for i in range(self.dimension):
                    for j in range(self.dimension):
                        res[(i+j)%self.dimension] = self.value[i]*other.value[j]    #loops back around
                return Cyclo(res,self.order)
            except:
                raise IndexError('dimensions on both side of the multiplication don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            res = [other*x for x in self.value]
            return Cyclo(res,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for *")

    def __add__(self,other):
        if isinstance(other, Cyclo):
            try:
                return Cyclo([self.value[k]+other.value[k] for k in range(self.dimension)],self.order)
            except:
                raise IndexError('dimensions on both side of the addition don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            sum = self.value.copy()
            sum[0] += other
            return Cyclo(sum,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for +")

    def __radd__(self,other):
        #right side
        if isinstance(other, Cyclo):
            try:
                return Cyclo([self.value[k]+other.value[k] for k in range(self.dimension)],self.order)
            except:
                raise IndexError('dimensions on both side of the addition don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            sum = self.value.copy()
            sum[0] += other
            return Cyclo(sum,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for +")

    def __sub__(self,other):
        if isinstance(other, Cyclo):
            try:
                return Cyclo([self.value[k]-other.value[k] for k in range(self.dimension)],self.order)
            except:
                raise IndexError('dimensions on both side of the addition don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            sum = self.value.copy()
            sum[0] -= other
            return Cyclo(sum,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for +")

    def __rsub__(self,other):
        if isinstance(other, Cyclo):
            try:
                return Cyclo([other.value[k]-self.value[k] for k in range(self.dimension)],self.order)
            except:
                raise IndexError('dimensions on both side of the addition don\'t match')
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, complex):
            sum = [(-1)*x for x in self.value]
            sum[0] += other
            return Cyclo(sum,self.order)
        else:
            raise TypeError("Unsupported operand type(s) for +")
        
    def __pow__(self,other):
        if isinstance(other, int) and i > 0:
            res, i = Cyclo(self.copy(),self.order), other
            while i > 1:
                res *= self
                i -= 1
            return res
        elif isinstance(other, int) and i == 0:
            return Cyclo([1]+[0]*(self.dimension-1),self.order)
        elif isinstance(other, int) and i < 0:
            #TODO probably implement negative powers using cyclotomic polynomials ?
            #/!\ IN CONSTRUCTION /!\
            pass
        else:
            raise TypeError("Unsupported operand type(s) for **")

    def Rshift(self,shift=1):
        if self.even and shift: #shift != 0
            h = self.minus(shift)
            norm_shift = (-1)*(shift%self.dimension)
            left_side = [(-h)*x for x in self.value[norm_shift:]]
            right_side = [h*x for x in self.value[:norm_shift]]
            res = left_side + right_side
            return Cyclo(res, self.order)
        elif self.even and shift == 0:
            return self
        else:
            norm_shift = (-1)*(shift%self.dimension)
            res = self.value[norm_shift:] + self.value[:norm_shift]
            return Cyclo(res, self.order)

    def Lshift(self,shift=1):
        #shouldnt really exist
        return self.Rshift(-shift)
    
    def firstcoord(self):
        #projection of number onto the first coordinate (the A-component of number z in A[w])
        return self.value[0]

    def copy(self):
        #deepcopy of self.value
        res = [x for x in self.value]
        return res
    
    def minus(self,index):
        #gets called when order is a multiple of 2
        if (index%self.order) < self.dimension:
            return 1
        else:
            return -1
        
def wpower(shift,order):
    #returns the first root of unity (2ipi/n) of some order to the power of shift
    dimension = order

    if order%2 == 0:
        dimension = int(order/2)

    new_i = shift%dimension
    h = 1
    if (shift%order) >= dimension:
        h = -1

    res = [0]*(new_i) + [h] + [0]*(dimension-new_i-1)
    return Cyclo(res,order)

def gcd(a, b):
    a, b = abs(a), abs(b)   #make sure it's positive
    while b:
        a, b = b, a%b
    return a

class Zeta:
    #takes exponent k and order n as arguments, defines Zeta(k,n) which is = \exp{2i\pi k/n} an n-th root of unity

    Basis = {}  #one basis as a list for each distinct order, syntax is -> order : basis

    Powers = {} #lists all other powers as linear combinations of the elements in the basis MAYBE NEEDS TO BE IN CYCLO ?

    Primes = [] #primes, just in case we need them
    Factor = {} #for each order, keeps track of its prime factorization, syntax is -> order : {prime : valuation}
    
    def __init__(self,exponent,order):

        self.exponent = exponent
        self.order = order

        self.minus = False  #whether the order is 2 mod n or not (if it's 2, then 2 as a prime doesn't add to the dimension, so, useless)
        if order%4 == 2:
            self.minus = True

        if self.Factor.get(order,None) is None:
            self.Factor[order] = {}
            primes = [p for p in prime_sieve(order+1)]
            self.Primes += primes[len(self.Primes):]
            for p in self.Primes:
                if self.order%p == 0:
                    k = 1
                    while self.order%(p**(k+1)) == 0:
                        k += 1
                    self.Factor[order][p] = k
            primes = []
            if self.minus:
                self.Factor[order].pop(2)

        if self.Basis.get(order,None) is None:
            self.Basis[order] = []  #makes sure that when other zetas are created in get_basis, no other init enters this block
            self.Basis[order] = self.get_basis() #change this function if you want to change the method of generating the basis
            #self.Basis[order] = self.get_basis2()  #PLUS comment the next line because get_basis2 already computes the powers
            self.get_powers()

    def __repr__(self):
        #how zetas appear in command
        return 'Zeta(%s,%s)' % (str(self.exponent), str(self.order))

    def get_basis(self):
        #class function designed to compute the canonical basis explicited in (Bosma, 1990)
        #notation is consistent with the article, eg. A represents each a, C[p] is c_p, etc

        basis_exponents = []
        n = self.order
        
        A,C,R = {}, {}, {}

        for p,k in self.Factor[n].items():
            if p != 2:
                A[p] = [a for a in range(p-1)]
                C[p] = [c for c in range(p**(k-1))]
                R[p] = (n//p, n//(p**k))
            else:
                A[2] = [0,1]
                C[2] = [c for c in range(2**(k-2))]
                R[p] = (n//4, n//(p**k))

        for tupa in product(*[a for a in A.values()]):
            for tupc in product(*[c for c in C.values()]):
                #here we create one element of the basis
                sum = 0
                p_index = 0
                
                for p in self.Factor[n].keys():
                    sum += tupa[p_index]*R[p][0] + tupc[p_index]*R[p][1]
                    p_index += 1
                
                basis_exponents.append(sum%n)
            
        return [Zeta(k,n) for k in basis_exponents]
    
    def get_powers(self):
        #evaluates each other power of zeta_n (powers ranging from 0 to n-1)
        #returns a polynomial in dictionary form to represent how this root of unity expresses itself as a linear combination of roots in the basis

        n = self.order
        if n not in self.Powers.keys():
            self.Powers[n] = {}

        basis_exponents = {b.exponent for b in self.Basis[n]}
        up_to_n = {k for k in range(n)}

        for exponent in up_to_n:
            if exponent in basis_exponents:
                self.Powers[n][exponent] = {exponent : 1}   #linear combination seen as a polynomial, note that the coefficient is 1
            else:
                A,C,R = {}, {}, {}
                self.Powers[n][exponent] = {}   #polynomial

                for p,k in self.Factor[n].items():
                    if p != 2:
                        R[p] = (n//p, n//(p**k))
                        if n%(p**2) == 0:
                            C[p] = (pow(R[p][1],-1,p)*exponent)%p
                            A[p] = (pow(R[p][1],-1,p)*((exponent - C[p]*R[p][1])/(p**(k-1))))%p
                        else:
                            C[p] = 0
                            A[p] = (pow(R[p][1],-1,p)*exponent)%p
                    else:
                        R[2] = (n//4, n//(2**k))
                        if n%8 == 0:
                            C[2] = (pow(R[2][1],-1,4)*exponent)%4
                            A[2] = (pow(R[2][1],-1,4)*((exponent - C[2]*R[2][1])/(2**(k-2))))%4
                        else:
                            C[2] = 0
                            A[2] = (pow(R[2][1],-1,4)*exponent)%4

                zeta_n_exponent = 0
                bad_exponents = []

                for p in self.Factor[n].keys():
                    if p != 2:
                        if (A[p] + 1)%p == 0:
                            bad_exponents.append(p)
                        else:
                            zeta_n_exponent += A[p]*R[p][0]
                        zeta_n_exponent += C[p]*R[p][1]
                    else:
                        if A[2]%4 == 2 or A[2]%4 == 3:
                            bad_exponents.append(2)
                        else:
                            zeta_n_exponent += A[2]*R[2][0]
                        zeta_n_exponent += C[2]*R[2][1]

                for prod in product(*[linearize(A[p],p) for p in bad_exponents]):    #every big product in the linear sum
                    mult = multiply(prod,n)
                    index = (mult[0]+zeta_n_exponent)%n
                    self.Powers[n][exponent][index] = self.Powers[n][exponent].get(index, 0) + mult[2]

    def get_basis2(self):
        #second way of computing a basis for Q(zeta_n); the simplest basis, the other powers are computed directly by doing successive polynomial divisions using the n-th cyclotomic polynomial
        n = self.order
        primes = self.Factor[n]

        if n not in self.Powers.keys():
            self.Powers[n] = {}

        numerator, denominator = {0 : 1}, {0 : 1}

        for repeat in range(1,len(primes)+1):
            for combo in combinations(primes, repeat):
                #multiplying each prime divisor of n with exponent 0 or 1, bc if 2 or more, then the mobius function is 0
                #here we compute the n-th cyclotomic polynomial as a quotient of two big products (numerator & denominator) applying the formula (use latex) :
                #$\Phi_n(x) = \prod_{d|n}\left(x^{n/d} - 1\right)^{\mu(d)}$ where \mu is the Möbius arithmetic function
                
                divisor = 1
                for p in combo:
                    divisor *= p

                if repeat % 2 == 1:
                    denominator = mult(denominator,{n // divisor : 1, 0 : -1}) #polynomial X^{n/d} - 1
                else:
                    numerator = mult(numerator,{n // divisor : 1, 0 : -1})

        #add the part where you don't choose any prime, so divisor = 1 and n // divisor = n
        numerator = mult(numerator,{n : 1, 0 : -1})

        #Phi is the n-th cyclotomic polynomial
        Phi = div(numerator,denominator)[0]
        degPhi = max([k for k in Phi.keys()]+[0])

        for power in range(n):
            if power < degPhi:
                self.Powers[n][power] = {power : 1}
            else:
                #this long expression is just euclidean division of Phi to get the next term in the sequences of powers of zeta
                #expressed as linear combinations in the basis zeta^0, zeta^1, ... , zeta^d-1 ; where d is the degree of Phi
                self.Powers[n][power] = {exponent : -self.Powers[n][power-1].get(degPhi-1,0)*Phi.get(exponent,0) + self.Powers[n][power-1].get(exponent-1,0) for exponent in range(degPhi)}
            
        return [Zeta(i,n) for i in range(degPhi)]

def linearize(k,p):
    #converts zeta(k,p) to its linear expression in the canonical basis, where p is prime, and k is generally p-1
    #used in get_basis
    if p != 2:
        if k == p-1:
            return [(i,p,-1) for i in range(p-1)]   #-1 is the minus sign of -e(i/p)
        else:
            return [(k,p,1)]
    else:
        if k == 2:
            return [(0,4,-1)]
        elif k == 3:
            return [(1,4,-1)]
        else:
            return [(k,p,1)]
    
def multiply(linear,n):
    #takes a bunch of zeta(k_p,p) for primes p in the factorization of n, and returns zeta(k,n) 
    #used in get_basis
    sum = 0
    sign = 1
    for zeta in linear:
        sum += zeta[0]*(n//zeta[1])
        sign *= zeta[2]
    return (sum%n, n, sign)

def prime_sieve(limit):
    #prime sieve function returning the next prime as an iterator object
    a = [True] * limit                          # Initialize the primality list
    a[0] = a[1] = False

    for (i, isprime) in enumerate(a):
        if isprime:
            yield i
            for n in range(i*i, limit, i):     # Mark factors non-prime
                a[n] = False

def mult(P,Q):
    #implemation of polynomial multiplication outside of the Polynomial class
    degP, degQ = max([k for k in P.keys()]+[0]), max([k for k in Q.keys()]+[0])
    Pi = {}
    for j in range(degP + degQ + 1):
        coeff = 0
        i = max(0,j-degQ)
        while i <= min(j,degP):   #making sure the index stays within each of the polynomials' indices
            coeff += P.get(i,0)*Q.get(j-i,0)
            i += 1
        Pi[j] = coeff
    return Pi

def div(A,B):
    #implemation of polynomial division outside of the Polynomial class
    degA, degB = max([k for k in A.keys()]+[0]), max([k for k in B.keys()]+[0])
    try:
        R,Q,n = A.copy(),{},degA
    except:
        print('the degree of A is less than the degree of B !')

    while n >= degB:
        Q[n-degB] = Q.get(n-degB,0) + R[n]/B[degB]
        for i in range(degB+1):
            R[n-degB+i] = R.get(n-degB+i,0)-B.get(i,0)*(R[n]/B[degB])
        for key in sorted(R.keys()):
            if not R[key]:
                del R[key]
        if len(R) == 0:
            R = {0:0}
        n = max([key for key in R.keys()]+[0])  #degree of R gets an update
    
    return (Q,R)

def copy(d):
    #deep copy of a dictionary
    dict = {k:x for (k,x) in d.items()}
    return dict

#START OF TESTS
###

#n is the order of your root of unity
#change between get_basis and get_basis2 to compare the time taken for each basis + powers to be computed

n = int(input('n = '))

start = time.time()
a = Zeta(1,n)    #equivalent to the element \exp{2i\pi/n}
print('time spent creating basis and powers = ' + str(time.time() - start))
print(len(a.Basis[n]))
