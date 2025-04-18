import cmath
import time
from itertools import product
from itertools import combinations
from rationals import *

#start of elements of the field Q(zeta) class Cyclo

class Cyclo:
    #implementation of the Field Q[w] (w being a root of unity)
    #TODO IMPLEMENT COMPATIBILITY BETWEEN DIFFERENT ORDERS USING W BOSMA

    Basis = {}  #one basis as a list for each distinct order, syntax is -> order : basis

    Powers = {} #lists all other powers as linear combinations of the elements in the basis, syntax is -> order : { exponent : polynomial }
    #polynomial is a dictionary of the form { exponent : coefficient }

    Cyclotomic_Polynomial = {}  #lists cyclotomic polynomials by order, syntax is -> order : polynomial

    def __init__(self, number, order):
        #takes a dictionary and an order as input, dic is a polynomial in w, and order is the order of the root of unity w

        self.value = number #this is a dictionary
        self.order = order  #this is an int
        self.zeta = Zeta(1,order)

        if self.order not in self.Basis:
            self.Basis[self.order] = self.zeta.Basis[self.order]
            self.Powers[self.order] = self.zeta.Powers[self.order]

            #optional
            self.Cyclotomic_Polynomial[self.order] = self.zeta.get_phi()
        
        self.dimension = len(self.Basis[self.order])

        #this part expresses the number in the canonical basis of Q[w] and removes 0 coefficients from the dictionary
        exponents_in_basis = {z.exponent for z in self.Basis[self.order]}
        buffer = self.value.copy()
        for key, value in self.value.items():
            if value == 0:
                del buffer[key]
            else:
                if key not in exponents_in_basis:
                    for exp, coeff in self.Powers[self.order][key%self.order].items():  #divide by order to allow for arbitrarily large powers
                        buffer[exp] = buffer.get(exp,0) + coeff*value
                    del buffer[key]
        self.value = buffer
        
    def __repr__(self):
        #how the number appears when called in command line (as seen inside a list or whatnot)
        return str(self.value)

    def __mul__(self,other):
        if isinstance(other, Cyclo) and self.order == other.order:
            res = {}
            for i, coeff_i in self.value.items():
                for j, coeff_j in other.value.items():
                    exponent = (i+j)%self.order
                    for bas_exp, coeff in self.Powers[self.order][exponent].items():
                        res[bas_exp] = res.get(bas_exp,0) + (coeff*coeff_i*coeff_j)
            return Cyclo(res,self.order)
        else:
            try:
                res = {key : value*other for key,value in self.value.items()}
                return Cyclo(res,self.order)
            except:
                error('*',self,other)
        
    def __rmul__(self,other):
        if isinstance(other, Cyclo) and self.order == other.order:
            res = {}
            for i, coeff_i in self.value.items():
                for j, coeff_j in other.value.items():
                    exponent = (i+j)%self.order
                    for bas_exp, coeff in self.Powers[self.order][exponent].items():
                        res[bas_exp] = res.get(bas_exp,0) + (coeff*coeff_j*coeff_i)
            return Cyclo(res,self.order)
        else:
            try:
                res = {key : other*value for key,value in self.value.items()}
                return Cyclo(res,self.order)
            except:
                error('*',self,other)

    def __add__(self,other):
        if isinstance(other, Cyclo):
            res = self.copy()
            for key, value in other.value.items():
                res[key] = res.get(key,0) + value
            return Cyclo(res,self.order)
        else:
            try:
                sum = self.value.copy()
                sum[0] = sum.get(0,0) + other
                return Cyclo(sum,self.order)
            except:
                error('+',self,other)          

    def __radd__(self,other):
        if isinstance(other, Cyclo):
            res = other.copy()
            for key, value in self.value.items():
                res[key] = res.get(key,0) + value
            return Cyclo(res,self.order)
        else:
            try:
                sum = self.value.copy()
                sum[0] = other + sum.get(0,0)
                return Cyclo(sum,self.order)
            except:
                error('+',self,other)

    def __sub__(self,other):
        if isinstance(other, Cyclo):
            res = self.copy()
            for key, value in other.value.items():
                res[key] = res.get(key,0) - value
            return Cyclo(res,self.order)
        else:
            try:
                sum = self.value.copy()
                sum[0] = sum.get(0,0) - other
                return Cyclo(sum,self.order)
            except:
                error('-',self,other)

    def __rsub__(self,other):
        if isinstance(other, Cyclo):
            res = other.copy()
            for key, value in self.value.items():
                res[key] = res.get(key,0) - value
            return Cyclo(res,self.order)
        else:
            try:
                sum = {key : -value for key,value in self.value}
                sum[0] = other + sum.get(0,0)
                return Cyclo(sum,self.order)
            except:
                error('-',self,other) 
            
    def __truediv__(self,other):
        if isinstance(other, Cyclo) and self.order == other.order:
            try:
                res = self * other.inversetimesnorm() / other.norm()
                return res
            except:
                raise ZeroDivisionError('division by zero')
        else:
            if not other:
                raise ZeroDivisionError('division by zero')
            else:
                try:
                    res = {key : value / other for key,value in self.value.items()}
                    return res
                except:
                    error('/',self,other)
                
    def __rtruediv__(self,other):
        if isinstance(other, Cyclo) and self.order == other.order:
            try:
                res = other * self.inversetimesnorm() / self.norm()
                return res
            except:
                raise ZeroDivisionError('division by zero')
        else:
            try:    #other could be some int or some float idk
                res = other * self.inversetimesnorm() / self.norm()
                return res
            except:
                error('/',self,other)
        
    def __pow__(self,other):
        if isinstance(other, int) and other > 0:
            res, i = Cyclo(self.copy(),self.order), other
            while i > 1:
                res *= self
                i -= 1
            return res
        elif isinstance(other, int) and other == 0:
            ord = self.order
            return Cyclo({0:1},ord)
        elif isinstance(other, int) and other < 0:
            negative, i = (1 / self), other
            while i < -1:
                res *= negative
                i += 1
            return res
        else:
            error('**',self,other)

    def Rshift(self,shift=1):
        #multiplies by w^shift, defaults to multiplication by w
        n = self.order
        shifted = {(key + shift)%n : value for key,value in self.value.items()}
        res = {}
        for power, value in shifted.items():
            for exp, coeff in self.Powers[n][power].items():
                res[exp] = res.get(exp,0) + coeff*value
        return Cyclo(res,n)

    def Lshift(self,shift=1):
        #doesnt really need to exist
        return self.Rshift(-shift)
    
    def firstcoord(self):
        #projection of number onto the first coordinate (the A-component of number z in Q[w])
        return self.value.get(0,'no coordinate in Q')

    def copy(self):
        #deepcopy of self.value
        res = {key : value for key,value in self.value.items()}
        return res
    
    def inversetimesnorm(self):
        #computes the inverse of self using conjugacy (in the extension) -> still need to divide by the norm
        coprimes, n, res = [], self.order, 1
        for k in range(2,n):
            if gcd(k,n) == 1:
                coprimes.append(k)  #coprimes except for 1
        for k in coprimes:
            new = Cyclo({k * key : value for key,value in self.value.items()},n)
            res *= new
        return res
    
    def norm(self):
        #computes the norm of self using the inverse method
        norm = self.inversetimesnorm() * self
        if len(norm.value) == 1:
            return norm.value[0]
        else:
            raise ValueError('norm isn\'t a real number')
        
#end of elements of the field Q(zeta) class Cyclo

def gcd(a,b):
    #a, b = abs(a), abs(b)  #BREAK GLASS IN CASE OF EMERGENCY
    while b:
        a, b = b, a%b
    return abs(a)

def lcm(a,b):
    #least common multiple
    return int((a*b)/gcd(a,b))

def gcd_extended(*args):
    #gcd for a tuple
    if len(args) > 1:
        return gcd_extended(gcd(args[0],args[1]),*args[2:])
    else:
        return args[0]
    
def lcm_extended(*args):
    #lcm for a tuple
    if len(args) > 1:
        return lcm_extended(lcm(args[0],args[1]),*args[2:])
    else:
        return args[0]

#start of roots of unity class Zeta

class Zeta:
    #TODO add some arithmetic

    Basis = {}  #one basis as a list for each distinct order, syntax is -> order : basis

    Powers = {} #lists all other powers as linear combinations of the elements in the basis, syntax is -> order : { exponent : polynomial }
    #polynomial is a dictionary of the form { exponent : coefficient }

    Primes = [] #primes, just in case we need them
    Factor = {} #for each order, keeps track of its prime factorization, syntax is -> order : {prime : valuation}
    
    def __init__(self,exponent,order):

        self.exponent = exponent
        self.order = order

        #self.primitive = False
        #if gcd(exponent,order) == 1:
        #    self.primitive = True
        #NOT NECESSARY ?

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
    
    def __str__(self):
        #used to be blank, comment to fix
        return 'Zeta(%s,%s)' % (str(self.exponent), str(self.order))
        #pass
    
    def __bool__(self):
        #boolean value of a root of unity, considers 1 as the neutral element, like 0, se returns False in this case
        return bool(self.exponent)

    def get_basis(self):

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
    
    def get_polycyclo(self,n,primes = [], divisors=[]):
        
        #TODO fix that shit
        if len(primes) == 0:
            primes = self.Primes
        if n in primes:
            return {k:1 for k in range(n)}
        if len(divisors) == 0:
            divisors = []
            for powers in product(*[range(k+1) for k in self.Factor[n].values()]):
                p_index = 0
                d = 1
                for p in self.Factor[n].keys():
                    d *= p**powers[p_index]
                    p_index += 1
                divisors.append(d)
        Phi = {n:1,0:-1}
        for d in divisors:
            new_div = []
            for e in divisors:
                if d%e == 0:
                    new_div.append(d//e)
            Phi = div(Phi,self.get_polycyclo(d,primes,new_div))[0]
        return Phi
    
    def get_phi(self):
        #returns the n-th cyclotomic polynomial, where n is the order of self
        #contained in get_basis2
        n = self.order
        primes = self.Factor[n]
        numerator, denominator = {0 : 1}, {0 : 1}

        for repeat in range(1,len(primes)+1):
            for combo in combinations(primes, repeat):
                #multiplying each prime divisor of n with exponent 0 or 1, bc if 2 or more, then the mobius function is 0
                divisor = 1
                for p in combo:
                    divisor *= p

                if repeat % 2 == 1:
                    denominator = mult(denominator,{n // divisor : 1, 0 : -1}) #polynomial X^{n/d} - 1
                else:
                    numerator = mult(numerator,{n // divisor : 1, 0 : -1})

        #add the part where you don't choose any prime, so divisor = 1 and n // divisor = n
        numerator = mult(numerator,{n : 1, 0 : -1})

        Phi = div(numerator,denominator)[0]
        return Phi

    def get_basis2(self):
        n = self.order
        primes = self.Factor[n]

        if n not in self.Powers.keys():
            self.Powers[n] = {}

        numerator, denominator = {0 : 1}, {0 : 1}

        for repeat in range(1,len(primes)+1):
            for combo in combinations(primes, repeat):
                #multiplying each prime divisor of n with exponent 0 or 1, bc if 2 or more, then the mobius function is 0
                divisor = 1
                for p in combo:
                    divisor *= p

                if repeat % 2 == 1:
                    denominator = mult(denominator,{n // divisor : 1, 0 : -1}) #polynomial X^{n/d} - 1
                else:
                    numerator = mult(numerator,{n // divisor : 1, 0 : -1})

        #add the part where you don't choose any prime, so divisor = 1 and n // divisor = n
        numerator = mult(numerator,{n : 1, 0 : -1})

        Phi = div(numerator,denominator)[0]
        degPhi = max([k for k in Phi.keys()]+[0])

        for power in range(n):
            #self.Powers[n][power] = {} -> maybe not useful since we update dictionaries entirely aftewards ?

            if power < degPhi:
                self.Powers[n][power] = {power : 1}
            else:
                #this long expression is just euclidean division of Phi to get the next term in the sequences of powers of zeta
                #expressed as linear combinations in the basis zeta^0, zeta^1, ... , zeta^d-1 ; where d is the degree of Phi
                self.Powers[n][power] = {exponent : -self.Powers[n][power-1].get(degPhi-1,0)*Phi.get(exponent,0) + self.Powers[n][power-1].get(exponent-1,0) for exponent in range(degPhi)}
            
        return [Zeta(i,n) for i in range(degPhi)]
    
#end of roots of unity class Zeta

#functions used in the creation of a basis in class Zeta

def linearize(k,p):
    #converts zeta(k,p) to its linear expression in the canonical basis, where p is prime, and k is generally p-1
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
    sum = 0
    sign = 1
    for zeta in linear:
        sum += zeta[0]*(n//zeta[1])
        sign *= zeta[2]
    return (sum%n, n, sign)

def prime_sieve(limit):
    a = [True] * limit                          # Initialize the primality list
    a[0] = a[1] = False

    for (i, isprime) in enumerate(a):
        if isprime:
            yield i
            for n in range(i*i, limit, i):     # Mark factors non-prime
                a[n] = False

def substract(P,Q):
    #substracts Q from P ie (P-Q)
    res = P.copy()
    for exp, coeff in Q.items():
        res[exp] = res.get(exp,0) - coeff
    
    #gets rid of zeroes
    buffer = res.copy()
    for exp, coeff in buffer.items():
        if coeff == 0:
            del res[exp]
    if len(res) == 0:
        res = {0:0}
    
    return res

def mult(P,Q):
    #returns a polynomial equal to the product PxQ with coefficients until the sum of their respective degrees
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
    #returns quotient and remainder of the euclidean division of A by B, should correspond to A = BQ + R, here, A = self
    #outside of poly class so A and B are dictionaries from the get go
    R, Q, n, b = A.copy(), {}, max([key for key in A.keys()]+[0]), max([key for key in B.keys()]+[0])
    bmax = max(b-1,0)
    last_coeff_B = B[b]

    while n > bmax:

        Q[n-b] = R[n]/last_coeff_B  #update to quotient Q

        for exp, coeff in B.items():    #update to remainder R
            R[n-b+exp] = R.get(n-b+exp,0)-coeff*(R[n]/last_coeff_B)
        
        for key in sorted(R.keys()):
            if not R[key]:
                del R[key]  #removes zeroes
            if not len(R):
                R = {0:0}   #makes sure R is not empty

        n = max([key for key in R.keys()]+[0])  #degree of R gets an update

    return (Q,R)

def copy(d):
    dict = {k:x for (k,x) in d.items()}
    return dict

# a = Cyclo({7:3},12)
# b = Cyclo({3:2},12)
# print(a / b)
# print(b.Basis)
# print(b**(-1))

#end of functions used in Zeta

#START OF TESTS
###

# n = int(input('n = '))

# start = time.time()
# a = Zeta(1,n)
# print('time spent creating basis and powers = ' + str(time.time() - start))
# print(len(a.Basis[n]))

# for k in range(n):
#     if len(a.Powers[n][k]) > 1:
#         print(k,len(a.Powers[n][k]))
#print(a.get_polycyclo(84))