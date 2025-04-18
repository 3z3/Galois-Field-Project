from finite_field_class import *
from rootofunity_class import *
import math
import random
import time

class Polynomial:

    def __init__(self,polynomial):
        #self.sequence is the dictionary form of the polynomial, with keys as powers and values as coefficients
        #self.deg is the degree of the polynomial, still 0 when P = 0 though
        if isinstance(polynomial,dict):
            self.sequence = polynomial  #expects only integer keys
        elif isinstance(polynomial,list):
            self.sequence = {k:polynomial[k] for k in range(len(polynomial))}   #converts a list to a dictionary
        elif isinstance(polynomial,str):
            self.sequence = self.interpret(polynomial)
        else:
            raise ValueError('what you entered cannot be interpreted as a polynomial')

        try:    #removes zeroes from the dictionary
            for key in sorted(self.sequence.keys()):
                if not bool(self.sequence[key]):
                    del self.sequence[key]
            if len(self.sequence) == 0:
                self.sequence = {0:0}
        except:
            self.sequence = {0:0}
            print('Given Polynomial was Empty!')
        self.deg = max([k for k in self.sequence.keys()]+[0])

    def __repr__(self):
        return str(self.sequence)
    
    def __str__(self):
        #how the polynomial is represented in print
        return self.display()
    
    def __bool__(self):
        #how the polynomial is evaluated in logic statements - ifs
        #returns False if poly is 0 and True if not
        return bool(not ((self.deg < 1) and (not bool(self.sequence.get(0,False)))))

    def __eq__(self,other):
        #equals ==
        try:
            return not bool(self - other)   #if P == Q, P-Q == 0, so P-Q is False, we want True, so we add not
        except:
            error('==',self,other)
    
    def __call__(self, x):
        #works with polynomials too, since add, mul, and pow are defined on Polynomial
        if isinstance(x, int) or isinstance(x, float) or isinstance(x, complex) or isinstance(x, Cyclo) or isinstance(x, GaloisFp) or isinstance(x, Polynomial):
            return self.evaluate(x)
        else:
            raise TypeError("%s object is not callable" % (type(x).__name__))
        
    def twovaluate(self, x):
        #weaker, larger-complexity, version of evaluate
        sum = 0
        for k,v in self.sequence.items():
            sum += v*(x**k)
        return sum

    def evaluate(self, x, keys = []):
        #evaluates the value of the polynomial at x with the horner recursive method
        #x might be a float, an int, or some other type that understands addition, multiplication and integer powers
        k = keys
        if len(k) == 0:
            k = [key for key in self.sequence.keys()]    #sorts keys coefficients in increasing order
            k.sort()
        
        if len(k) == 1:
            return self.sequence[k[0]]
        
        return (self.sequence[k[0]]+(x**(k[1]-k[0]))*self.evaluate(x,k[1:]))

    def __add__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) + other.sequence.get(key,0)
            return Polynomial(sum_poly)
        else:
            try:
                sum_poly = self.copy()
                sum_poly[0] = sum_poly.get(0,0) + other
                return Polynomial(sum_poly)
            except:
                raise TypeError("Unsupported operand type(s) for +: %s and %s" % (type(self).__name__,type(other).__name__))
            
    def __radd__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) + other.sequence.get(key,0)
            return Polynomial(sum_poly)
        else:
            try:
                sum_poly = self.copy()
                sum_poly[0] = other + sum_poly.get(0,0)
                return Polynomial(sum_poly)
            except:
                raise TypeError("Unsupported operand type(s) for +: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __sub__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) - other.sequence.get(key,0)
            return Polynomial(sum_poly)
        else:
            try:
                sum_poly = self.copy()
                sum_poly[0] = sum_poly.get(0,0) - other
                return Polynomial(sum_poly)
            except:
                raise TypeError("Unsupported operand type(s) for -: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __rsub__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) - other.sequence.get(key,0)
            return Polynomial(sum_poly)
        else:
            try:
                sum_poly = {exp : -coeff for exp,coeff in self.sequence.items()}
                sum_poly[0] = other + sum_poly.get(0,0)
                return Polynomial(sum_poly)
            except:
                print(self,other)
                raise TypeError("Unsupported operand type(s) for -: %s and %s" % (type(self).__name__,type(other).__name__))
    
    def __truediv__(self, other):
        if isinstance(other,Polynomial):
            raise TypeError('a polynomial isnt necessarily divisible by another')
        try:
            return Polynomial({k:x/other for (k,x) in self.sequence.items()})
        except:
            raise TypeError("Unsupported operand type(s) for /: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __mul__(self, other):
        if isinstance(other, Polynomial):   #poly multiplication
            return self.mult(other)
        else:
            try:
                return Polynomial({k:x*other for (k,x) in self.sequence.items()})
            except:
                stop(self,other)
                raise TypeError("Unsupported operand type(s) for *: %s and %s" % (type(self).__name__,type(other).__name__))

    def __rmul__(self, other):
        if isinstance(other, Polynomial):   #poly mult to the right
            return other.mult(self)
        else:
            try:
                return Polynomial({k:x*other for (k,x) in self.sequence.items()})
            except:
                print(self,other)
                raise TypeError("Unsupported operand type(s) for *: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __pow__(self, other):
        #integer powers of a polynomial
        if isinstance(other, int):
            if other != 0:
                res = Polynomial(self.copy())
                for i in range(1,other):
                    res *= self
                return res
            else:
                return Polynomial({0:1})
        else:
            raise TypeError("Unsupported operand type(s) for **: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __mod__(self, other):
        if isinstance(other, Polynomial):   #remainder of polynomial euclidean division
            return self.div(other)[1]
        elif isinstance(other, int):
            return Polynomial({k:x%other for (k,x) in self.sequence.items()})
        elif isinstance(other, GaloisFp):
            return Polynomial({k:x%other.value for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for %: %s and %s" % (type(self).__name__,type(other).__name__))
        
    def __floordiv__(self, other):
        if isinstance(other, Polynomial):   #quotient of polynomial euclidean division
            return self.div(other)[0]
        elif isinstance(other, int):
            return Polynomial({k:x//other for (k,x) in self.sequence.items()})
        elif isinstance(other, GaloisFp):
            return Polynomial({k:x//other.value for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for //: %s and %s" % (type(self).__name__,type(other).__name__))

    def copy(self):
        #copies a dictionary using another ref as to not ref the first dictionary when modifying the elements
        dict = {k:x for (k,x) in self.sequence.items()}
        return dict
    
    def listform(self, deg = None):
        #returns the polynomial in the form of a list
        if deg is None:
            deg = self.deg
        return [self.sequence.get(k,0) for k in range(deg+1)]

    def display(self):
        #returns an algebraic representation of the polynomial object as a string
        string = ''
        for key,coeff in sorted(self.sequence.items()):
            #coeff string
            if not key:
                coeff_string = str(coeff)   #case key == 0
                X_string = ''
            elif bool(coeff + 1) and bool(coeff -1):
                coeff_string = str(coeff)   #case key > 0 and abs(coeff) != 1
            elif not bool(coeff + 1):
                coeff_string = '-'
            else:
                coeff_string = ''   #case abs(coeff) == 1

            #X string
            if key > 1:
                X_string = 'X^%s' % key
            elif key == 1:
                X_string = 'X'

            #small tweaks, sign and xxx.0 floating points
            if coeff > 0:
                coeff_string = '+' + coeff_string
            if coeff_string[-2:] == '.0':
                coeff_string = coeff_string[:-2]

            #uncomment if you want the display to be less bulky
            
            if isinstance(coeff,float) or isinstance(coeff,Ratio):
                coeff_string = coeff_string[:coeff_string.index('.')+3] #2 decimal places

            #updating the string
            string += coeff_string + X_string
        
        if string[0] == '+':
            string = string[1:]

        return string
    
    def interpret(self,string):
        #interprets a string of polynomial terms as a dictionary corresponding to the wanted polynomial /!\ WORKS
        def convert(x_str):
            #converts a string to a number (integer or float)
            x = float(x_str)
            if int(x) == x:
                x = int(x)
            return x
        
        pows, coeffs, poly = [], [], {}
        splitplus = string.split('+')
        #meant to split the polynomial into + and - then adapt (using sign) the sign of the coefficient
        for plusdivided in splitplus:
            sign = -1
            for term in plusdivided.split('-'):
                if 'X' in term: #this means pow > 0
                    splitted = term.split('X')  #splits the term in coefficient , power along 'X'
                    splitted_pow = splitted[1][int('^' in term):]   #begins at 1 if True, 0 if False, skips ^ either way
                    splitted_coeff = splitted[0]
                    if len(splitted_pow):
                        pows.append(convert(splitted_pow))
                    else:
                        pows.append(1)  #if empty, that means power is 1
                    if len(splitted_coeff): 
                        coeffs.append(-sign * convert(splitted_coeff))
                    else:
                        coeffs.append(-sign)    #if empty that means coefficient is 1 or -1 depending on the sign
                else:   #this means pow == 0
                    pows.append(0)
                    coeffs.append(-sign * convert(term))
                sign = abs(sign)

        for pow,coeff in zip(pows,coeffs):  #updates the polynomial (dictionary) with all the info collected
            poly[pow] = poly.get(pow,0) + coeff
        return poly
    
    def mult(self,Q):
        #returns a polynomial equal to the product PxQ with coefficients until the sum of their respective degrees
        #used in the mult magic method
        Pi = []
        for j in range(self.deg+Q.deg+1):
            coeff = 0
            i = max(0,j-Q.deg)
            while i <= min(j,self.deg):   #making sure the index stays within each of the polynomials' indices
                coeff += self.sequence.get(i,0)*Q.sequence.get(j-i,0)
                i += 1
            Pi.append(coeff)
        return Polynomial(Pi)
    
    def div(self,B):
        #returns quotient and remainder of the euclidean division of A by B, should correspond to A = BQ + R, here, A = self
        b = B.deg
        if B.deg == 0:
            const = B.sequence[0]
            return (self / const, 0)
        R, Q, n = self.copy(), {}, self.deg
        last_coeff_B = B.sequence[b]
        
        while n >= b:

            last_coeff_R = R[n]

            Q[n-b] = last_coeff_R/last_coeff_B  #update to quotient Q

            for exp, coeff in B.sequence.items():    #update to remainder R
                R[n-b+exp] = R.get(n-b+exp,0)-coeff*(last_coeff_R/last_coeff_B)
            
            for key in sorted(R.keys()):
                if not R[key]:
                    del R[key]  #removes zeroes
                if not len(R):
                    R = {0:0}   #makes sure R is not empty

            n = max([key for key in R.keys()]+[0])  #degree of R gets an update

        return (Polynomial(Q),Polynomial(R))
    
    def Qdiv(self,B):
        #returns quotient and remainder of the euclidean division of A by B, should correspond to A = BQ + R, here, A = self
        #adapted to the Ratio class in order to divide in Q exactly
        if B.deg == 0:
            const = B.sequence[0]
            return (self / const, 0)
        R, Q = Polynomial(self.copy()), 0
        b, n = B.deg, self.deg
        last_coeff_B = B.sequence[b]
        
        while n >= b:
            delta = n-b
            last_coeff_R = R.sequence[n]
            ratio_check = bool(isinstance(last_coeff_R,Ratio) or isinstance(last_coeff_B,Ratio))
            #do we really need ratio check ?
            if ratio_check:
                factor = Polynomial({delta:last_coeff_R/last_coeff_B})
                Q += factor
                R -= B * factor
            else:
                factor = Polynomial({delta:Ratio(last_coeff_R,last_coeff_B)})
                Q += factor
                R -= B * factor

            n = R.deg
        
        return (Q,R)
    
    def derivative(self):
        #outputs the algebraic derivative of the polynomial
        pprime = {}
        for k,v in self.sequence.items() and k != 0:
            pprime[k-1] = k*v
        
        if self.deg == 0:
            pprime = {0:0}
        
        return Polynomial(pprime)
    
    def Bezout(A,B):
        #returns Bezout coefficients and remainder of polynomials A and B
        #sketchy, dont look at it too much
        R, Rplus = A, B
        U, Uplus = 1, 0
        V, Vplus = 0, 1

        while Rplus:
            #making R and Rplus have 1 leading coefficients
            leadingR, leadingRplus = R.sequence[R.deg], Rplus.sequence[Rplus.deg]
            check_Rplus = isinstance(leadingRplus,Ratio)
            #check_R = isinstance(leadingR,Ratio)

            # if check_R:
            #     R,U,V = R * (1/leadingR), U * (1/leadingR), V * (1/leadingR)
            # else:
            #     R,U,V = R * Ratio(1,leadingR), U * Ratio(1,leadingR), V * Ratio(1,leadingR)
            if check_Rplus:
                Rplus,Uplus,Vplus = Rplus * (1/leadingRplus), Uplus * (1/leadingRplus), Vplus * (1/leadingRplus)
            else:
                Rplus,Uplus,Vplus = Rplus * Ratio(1,leadingRplus), Uplus * Ratio(1,leadingRplus), Vplus * Ratio(1,leadingRplus)

            print('this is a comparison of AU + BV and R (should be equal) :')
            print(A*U + B*V, R, A*U + B*V == R)
            print(' and this is a comparison of AUplus + BVplus and Rplus (should be equal) :')
            stop(A*Uplus + B*Vplus, Rplus, A*Uplus + B*Vplus == Rplus)

            quotient, remainder = R.Qdiv(Rplus)
            R, Rplus = Rplus, remainder
            U, Uplus = Uplus, U - (quotient * Uplus)
            V, Vplus = Vplus, V - (quotient * Vplus)

        #/!\ this commented part makes sure the resulting polynomials are Z-valued
        #
        # denominators = []
        # prod = 1
        # #/!\ sometimes U or V arent polynomials
        # if isinstance(U,Polynomial):
        #     Udeno = U.sequence.values()
        # else:
        #     Udeno = U
        # if isinstance(V,Polynomial):
        #     Vdeno = V.sequence.values()
        # else:
        #     Vdeno = V
        # values = set(Udeno) | set(Vdeno)    #union of sets
        # for value in values:
        #     if isinstance(value, Ratio):
        #         q = value.denominator
        #         if q != 1 and q != -1:
        #             denominators.append(q)
        #             prod *= q
        #     elif isinstance(value, float):
        #         stop(value)

        # if len(denominators):
        #     lcmext = lcm_extended(*denominators)
        #     U *= lcmext; V *= lcmext; R *= lcmext

        return (U,V,R)    #in the following order : 1st and 2nd Bezout coefficients of A and B, then gcd of A and B (ie AU + BV = R)
    
#end of polynomial class
    
def ComplexFFT(coeff, zeta):
    #fast fourier transform algorithm
    #coeff must be of size a power of 2
    #zeta is a root of unity
    #CORRECT FFT ALGO /!\
    n = len(coeff)
    if n == 1:
        #end of recursion
        res = [coeff[0]]
        return res
    else:
        #this is where the stuff happens
        midpoint = n//2
        values = [0]*n

        zeta_squared = zeta*zeta
        evens = ComplexFFT([coeff[2*k] for k in range(midpoint)], zeta_squared) #divide & conquer process, from 0 up to N/2-1
        odds = ComplexFFT([coeff[2*k+1] for k in range(midpoint)], zeta_squared)  #recursive process

        zetak = 1
        for k in range(midpoint):
            #making the most out of half the info; N evaluations of the poly P are obtained with only N/2 eval of poly E and O (evens and odds)
            values[k] = evens[k] + (zetak * odds[k])  #P(w^k) = A(w^2k) + (w^k)*B(w^2k), omega multiplication being the Rshift method
            values[k+midpoint] = evens[k] - (zetak * odds[k])
            zetak *= zeta

        return values
    
def FFTMult(A,B):
    #multiplication of two polynomials with the FFT algorithm
    #does not work with coefficients up to 5000 coupled with degree 50000 polynomials
    m = max(A.deg,B.deg)

    #setting up the FFT degree of C = A x B
    pow = len(bin(2*m))-2
    n = 2**pow

    #the part is for computing the infinity-norm of A and B
    normA, normB = 0, 0
    for a in A.sequence.values():
        aplus = abs(a)
        if aplus > normA:
            normA = aplus
    for b in B.sequence.values():
        bplus = abs(b)
        if bplus > normB:
            normB = bplus

    #calculating constants in the function f(x) = x(logx + beta) - gamma -> needed to approximate alpha
    gamma = math.log((n/2)*normA*normB) + 3 * math.log(n) + math.log(2) #maybe dont need to compute the log in its entirety ?
    beta = math.log(n/(2*math.pi)) - 1

    def newtonalpha(x): #newton's method for f
        res = (x + gamma) / (math.log(x) + beta + 1)
        return res
    
    a = 2 * gamma
    aplus = newtonalpha(a)
    delta = abs(a - aplus)
    
    while delta > 1:    #repeating newton's method until the difference between consecutive terms is less than 1
        a, aplus = aplus, newtonalpha(aplus)
        delta = abs(a - aplus)
        
    order = math.floor(aplus)    #alpha = order must be an integer, so taking the next integer after root of f that still verifies the boundary, and substracting 1 (ceil - 1 = floor)

    zeta = exp(2*math.pi*complex('j')/n, order-1)
    zeta_minus = exp(-2*math.pi*complex('j')/n, order-1)

    newA, newB = A.listform() + [0]*(n - A.deg - 1), B.listform() + [0]*(n - B.deg - 1) #so the length of coefficients fit the input

    Aeval, Beval = ComplexFFT(newA, zeta), ComplexFFT(newB, zeta)   #FFT takes lists as input, not dictionaries let alone polys
    Ceval = [x*y for x, y in zip(Aeval, Beval)]

    C = ComplexFFT(Ceval, zeta_minus)
    C = [c/n for c in C]
    for k in range(len(C)):
        #CAREFUL /!\ THIS ROUNDS COEFFICIENT OF THE PRODUCT TO INTEGERS !!! DO NOT MULTIPLY FLOAT POLYNOMIALS WITH THAT
        try:
            x = round(C[k].real)
        except:
            print(C[k])
            print(str(input('continue ? ')))
        y = round(C[k].imag)
        if y == 0:
            new_z = x
        else:
            new_z = complex(x,y)
        C[k] = new_z
    return Polynomial(C)

def exp(x,order):
    #taylor series of exponential evaluted at x of order 'order'
    #can take complex numbers as input
    #uses horner method of evaluation
    k = order
    res = 1 + (x/k)
    while k > 1:
        k -= 1
        res = 1 + (x/k)*res
    return res

#random tests to evaluate the speed of the FFT multiplication algorithm

A,B = Polynomial([random.randint(-50,50) for x in range(500000)]+[1]), Polynomial([random.randint(-50,50) for x in range(500000)]+[1])
start1 = time.time()
#print(FFTMult(A,B))
C = FFTMult(A,B)
print('computing A x B took ' + str(time.time()-start1) + ' seconds')

#end of tests

#start of order q Galois Field class

class GaloisFq(GaloisFp):

    PolyExists = {}
    split = {}
    big = {}

    def __init__(self, number, p, n=1): #TODO add an irreducible polynomial as variable input
        super().__init__(number, p, n)
        self.q = p**n

        #TODO set the polynomial to be the splitting poly of the field F_q (in PolyExists dictionary)

        if self.q not in self.PolyExists:
            self.PolyExists[self.q] = False
            self.big[self.q] = Polynomial({1:GaloisFp(-1,self.p),self.q:GaloisFp(1,self.p)})
        #self.powers.append(self.power)
        #self.orders.append(self.q)
        
        tested = 0
        while not self.PolyExists[self.q]:
            tested += 1
            split_test = []
            split_test.append(GaloisFp(random.randint(1,self.p-1),self.p))
            i = 1
            while i < n:
                split_test.append(GaloisFp(random.randint(0,self.p-1),self.p))
                i += 1
            split_test.append(1)
            split_test = Polynomial(split_test)
            if not self.big[self.q].div(split_test)[1]:
                self.PolyExists[self.q] = True
                self.split[self.q] = split_test

        print('tested ' + str(tested) + ' times')

        #TODO after splitting poly has been eventually set, write number in the form of a polynomial with variable a root of the splitting poly
        #for example if i is a root of P = X^2+1 in Z/pZ (which is irreducible), then, among others, 2+3i is an element of F_49

    def split_show(self):
        print(self.split[self.q])

    def display(self):
        pass

    #TODO write arithmetic rules in the field F_q (with magic methods)

#end of order q Galois field class
