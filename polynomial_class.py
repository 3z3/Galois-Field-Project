from finite_field_class import *
from rootofunity_class import *
import random
import time

class Polynomial:

    def __init__(self,polynomial,modulo=None):
        if isinstance(polynomial,list):
            self.sequence = {k:polynomial[k] for k in range(len(polynomial))}   #converts a list to a dictionary
        else:
            self.sequence = polynomial  #expects a dictionary
        if not (modulo is None):
            try:
                self.sequence = {k:x%modulo for (k,x) in self.sequence.items()}
            except:
                if isinstance(modulo,int):
                    raise TypeError('Given Polynomial has non-number(s) coefficient(s)')
                else:
                    raise TypeError('Modulo is not an int')
        try:    #removes zeroes from the dictionary
            for key in sorted(self.sequence.keys()):
                if not self.sequence[key]:
                    del self.sequence[key]
            if len(self.sequence) == 0:
                self.sequence = {0:0}
        except:
            self.sequence = {0:0}
            print('Given Polynomial was Empty!')
        self.deg = max([k for k in self.sequence.keys()]+[0])
        #if self.deg == [0]:
        #    self.deg = -1 #minus infinity

    def __repr__(self):
        return str(self.sequence)
    
    def __str__(self):
        #how the polynomial is represented in print
        return self.display()
    
    def __bool__(self):
        #how the polynomial is evaluated in logic statements - ifs
        #returns False if poly is 0 and True if not
        return bool(not ((self.deg < 1) and (not self.sequence.get(0,False))))

    def __add__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) + other.sequence.get(key,0)
            return Polynomial(sum_poly)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):
            sum_poly = self.copy()
            sum_poly[0] += other
            return Polynomial()
        else:
            raise TypeError("Unsupported operand type(s) for +")
        
    def __sub__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = self.copy()
            for key in other.sequence.keys():
                sum_poly[key] = sum_poly.get(key,0) - other.sequence.get(key,0)
            return Polynomial(sum_poly)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):
            sum_poly = self.copy()
            sum_poly[0] -= other
            return Polynomial()
        else:
            raise TypeError("Unsupported operand type(s) for -")
    
    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):   #division of a polynomial by a scalar
            return Polynomial({k:x/other for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for /")
        
    def __mul__(self, other):
        if isinstance(other, Polynomial):   #poly multiplication
            return self.mult(other)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):    #external scalar multiplication
            return Polynomial({k:x*other for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for *")

    def __rmul__(self, other):
        if isinstance(other, Polynomial):   #poly mult to the right
            return other.mult(self)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):    #external scalar multiplication to the right
            return Polynomial({k:other*x for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for *")
        
    def __mod__(self, other):
        if isinstance(other, Polynomial):   #remainder of polynomial euclidean division
            return self.div(other)[1]
        elif isinstance(other, int):
            return Polynomial({k:x%other for (k,x) in self.sequence.items()})
        elif isinstance(other, GaloisFp):
            return Polynomial({k:x%other.value for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for %")
        
    def __floordiv__(self, other):
        if isinstance(other, Polynomial):   #quotient of polynomial euclidean division
            return self.div(other)[0]
        elif isinstance(other, int):
            return Polynomial({k:x//other for (k,x) in self.sequence.items()})
        elif isinstance(other, GaloisFp):
            return Polynomial({k:x//other.value for (k,x) in self.sequence.items()})
        else:
            raise TypeError("Unsupported operand type(s) for //")

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
        poly_str = ''
        for key,coeff in self.sequence.items():
            term = ''
            if coeff:
                if len(poly_str) != 0:
                    term += '+'
                if key != 0:
                    if (coeff - 1):
                        if key != 1:
                            term += str(coeff) + 'X^%d' % (key)
                        else:
                            term += str(coeff) + 'X'
                    else:
                        if key != 1:
                            term += 'X^%d' % (key)
                        else:
                            term += 'X'
                else:
                    term += str(coeff)
            poly_str += term
        
        #gets rid of +-
        k=0
        while k < len(poly_str):
            if poly_str[k] == '-' and poly_str[k-1] == '+':
                poly_str = poly_str[:k-1] + poly_str[k:]
            else:
                k += 1

        if len(poly_str) == 0:
            poly_str += '0'

        return poly_str
    
    def mult(self,Q):
        #returns a polynomial equal to the product PxQ with coefficients until the sum of their respective degrees
        #used in the mult magic method

        #TODO NEEDS FFT
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
        try:
            R,Q,n = self.copy(),{},self.deg
        except:
            print('the degree of A is less than the degree of B !')

        while n >= B.deg:
            Q[n-B.deg] = Q.get(n-B.deg,0) + R[n]/B.sequence[B.deg]
            for i in range(B.deg+1):
                R[n-B.deg+i] = R.get(n-B.deg+i,0)-B.sequence.get(i,0)*(R[n]/B.sequence[B.deg])
            try:
                for key in sorted(R.keys()):
                    if not R[key]:
                        del R[key]
                if len(R) == 0:
                    R = {0:0}
            except:
                print(R)
                print('hello it\'s me')
            n = max([key for key in R.keys()]+[0])  #degree of R gets an update
        
        return (Polynomial(Q),Polynomial(R))
    
#end of polynomial class

def FFT(coeff, shift=1, order=1, not_power=True):
    #fast fourier transform algorithm
    #shift must be a nonzero integer, shift = power to which w is raised
    #order is the order of w (as a root of unity)
    n = len(coeff)
    if n == 1:
        #end of recursion
        if isinstance(coeff[0], Cyclo):
            res = [coeff[0]]
        else:    
            res = [Cyclo([coeff[0]] + [0]*(int(order/2)-1), order)]
        return res
    
    else:
        #this is where the stuff happens
        N, new_coeff = int(order/abs(shift)), coeff

        if not_power:
            #initializing which root of unity w is (rounding up to a power of 2)
            pow = len(bin(n))-2
            N = 2**pow
            if 2*n == N:
                N = n
                order = n
            else:
                new_coeff = coeff + [0]*(N-n)
                order = N

        if N != n:
            new_coeff = coeff + [0]*(N-n)   #adding padding if coefficient list is too short for the algorithm
        midpoint = int(N/2)
        values = [0]*N

        evens = FFT([new_coeff[2*k] for k in range(midpoint)], shift*2, order, False) #divide & conquer process, from 0 up to N/2-1
        odds = FFT([new_coeff[2*k+1] for k in range(midpoint)], shift*2, order, False)  #recursive process

        alpha_shift = 0
        for k in range(midpoint):
            #making the most out of half the info; N evaluations of the poly P are obtained with only N/2 eval of poly E and O (evens and odds)
            values[k] = evens[k] + (odds[k].Rshift(alpha_shift))  #P(w^k) = A(w^2k) + (w^k)*B(w^2k), omega multiplication being the Rshift method
            values[k+midpoint] = evens[k] - (odds[k].Rshift(alpha_shift))
            alpha_shift += shift

        return values
    
def FourierMult(A,B):
    #multiplication of two polynomials with the FFT algorithm
    m = max(A.deg,B.deg)

    #setting up the FFT degree of C = A x B
    pow = len(bin(2*m))-2
    n = 2**pow

    Aeval, Beval = FFT(A.listform(), 1, n, False), FFT(B.listform(), 1, n, False)   #FFT takes lists as input, not dictionaries let alone polys
    Ceval = [x*y for x,y in zip(Aeval,Beval)]   #multiply evaluations, eval by eval: C(1) = A(1)xB(1), C(w) = A(w)xB(w), etc

    coeffw = FFT(Ceval, -1, n, False)   #apply inverse FFT to evaluations of C at points 1, w, w^2, ... , w^{n-1}
    return Polynomial([x.firstcoord()/n for x in coeffw])


a = FFT([1,3,1,1,-1,1,2,1])

print(a)
print('this is the fourier transform of your polynomial with coordinates in Z[w]')

b = FFT(a, -1)

print(b)
print('and this is the inverse fourier transform of the fourier transform')


A,B = Polynomial([1,1,1]), Polynomial([1,2,3])
print(FourierMult(A,B))
print('this is the result of A x B')

C = Cyclo([-2,-3,6,-5], 8)
print(C.Rshift(0))

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

#random polynomial test zone /!\
'''
start = time.time()
for i in range(10000):
    deg1 = random.randint(1,50)
    deg2 = random.randint(1,50)
    P,Q,R = {},{},Polynomial({0:0})
    for k in range(deg1+1):
        P[k] = random.randint(0,100)
    for k in range(deg2+1):
        Q[k] = random.randint(0,100)
    P = Polynomial(P)
    Q = Polynomial(Q)
    R = P+Q
end = time.time()
print(end-start)
'''

'''
start = time.time()
P = Polynomial([-1,-2,1,1,0,1])

Q = Polynomial([-1,0,1])

R = Polynomial([1,2,0,1])

print(P.div(Q))

print(P*Q)

print(P-Q*R)    #supposed to be zero

x = GaloisFp(19,23)

y = GaloisFp(12,23)

print(P*x+R*y)

#print((1/x)*R) ------> THIS DOESNT WORK

print((1/x)*P)

end = time.time()
print(end-start)
'''


'''
p,n = int(input('give me p: ')), int(input('give me n: '))

start = time.time()

x = GaloisFq(2,p,n)

y = GaloisFq(5,p,n)

x.split_show()

y.split_show()

end = time.time()
print('time taken = ' + str(end-start))
'''