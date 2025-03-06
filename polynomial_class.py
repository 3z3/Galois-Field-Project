from finite_field_class import *
import random
import time

class Polynomial:

    def __init__(self,polynomial,modulo=None):
        self.sequence = polynomial
        if not (modulo is None):
            try:
                self.sequence = [x%modulo for x in self.sequence]
            except:
                if isinstance(modulo,int):
                    raise TypeError('Given Polynomial has non-number(s) coefficient(s)')
                else:
                    raise TypeError('Modulo is not an int')
        try:
            while self.sequence[-1] == 0 and len(self.sequence) > 1:
                self.sequence.pop(-1)
        except:
            self.sequence = [0]
            print('Given Polynomial was Empty!')
        self.deg = len(self.sequence)-1
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
        return bool(not ((self.deg < 1) and (0 in self.sequence)))

    def __add__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = [x + y for x,y in zip(self.sequence+[0]*max(0,other.deg-self.deg),other.sequence+[0]*max(0,self.deg-other.deg))] #terms are zero if beyond each other's degree
            return Polynomial(sum_poly)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):
            return Polynomial([self.sequence[0]+other]+self.sequence[1:])
        else:
            raise TypeError("Unsupported operand type(s) for +")
        
    def __sub__(self, other):
        if isinstance(other, Polynomial):
            sum_poly = [x - y for x,y in zip(self.sequence+[0]*max(0,other.deg-self.deg),other.sequence+[0]*max(0,self.deg-other.deg))] #terms are zero if beyond each other's degree
            return Polynomial(sum_poly)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):
            return Polynomial([self.sequence[0]-other]+self.sequence[1:])
        else:
            raise TypeError("Unsupported operand type(s) for -")
    
    def __truediv__(self, other):
        if isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):   #division of a polynomial by a scalar
            return Polynomial([x/other for x in self.sequence])
        else:
            raise TypeError("Unsupported operand type(s) for /")
        
    def __mul__(self, other):
        if isinstance(other, Polynomial):   #poly multiplication
            return self.mult(other)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):    #external scalar multiplication
            return Polynomial([x*other for x in self.sequence])
        else:
            raise TypeError("Unsupported operand type(s) for *")

    def __rmul__(self, other):
        if isinstance(other, Polynomial):   #poly mult to the right
            return other.mult(self)
        elif isinstance(other, int) or isinstance(other, float) or isinstance(other, GaloisFp):    #external scalar multiplication to the right
            return Polynomial([other*x for x in self.sequence])
        else:
            raise TypeError("Unsupported operand type(s) for *")
        
    def __mod__(self, other):
        if isinstance(other, Polynomial):   #remainder of polynomial euclidean division
            return self.div(other)[1]
        elif isinstance(other, int):
            return Polynomial([x%other for x in self.sequence])
        elif isinstance(other, GaloisFp):
            return Polynomial([x%other.value for x in self.sequence])
        else:
            raise TypeError("Unsupported operand type(s) for %")
        
    def __floordiv__(self, other):
        if isinstance(other, Polynomial):   #quotient of polynomial euclidean division
            return self.div(other)[0]
        elif isinstance(other, int):
            return Polynomial([x//other for x in self.sequence])
        elif isinstance(other, GaloisFp):
            return Polynomial([x//other.value for x in self.sequence])
        else:
            raise TypeError("Unsupported operand type(s) for //")

    def copy(self):
        #copies a list using another ref as to not ref the first list when modifying the elements
        l = []
        for i in self.sequence:
            l.append(i)
        return l

    def display(self):
        #returns an algebraic representation of the polynomial object as a string
        poly_str = ''
        i = 0
        for coeff in self.sequence:
            term = ''
            if coeff:
                if len(poly_str) != 0:
                    term += '+'
                if i != 0:
                    if (coeff - 1):
                        if i != 1:
                            term += str(coeff) + 'X^%d' % (i)
                        else:
                            term += str(coeff) + 'X'
                    else:
                        if i != 1:
                            term += 'X^%d' % (i)
                        else:
                            term += 'X'
                else:
                    term += str(coeff)
            poly_str += term
            i += 1
        
        #gets rid of +-
        k=0
        while k < len(poly_str):
            if poly_str[k] == '-' and poly_str[k-1] == '+':
                poly_str = poly_str[:k-1] + poly_str[k:]
            else:
                k += 1
        '''
        k = 0
        while k < len(poly_str):
            if poly_str[k] == '0':
                if k == 0:
                    poly_str = poly_str[1:]
                    if poly_str[0] == 'X':
                        poly_str = poly_str[1:]
                        if poly_str[0] == '^':
                            poly_str = poly_str[3:]
                        else:
                            poly_str = poly_str[1:]
                elif poly_str[k-1] not in [1,2,3,4,5,6,7,8,9]:
                    poly_str = poly_str[:k] + poly_str[k+1:]
                    if poly_str[k] == 'X':
                        poly_str = poly_str[:k] + poly_str[k+1:]
                        if poly_str[k] == '^':
                            poly_str = poly_str[:k] + poly_str[k+3:]
                        else:
                            poly_str = poly_str[:k] + poly_str[k+1:]
                else:
                    k += 1
            elif poly_str[k] == '-' and poly_str[k-1] == '+':
                poly_str = poly_str[:k-1] + poly_str[k:]
            else:
                k += 1
        '''

        if len(poly_str) == 0:
            poly_str += '0'

        return poly_str
    
    def mult(self,Q):
        #returns a polynomial equal to the product PxQ with coefficients until the sum of their respective degrees
        #used in the mult magic method
        Pi = []
        for j in range(self.deg+Q.deg+1):
            coeff = 0
            i = max(0,j-Q.deg)
            while i <= min(j,self.deg):   #making sure the index stays within each of the polynomials' indices
                coeff += self.sequence[i]*Q.sequence[j-i]
                i += 1
            Pi.append(coeff)
        return Polynomial(Pi)
    
    def div(self,B):
        #returns quotient and remainder of the euclidean division of A by B, should correspond to A = BQ + R, here, A = self
        try:
            R,Q,n = self.copy(),[0]*(self.deg-B.deg+1),self.deg
        except:
            print('the degree of A is less than the degree of B !')

        while n >= B.deg:
            Q[n-B.deg] = Q[n-B.deg]+(R[-1]/B.sequence[-1])
            #assert isinstance(B,list), 'polynomial B is ' + str(B)
            for i in range(B.deg+1):
                R[n-B.deg+i] = R[n-B.deg+i]-B.sequence[i]*(R[-1]/B.sequence[-1])
            try:
                while R[-1] == 0 and R != [0]:
                    R.pop(-1)
            except:
                print(R)
                #print('hello it\'s me')
            n = len(R)-1
        
        return (Polynomial(Q),Polynomial(R))
    
#end of polynomial class

class GaloisFq(GaloisFp):

    PolyExists = {}
    split = {}
    big = {}

    def __init__(self, number, p, n=1):
        super().__init__(number, p, n)
        self.q = p**n
        if self.p not in self.PolyExists:
            self.PolyExists[self.p] = False
            self.big[self.p] = Polynomial([0,-1]+([0]*(-2+self.q))+[1])
        #self.powers.append(self.power)
        #self.orders.append(self.q)
        
        while not self.PolyExists[self.p]:
            split_test = []
            split_test.append(random.randint(1,self.p-1))
            i = 1
            while i < n:
                split_test.append(random.randint(0,self.p-1))
                i += 1
            split_test.append(1)
            split_test = Polynomial(split_test)
            if not self.big[self.p].div(split_test)[1]:
                self.PolyExists[self.p] = True
                self.split[self.p] = split_test

    def split_show(self):
        print(self.split[self.p])

    print('hello')

#end of order q Galois field class
        
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

x = GaloisFq(2,23,5)

y = GaloisFq(5,23,5)

x.split_show()

y.split_show()