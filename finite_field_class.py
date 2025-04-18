class GaloisFp:

    def __init__(self,number,p,n=1):
        self.p = p
        self.value = number%(self.p)
        self.power = n
        #self.q = p**n
        if not isinstance(self.p,int) or not isinstance(self.power,int):
            raise TypeError('p or n are not integers')
        
    def __repr__(self):
        #how the number appears when called in command line
        return str(self.value) #+ ' (mod %d)' (self.p)
    
    def __str__(self):
        #how the number appears in print
        return str(self.value)
    
    def __bool__(self):
        #how the number is evaluated in logic statements - ifs
        return bool(self.value)
        
    def __add__(self,other):
        #addition over a finite field
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((self.value+other.value)%self.p,self.p)
            else:
                raise TypeError('Tried adding integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((self.value+other)%self.p,self.p)    #external addition with an integer
        else:
            raise TypeError("Unsupported operand type(s) for +")
        
    def __radd__(self,other):
        #addition to the right
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((other.value+self.value)%self.p,self.p)
            else:
                raise TypeError('Tried adding integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((other+self.value)%self.p,self.p)    #external addition with an integer
        else:
            raise TypeError("Unsupported operand type(s) for +")
    
    def __sub__(self,other):
        #substraction over a finite field
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((self.value-other.value)%self.p,self.p)
            else:
                raise TypeError('Tried substracting integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((self.value-other)%self.p,self.p)
        else:
            raise TypeError("Unsupported operand type(s) for -")
        
    def __rsub__(self,other):
        #substraction to the right
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((other.value-self.value)%self.p,self.p)
            else:
                raise TypeError('Tried substracting integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((other-self.value)%self.p,self.p)
        else:
            raise TypeError("Unsupported operand type(s) for -")

    def __mul__(self,other):
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((self.value*other.value)%self.p,self.p)
            else:
                raise TypeError('Tried multiplying integers over distinct fields')
        elif isinstance(other, int) or isinstance(other, float):
            return GaloisFp((self.value*int(other))%self.p,self.p)
        else:
            #WHY DO I NEED TO DO THAT ???
            try:
                return other.__rmul__(self)
            except:
                raise TypeError("Unsupported operand type(s) for *")

    def __rmul__(self,other):
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((other.value*self.value)%self.p,self.p)
            else:
                raise TypeError('Tried multiplying integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((other*self.value)%self.p,self.p)
        else:
            raise TypeError("Unsupported operand type(s) for *")

    def __truediv__(self,other):
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p and other.value != 0:
                return GaloisFp((self.value*inv(other.value,self.p))%self.p,self.p)
            else:
                raise TypeError('Tried dividing integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((self.value*inv(other,self.p))%self.p,self.p)
        else:
            raise TypeError("Unsupported operand type(s) for /")

    def __rtruediv__(self,other):
        if isinstance(other, GaloisFp) and self.power == 1:
            if self.p == other.p:
                return GaloisFp((inv(self.value,self.p)*other.value)%self.p,self.p)
            else:
                raise TypeError('Tried dividing integers over distinct fields')
        elif isinstance(other, int):
            return GaloisFp((inv(self.value,self.p)*other)%self.p,self.p)
        else:
            raise TypeError("Unsupported operand type(s) for /")

    def __mod__(self,other):
        pass

def xgcd(a,b):
    #returns gcd and bezout coefficients of a,b integers, in particular, x is the mult inverse of a mod b
    prevx, x = 1, 0; prevy, y = 0, 1
    while b:
        q = a//b
        x, prevx = prevx - q*x, x
        y, prevy = prevy - q*y, y
        a, b = b, a % b
    return a, prevx, prevy

def inv(a, m):
    #calls the xgcd function to return the modular inverse of a mod m
    g, x, y = xgcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
    
#x = GaloisFp(13,127)
#y =  GaloisFp(53,127)
#print((1/y)*3)
