def error(operation,*variables):
    #easy way to write type errors
    #operation must be a string
    endstring = ''
    for var in variables[:-1]:
        endstring += type(var).__name__ + ' and '
    endstring += type(variables[-1]).__name__
    mainstring = 'Unsupported operand type(s) for ' + operation + ': ' + endstring
    raise TypeError(mainstring)

def stop(*variables):
    #function that prints a bunch of variables then stops while printing a message
    print(*variables)
    a = str(input('do you want to continue ? '))
    pass

class Ratio:
    #class representing the rational numbers as couples of integers with a nonzero denominator

    def __init__(self,p,q):
        if isinstance(p,int) and isinstance(q,int) and q:
            a,b = p,q
            while b:
                a, b = b, a%b
            if q < 0:
                p *= -1
                q *= -1
            if abs(a) == 1:
                self.numerator = p
                self.denominator = q
            else:
                self.numerator = int(p/a)
                self.denominator = int(q/a)
            self.float = self.numerator/self.denominator
        elif isinstance(p,float) and isinstance(q,float) and int(p) == p and int(q) == q:
            temp = Ratio(int(p),int(q))
            self.numerator = temp.numerator
            self.denominator = temp.denominator
            self.float = temp.float
        elif isinstance(q,int) and isinstance(p,float) and int(p) == p:
            temp = Ratio(int(p),q)
            self.numerator = temp.numerator
            self.denominator = temp.denominator
            self.float = temp.float
        elif isinstance(p,int) and isinstance(q,float) and int(q) == q:
            temp = Ratio(p,int(q))
            self.numerator = temp.numerator
            self.denominator = temp.denominator
            self.float = temp.float
        elif isinstance(p,int) and isinstance(q,int) and not q:
            raise ValueError('tried to divide by zero')
        else:
            raise TypeError('numerator or denominator are not integers')
        
    def __repr__(self):
        return str(self.float)
    
    def __str__(self):
        return str(self.float)
    
    def __bool__(self):
        return bool(self.numerator)
    
    def __neg__(self):
        #fraction times minus one
        return Ratio(-self.numerator,self.denominator)
    
    def __lt__(self,other):
        #less than <
        if isinstance(other,int) or isinstance(other,float):
            return bool(self.numerator < other * self.denominator)
        elif isinstance(other,Ratio):
            return bool(self.numerator * other.denominator < other.numerator * self.denominator)
        else:
            error('<',self,other)

    def __gt__(self,other):
        #greater than >
        if isinstance(other,int) or isinstance(other,float):
            return bool(self.numerator > other * self.denominator)
        elif isinstance(other,Ratio):
            return bool(self.numerator * other.denominator > other.numerator * self.denominator)
        else:
            error('>',self,other)

    def __add__(self,other):
        if isinstance(other,Ratio):
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            if b == d:
                return Ratio(a + c, b)
            elif b%d == 0:
                return Ratio(a + c * (b//d), b)
            elif d%b == 0:
                return Ratio(a * (d//b) + c, d)
            else:
                return Ratio(a * d + b * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(a + other * b, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a + int(other) * b, b)
        else:
            try:
                return other.__radd__(self)
            except:
                error('+',self,other)

    def __radd__(self,other):
        if isinstance(other,Ratio):
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            if b == d:
                return Ratio(a + c, b)
            elif b%d == 0:
                return Ratio(a + c * (b//d), b)
            elif d%b == 0:
                return Ratio(a * (d//b) + c, d)
            else:
                return Ratio(a * d + b * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(a + other * b, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a + int(other) * b, b)
        else:
            error('+',self,other)

    def __sub__(self,other):    
        if isinstance(other,Ratio):
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            if b == d:
                return Ratio(a - c, b)
            elif b%d == 0:
                return Ratio(a - c * (b//d), b)
            elif d%b == 0:
                return Ratio(a * (d//b) - c, d)
            else:
                return Ratio(a * d - b * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(a - other * b, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a - int(other) * b, b)
        else:
            try:
                return other.__rsub__(self)
            except:
                error('-',self,other)

    def __rsub__(self,other):    
        if isinstance(other,Ratio):
            c,d,a,b = self.numerator,self.denominator,other.numerator,other.denominator
            if b == d:
                return Ratio(a - c, b)
            elif b%d == 0:
                return Ratio(a - c * (b//d), b)
            elif d%b == 0:
                return Ratio(a * (d//b) - c, d)
            else:
                return Ratio(a * d - b * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(other * b - a, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(int(other) * b - a, b)
        else:
            error('-',self,other)

    def __mul__(self,other):    
        if isinstance(other,Ratio):
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            return Ratio(a * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(a * other, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a * int(other), b)
        else:
            try:
                return other.__rmul__(self)
            except:
                error('*',self,other)

    def __rmul__(self,other):    
        if isinstance(other,Ratio):
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            return Ratio(a * c, b * d)
        elif isinstance(other,int):
            a,b = self.numerator,self.denominator
            return Ratio(a * other, b)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a * int(other), b)
        else:
            error('*',self,other)

    def __truediv__(self,other):
        check_div = bool(other)
        if isinstance(other,Ratio) and check_div:
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            return Ratio(a * d, b * c)
        elif isinstance(other,int) and check_div:
            a,b = self.numerator,self.denominator
            return Ratio(a, b * other)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(a, b * int(other))
        else:
            try:
                return other.__rtruediv__(self)
            except:
                error('/',self,other)

    def __rtruediv__(self,other):
        check_div = bool(self)
        if isinstance(other,Ratio) and check_div:
            a,b,c,d = self.numerator,self.denominator,other.numerator,other.denominator
            return Ratio(c * b, d * a)
        elif isinstance(other,int) and check_div:
            a,b = self.numerator,self.denominator
            return Ratio(other * b, a)
        elif isinstance(other,float) and int(other) == other:
            a,b = self.numerator,self.denominator
            return Ratio(b * int(other), a)
        elif not check_div:
            raise ZeroDivisionError('tried to divide by zero')
        else:
            error('/',self,other)

    def __pow__(self,other):
        if isinstance(other,int):
            a,b = self.numerator, self.denominator
            return Ratio(a**other, b**other)
        if isinstance(other,float) and int(other) == other:
            a,b = self.numerator, self.denominator
            return Ratio(a**other, b**other)
        else:
            error('**',self,other)
