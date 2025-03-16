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
        # /!\ THERES A PROBLEM WITH MUL /!\
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
        elif isinstance(other, int) or isinstance(other, float):
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
        elif isinstance(other, int) or isinstance(other, float):
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
        elif isinstance(other, int) or isinstance(other, float):
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
        elif isinstance(other, int) or isinstance(other, float):
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
        elif isinstance(other, int) or isinstance(other, float):
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
        elif isinstance(other, int) or isinstance(other, float):
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