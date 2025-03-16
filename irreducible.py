from polynomial_class import *

q = { 1: [[1]] }

def decompose(n):
    #executes quickly approx up to 60
    try:
        return q[n]
        #for part in q[n]:
        #    yield part
    except:
        pass

    result = [[n]]
    #yield [n]

    for i in range(1, n):
        a = n-i
        R = decompose(i)
        for r in R:
            if r[0] <= a:
                result.append([a] + r)
                #yield ([a] + r)

    q[n] = result
    #result.remove([n])
    return result

'''
def decompose_iter(n):
    #returns an iterable
    # n = a + b
    if n == 1:
        yield [1]

    for i in range(1,n):
        a = n-i
        for b in decompose_iter(i):
            if b[0] <= a:
                yield [a] + b

    #MUCH SLOWER THAN ITS NON-ITERABLE COUNTERPART
'''

def exhaust(part):
    #takes a partition of an integer as input, and returns a list of every product of irreducible polys with degree corresponding to the partition
    #eg. if part = [3,2] exhaust will return a list of every product of irreducible polys of degree 3 and degree 2

    if len(part) == 1:
        try:
            return Irred[part[0]]
        except:
            print(Irred)
            print(part[0])
            b = str(input('continue ? '))
    else:
        res = []
        products_before = exhaust(part[1:])

        for P in Irred[part[0]]:
            res += [P*Q for Q in products_before]

        #print(res)
        return res
    
def PolyGen(n,mod):
    #generates every monic polynomial in Z/pZ from degree 2 to degree n, here p = mod
    #sorted by degree (in a dictionary -> res)
    deg = 2
    res = {2:[Polynomial({2:GaloisFp(1,mod)})]}
    #print(res)
    #b = str(input('continue ? '))
    while deg <= n:
        curr = res[deg][-1].copy()
        not_done = True
        i = 0
        while not_done and i < deg:
            curr[i] = curr.get(i,0) + GaloisFp(1,mod)
            if not curr[i]:
                i += 1
            else:
                not_done = False
        else:
            if i == deg:
                deg += 1
                res[deg] = []
                curr = {deg:GaloisFp(1,mod)}
                #print(curr)
                #b = str(input('continue ? '))
        res[deg].append(Polynomial(curr))
        #print(res)
        #b = str(input('continue ? '))

    res[deg].pop(-1) #pops the polynomial X^{n+1} that eventually gets added

    return res

p = int(input('p = '))
n = int(input('max degree n = '))

Irred = {1:[Polynomial([GaloisFp(x,p),GaloisFp(1,p)]) for x in range(p)]}

Polys = PolyGen(n,p)

start1 = time.time()

for deg in range(2,n+1):
    Reducible = []
    for part in decompose(deg):
        if part != [deg]:
            try:
                Reducible += exhaust(part)
            except:
                print(Irred)
                print('irred')
                print(part)
                print('partition')
                print(str(input('continue ? ')))
        '''
        for prod in exhaust(part):
            try:
                Reducible = Reducible + prod
            except:
                print(type(Reducible))
                print(type(prod))
                print(str(input('continue ? ')))
        '''

    for P in Reducible:
        for Q in Polys[deg]:
            if (not (P-Q)):
                Polys[deg].remove(Q)

    #print(Reducible)
    #print('these are reducible polys')

    Irred[deg] = [P for P in Polys[deg]]

    #print(Irred[deg])
    #print(str(input('continue ? ')))

end1 = time.time()

querry = str(input('continue ? '))
print({key:len(Irred[key]) for key in Irred.keys()})
print(str(end1 - start1) + ' seconds taken')

#print(len(PolyGen(7,p)))

#TODO write another algorithm which takes every product and instead of doing every unique mult, multiply each irred by previously obtained polynomials (may be composite)

Irred = {1:[Polynomial([GaloisFp(x,p),GaloisFp(1,p)]) for x in range(p)]}

Polys = PolyGen(n,p)
Polys[1] = Irred[1]
SecondPolys = PolyGen(n,p)

start2 = time.time()

Reducible = []

for deg in range(2,n+1):

    for k in Irred.keys():
        for P in Irred[k]:
            for Q in Polys[deg-k]:
                Reducible.append(P*Q)

    for P in Reducible:
        for Q in SecondPolys[deg]:
            if (not (P-Q)):
                SecondPolys[deg].remove(Q)

    Irred[deg] = [P for P in SecondPolys[deg]]

end2 = time.time()

querry = str(input('continue ? '))
print({key:len(Irred[key]) for key in Irred.keys()})
print(str(end2 - start2) + ' seconds taken')