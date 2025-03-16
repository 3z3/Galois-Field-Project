from polynomial_class import *

q = { 1: [[1]] }

def decompose(n):
    #returns integer partitions of n, eg. 4 -> [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]
    #executes quickly approx up to n = 60
    try:
        return q[n]
    except:
        pass

    result = [[n]]

    for i in range(1, n):
        a = n-i
        R = decompose(i)
        for r in R:
            if r[0] <= a:
                result.append([a] + r)

    q[n] = result
    return result

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

        return res
    
def PolyGen(n,mod):
    #generates every monic polynomial in Z/pZ from degree 2 to degree n, here p = mod
    #sorted by degree (in a dictionary -> res)
    deg = 2
    res = {2:[Polynomial({2:GaloisFp(1,mod)})]}
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
        res[deg].append(Polynomial(curr))

    res[deg].pop(-1) #pops the polynomial X^{n+1} that eventually gets added

    return res

p = int(input('p = '))
n = int(input('max degree n = '))

Irred = {1:[Polynomial([GaloisFp(x,p),GaloisFp(1,p)]) for x in range(p)]}

Polys = PolyGen(n,p)

start1 = time.time()

for deg in range(2,n+1):
    #sieve method
    #at the end, shows irreducible polynomials of degree up to n, in Z/pZ
    Reducible = []
    for part in decompose(deg):
        if part != [deg]:
            try:
                Reducible += exhaust(part)

    for P in Reducible:
        for Q in Polys[deg]:
            if (not (P-Q)):
                Polys[deg].remove(Q)

    Irred[deg] = [P for P in Polys[deg]]

end1 = time.time()

querry = str(input('continue ? '))
print({key:len(Irred[key]) for key in Irred.keys()})
print(str(end1 - start1) + ' seconds taken')

Irred = {1:[Polynomial([GaloisFp(x,p),GaloisFp(1,p)]) for x in range(p)]}

Polys = PolyGen(n,p)
Polys[1] = Irred[1]
SecondPolys = PolyGen(n,p)

start2 = time.time()

Reducible = []

for deg in range(2,n+1):
    #other method where instead of sieving according to every unique product we do every multiplication (might encounter multiples, for example (X-1)*(X^2 + X) and (X^2 - 1)*X alike)
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

#at the end both methods return the same irreducible polynomials, and the time taken with each method is compared
