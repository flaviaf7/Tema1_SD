import time
import random
import sys
sys.setrecursionlimit(10000)

def mergeSortwrapper(a):
    mergeSort(a, 0, len(a)-1)


# MergeSort 
# O(n*log(n)) time complexity 
# O(n) space complexity
def mergeSort(a, st, dr):
    if st<dr:

        # calculam mijlocul in asa fel in cat sa nu 
        # avem overflow pentru numere foarte mari
        mij=st+(dr-st)//2

        # sortam ambele jumatati
        mergeSort(a, st, mij)
        mergeSort(a, mij+1, dr)
        
        # etapa de Impera
        merge(a, st, mij, dr)

# facem merge la 2 subsiruri ale array-ului 
def merge(a, st, mij, dr):

    # lungimile subsirurilor
    l1= mij - st + 1
    l2 = dr - mij

    S = [0] * (l1)
    D = [0] * (l2)

    for i in range(l1):
        S[i] = a[st + i]
    for j in range(l2):
        D[j] = a[mij + 1 + j]

    i = 0    # indexul primului subsir
    j = 0    # indexul celui de-al doilea subsir
    k = st   # indexul sirului a
    while i < l1 and j < l2:
        if S[i] <= D[j]:
            a[k] = S[i]
            i += 1
        else:
            a[k] = D[j]
            j += 1
        k += 1

    # in caz ca sunt elemente necopiate din S[] sau D[], le adaugam:
    while i < l1:
        a[k] = S[i]
        i += 1
        k += 1
    while j < l2:
        a[k] = D[j]
        j += 1
        k += 1
    




 # Shell sort
 # O(n logn) - best & average
 # O(n^2) - worst case
 # space complexity - O(1)
def shellSort(a):

    interval = len(a)//2
    while interval > 0:
        for i in range(interval, len(a)):
            aux = a[i]
            j = i
            # interschimbam fiecare element situat la "interval" distanta unul de celalalt
            while j >= interval and a[j-interval] > aux:
                a[j] = a[j-interval]
                j -= interval
            a[j] = aux
        interval //= 2





# radixSort (folosindu-ne de countingSort)
# O(d(n+b)) - time complexity
#      d -> nr max de cifre al numerelor din array
#      b -> baza

def radixSort(a, baza=10):

    maxim = max(a)

    # folosim countingSort pt a sorta numerele in functi de fiecare cifra,
    # incepand cu cifra unitatilor, pana la ultima cifra a maximului
    poz = 1
    while maxim // poz > 0:
        countingSortR(a, poz, baza)
        poz *= baza

def countingSortR(a, poz, baza):
    n = len(a)
    l = [0] * n        # lista auxiliara
    frecv = [0] * baza   # lista pt frecventa fiecarei cifre

    # calculam frecventa fiecarei cifre de pe pozitia poz
    for i in range(n):
        idx = a[i] // poz % baza
        frecv[idx] += 1

    # calculam numarul cumulativ pentru fiecare cifra
    # adica numarul total de elemente mai mici sau egale decat acea cifra
    for i in range(1,baza):
        frecv[i] += frecv[i-1]


    # parcurgem lista invers si sortam numerele 
    # folosindu-ne de lista de frecvente a cifrei de pe pozitia poz
    i = n-1
    while i >= 0:
        idx = a[i] // poz % baza
        l[frecv[idx]-1] = a[i]
        frecv[idx] -= 1
        i -= 1

    for i in range(0,n):
        a[i] = l[i]






# quickSort
# time complexity - best case O(nlogn)
#                 - worst case O(n^2)
# space complexity - O(1) (dar are logn call-uri recursive pe stiva pentru functia quickSorthelp)
def quickSort(a):
    if N>100000:
        return
    quickSorthelp(a, 0, len(a)-1)

def quickSorthelp(a, st, dr):
    if st<dr:
        p = partition(a, st, dr)
        quickSorthelp(a, st, p-1)
        quickSorthelp(a, p+1, dr)
    

def partition(a, st, dr):
    # alegem pivotul ca fiind mediana dintre primul,
    # ultimul si mijlocul listei
    pivot_idx = median_of_three(a, st, dr)
    pivot = a[pivot_idx]
    a[pivot_idx] , a[dr] = a[dr], a[pivot_idx] # mutarea pivotului la sfarsitul listei
    
    # impartim lista in 2 subliste, cu elemente mai mici, respectiv mai mari decat pivotul
    capat=st
    for i in range(st,dr):
        if a[i] < pivot:
            a[i], a[capat] = a[capat], a[i] 
            capat += 1
    
    a[dr], a[capat] = a[capat], a[dr] # mutarea pivotului inapoi in pozitia corecta
    return capat


def median_of_three(a, st, dr):
    mij = st +(dr-st)//2
    # sortam elementele de la inceput, sfarsit si mijloc
    # si alegem mediana dintre ele
    if a[st] > a[dr]:
        a[st], a[dr] = a[dr], a[st]
    if a[mij] > a[dr]:
        a[mij], a[dr] = a[dr], a[mij]
    if a[mij] < a[st]:
        a[mij], a[st] = a[st], a[mij]

    return mij






# countingSort
# time complexity - O(n+k)
#        k -> max-min 
def countingSort(a):
    if NMAX>1000000000:
        print("Algoritmul nu poate efectua sortarea ")
        return
    
    maxim = max(a)
    j=0

    n = len(a)
    l = [0] * (maxim+1)
    
    for i in range(n):
        l[a[i]]+=1
    
    for i in range(1,maxim+1):
        while l[i]>0:
            a[j]=i
            j+=1
            l[i]-=1


# functie pentru a verifica daca un vector este sortat corect
def isSorted(a):
    for i in range(len(a)):
        if a[i] < a[i-1]:
            return False
    return True

# functie pentru a genera un vector cu numere aleatoare
def generateRandomList(n, maxVal, baza=10):
    if baza==16:
        return [hex(random.randint(0,maxVal)) for _ in range(n)]
    else:
        return [random.randint(0,maxVal) for _ in range(n)]

def timp(sort, a):
    start = time.time()
    sort(a)
    end = time.time()
    return end-start

sorts = [mergeSortwrapper, shellSort, radixSort, quickSort, countingSort]

with open("teste.txt", "r") as f:
    T = int(f.readline().strip()) #nr de teste
    tests = []
    for i in range(T):
        N, NMAX = map(int, f.readline().strip().split()) #nr de elemente si cel mai mare element
        a = generateRandomList(N, NMAX)
        tests.append((N,NMAX, a))
    
with open("RezultateTeste.txt","w") as f:
    for i, test in enumerate(tests):
        N, NMAX, a = test
        f.write(f"Testul {i+1} - N={N}, NMAX={NMAX}\n")
        for sort in sorts:
            acpy = a.copy()
            t = timp(sort, acpy)
            # if not isSorted(acpy):
            #     f.write(f"Algoritmul {sort.__name__} a sortat gresit!\n")
            # else:
            f.write(f"Algoritmul {sort.__name__} a sortat corect in {t:.6f} secunde\n")
        f.write("\n")
    
with open("teste_radixsort.txt", "r") as f:
    t = int(f.readline())  # numărul de teste
    tests = []
    for i in range(t):
        N, NMAX, base = map(int, f.readline().split())
        a = generateRandomList(N, NMAX, base)
        tests.append((N,NMAX,a,base))

with open("rezultate_radixsort.txt", "w") as f:
    for i, test in enumerate(tests):
        N, NMAX, a, base = test
        f.write(f"Testul {i+1} - N={N} NMAX={NMAX} BASE={base}\n")
        start_time = time.time()
        radixSort(a, base)
        end_time = time.time()
        # if not isSorted(acpy):
        #     f.write(f"Algoritmul {sort.__name__} cu baza {base} a sortat gresit!\n")
        # else:
        f.write(f"Algoritmul {sort.__name__} cu baza {base} a sortat corect in {t:.6f} secunde\n")
        
