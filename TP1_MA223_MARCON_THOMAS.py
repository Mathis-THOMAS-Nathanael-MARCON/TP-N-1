""" MARCON Nathanael 
    THOMAS Mathis
    2PF2
    TP1 génie mathématique """

########## IMPORTS ##########
import numpy as np
import time
import matplotlib.pyplot as pp

########## FONCTIONS ##########

##### Gauss #####
def ReductionGauss(A):
   for i in range(0, len(A)):
        for k in range(i,len(A)-1):
            g = A[k+1][i]/A[i][i]
            A[k+1] =A[k+1] - g*A[i]
   return A

def ResolutionSystTrisup(Taug):
    
    taille = np.shape(Taug)
    ligne = taille[0]
    
    sol = [0]*ligne
    
    for i in range(ligne-1,-1,-1):
        sol[i] =  Taug[i][ligne]
        
        for k in range(i+1,ligne):
            sol[i] =  sol[i] - Taug[i][k]*sol[k]
        sol[i] = sol[i]/Taug[i][i]
    return sol

def Gauss(A, B):
    tailleA = np.shape(A)
    ligneA = tailleA[0]

    Aaug = np.c_[A, B]
    Taug = ReductionGauss(Aaug)
    sol = ResolutionSystTrisup(Taug)
    return sol

##### LU #####
def DecompositionLU(A):
    matrice_L = np.eye(len(A))
    for i in range(0, len(A)):
        for k in range(i+1,len(A)):
            g = A[k][i]/A[i][i]
            A[k] =A[k] - g*A[i]
            matrice_L[k,i] = g
    return [matrice_L,A]

def ResolutionLU(L, U, B):
    matrice_LB = np.c_[L, B]
    taille = np.shape(matrice_LB)
    ligne = taille[0]
    
    matrice_y = [0]*ligne
    
    for i in range(0,ligne):
        matrice_y[i] =  matrice_LB[i][ligne]/matrice_LB[i][i]
        
        for k in range(0,i):
            matrice_y[i] =  matrice_y[i] - matrice_LB[i][k]*matrice_y[k]
            
    matrice_Ux = np.c_[U, matrice_y]
    sol = ResolutionSystTrisup(matrice_Ux)
    return sol

def MethodeLU(A, B):
    LU = DecompositionLU(A)
    L = LU[0]
    U = LU[1]
    sol = ResolutionLU(L, U, B)
    return sol

##### Pivot Partiel #####
def trig_pivot_partiel(A, B):
    A = np.c_[A, B]

    for i in range(0, len(A)):
        max = abs(A[i][i])
        indiceMax = i
        for j in range(i, len(A)):
            if abs(A[j][i]) > max:
                max = A[j][i]
                indiceMax = j
                
                X = np.copy(A[i])
                A[i] = A[j]
                A[j] = X

        for k in range(i,len(A)-1):
            g = A[k+1][i]/A[i][i]
            A[k+1] =A[k+1] - g*A[i]
    return A

def pivot_partiel(A,B):
    Taug = trig_pivot_partiel(A, B)
    sol = ResolutionSystTrisup(Taug)
    return sol

##### Pivot Partiel #####

def trig_pivot_total(A, B):
    A = np.c_[A, B]
    global Li
    Li= []
    global Lj
    Lj= []

    for i in range(0, len(A)):
        max = abs(A[i][i])
        indiceMax = i
        for j in range(i, len(A)):
            if abs(A[i][j]) > max:
                max = A[i][j]
                indiceMax = j

                X = np.copy(A[:,i])
                A[:,i] = A[:,j]
                A[:,j] = X
                
                Li.append(i)
                Lj.append(j)
        for k in range(i,len(A)-1):
            g = A[k+1][i]/A[i][i]
            A[k+1] =A[k+1] - g*A[i]

    return A

def Resolution_pivot_total(Taug):
    taille = np.shape(Taug)
    ligne = taille[0]
    
    sol = [0]*ligne

    for ldown in range(ligne-1, -1, -1):
        sol[ldown] = Taug[ldown, -1]
        
        for lup in range(ldown+1, ligne):          
            sol[ldown] = sol[ldown] - Taug[ldown, lup]*sol[lup]
     
        sol[ldown] = sol[ldown]/Taug[ldown, ldown]

    for o in range(len(Li)-1, -1, -1):
        X = sol[Li[o]]
        sol[Li[o]] = sol[Lj[o]]
        sol[Lj[o]] = X
    return sol

def pivot_total(A,B):
    A = trig_pivot_total(A, B)
    sol = Resolution_pivot_total(A)
    return sol
    
##### Graphs #####    



def Graph_temps():
    x = []
    y1 = []
    y2 =[]
    y3 = []
    y4 = []
    y5 = []
    for i in range(250,1050,50):
        print(i)
        matrice_A = np.random.rand(i,i)
        matrice_B = np.random.rand(i,1)
        matrice_C = np.copy(matrice_A)
        matrice_D = np.copy(matrice_A)
        matrice_E = np.copy(matrice_A)
        
        time1_1 = time.time()
        Gauss(matrice_A, matrice_B)
        time2_1 = time.time()
        durée1 = time2_1 - time1_1
        #print(durée1)
        
        time1_2 = time.time()
        MethodeLU(matrice_A, matrice_B)
        time2_2 = time.time()
        durée2 = time2_2 - time1_2
       #print(durée2)
        
        time1_3 = time.time()
        pivot_partiel(matrice_C, matrice_B)
        time2_3 = time.time()
        durée3 = time2_3 - time1_3
        #print(durée3)

        time1_4 = time.time()
        pivot_total(matrice_D, matrice_B)
        time2_4 = time.time()
        durée4 = time2_4 - time1_4
        #print(durée4)

        time1_5 = time.time()
        np.linalg.solve(matrice_E, matrice_B)
        time2_5 = time.time()
        durée5 = time2_5 - time1_5
        #print(durée5)
        
        x.append(i)
        y1.append(durée1)
        y2.append(durée2)
        y3.append(durée3)
        y4.append(durée4)
        y5.append(durée5)

    pp.plot (x, y1, "o--", color = 'r', label = 'Gauss')
    pp.plot (x, y2, "o--", color = 'b', label = 'LU')
    pp.plot (x, y3, "o--", color = 'g', label = 'partiel')
    pp.plot (x, y4, "o--", color = 'y', label = 'total')
    pp.plot (x, y5, "o--", color = 'm', label = 'linalg')
    pp.yscale('log')
    pp.xscale('log')
    pp.title("temps d'execution en fonction de la taille de la matrice")
    pp.xlabel("taille de la matrice")
    pp.ylabel("temps d'execution (s)")
    pp.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', title="Legend")                                                             #afficher la grille
    pp.show()
    pp.savefig('figure courbe tp1')

def Graph_erreur():
    x = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    for i in range(250,1050,50):
        print(i)
        
        matrice_A = np.random.rand(i,i)
        matrice_B = np.random.rand(i,1)
        matrice_C = np.copy(matrice_A)
        matrice_D = np.copy(matrice_A)
        matrice_E = np.copy(matrice_A)
        matrice_F = np.copy(matrice_A)
        

        X1 = Gauss(matrice_A, matrice_B)
        Xabs1 = np.linalg.norm(np.dot(matrice_A, X1)-np.ravel(matrice_B))


        X2 = MethodeLU(matrice_A, matrice_B)
        Xabs2 = np.linalg.norm(np.dot(matrice_C, X2)-np.ravel(matrice_B))


        X3 = pivot_partiel(matrice_C, matrice_B)
        Xabs3 = np.linalg.norm(np.dot(matrice_C, X3)-np.ravel(matrice_B))


        X4 = pivot_total(matrice_D, matrice_B)
        Xabs4 = np.linalg.norm(np.dot(matrice_D, X4)-np.ravel(matrice_B))


        X5 = np.linalg.solve(matrice_E, matrice_B)
        Xabs5 = np.linalg.norm(np.dot(matrice_E, X5)-matrice_B)

        
        x.append(i)
        y1.append(Xabs1)
        y2.append(Xabs2)
        y3.append(Xabs3)
        y4.append(Xabs4)
        y5.append(Xabs5) 
    
    pp.plot (x, y1, "o--", color = 'r', label = 'Gauss')
    pp.plot (x, y2, "o--", color = 'b', label = 'LU')
    pp.plot (x, y3, "o--", color = 'g', label = 'partiel')
    pp.plot (x, y4, "o--", color = 'y', label = 'total')
    pp.plot (x, y5, "o--", color = 'm', label = 'linalg')
    pp.yscale('log')
    #pp.xscale('log')
    pp.title("Erreur en fonction de la taille de la matrice")
    pp.xlabel("taille de la matrice")
    pp.ylabel("Erreur")
    pp.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', title="Legend")                                                             #afficher la grille
    pp.show()
    pp.savefig('figure courbe 2 tp1')


'''
Appel fonctions graphiques :
#Graph_erreur()
#Graph_temps()

matrices de test :
A = np.random.rand(4,4)
B = np.random.rand(4,1)
#A = np.array([[3, 2, -1, 4], [-3, -4, 4, -2], [6, 2, 2, 7], [9, 4, 2, 18]])
C = np.copy(A)
#B = np.array([4, -5, -2, 13])

x = choix_pivot(A,B)
print(np.dot(A,x))
print(B)

T = trig_choix_pivot(A)
print(T)

A = np.random.rand(4,4)
B = np.random.rand(4,1)


x = MethodeLU(A, B)
AX = np.dot(C, x)

print(A)
print(x)
print(B)
print(AX)


LU = DecompositionLU(matrice_testA)
print(LU)
ResolutionLU(LU[0], LU[1], matrice_testB)

########## TEST DES FONCTIONS ##########

if __name__ == "__main__":
    
    print("Tests :", '\n')
    matrice_test = np.array([[1, 1, 1, 1, 1], [2, 4, -3, 2, 1], [-1, -1, 0, -3, 2], [1, -1, 4, 9, -8]])
    print("matrice A augmentée :", '\n')
    print(matrice_test, '\n')
    testGauss = ReductionGauss(matrice_test)
    print("matrice A augmentée triangulaire :", '\n')
    print(testGauss)
   
    print("=======================================")
     
    print("solutions du système :", '\n')
    Remontée = ResolutionSystTrisup(testGauss)
    print(Remontée)
    
    print("=======================================")
    
    matrice_testA = np.array([[1, 1, 1, 1], [2, 4, -3, 2], [-1, -1, 0, -3], [1, -1, 4, 9]])
    #matrice_testA = np.random.rand(500,500)
    print(matrice_testA)
    matrice_testB = np.array([1, 1, 2, -8])
    #matrice_testB = np.random.rand(500,1)
    print(matrice_testB)
    
    print("solutions du système (complet):", '\n')

    C = Gauss(matrice_testA, matrice_testB)
    print(C)

    print("Test 2 :", '\n')
    matrice_test = np.array([[3, 2, -1, 4, 4], [-3, -4, 4, -2, -5], [6, 2, 2, 7, -2], [9, 4, 2, 18, 13]])
    print("matrice A augmentée :", '\n')
    print(matrice_test, '\n')
    testGauss = ReductionGauss(matrice_test)
    print("matrice A augmentée triangulaire :", '\n')
    print(testGauss)
   
    print("=======================================")
     
    print("solutions du système :", '\n')
    Remontée = ResolutionSystTrisup(testGauss)
    print(Remontée)
    
    print("=======================================")
    
    matrice_testA = np.array([[3, 2, -1, 4], [-3, -4, 4, -2], [6, 2, 2, 7], [9, 4, 2, 18]])
    #matrice_testA = np.random.rand(500,500)
    print(matrice_testA)
    matrice_testB = np.array([4, -5, -2, 13])
    #matrice_testB = np.random.rand(500,1)
    print(matrice_testB)
    
    print("solutions du système (complet):", '\n')
    C = Gauss(matrice_testA, matrice_testB)
    print(C)
'''
