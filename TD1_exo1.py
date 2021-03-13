""" MARCON Nathanael 
    THOMAS Mathis
    2PF2
    TP1 génie mathématique """

import numpy as np
import time
import matplotlib.pyplot as pp

def ReductionGauss(Aaug):
    
    taille = np.shape(Aaug)
    ligne = taille[0]
    colonne = taille[1]
    
    for lup in range(0,ligne-1):
        
        for ldown in range(lup+1, ligne):

            g = Aaug[ldown,lup]/Aaug[lup,lup]
            
            for c in range(lup, colonne):
                
                Aaug[ldown,c] = Aaug[ldown,c] - g*Aaug[lup,c]
    return Aaug

def ResolutionSystTrisup(Taug):
    
    taille = np.shape(Taug)
    ligne = taille[0]
    
    sol = [0]*ligne
    sol[ligne-1] = Taug[-1,-1]/Taug[-1,-2]

    for ldown in range(ligne-2, -1, -1):
        print(sol)
        sol[ldown] = Taug[ldown, -1]
        
        for lup in range(ldown+1, ligne):
            
            sol[ldown] = sol[ldown] - Taug[ldown, lup]*sol[lup]
     
        sol[ldown] = sol[ldown]/Taug[ldown, ldown]
    
    return sol

def Gauss(A, B):
    tailleA = np.shape(A)
    ligneA = tailleA[0]

    Aaug = np.c_[A, B]
    Taug = ReductionGauss(Aaug)
    sol = ResolutionSystTrisup(Taug)
    
    return sol

def DecompositionLU(A):
    
    taille = np.shape(A)
    ligne = taille[0]
    colonne = taille[1]
    matrice_L = np.eye(ligne)
    
    for lup in range(0,ligne-1):
        
        for ldown in range(lup+1, ligne):

            g = A[ldown,lup]/A[lup,lup]
            
            for c in range(lup, colonne):
                
                A[ldown,c] = A[ldown,c] - g*A[lup,c]
                matrice_L[ldown,lup] = g 
    return [matrice_L,A]


    

def Graph_temps():
    x = []
    y = []
    for i in range(50,1500,50):
        matrice_A = np.random.rand(i,i)
        matrice_B = np.random.rand(i,1)
        time1 = time.time()
        print(Gauss(matrice_A, matrice_B))
        time2 = time.time()
        durée = time2 - time1
        print(durée)
        x.append(i)
        y.append(durée)
    pp.plot (x, y, "o--", color = 'r')
    pp.title("temps d'execution en fonction de la taille de la matrice")
    pp.xlabel("taille de la matrice")
    pp.ylabel("temps d'execution")
    pp.legend()                                                               #afficher la grille
    pp.show()
    pp.savefig('figure courbe tp1')

def Graph_erreur():
    x = []
    y = []
    for i in range(2,80,1):
        matrice_A = np.random.rand(i,i)
        matrice_B = np.random.rand(i,1)
        X = Gauss(matrice_A, matrice_B)
        print("SOL",X)
        print('\n')
        print(np.ravel(matrice_B))
        print('\n')
        print(matrice_A)
        print('\n')
        print("PRODUIT",matrice_A.dot(X))
        Xabs = np.linalg.norm(np.dot(matrice_A, X)-np.ravel(matrice_B))
        x.append(i)
        y.append(Xabs)
    pp.plot (x, y, "o--", color = 'r')
    pp.title("Erreur en fonction de la taille de la matrice")
    pp.xlabel("taille de la matrice")
    pp.ylabel("Erreur")
    pp.legend()                                                               #afficher la grille
    pp.show()
    pp.savefig('figure courbe tp1')

Graph_temps()

#matrice_testA = np.array([[3, 2, -1, 4], [-3, -4, 4, -2], [6, 2, 2, 7], [9, 4, 2, 18]])
#LU = DecompositionLU(matrice_testA)
#print(LU)
########## TEST DES FONCTIONS ##########
'''
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