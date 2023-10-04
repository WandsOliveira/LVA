"""
CONJUNTO DE FUNÇÕES UTILIZADAS..
Autor : Wanderson Vinicius de Oliveira MOnteiro

Obs1: Realizar documentação
Obs2: Mudar nome arquivo

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy
from tqdm import tqdm
class Rod:

    def __init__(self) -> None:
        pass
    
    class FEM:
        def __init__(self, E, rho, eta, b, h, L, Nc, Frequency, Nos) -> None:
            self.b = b
            self.h = h
            self.L = L
            self.Nc = Nc
            self.Nos = Nos
            self.Frequency = np.arange(1, Frequency + 1)
            self.Rod_Number = Nc * len(E)
            self.List_GL_Materials = list(range(0, self.Rod_Number  ))
            self.E = np.tile(E*(1 + 1j*eta), int(len(self.List_GL_Materials)/len(E)))
            self.rho = np.tile(rho, int(len(self.List_GL_Materials)/len(rho)))
            self.A = self.b * self.h  # Área
            self.I = self.b * self.h ** 3 / 12  # Momento de Inércia
            self.Nos_FEM = self.Nos
            self.Elementos_FEM = self.Nos_FEM - 1
            self.Gdl = 1
            self.Nos_Principais = len(E) + 1 + ((len(E) * (self.Nc -1 )))
            self.Global_FEM = self.Gdl * (self.Nos_Principais + (self.Elementos_FEM - 1) * self.Rod_Number)
            self.GL_Elemento = np.zeros((self.Rod_Number, self.Gdl * (2 + (self.Elementos_FEM - 1))))
            self.L_FEM = self.L / (2 * self.Elementos_FEM)
            self.GDL_Principais = np.arange(0, self.Nos_Principais)
            self.u_max = self.GDL_Principais[-1]
            # Preenchimento dos nós principais
            vv = 0
            for ii in range(self.Rod_Number):
                self.GL_Elemento[ii, 0] = self.GDL_Principais[vv]
                self.GL_Elemento[ii, -1] = self.GDL_Principais[vv + 1]
                vv = vv + 1
            ## Preenchimento dos nós internos
            M_aux = np.zeros((self.Rod_Number, self.Gdl * (2 + (self.Elementos_FEM - 1) )))
            ini =  np.amax(self.GL_Elemento)
            seq = np.arange(ini + 1, ini + self.Global_FEM)
            zz = 0

            for i in np.arange(len(M_aux)):
                for j in np.arange(1,len(M_aux[0]) - 1):
                    M_aux[i, j] = seq[zz]
                    zz += 1 

            self.GL_Elemento = self.GL_Elemento + M_aux

            self.KG, self.MG = self.Barra_Matriz_Massa_Rigidez_Global(self.E, self.rho, self.L_FEM, self.A, self.Rod_Number, self.GL_Elemento, self.Global_FEM)
        def Barra_Matriz_Massa_Rigidez_Global(self, E, rho, L_FEM, A, Beam_Number, GL_Elemento, Global):
            ## Alocação de memoria matrizes de Massa e Rigidez
            KG = np.zeros((Global, Global), dtype= complex) # Rigidez
            MG = np.zeros((Global, Global), dtype= complex) # Massa

            # Alocando materiais

            for i in np.arange(Beam_Number):

                row_aux = GL_Elemento[i,:]
                for ww in range(0, len(row_aux) - 1, 1):
                    row_columns = row_aux[ ww: ww + 2]
                    # Matrizes locais
                    Ke = (E[i] * A / L_FEM) * np.array([[1/2, -1/2], [-1/2, 1/2]])  # Matriz de rigidez local
                    Me = rho[i] * A * L_FEM * np.array([[2/3, 1/3], [1/3, 2/3]])  # Matriz de Massa loca
                    KGe =   Ke  
                    MGe =   Me 
                    # Montagem da matriz de massa e rigidez global
                    Matrix_Aux = np.zeros((3, 3), dtype= complex)    # Matriz auxiliar da rigidez
                    Matrix_Aux_m = np.zeros((3, 3), dtype= complex)  # Matriz auxiliar da massa
                    Matrix_Aux[1 : 3, 1 : 3] =  KGe[:,:]
                    Matrix_Aux_m[1 : 3, 1 : 3] =  MGe[:,:]
                    Matrix_Aux[0, 1 : 3] =  row_columns[0:2]
                    Matrix_Aux[1 : 3, 0] =  row_columns[0:2]         
                    Matrix_Aux_m[0, 1 : 3] =  row_columns[0:2]
                    Matrix_Aux_m[1 : 3, 0] =  row_columns[0:2]  
                    for j in range(2):
                        Columns_aux = int(np.real(Matrix_Aux[0, j + 1]))
                        Columns_aux_m = int(np.real(Matrix_Aux_m[0, j + 1]))
                        for k in range(2):
                            Line_aux = int(np.real(Matrix_Aux[k + 1, 0]))
                            Line_aux_m = int(np.real(Matrix_Aux_m[k + 1, 0]))
                            KG[Line_aux, Columns_aux] += Matrix_Aux[k + 1, j + 1]           # Matriz rigidez
                            MG[Line_aux_m, Columns_aux_m] += Matrix_Aux_m[k + 1, j + 1]     # Matriz massa
             
            return KG, MG

        # Função para o spy
        def Spy_matriz(self, Matriz):
            non_zeros_indices = np.where(np.abs(Matriz) > 0)   # Valores não nulos

            row_indices, col_indices = non_zeros_indices
            plt.scatter(col_indices, row_indices, marker='s', color='b')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.grid(True)
            plt.title("Matriz de Rigidez")
            plt.gca().invert_yaxis()  # Convenção do Matlab
            plt.show()

        def Barra_Condicao_Contorno(self, Matriz, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento):
            F = np.zeros((len(Matriz)))               # Alocação de Memória : Força 
            F[GDL_Forca] = Intensidade_Forca        # Contorno da força

            Matriz_Nova = None 
            F_Novo = None

            try:
                if Contorno_Deslocamento == 'Livre - Livre':
                    Matriz_Nova = np.copy(Matriz)
                    F_Novo = np.copy(F)
                elif Contorno_Deslocamento == 'Engastado - Livre':
                    Matriz_Nova = Matriz[1:, 1:]
                    F_Novo  = F[1:]
                elif Contorno_Deslocamento == 'Livre - Engastado':
                    Matriz_Nova = Matriz[:-1, :-1]
                    F_Novo = F[:-1]
                else:
                    raise ValueError("Condição de deslocamento inválida")
            except ValueError as e:
                print(f"Erro: {e}")

            return F_Novo , Matriz_Nova 


        def Barra_Final_Response(self, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento, Type):

            u = np.zeros((1 , len(self.Frequency)), dtype=complex)  # Alocação de Memória: Deslocamento

            for i, omega in tqdm(enumerate(self.Frequency), total=len(self.Frequency), desc="Calculando FEM:" ):
                omega = 2 * np.pi * self.Frequency[i-1]
                U = np.zeros((self.Global_FEM))
                Dinamica = self.KG - (omega ** 2) * self.MG  # Matriz de rigidez dinâmica
                F, D = self.Barra_Condicao_Contorno(Dinamica, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento)
                U = np.linalg.solve(D, F)
                try: 
                    if Type == 'Receptancia':
                        u[0, i -1] = U[self.u_max]
                    elif Type == 'Transmitancia' :
                        u[0, i -1] = U[self.u_max] / U[0]
                    
                    elif Type == 'Acelerancia' :
                        u[0, i -1] = - omega ** 2 * U[self.u_max] 
                    else:
                        raise ValueError("Tipo de plot inválido!")
                except ValueError as e:
                        print(f"Erro: {e}")
            return u
    class SEM:
         
        def __init__(self, E, rho, eta, b, h, L, Nc, Frequency) -> None:

            self.E = E * (1 + 1j * eta)
            self.rho = rho
            self.b = b
            self.h = h
            self.L = L
            self.Nc = Nc
            self.Frequency = np.arange(1, Frequency + 1)
            self.Rod_Number = Nc * len(E)
            self.List_GL_Materials = list(range(0, self.Rod_Number  ))
            self.E = np.tile(E * (1 + 1j*eta), int(len(self.List_GL_Materials)/len(E)))
            self.rho = np.tile(rho, int(len(self.List_GL_Materials)/len(rho)))
            self.A = self.b * self.h  # Área
            self.Nos_SEM = 2
            self.Elementos_SEM = self.Nos_SEM - 1
            self.Gdl = 1
            self.Nos_Principais_SEM =  len(E) + 1 + (len(E) * (self.Nc - 1))
            self.Global_SEM= self.Gdl * (self.Nos_Principais_SEM + (self.Elementos_SEM - 1) * self.Rod_Number)
            self.GL_Elemento_SEM = np.zeros((self.Rod_Number, self.Gdl * (2 + (self.Elementos_SEM - 1))))
            self.GDL_Principais_SEM = np.arange(0, self.Nos_Principais_SEM)
            self.u_max_SEM = self.GDL_Principais_SEM[-1]
            # Preenchimento dos nós principais
            vv = 0
            for ii in range(self.Rod_Number):
                self.GL_Elemento_SEM[ii, 0] = self.GDL_Principais_SEM[vv]
                self.GL_Elemento_SEM[ii, -1] = self.GDL_Principais_SEM[vv + 1]
                vv = vv + 1




        def Spectral_Dynamics_Matrix(self, E, rho, L, A, Rod_Number, GL_Elemento_SEM, Global_SEM, omega):
            ## Alocação de memoria matrizes de Massa e Rigidez
            S = np.zeros((Global_SEM, Global_SEM), dtype= complex) #  Spectral Dynamics Matrix

            for i in np.arange(Rod_Number):

                row_aux = GL_Elemento_SEM[i,:]
                for ww in range(0, len(row_aux) - 1, 1):
                    row_columns = row_aux[ ww: ww + 2]
                    # Matrizes locais


                    K_L = omega * np.sqrt(rho[i]/E[i])  # Número de Onda da barra

                    Se = (E[i] * A / L) * np.array([[K_L * L * (1/np.tan(K_L * L)), - K_L * L * (1/np.sin(K_L * L))], [- K_L * L * (1/np.sin(K_L * L)), K_L * L * (1/np.tan(K_L * L))]])
                    
                    # Montagem da matriz de massa e rigidez global
                    Matrix_Aux = np.zeros((3, 3), dtype= complex)    # Matriz auxiliar da rigidez
                    Matrix_Aux[1 : 3, 1 : 3] =  Se[:,:]
                    Matrix_Aux[0, 1 : 3] =  row_columns[0:2]
                    Matrix_Aux[1 : 3, 0] =  row_columns[0:2]         
 
                    for j in range(2):
                        Columns_aux = int(np.real(Matrix_Aux[0, j + 1]))

                        for k in range(2):
                            Line_aux = int(np.real(Matrix_Aux[k + 1, 0]))

                            S[Line_aux, Columns_aux] += Matrix_Aux[k + 1, j + 1]           # Matriz rigidez

             
            return S
        
        def Barra_Condicao_Contorno(self, Matriz, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento):
            F = np.zeros((len(Matriz)))               # Alocação de Memória : Força 
            F[GDL_Forca] = Intensidade_Forca          # Contorno da força

            Matriz_Nova = None 
            F_Novo = None

            try:
                if Contorno_Deslocamento == 'Livre - Livre':
                    Matriz_Nova = np.copy(Matriz)
                    F_Novo = np.copy(F)
                elif Contorno_Deslocamento == 'Engastado - Livre':
                    Matriz_Nova = Matriz[1:, 1:]
                    F_Novo  = F[1:]
                elif Contorno_Deslocamento == 'Livre - Engastado':
                    Matriz_Nova = Matriz[:-1, :-1]
                    F_Novo = F[:-1]
                else:
                    raise ValueError("Condição de deslocamento inválida")
            except ValueError as e:
                print(f"Erro: {e}")

            return F_Novo , Matriz_Nova 


        def Barra_Final_Response(self, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento, Type):

            u = np.zeros((self.Global_SEM, len(self.Frequency)), dtype=complex)  # Alocação de Memória: Deslocamento

            for i, omega in tqdm(enumerate(self.Frequency), total=len(self.Frequency), desc="Calculando SEM:" ):
                omega = 2 * np.pi * self.Frequency[i-1]
                self.S = self.Spectral_Dynamics_Matrix(self.E, self.rho, self.L, self.A, self.Rod_Number, self.GL_Elemento_SEM, self.Global_SEM, omega)
                U = np.zeros((self.Global_SEM))
                F, D = self.Barra_Condicao_Contorno(self.S, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento)
                U = np.linalg.solve(D, F)
                try: 
                        if Type == 'Receptancia':
                            u[0, i -1] = U[self.u_max_SEM]
                        elif Type == 'Transmitancia' :
                            u[0, i -1] = U[self.u_max_SEM] / U[0]
                        elif Type == 'Acelerancia' :
                            u[0, i -1] = - omega ** 2 * U[self.u_max_SEM] 
                        else:
                            raise ValueError("Tipo de plot inválido!")
                except ValueError as e:
                        print(f"Erro: {e}")
            
            return u
    class WFE: 
        def __init__(self, E, rho, eta, b, h, L, Nc, Frequency,  Nos) -> None:
            self.b = b
            self.h = h
            self.L = L
            self.Nc = Nc
            self.Nos = Nos
            self.Frequency = np.arange(1, Frequency + 1)
            self.Rod_Number_WFE = len(E)
            self.List_GL_Materials = list(range(0, self.Rod_Number_WFE  ))
            self.E = np.tile(E * (1 + 1j*eta), int(len(self.List_GL_Materials)/len(E)))
            self.rho = np.tile(rho, int(len(self.List_GL_Materials)/len(rho)))
            self.A = self.b * self.h  # Área
            self.I = self.b * self.h ** 3 / 12  # Momento de Inércia
            self.Nos_WFE = self.Nos
            self.Elementos_WFE = self.Nos_WFE - 1
            self.Gdl = 1
            self.Nos_Principais_WFE = len(E) + 1 
            self.Global_WFE = self.Gdl * (self.Nos_Principais_WFE + (self.Elementos_WFE - 1) * self.Rod_Number_WFE)
            self.GL_Elemento_WFE = np.zeros((self.Rod_Number_WFE, self.Gdl * (2 + (self.Elementos_WFE - 1))))
            self.L_WFE = self.L / (2 * self.Elementos_WFE)
            self.GDL_Principais_WFE = np.arange(0, self.Nos_Principais_WFE)
            # Preenchimento dos nós principais
            vv = 0
            for ii in range(self.Rod_Number_WFE):
                self.GL_Elemento_WFE[ii, 0] = self.GDL_Principais_WFE[vv]
                self.GL_Elemento_WFE[ii, -1] = self.GDL_Principais_WFE[vv + 1]
                vv = vv + 1
            ## Preenchimento dos nós internos
            M_aux = np.zeros((self.Rod_Number_WFE, self.Gdl * (2 + (self.Elementos_WFE - 1) )))
            ini =  np.amax(self.GL_Elemento_WFE)
            seq = np.arange(ini + 1, ini + self.Global_WFE)
            zz = 0

            for i in np.arange(len(M_aux)):
                for j in np.arange(1,len(M_aux[0]) - 1):
                    M_aux[i, j] = seq[zz]
                    zz += 1 

            self.GL_Elemento_WFE = self.GL_Elemento_WFE + M_aux

            self.nl = np.array(self.GL_Elemento_WFE[0,0])
            self.nl = self.nl.astype(int)
            self.nr = np.array(self.GL_Elemento_WFE[-1, -1])
            self.nr = self.nr.astype(int)
            self.ni = np.array([self.GL_Elemento_WFE[1:-1, 1:-1]])
            self.ni = self.GL_Elemento_WFE.flatten()
            self.ni = self.ni[self.ni != self.nl]  
            self.ni = self.ni[self.ni != self.nr] 
            self.ni = self.ni.astype(int)
            self.ni = np.unique(self.ni)



            self.KG_WFE, self.MG_WFE = self.Barra_Matriz_Massa_Rigidez_Global(self.E, self.rho, self.L_WFE, self.A, self.Rod_Number_WFE, self.GL_Elemento_WFE, self.Global_WFE)

        def Barra_Matriz_Massa_Rigidez_Global(self, E, rho, L_FEM, A, Beam_Number, GL_Elemento, Global):
            ## Alocação de memoria matrizes de Massa e Rigidez
            KG = np.zeros((Global, Global), dtype= complex) # Rigidez
            MG = np.zeros((Global, Global), dtype= complex) # Massa

            for i in np.arange(Beam_Number):

                row_aux = GL_Elemento[i,:]
                for ww in range(0, len(row_aux) - 1, 1):
                    row_columns = row_aux[ ww: ww + 2]
                    # Matrizes locais
                    Ke = (E[i] * A / L_FEM) * np.array([[1/2, -1/2], [-1/2, 1/2]])  # Matriz de rigidez local
                    Me = rho[i] * A * L_FEM * np.array([[2/3, 1/3], [1/3, 2/3]])  # Matriz de Massa loca

                    KGe =   Ke  
                    MGe =   Me 
                    # Montagem da matriz de massa e rigidez global
                    Matrix_Aux = np.zeros((3, 3), dtype= complex)    # Matriz auxiliar da rigidez
                    Matrix_Aux_m = np.zeros((3, 3), dtype= complex)  # Matriz auxiliar da massa
                    Matrix_Aux[1 : 3, 1 : 3] =  KGe[:,:]
                    Matrix_Aux_m[1 : 3, 1 : 3] =  MGe[:,:]
                    Matrix_Aux[0, 1 : 3] =  row_columns[0:2]
                    Matrix_Aux[1 : 3, 0] =  row_columns[0:2]         
                    Matrix_Aux_m[0, 1 : 3] =  row_columns[0:2]
                    Matrix_Aux_m[1 : 3, 0] =  row_columns[0:2]  
                    for j in range(2):
                        Columns_aux = int(np.real(Matrix_Aux[0, j + 1]))
                        Columns_aux_m = int(np.real(Matrix_Aux_m[0, j + 1]))
                        for k in range(2):
                            Line_aux = int(np.real(Matrix_Aux[k + 1, 0]))
                            Line_aux_m = int(np.real(Matrix_Aux_m[k + 1, 0]))
                            KG[Line_aux, Columns_aux] += Matrix_Aux[k + 1, j + 1]           # Matriz rigidez
                            MG[Line_aux_m, Columns_aux_m] += Matrix_Aux_m[k + 1, j + 1]     # Matriz massa
                
            return KG, MG

        # Função para o spy
        def Spy_matriz(self, Matriz):
            non_zeros_indices = np.where(np.abs(Matriz) > 0)   # Valores não nulos

            row_indices, col_indices = non_zeros_indices
            plt.scatter(col_indices, row_indices, marker='s', color='b')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.grid(True)
            plt.title("Matriz de Rigidez")
            plt.gca().invert_yaxis()  # Convenção do Matlab
            plt.show()
        def Barra_Condicao_Contorno(self, Matriz, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento):
            F = np.zeros((len(Matriz)))               # Alocação de Memória : Força 
            F[GDL_Forca] = Intensidade_Forca        # Contorno da força

            Matriz_Nova = None 
            F_Novo = None

            try:
                if Contorno_Deslocamento == 'Livre - Livre':
                    Matriz_Nova = np.copy(Matriz)
                    F_Novo = np.copy(F)
                elif Contorno_Deslocamento == 'Engastado - Livre':
                    Matriz_Nova = Matriz[1:, 1:]
                    F_Novo  = F[1:]
                elif Contorno_Deslocamento == 'Livre - Engastado':
                    Matriz_Nova = Matriz[:-1, :-1]
                    F_Novo = F[:-1]
                else:
                    raise ValueError("Condição de deslocamento inválida")
            except ValueError as e:
                print(f"Erro: {e}")

            return F_Novo , Matriz_Nova 

        def Matriz_Rigidez_Dinamica_Condensada(self, D, nl, nr, ni, n):


            DLL = D[nl,nl]


            DRR = D[nr,nr]

            DII = np.zeros((len(ni),len(ni)), dtype= complex)
            
            for j in range(len(ni)):
                for i in range(len(ni)):
                    DII[i,j] =  D[ni[j]][ni[i]]

            DIL = np.zeros((len(ni), 1), dtype = complex)


            for j in range(len(ni)):
                DIL[j, nl] =  D[ni[j]][nl]


            DIR = np.zeros((len(ni), 1), dtype = complex)

            for j in range(len(ni)):
                DIR[j,0] =  D[ni[j]][nr]

            DLI = np.zeros((1, len(ni)), dtype = complex)

            for i in range(len(ni)):
                DLI[0,i] =  D[nl][ni[i]]


            DLR = np.zeros((1, 1))

            DLR=  D[nl,nr]

            DRI = np.zeros((1, len(ni)), dtype = complex)

            for i in range(len(ni)):
                DRI[0,i] =  D[nr][ni[i]]

 
            DRL = np.zeros((1, 1))
            DRL=  D[nr][nl]    
            
            inv_DII = np.linalg.inv(DII)
            DLL = DLL - np.dot(DLI, np.dot(inv_DII , DIL))
            DRL = DRL - np.dot(DRI, np.dot(inv_DII , DIL))
            DLR = DLR - np.dot(DLI, np.dot(inv_DII , DIR))
            DRR = DRR - np.dot(DRI, np.dot(inv_DII , DIR))

            
            N = np.block([[np.zeros((1), dtype=complex), np.eye((1), dtype = complex)], [- DRL, -DRR]])
            
            L = np.block([[np.eye((1), dtype = complex), np.zeros((1), dtype= complex)], [DLL, DLR]])

    
            eigenvalue, eigenvector = scipy.linalg.eig(N, L)
            ind = np.argsort(np.abs(eigenvalue))
            ind = ind[1]
            eigenvalue = eigenvalue[ind]
            eigenvector = N * eigenvector
            eigenvector = eigenvector[:,ind]
            eigenvector = eigenvector.reshape(-1, 1)


            return eigenvalue, eigenvector

        def FRF(self, eigenvector, eigenvalue, gdl, pos, Nc, F0, FL):
            Ra = np.ones(1, dtype = complex)
            Ra[pos] = -1
            Ra = np.diag(Ra)
            R = np.block([[Ra, np.zeros(gdl, dtype = complex)], [np.zeros(gdl, dtype = complex), -1 * Ra]])
            phi_q = eigenvector[0, :]
            phi_q = np.array([phi_q])
            phi_F = eigenvector[1, :]
            phi_F = np.array([phi_F])
            phi_star = np.dot(R, eigenvector)
            phi_star_q = phi_star[0,:]
            phi_star_q = np.array([phi_star_q])
            phi_star_F = phi_star[1,:]
            phi_star_F = np.array([phi_star_F])

            inv_phi_F = np.linalg.inv(phi_F)
            inv_phi_star_F = np.linalg.inv(phi_star_F)
            A =  np.block([[np.eye((1), dtype = complex), np.dot(np.dot(inv_phi_F, phi_star_F), eigenvalue ** Nc)], [np.dot(np.dot(inv_phi_star_F, phi_F), eigenvalue ** Nc), np.eye((1), dtype= complex)]])
            Fo = np.block([[- np.dot(inv_phi_F, F0)],[ np.dot(inv_phi_star_F, FL)]])
            inv_A = np.linalg.inv(A)
            Qe = np.dot(inv_A, Fo)
            Q = Qe[0,:]
            Q_star = Qe[1,:]

            qe = np.zeros(([1, Nc + 1]), dtype = complex)
            Fe = np.zeros(([1, Nc + 1]), dtype = complex)
            for ke in range(1, Nc+2):
                qe[:,ke - 1] = np.dot(np.dot(phi_q, eigenvalue ** (ke -1)), Q ) + np.dot(np.dot(phi_star_q, eigenvalue ** (Nc + 1 -ke)), Q_star ) 
                Fe[:,ke - 1] = np.dot(np.dot(phi_F, eigenvalue ** (ke -1)), Q ) + np.dot(np.dot(phi_star_F, eigenvalue ** (Nc + 1 -ke)), Q_star ) 

            return qe, Fe
        def Barra_Final_Response(self, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento, Type, F0, FL):
            F0 = np.array([F0])
            FL = np.array([FL])
            u = np.zeros((1, len(self.Frequency)), dtype=complex)  # Alocação de Memória: Deslocamento
            beta = np.zeros((len(self.Frequency), self.Gdl), dtype= complex)
            for i, omega in tqdm(enumerate(self.Frequency), total=len(self.Frequency),desc="Calculando WFE:" ):
                omega = 2 * np.pi * self.Frequency[i-1]
                Dinamica = self.KG_WFE - (omega ** 2) * self.MG_WFE  # Matriz de rigidez dinâmica
                F, D = self.Barra_Condicao_Contorno(Dinamica, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento)
                Lambda, Phi = self.Matriz_Rigidez_Dinamica_Condensada(D,self.nl, self.nr, self.ni, int(self.Gdl))

                # # Teorema de Floquet- Bloch
                mu = np.log(Lambda) / (-1j)
                beta[i-1, :] = mu

                qe, Fe = self.FRF(Phi, Lambda, self.Gdl, 0, self.Nc,F0, FL )

                try: 
                    if Type == 'Receptancia':
                        u[0, i -1] = qe[0,-1]
                    elif Type == 'Transmitancia':
                        u[0, i -1] = qe[0,-1] / qe[0,0]
                    else:
                        raise ValueError("Tipo de plot inválido!")
                except ValueError as e:
                 print(f"Erro: {e}")

            return  u, beta
                