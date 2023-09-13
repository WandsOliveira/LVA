"""
CONJUNTO DE FUNÇÕES UTILIZADAS..
Autor : Wanderson Vinicius de Oliveira MOnteiro

Obs1: Realizar documentação
Obs2: Mudar nome arquivo

"""

import matplotlib.pyplot as plt

import numpy as np

# Importações necessárias
import matplotlib.pyplot as plt
import numpy as np

class Rod:

    def __init__(self) -> None:
        pass
    
    class FEM:
        def __init__(self, E, rho, eta, b, h, L, Nc, Frequency, Rod_Number, Nos) -> None:
            self.E = E * (1 + 1j * eta)
            self.rho = rho
            self.b = b
            self.h = h
            self.L = L
            self.Nc = Nc
            self.Nos = Nos
            self.Frequency = np.arange(1, Frequency + 1)
            self.Rod_Number = Rod_Number
            self.A = self.b * self.h  # Área
            self.I = self.b * self.h ** 3 / 12  # Momento de Inércia
            self.Nos_FEM = self.Nos
            self.Elementos_FEM = self.Nos_FEM - 1
            self.Gdl = 1
            self.Nos_Principais = 2 + (self.Nc - 1)
            self.Global_FEM = self.Gdl * (self.Nos_Principais + (self.Elementos_FEM - 1) * self.Rod_Number)
            self.GL_Elemento = np.zeros((self.Rod_Number, self.Gdl * (2 + (self.Elementos_FEM - 1))))
            self.L_FEM = self.L / (2 * self.Elementos_FEM)
            self.GDL_Principais = np.arange(0, self.Nos_Principais)
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

            self.KG, self.MG = self.Barra_Matriz_Massa_Rigidez_Global(E, rho, self.L_FEM, self.A, self.Rod_Number, self.GL_Elemento, self.Global_FEM)

        def Barra_Matriz_Massa_Rigidez_Global(self, E, rho, L_FEM, A, Beam_Number, GL_Elemento, Global):
            ## Alocação de memoria matrizes de Massa e Rigidez
            KG = np.zeros((Global, Global), dtype= complex) # Rigidez
            MG = np.zeros((Global, Global), dtype= complex) # Massa

            for i in np.arange(Beam_Number):

                row_aux = GL_Elemento[i,:]
                for ww in range(0, len(row_aux) - 1, 1):
                    row_columns = row_aux[ ww: ww + 2]
                    # Matrizes locais
                    Ke = (E[0] * A / L_FEM) * np.array([[1/2, -1/2], [-1/2, 1/2]])  # Matriz de rigidez local
                    Me = rho[0] * A * L_FEM * np.array([[2/3, 1/3], [1/3, 2/3]])  # Matriz de Massa loca

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

        def Plot_FRF(self, Frequencia, u, Type):
            plt.figure( figsize= (12, 4))
            try: 
                if Type == 'Receptancia':
                    plt.plot(Frequencia, 20* np.log10(np.abs(u[1, :])), color='red', linestyle='-', label = ' FEM')
                    plt.ylabel('Receptância dB ')
                    plt.xlabel('Frequência (Hz)')
                    plt.grid(True)
                    plt.legend()
                    plt.show()
                elif Type == 'Transmitancia':
                    plt.plot(Frequencia, 20* np.log10(np.abs(u[1, :]) / np.abs(u[0, :])), color='red', linestyle='-', label = ' FEM')
                    plt.ylabel('Transmitancia dB ')
                    plt.xlabel('Frequência (Hz)')
                    plt.grid(True)
                    plt.legend()
                    plt.show()
                else:
                    raise ValueError("Tipo de plot inválido!")
            except ValueError as e:
                print(f"Erro: {e}")

        def Barra_Final_Response(self, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento, Type):

            u = np.zeros((len(self.KG), len(self.Frequency)), dtype=complex)  # Alocação de Memória: Deslocamento

            for i, omega in enumerate(self.Frequency):
                omega = 2 * np.pi * self.Frequency[i-1]
                Dinamica = self.KG - (omega ** 2) * self.MG  # Matriz de rigidez dinâmica

                F, D = self.Barra_Condicao_Contorno(Dinamica, GDL_Forca, Intensidade_Forca, Contorno_Deslocamento)
                u[:, i-1] = np.linalg.solve(D, F)

            return  self.Plot_FRF(self.Frequency, u, Type)  # Chama a função Plot_FRF para plotar a receptância



        