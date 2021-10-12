#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 09:04:11 2021

@author: marcio
"""

from numpy import pi, exp, sqrt, arccos
import readline
import sys
readline.parse_and_bind('tab: complete')
atomos = [] ## lista para adicionar os tipos de átomos
PARAMETERS = [] ## lista para os cell_parameters
Simb = [] ## lista para os simbolos dos atomos do Atomic_positions
X = [] ## coordenada X do Atomic_positions
Y = [] ## coordenada Y do Atomic_positions
Z = [] ## coordenada Z  do Atomic_positions
nat = 0 ## número de átomos na célula ## variável global (tem o valor modificado)
ntp = 0 ## número de tipos de átomos ## variável global (tem o valor modificado)
#-------------------início da leitura do arquivo vcrelax.out------------------
try:
    with open(sys.argv[1], 'r') as vcrelax:
        for lines in vcrelax:
            if 'number of atoms/cell' in lines:
                nat += int(lines[40:45]) ## modifica o valor da variável nat
                ntp += int((vcrelax.readline().rstrip().split('='))[1]) # modifica ntp
            if 'Begin final coordinates' in lines:
                vcrelax.readline(), vcrelax.readline(), vcrelax.readline(), vcrelax.readline() ## remove 4 linhas abaixo de 'Begin final coordinates'
                for i in range(0,3):
                    PARAMETERS.append(vcrelax.readline().split()) # matriz de vetores
                break
        for lines in vcrelax: # lê o arquivo a partir da linha 'Begin final coordinates'
            if 'ATOMIC_POSITIONS' in lines:
                break
        try:
            for lines in vcrelax:
                dados = lines.split()
                Simb.append(dados[0]) # armazena os símbolos dos átomos
                X.append(float(dados[1])) # armazena as coordenadas X
                Y.append(float(dados[2])) # armazena as coordenadas Y
                Z.append(float(dados[3])) # armazena as coordenadas Z
        except:
            Simb.pop() #elimina caracter indesejado (variável de controle)
            pass # interrompe o laço for e dá sequência no código
        for lines in vcrelax:
            if 'valence' in lines: ## busca os tipos de átomos
                for i in range(0, ntp):
                    atomos.append(str((vcrelax.readline().rstrip().split())[0]))
##--------------término da leitura do arquivo vcrelax.out---------------------
#---------------contar o número de repetições dos atomos----------------------
    repeticao = [] ## lista para armazenar o número de repetição de cada tipo de átomo
    for i in atomos:
        repeticao.append(Simb.count(i))
    atoms = " ".join(map(str,atomos)) ## concatena a lista atomos
    cont ="  ".join(map(str,repeticao)) ## concatena a lista repetição

##----------------matriz de vetores da célula---------------------------------
    vx1, vy1, vz1 = float((PARAMETERS[0])[0]), float((PARAMETERS[0])[1]), float((PARAMETERS[0])[2])
    vx2, vy2, vz2 = float((PARAMETERS[1])[0]), float((PARAMETERS[1])[1]), float((PARAMETERS[1])[2])
    vx3, vy3, vz3 = float((PARAMETERS[2])[0]), float((PARAMETERS[2])[1]), float((PARAMETERS[2])[2])
##-----------------------determinante da matriz-------------------------------
    det = (vx1*vy2*vz3 + vx2*vy3*vz1 + vx3*vy1*vz2 - vz1*vy2*vx3 - vz2*vy3*vx1 - vz3*vy1*vx2)
##----------------------paramentros de rede-----------------------------------
    r1 = sqrt(vx1**2 + vy1**2 + vz1**2)
    r2 =  sqrt(vx2**2 + vy2**2 + vz2**2)
    r3 =  sqrt(vx3**2 + vy3**2 + vz3**2)
##-------------------------angulos--------------------------------------------
    alpha = (180/pi)*arccos((vx2*vx3+vy2*vy3+vz2*vz3)/(r2*r3))
    beta = (180/pi)*arccos((vx1*vx3+vy1*vy3+vz1*vz3)/(r1*r3))
    gamma = (180/pi)*arccos((vx1*vx2+vy1*vy2+vz1*vz2)/(r1*r2))
#-----------------------Calculando os cofatores-------------------------------
    t11 = (vy2*vz3 - vz2*vy3)/det
    t21 = -(vy1*vz3 - vz1*vy3)/det
    t31 = (vy1*vz2 - vz1*vy2)/det
    t12 = -(vx2*vz3 - vz2*vx3)/det
    t22 = (vx1*vz3 - vz1*vx3)/det
    t32 = -(vx1*vz2 - vz1*vx2)/det
    t13 = (vx2*vy3 - vy2*vx3)/det
    t23 = -(vx1*vy3 - vy1*vx3)/det
    t33 = (vx1*vy2 - vy1*vx2)/det
##------------------calcular as novas coordenadas ----------------------------
    n_coord = [] ## lista para armazenar as novas coordenadas
    for at, i, j, k in zip (Simb, X, Y, Z):
        a = (t11*i + t12*j + t13*k)
        b = (t21*i + t22*j + t23*k)
        c = (t31*i + t32*j + t33*k)
        texto = str(f'{at:<2} {(atomos.index(at)+1):^10d} {a:^10.6f} {b:^10.6f} {c:^10.6f}'+'\n')
        n_coord.append(texto)
    positions ="".join(map(str,n_coord)) # cria uma string com o Simbolo, ID e cordenadas X Y Z
#-----------------------------------------------------------------------------
    with open("cryst.res", 'w') as saida: ## cria um arquivo .res
        saida.write("TITL insert-title\n")
        saida.write(f"CELL  0.71073 {r1:.6f} {r2:.6f} {r3:.6f} {alpha:.6f}  {beta:.6f} {gamma:.6f}\n")
        saida.write("LATT  -1\n")
        saida.write(f"SFAC {atoms}\n")
        saida.write(f"UNIT\t{cont}\n")
        saida.write("MERG 2\nOMIT 0\nFMAP 2\nPLAN 5\nBOND\nL.S. 10\nWGHT 0.0555  0.0000\n")
        saida.write(f"{positions}\n")
        saida.write("HKLF 4\nEND")
    print(50*'=')
    print("Dados salvos como: cryst.res")
except:
    print("""Cryst-1.0:
USO:
    [cryst.py] [arquivo .out do cálculo vc-relax]
    """)
