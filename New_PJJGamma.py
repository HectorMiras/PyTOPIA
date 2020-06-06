# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 08:32:39 2020

@author: josem
"""

import tkinter as tk
from tkinter import messagebox
from tkinter.filedialog import askopenfilename, asksaveasfilename
from tkinter.simpledialog import askstring
import numpy as np
import cv2
from PIL import Image, ImageTk
import sys
import math
import os.path as path
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import rotate
import FilmDoseClass

fichero_dosimetria = "Dosimetria.txt"                       # Fichero con las curvas de calibracion
modelos_ajuste_ids = ["Polinomial_3Grado" , "Percolacion"]  # Ids para los modelos de ajuste dosis=f(DoN). Sólo tif
foco = -1                                                   # Variable para controlar el foco en la ventena principal
flag_velo_1 = False                                         # Variable para establecer el velo en la imagen izquierda. Sólo tif
flag_velo_2 = False                                         # Variable para establecer el velo en la imagen derecha. Sólo tif
flag_dosis_conocida_1 = False                               # Variable para establecer la dosis conocida en pelculas tif
flag_dosis_conocida_2 = False                               # Variable para establecer la dosis conocida en pelculas tif en canvas derecha
rectangulo_velo_200 = None                                  # Rectangulo a dibujar y/o borrar en el establecimiento del velo o de la dosis conocida
flag_perfil = False                                         # Variable para controlar los perfiles de la imagen seleccionada
flag_recorte = False                                        # Variable para controlar recortes d la imagen seleccionada
_start = -1                                                 # Variable para controlar posición inicial de recorte
_end = -1                                                   # Variable para controlar posición final de recorte
flag_recorte_ambas = False                                  # Variable para controlar recorte en ambas imágenes si están calibradas

#
# Clase para pasar de Do neta a valores de pixel
#
class DoN_to_Pixel:
    def __init__(self, DoN, velo):
        self.pixel = np.zeros((3, len(DoN[0]), len(DoN[0][0])))
        for i in range(3):
            self.pixel[i] = velo[i]*(np.power(10.0, DoN[i]))

#
# Obtiene la matriz gamma de dos matrices 2D.
# Ambas deben ser del mismo tamaño y tener el mismo factor_pixel_mm
# colores y rango deben tener el mismo numero de elementos y estar ordenados de menor a mayor
#
class Get_gamma:
    def __init__(self, matrizA, matrizB, dosis_tol, radio_pixel_tol, colores, rango, nombre_color):
        width = len(matrizA[0])
        height = len(matrizA)
    
        self.pasa_tantos_por_cien = 0
        self.no_pasa_tantos_por_cien = 0
        
        self.valores_gamma = np.zeros((height, width))
        for i in range(0, height):
            for j in range(0, width):
                x1 = j - radio_pixel_tol
                x2 = j + radio_pixel_tol
                y1 = i - radio_pixel_tol
                y2 = i + radio_pixel_tol
                x1 = max(0, x1)
                x2 = min(x2, width)
                y1 = max(0, y1)
                y2 = min(y2, height)
                index_gamma = -1                
                index_gamma = self.Get_gamma_recorte(i, j, matrizA[i][j], matrizB, x1, y1, x2, y2, dosis_tol, 
                                                     radio_pixel_tol, rango)
                self.valores_gamma[i][j] = index_gamma
                  
        self.imagen_gamma = np.zeros((len(self.valores_gamma), len(self.valores_gamma[0]), 3))
        contador_pasa = 0
        contador_no_pasa = 0
        for i in range(len(self.valores_gamma)):
            for j in range(len(self.valores_gamma[0])):
                indice = self.where_is(self.valores_gamma[i][j], rango)                
                if self.valores_gamma[i][j] <= 1.0:
                    contador_pasa +=1
                else:
                    contador_no_pasa +=1
                
                self.imagen_gamma[i][j][0] = colores[indice][0]
                self.imagen_gamma[i][j][1] = colores[indice][1]
                self.imagen_gamma[i][j][2] = colores[indice][2]
        
        self.pasa_tantos_por_cien = 100.0 * contador_pasa/(contador_pasa + contador_no_pasa)
        self.no_pasa_tantos_por_cien = 100.0 * contador_no_pasa/(contador_pasa + contador_no_pasa)
    
    def where_is(self, valor_gamma, rango):
        if valor_gamma <= rango[0]:
            return 0
        if valor_gamma >= rango[len(rango)-1]:
            return len(rango)
        
        for i in range(len(rango)):
            if valor_gamma <= rango[i]:
                where = i
                break
        where = i
        return where

    def distancia_puntos(self, x1, y1, x2, y2):
        x = x2 - x1
        y = y2 - y1
        return abs(np.sqrt(np.power(x,2) + np.power(y, 2)))
    
    def Get_gamma_recorte(self, i, j, valor_ref, matrizB, x1, y1, x2, y2, dosis_tol, radio_pixel_tol, rango):
        valores_gamma = []
        for n in range (y1, y2):
            for m in range(x1, x2):
                dst = self.distancia_puntos(i, j, n, m)
                if (dst <= radio_pixel_tol):
                    dif_d = abs(valor_ref - matrizB[n][m])
                    valor = math.sqrt(pow(dst/radio_pixel_tol, 2) + pow(dif_d/dosis_tol, 2))
                    valores_gamma.append(valor)
                    #if valor <= rango[0]:
                    if valor <= 1.0: # Más rápido, menos detalle
                        return valor
        indice_gamma = min(valores_gamma)
        return indice_gamma

#
# Clase que obtiene la matriz diferencia de dos matrices 2D de datos
# Considera el signo tomando como referencia la matrizA
# Ambas deben ser del mismo tamaño y tener el mismo factor_pixel_mm
#
class Get_diferencia:
    def __init__(self, matrizA, matrizB):
        self.matriz_diferencia = matrizB - matrizA

#
# Clase que obtiene los perfiles para la ecualizacion
#
class Perfil_ecualizacion:
    def __init__(self, ccx, ccy, desplazamiento_lateral, desplazamiento_vertical, angulo_giro, objeto_ppal_1, objeto_ppal_2, canal, factor_dosis_2):
        perfil1 = objeto_ppal_1.get_dosis(3)
        
        if ccx < 0:
            ccx = 0
        if ccy < 0:
            ccy = 0
        if ccx >= objeto_ppal_2.get_width():
            ccx = objeto_ppal_2.get_width()-1
        if ccy >= objeto_ppal_2.get_height():
            ccy = objeto_ppal_2.get_height()-1
            
        buffer_y1 = perfil1[canal, int(ccy), :]
        buffer_x1 = perfil1[canal, :, int(ccx)]
        
        perfil2_dust = objeto_ppal_2.get_dosis(3)
        # Aplica las transformaciones al canal de la matriz de dosis del objeto ppal 2
        M = cv2.getRotationMatrix2D((int(objeto_ppal_2.get_width()/2), int(objeto_ppal_2.get_height()/2)), float(angulo_giro), 1.0)
        buffer = perfil2_dust[canal]
        rows = len(buffer)
        cols = len(buffer[0])
        dstg = cv2.warpAffine(buffer,M,(cols,rows))
      
        rows = len(dstg)
        cols = len(dstg[0])
        M = np.float32([[1,0,int(desplazamiento_lateral)],[0,1,int(desplazamiento_vertical)]])
        
        dstt = cv2.warpAffine(dstg, M, (cols,rows))
        
        perfil2 = dstt
        perfil2 *= float(factor_dosis_2)
        buffer_y2 = perfil2[int(ccy), :]      
        buffer_x2 = perfil2[:, int(ccx)]
        self.px_1= buffer_x1
        self.py_1= buffer_y1
        self.px_2= buffer_x2
        self.py_2= buffer_y2

#
# Clase para rotar una matriz 3D tipo imagen un angulo respecto de un cualquier punto
# Calcula una matriz 3D tipo datos
#
class Rotar_respecto_any_point:
    def __init__(self, matriz3D_im, cx, cy, angulo_deg, factor_escala):
        rows = len(matriz3D_im)
        cols = len(matriz3D_im[0])
        M = cv2.getRotationMatrix2D((cx,cy), angulo_deg, factor_escala)
        buffer = matriz3D_im[:, :, 0]
        self.dst = np.zeros((3, len(buffer), len(buffer[0])))
        dst_ = cv2.warpAffine(buffer,M,(cols,rows))
        self.dst[0] = dst_
        for i in range(1, 3):
            buffer = matriz3D_im[:, :, i]
            dst_ = cv2.warpAffine(buffer,M,(cols,rows))
            self.dst[i] = dst_
#
# Clase para superponer visualmente dos imagenes
#
class Alfa:
    def __init__(self, matriz_3D_1, matriz_3D_2, alfa):
        beta = 1.0 - alfa
        self.dst = cv2.addWeighted(matriz_3D_1, alfa, matriz_3D_2, beta, 0, dtype = cv2.CV_32F)

#
# Clase para dibujar una isodosis.
# matriz_datos es una matriz3D de datos
# canal es el canal de color a evaluar (0, 1 ó 2 -> R, G o B)
# valor es la isodosis
# ventana es el margen de diferencia tolerable para dibujarla
# color es el color de la isodosis
# fondo permite escrbir sobre blanco o utilizar la imagen como fondo
# 
# imagen_isodosis es una matriz 3D tipo imagen para presentar
class Get_isodosis:
    def __init__(self, matriz_datos, canal, valor, ventana, color, fondo):
        height = len(matriz_datos[0])
        width = len(matriz_datos[0][0])
        buffer = np.zeros((height, width, 3))
            
        if fondo == True:
            for i in range(0, height):
                for j in range(0, width):
                    dif = abs(matriz_datos[canal][i][j] - valor)
                    if dif <= ventana:
                        buffer[i][j][2] = color[2]
                        buffer[i][j][1] = color[1]
                        buffer[i][j][0] = color[0]
                    else:
                        buffer[i][j][2] = matriz_datos[0][i][j]
                        buffer[i][j][1] = matriz_datos[1][i][j]
                        buffer[i][j][0] = matriz_datos[2][i][j]
            self.imagen_isodosis = buffer
        else:
            for i in range(0, height):
                for j in range(0, width):
                    dif = abs(matriz_datos[canal][i][j] - valor)
                    if dif <= ventana:
                        buffer[i][j][2] = color[2]
                        buffer[i][j][1] = color[1]
                        buffer[i][j][0] = color[0]
                    else:
                        buffer[i][j][2] = 255
                        buffer[i][j][1] = 255
                        buffer[i][j][0] = 255
            self.imagen_isodosis = buffer
            
#
# Clase para dibujar una isodosis.
# matriz_datos es una matriz2D de datos
# valor es la isodosis
# ventana es el margen de diferencia tolerable para dibujarla
# color es el color de la isodosis
# fondo permite escrbir sobre blanco o utilizar la imagen como fondo
# 
# imagen_isodosis es una matriz 3D tipo imagen para presentar
class Get_isodosis_2D:
    def __init__(self, matriz_datos, valor, ventana, color, fondo):
        height = len(matriz_datos)
        width = len(matriz_datos[0])
        buffer = np.zeros((height, width, 3))
            
        if fondo == True:
            for i in range(0, height):
                for j in range(0, width):
                    dif = abs(matriz_datos[i][j] - valor)
                    if dif <= ventana:
                        buffer[i][j][2] = color[2]
                        buffer[i][j][1] = color[1]
                        buffer[i][j][0] = color[0]
                    else:
                        buffer[i][j][2] = matriz_datos[i][j]
                        buffer[i][j][1] = matriz_datos[i][j]
                        buffer[i][j][0] = matriz_datos[i][j]
            self.imagen_isodosis = buffer
        else:
            for i in range(0, height):
                for j in range(0, width):
                    dif = abs(matriz_datos[i][j] - valor)
                    if dif <= ventana:
                        buffer[i][j][2] = color[2]
                        buffer[i][j][1] = color[1]
                        buffer[i][j][0] = color[0]
                    else:
                        buffer[i][j][2] = 255
                        buffer[i][j][1] = 255
                        buffer[i][j][0] = 255
            self.imagen_isodosis = buffer

#
# Clase para redimensionar una matriz 3D de datos a new_height x new_width
#
class Resize:
    def __init__(self, matriz, new_height, new_width,):
        self.dst = np.zeros((3, int(new_height), int(new_width)))
        for i in range(3):
            self.dst[i] = cv2.resize(matriz[i],(int(new_width), int(new_height)), interpolation = cv2.INTER_CUBIC)

#
# Clase para crear los lienzos de trabajo previo a la ecualizacion de ambas imagenes.
#
# El lienzo es una imagen del tamaño y factor_pixel_mm necesario para que ambas imagenes presenten el mismo numero de pixeles.
# La anchura del lienzo será la mayor en mm de las dos imagenes.
# La altura del lienzo será la mayor en mm de las dos imagenes.
# El factor_pixel_mm será el que presente mayor resolución.
#
class Lienzo:
    def __init__(self, objeto_ppal_1, objeto_ppal_2):
        lw1_pixel = objeto_ppal_1.get_width()
        lh1_pixel = objeto_ppal_1.get_height()
        lw2_pixel = objeto_ppal_2.get_width()
        lh2_pixel = objeto_ppal_2.get_height()
        f1 = objeto_ppal_1.get_factor_pixel_mm()
        f2 = objeto_ppal_2.get_factor_pixel_mm()
        
        lw1_mm = lw1_pixel/f1
        lh1_mm = lh1_pixel/f1
        lw2_mm = lw2_pixel/f2
        lh2_mm = lh2_pixel/f2
        
        lw = max(lw1_mm, lw2_mm)
        lh = max(lh1_mm, lh2_mm)
        self.f = max(f1, f2)
        
        self.lw_pixel_new = lw*self.f
        self.lh_pixel_new = lh*self.f
        
        # Ambas imagenes deben tener informacion dosimetrica
        if objeto_ppal_1.get_tipo() == "pnc":
            dosis_1 = objeto_ppal_1.get_dosis(3)
            if objeto_ppal_1.get_width() == self.lw_pixel_new and objeto_ppal_1.get_factor_pixel_mm == self.f and objeto_ppal_1.get_height() == self.lh_pixel_new:
                self.dosis_1_new = dosis_1
            else:
                self.dosis_1_new = Resize(dosis_1, self.lh_pixel_new, self.lw_pixel_new).dst
        else:
            pixel_1 = objeto_ppal_1.get_pixel(3)
            don_1 = objeto_ppal_1.get_don(3)
            dosis_1 = objeto_ppal_1.get_dosis(3)
            if objeto_ppal_1.get_width() == self.lw_pixel_new and objeto_ppal_1.get_factor_pixel_mm == self.f and objeto_ppal_1.get_height() == self.lh_pixel_new:
                self.pixel_1_new = pixel_1
                self.don_1_new = don_1
                self.dosis_1_new = dosis_1
            else:
                self.pixel_1_new = Resize(pixel_1, self.lh_pixel_new, self.lw_pixel_new).dst
                self.don_1_new = Resize(don_1, self.lh_pixel_new, self.lw_pixel_new).dst
                self.dosis_1_new = Resize(dosis_1, self.lh_pixel_new, self.lw_pixel_new).dst
            
        if objeto_ppal_2.get_tipo() == "pnc":
            dosis_2 = objeto_ppal_2.get_dosis(3)
            if objeto_ppal_2.get_width() == self.lw_pixel_new and objeto_ppal_2.get_factor_pixel_mm == self.f and objeto_ppal_2.get_height() == self.lh_pixel_new:
                self.dosis_2_new = dosis_2
            else:
                self.dosis_2_new = Resize(dosis_2, self.lh_pixel_new, self.lw_pixel_new).dst
        else:
            pixel_2 = objeto_ppal_2.get_pixel(3)
            don_2 = objeto_ppal_2.get_don(3)
            dosis_2 = objeto_ppal_2.get_dosis(3)
            if objeto_ppal_2.get_width() == self.lw_pixel_new and objeto_ppal_2.get_factor_pixel_mm == self.f and objeto_ppal_2.get_height() == self.lh_pixel_new:
                self.pixel_2_new = pixel_2
                self.don_2_new = don_2
                self.dosis_2_new = dosis_2
            else:
                self.pixel_2_new = Resize(pixel_2, self.lh_pixel_new, self.lw_pixel_new).dst
                self.don_2_new = Resize(don_2, self.lh_pixel_new, self.lw_pixel_new).dst
                self.dosis_2_new = Resize(dosis_2, self.lh_pixel_new, self.lw_pixel_new).dst

#
# Clase para obtener perfiles de la imagen seleccionada.
# matriz es una matriz 2D de datos.
# index es la coordenada donde tomar el perfil, "x" en los verticales e "y" en los horizontales
# tipo es el tipo de perfil, 0 para los horizontales y 1 para verticales
# exclusion para dejar un margen en el perfil, 0 toma el perfil completo
#
class Get_profile:
    def __init__(self, matriz, index, tipo, exclusion):
            self.profile_y = []
            self.profile_x = []
            dim_1 = len(matriz) # Numero de filas
            dim_2 = len(matriz[0]) # Numero de columnas
    
            if tipo == 0:   # Perfil Horizontal
                zona_i = exclusion
                zona_f = dim_2 - exclusion
                for i in range(zona_i, zona_f):
                    self.profile_y.append(matriz[index][i])
                    self.profile_x.append(i)
            else:           # Perfil Vertical
                zona_i = exclusion
                zona_f = dim_1 - exclusion
                for i in range(zona_i, zona_f):
                    self.profile_y.append(matriz[i][index])
                    self.profile_x.append(i)

#
# Clase que escribe a partir de una matriz2D, un fichero en formato pnc
# factor_pixel_mm es el numero de pixeles por cada milimetro. En funcion 
# del valor maximo de la matriz leerá 8 ó 16 bits por canal para representar.
# 
class Escribe_PNC:
    def __init__(self, objeto_ppal, canal, fichero_salida):
        extension = Get_extension(fichero_salida).extension
        if extension != "pnc":
            fichero_salida += ".pnc"
        factor_pixel_mm = objeto_ppal.get_factor_pixel_mm()
        matriz2D = objeto_ppal.get_dosis(int(canal))
        width =len(matriz2D[0])
        height = len(matriz2D)
        fichero_out = open(fichero_salida, "w")
        fichero_out.write(fichero_salida + "\n")
        fichero_out.write("Version:,0.0\n\n")

        for i in range(0, width):
            fichero_out.write("," + str(i/factor_pixel_mm/10))
        fichero_out.write(",\n")
        for i in range(0, height):
            fichero_out.write(str(10*i/factor_pixel_mm/10) + ",")
            for j in range(0, width):
                fichero_out.write(str(matriz2D[i][j]) + ",")
            fichero_out.write("\n")   
        fichero_out.close()

#
# Clase que corta una matriz3D de datos
#
class Cut:
    def __init__(self, matriz, si_x, si_y, id_x, id_y):
        buffer = matriz[0]       
        dst_ = buffer[si_y:id_y, si_x:id_x]
        self.dst = np.zeros((3, len(dst_), len(dst_[0])))
        self.dst[0] = dst_
        for i in range(1, 3):
            buffer = matriz[i]       
            dst_ = buffer[si_y:id_y, si_x:id_x]
            self.dst[i] = dst_

#
# Clase que gira una matriz 3D de datos respacto del centro de la imagen sin perdida de informacion por la
# periferia. Se utiliza en transformaciones individuales que no impliquen a ambas imagenes. En caso contrario
# seria necesario volver a ecualizar.
# angulo_def >0 -> ccw
# angulo_deg <0 -> cw
class Rotar_centro:
    def __init__(self, matriz, angulo_deg):
        buffer = rotate(matriz[0], angulo_deg)
        self.dst = np.zeros((3, len(buffer), len(buffer[0])))
        for i in range(1, 3):
            self.dst[i] = rotate(matriz[i], angulo_deg)
        self.dst[0] = buffer

#
# Clase para desplazar horizontal y/o verticalmente una cantidad de pixeles una matriz de datos 3D
# dx -> desplazamiento horizontal
# dy -> desplazamiento vertical
#
class Traslacion:
    def __init__(self, matriz, dx, dy):
        rows = len(matriz[0])
        cols = len(matriz[0][0])
        self.dst = np.zeros((3, rows, cols))
        M = np.float32([[1,0,dx],[0,1,dy]])
        
        for i in range(3):
            dst = cv2.warpAffine(matriz[i], M, (cols,rows))
            self.dst[i] = dst

#
# Clase que rota 90 o 270 una matriz de datos 3D, respecto de su centro
# sentido -> 1 CCW
# sentido -> -1 CW
#
class Rota_90:
    def __init__(self, matriz, sentido):
        self.dst = np.zeros((3, len(matriz[0][0]), len(matriz[0])))
        for i in range(3):
            self.dst[i] = np.rot90(matriz[i], sentido)

#
# Clase que hace un flip de una matriz de datos 3D
#
# eje -> 0 Flip respecto de un eje horizontal
# eje -> 1 Flip respecto de un eje vertical
#
class Flip:
    def __init__(self, matriz, eje):
        self.dst = np.zeros((len(matriz), len(matriz[0]), len(matriz[0][0])))
        for i in range(3):
            self.dst[i] = np.flip(matriz[i], eje)

#
# Clase que obtiene los valores de DoN para los tres canales a partir de un unico valor de pixel de tres canales
# pixel es un array de tres valores correspondientes a cada canal
# pixel_velo tiene dimension 3
#
class Pixel3canales_to_DoN:
    def __init__(self, pixel_velo, pixel):
        self.DoN = np.zeros((3))
        for k in range(0, 3):
            if pixel[k] <= 0: # Para salvar la discontinuidad del logaritmo
                pixel[k] = 1
            buffer = -math.log10(pixel[k]) + math.log10(pixel_velo[k])
            if buffer < 0: # Para valores de pixeles menores al del velo
                buffer = 0
            self.DoN[k] = buffer

#
# Clase que permite obtener la densidad optica neta asociada con un valor de dosis
#
class Dosis_to_DoN:
    def __init__(self, lc, dosis):
        incremento_don = 0.0001
        self.DoN = np.zeros((3))
        
        # Ver que la dosis se encuentre en la curva de calibracion
        for k in range(0, 3):
            minima_dosis = min(lc.dosis[k])
            maxima_dosis = max(lc.dosis[k])
            #print("Canal: " + str(k))
            #print("Rango dosis: " + str(minima_dosis) + " <--> " + str(maxima_dosis))
            #print("Buscando dosis: " + str(dosis))
            if dosis >= minima_dosis and dosis <= maxima_dosis:
                #messagebox.showinfo(message="Rango de dosis para el canal " + str(k) + ": " + str(minima_dosis) + " <--> " + str(maxima_dosis) + "!!!", title="Aviso!!!")
                # Busca la DoN para esa dosis
                if lc.modelo == 0:
                    DoNi = 0
                    while True:
                        dosisi = lc.aj.Pol_3(DoNi, lc.aj.A[k], lc.aj.B[k], lc.aj.C[k], lc.aj.D[k])
                        if dosisi < dosis:
                            DoNi += incremento_don
                        else:
                            self.DoN[k] = DoNi
                            break
                if lc.modelo == 1:
                    DoNi = 0
                    while True:
                        dosisi = lc.aj.Percolacion(DoNi, lc.aj.A[k], lc.aj.B[k], lc.aj.C[k], lc.aj.D[k])
                        if dosisi < dosis:
                            DoNi += incremento_don
                        else:
                            self.DoN[k] = DoNi
                            break
            else:
                messagebox.showinfo(message=str(dosis) + " fuera del rango dosimetrico de la curva de calibracion " + str(lc.id[k]) + "!!!", title="Aviso!!!")
                sys.exit()
#
# Clase para linealizar la matriz de DoN de la imagen a la curva de calibración
#
class Get_K:
    def __init__(self, objeto_ppal):
        self.K = np.zeros((3))
        if objeto_ppal.get_tipo != "pnc" and objeto_ppal.get_cal(0) != -1 and objeto_ppal.get_cal(1) != -1 and objeto_ppal.get_cal(2) != -1 and objeto_ppal.get_modelo_ajuste() != -1: # tif calibrada
            if objeto_ppal.get_pixel_velo(0) != -1 and objeto_ppal.get_pixel_dosis_conocida(0) != -1 and objeto_ppal.get_pixel_velo(1) != -1 and objeto_ppal.get_pixel_dosis_conocida(1) != -1 and objeto_ppal.get_pixel_velo(2) != -1 and objeto_ppal.get_pixel_dosis_conocida(2) != -1 : # Velo y dosis conocida
                lc = Lee_curva_calibracion(objeto_ppal.get_cal(0), objeto_ppal.get_modelo_ajuste(), False, False, False)
                # Busca en el ajuste de la curva de calibracion los valores de pixel de la dosis conocida y del velo
                Don_dosis_conocida_ajuste_curva = Dosis_to_DoN(lc, objeto_ppal.get_dosis_linealizacion()).DoN
                pixel_velo_imagen = objeto_ppal.get_pixel_velo(3)
                pixel_dosis_conocida_imagen = objeto_ppal.get_pixel_dosis_conocida(3)
                Don_dosis_conocida_imagen = Pixel3canales_to_DoN(objeto_ppal.get_pixel_velo(3), pixel_dosis_conocida_imagen).DoN
                for k in range(0, 3):
                    self.K[k] = Don_dosis_conocida_ajuste_curva[k]/Don_dosis_conocida_imagen[k]        
            else:
                messagebox.showinfo(message="Establecer pixel de velo y dosis conocida!!!!", title="Aviso!!!")
        else:
             messagebox.showinfo(message="Solo tifs calibradas!!!!", title="Aviso!!!")

#
# Clase para obtener la matriz 2D de dosis a partir de la matriz 2D de DoN con el modelo de ajuste y los datos de la curva asociada
# al canal correspondiente
class DoN2D_to_Dosis:
    def __init__(self, canal, DoN, lc):
        
        self.dosis = np.zeros((len(DoN), len(DoN[0])))
        if lc.modelo == 0:
            for i in range(len(DoN)):
                for j in range(len(DoN[0])):
                        self.dosis[i][j] = lc.aj.Pol_3(DoN[i][j], lc.aj.A[canal], lc.aj.B[canal], lc.aj.C[canal], lc.aj.D[canal])
        if lc.modelo == 1:
            for i in range(len(DoN)):
                for j in range(len(DoN[0])):
                        self.dosis[i][j] = lc.aj.Percolacion(DoN[i][j], lc.aj.A[canal], lc.aj.B[canal], lc.aj.C[canal], lc.aj.D[canal])

#
# Clase para obtener la matriz 3D de dosis a apartir de la matriz 3D de DoN obtenida del objeto principal
# Como entrada necesita el objeto lee_curva_calibracion para obtener los coeficientes del ajuste y modelo para
# conocer que funcion de ajuste utilizar
#
class DoN3D_to_Dosis:
    def __init__(self, lc, objeto_ppal):
        self.Dosis = np.zeros((3, len(objeto_ppal.get_don(3)[0]), len(objeto_ppal.get_don(3)[0][0])))
        for k in range(3):
            self.Dosis[k] = DoN2D_to_Dosis(k, objeto_ppal.get_don(k), lc).dosis

#
# Clase de ajuste de la curva de puntos (DoN, dosis) al modelo elegido
# Ajusta los tres canales
#
class Ajuste:
    def __init__(self, DoN, dosis, modelo):
        if modelo == 0:
            self.A = np.zeros((3))
            self.B = np.zeros((3))
            self.C = np.zeros((3))
            self.D = np.zeros((3))
            for k in range(3):
                (self.A[k], self.B[k], self.C[k], self.D[k]), _ = curve_fit(self.Pol_3, DoN[k], dosis[k])
                #messagebox.showinfo(message=str(round(self.A[k], 3)) + " " + str(round(self.B[k], 3)) + " " + str(round(self.C[k], 3)) + " " + str(round(self.D[k], 3)), title="Aviso!!!")
                
        if modelo == 1:
            self.A = np.zeros((3))
            self.B = np.zeros((3))
            self.C = np.zeros((3))
            self.D = np.zeros((3))
            for k in range(3):
                (self.A[k], self.B[k], self.C[k], self.D[k]), _ = curve_fit(self.Percolacion, DoN[k], dosis[k], 
                                                                            bounds =((0, 5, 100, 0), (5, 20, 1000, 1)))
                #messagebox.showinfo(message=str(round(self.A[k], 3)) + " " + str(round(self.B[k], 3)) + " " + str(round(self.C[k], 3)) + " " + str(round(self.D[k], 3)), title="Aviso!!!")
                
    def Pol_3(self, x, A, B, C, D):
        return A + B*x + C*x**2 + D*x**3
    
    def Percolacion(self, x, A, B, C, D):
        return abs(((B/(A-x)))**(1/D))-C

#
# Clase para leer el fichero de dosimetria y buscar las curvas de calibracion del objeto_ppal.
# El fichero de dosimetria tiene un formato rigido, existiendo una linea por cada curva de calibracion
# En cada linea de calibracion, el primer punto (dosis, pixel) debe corresponder al velo
# El argumento indexr es el indice de la curva correspondiente asociada al canal rojo.
# Las asociadas al canal verde y al azul son las siguientes en orden
# Calcula la matriz de Densidad optica neta de los tres canales (DoN)
# Ajusta la curva (Dosis, DoN) segun el modelo solicitado y devuelve los coeficientes del ajuste
#   0 -> Polinomio de grado tres para cada canal
#   1 -> Percolacion
# Si flag_r == True, muestra la curva del canal rojo
# Si flag_g == True, muestra la curva del canal verde
# Si flag_b == True, muestra la curva del canal azul
#
class Lee_curva_calibracion:
    def __init__(self, indexr, modelo, flag_r, flag_g, flag_b):
        global fichero_dosimetria
        global modelos_ajuste_ids
        
        self.modelo = modelo
        # Busca la curva de calibracio "indexr"
        f = open(fichero_dosimetria)
        contador = 0
        while True:
            line = f.readline()
            if contador == indexr or line =="":
                break
            else:
                contador +=1
            
        self.id = []
        
        for k in range(3):
            data = line.split("#")
            self.id.append(data[1])
            if k == 0:
                 num_puntos = int(data[3])
                 fondo_cal = np.zeros((3))
                 self.dosis = np.zeros((3, num_puntos))
                 self.pixel = np.zeros((3, num_puntos))
                 self.pixel_menos_fondo_cal = np.zeros((3, num_puntos))
                 self.DoN = np.zeros((3, num_puntos))
            j = 0
            for i in range(4, len(data), 2):
                self.dosis[k][j] = data[i]
                self.pixel[k][j] = data[i+1]
                if j == 0:
                    fondo_cal[k] = self.pixel[k][0]
                self.pixel_menos_fondo_cal[k][j] = self.pixel[k][j] - fondo_cal[k]
                j += 1
            self.DoN[k] = Pixel_to_DoN_1(fondo_cal[k], self.pixel[k]).DoN
            if k != 2:
                line = f.readline()
        
        f.close()
        DoN_minimo = np.zeros((3))
        DoN_maximo = np.zeros((3))
        paso = np.zeros((3))
        num_puntos_plot = 50
        for i in range(3):
            DoN_minimo[i] = min(self.DoN[i])
            DoN_maximo[i] = max(self.DoN[i])
            paso[i] = (DoN_maximo[i] - DoN_minimo[i])/num_puntos_plot
        DoN_ajustada = np.zeros((3, num_puntos_plot))
        dosis_ajustada = np.zeros((3, num_puntos_plot))
        
        if self.modelo == 0:
            # Ajuste (DoN, dosis) a polinomio de tercer grado para cada canal
            self.aj = Ajuste(self.DoN, self.dosis, modelo)
            
            for k in range(3):
                j = 0
                i = DoN_minimo[k]
                while j < num_puntos_plot:
                    DoN_ajustada[k][j] = i
                    dosis_ajustada[k][j] = self.aj.Pol_3(i, self.aj.A[k], self.aj.B[k], self.aj.C[k], self.aj.D[k])
                    j += 1
                    i += paso[k]
            if flag_r == True:
                self.Plot(self.DoN[0], self.dosis[0], DoN_ajustada[0], dosis_ajustada[0], "Canal Rojo. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'red')
            if flag_g == True:
                self.Plot(self.DoN[1], self.dosis[1], DoN_ajustada[1], dosis_ajustada[1], "Canal Verde. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'green')
            if flag_g == True:
                self.Plot(self.DoN[2], self.dosis[2], DoN_ajustada[2], dosis_ajustada[2], "Canal Azul. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'blue')
        
        if self.modelo == 1:
            # Ajuste (DoN, dosis) a percolacion
            self.aj = Ajuste(self.DoN, self.dosis, modelo)
            
            for k in range(3):
                j = 0
                i = DoN_minimo[k]
                while j < num_puntos_plot:
                    DoN_ajustada[k][j] = i
                    dosis_ajustada[k][j] = self.aj.Percolacion(i, self.aj.A[k], self.aj.B[k], self.aj.C[k], self.aj.D[k])
                    j += 1
                    i += paso[k]
            if flag_r == True:
                self.Plot(self.DoN[0], self.dosis[0], DoN_ajustada[0], dosis_ajustada[0], "Canal Rojo. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'red')
            if flag_g == True:
                self.Plot(self.DoN[1], self.dosis[1], DoN_ajustada[1], dosis_ajustada[1], "Canal Verde. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'green')
            if flag_g == True:
                self.Plot(self.DoN[2], self.dosis[2], DoN_ajustada[2], dosis_ajustada[2], "Canal Azul. " + modelos_ajuste_ids[modelo], "DoN", "Dosis", 'black', 'blue')
               
    def Plot(self, x1, y1, x2, y2, titulo, x_label, y_label, color1, color2):
        plt.scatter(x1, y1, s = 7, label = "Medidos", color = color1)
        plt.scatter(x2, y2, s = 2, label = "Ajustados", color = color2)
        plt.legend(loc="upper left")
        plt.title(titulo)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.show()

#
# Clase que actualiza el objeto principal con la elección del modelo de ajuste. Si la imagen tiene ya asignadas curvas de calibración
# calcula la matriz de dosis y actualiza el objeto principal
#
class Set_modelo:
    
    def __init__(self, objeto_ppal, selection, ids):
       
        objeto_ppal.set_modelo_ajuste(selection)
        objeto_ppal.set_modelo_ajuste_id(ids[selection])

        messagebox.showinfo(message="Modelo de ajuste registrado!!!", title="Aviso!!!")
        
        if objeto_ppal.get_cal(0)!= -1:
            messagebox.showinfo(message="Se va a generar la matriz de dosis!!!", title="Aviso!!!")
            lc = Lee_curva_calibracion(objeto_ppal.get_cal(0), objeto_ppal.get_modelo_ajuste(), True, True, True)
            dosis = DoN3D_to_Dosis(lc, objeto_ppal).Dosis
            objeto_ppal.set_dosis(dosis)
            messagebox.showinfo(message="Imagen calibrada!!!", title="Aviso!!!")


#
# Clase que verifica que exista el fichero de dosimetria, lo lee y presenta una lista para seleccionar la curva
# correspondiente al canal rojo. Supone que las consecutivas corresponden al canal verde y azul
# Finalmente actualiza el objeto principal
#
class Selecciona_modelo:
    def __init__(self, foco, canvas, objeto_ppal):
        global modelos_ajuste_ids
        
        self.foco = foco
        self.canvas = canvas
        self.objeto_ppal = objeto_ppal
        
        self.ids = []
        self.ids = modelos_ajuste_ids
        
        if foco == 1:
            pos_x = int(canvas.winfo_width()/2)
            pos_y = int(canvas.winfo_height()/2)
        else:
            pos_x = int(3 * canvas.winfo_width()/2)
            pos_y = int(canvas.winfo_height()/2)
            
        self.lista_modelos = tk.Listbox(width = 30, height = len(self.ids), exportselection=False, selectmode = 'SINGLE')
        self.lista_modelos.place(x = pos_x, y = pos_y)
        self.lista_modelos.insert(0, *self.ids)
        self.lista_modelos.configure(state='normal')
        self.lista_modelos.bind('<<ListboxSelect>>', self.onselect)
        scroll1 = tk.Scrollbar(orient=tk.VERTICAL, width = 10)
        scroll1.place(in_=self.lista_modelos, relx=1, relheight=1, bordermode="outside")
        self.lista_modelos['yscrollcommand'] = scroll1.set
        self.lista_modelos.configure(state='normal')      
        scroll1.configure(command= self.lista_modelos.yview)
        self.label_modelos = tk.Label()
        x_label = max(0, pos_x - 25)
        y_label = max(0, pos_y - 25)
        self.label_modelos.place(x = x_label, y = y_label)
        self.label_modelos['text'] = "Modelos de ajuste dosis=f(DoNeta)"
    
    def onselect(self, event):
         widget = event.widget
         selection = widget.curselection()
         self.selection = selection[0]
         self.lista_modelos.place_forget()
         self.label_modelos.place_forget()
         Set_modelo(self.objeto_ppal, self.selection, self.ids)

#
# Clase que actualiza el objeto principal con la elección de la curva de calibración. Si la imagen tiene asignado
# modelo de ajuste, calcula la matriz de dosis y actualiza el objeto principal
#
class Set_curva_dosimetria:
    
    def __init__(self, objeto_ppal, selection, ids):
       
        objeto_ppal.set_cal(0, selection*3)
        objeto_ppal.set_cal(1, selection*3+1)
        objeto_ppal.set_cal(2, selection*3+2)
        objeto_ppal.set_cal_id(0, str(ids[selection*3]))
        objeto_ppal.set_cal_id(1, str(ids[selection*3+1]))
        objeto_ppal.set_cal_id(2, str(ids[selection*3+2]))

        messagebox.showinfo(message="Curvas de calibración registradas!!!", title="Aviso!!!")
        
        if objeto_ppal.get_modelo_ajuste()!= -1:
            messagebox.showinfo(message="Se va a generar la matriz de dosis!!!", title="Aviso!!!")
            lc = Lee_curva_calibracion(objeto_ppal.get_cal(0), objeto_ppal.get_modelo_ajuste(), True, True, True)
            dosis = DoN3D_to_Dosis(lc, objeto_ppal).Dosis
            objeto_ppal.set_dosis(dosis)
            messagebox.showinfo(message="Imagen calibrada!!!", title="Aviso!!!")

#
# Clase que verifica que exista el fichero de dosimetria, lo lee y presenta una lista para seleccionar la curva
# correspondiente al canal rojo. Supone que las consecutivas corresponden al canal verde y azul
# Finalmente actualiza el objeto principal
#
class Selecciona_curva:
    def __init__(self, foco, canvas, objeto_ppal):
        global fichero_dosimetria
        
        self.foco = foco
        self.canvas = canvas
        self.objeto_ppal = objeto_ppal
        
        if path.exists(fichero_dosimetria):
            pass
        else:
            messagebox.showinfo(message="No hay acceso al fichero de dosismetria¡¡¡\'" + str(fichero_dosimetria) + "\'", title="Aviso!!!")
            sys.exit()
        
        self.ids = []
        self.curvas = []
        f = open(fichero_dosimetria)
        contador = 0
        while True:
            line = f.readline()
            if line =="":
                break
            else:
                data = line.split("#")
                self.ids.append(data[1])
                if contador%3 ==0:
                    self.curvas.append(data[1])
                contador +=1
        f.close()
        
        if foco == 1:
            pos_x = int(canvas.winfo_width()/2)
            pos_y = int(canvas.winfo_height()/2)
        else:
            pos_x = int(3 * canvas.winfo_width()/2)
            pos_y = int(canvas.winfo_height()/2)
            
        self.lista_curvas = tk.Listbox(width = 30, height = len(self.curvas), exportselection=False, selectmode = 'SINGLE')
        self.lista_curvas.place(x = pos_x, y = pos_y)
        self.lista_curvas.insert(0, *self.curvas)
        self.lista_curvas.configure(state='normal')
        self.lista_curvas.bind('<<ListboxSelect>>', self.onselect)
        scroll1 = tk.Scrollbar(orient=tk.VERTICAL, width = 10)
        scroll1.place(in_=self.lista_curvas, relx=1, relheight=1, bordermode="outside")
        self.lista_curvas['yscrollcommand'] = scroll1.set
        self.lista_curvas.configure(state='normal')      
        scroll1.configure(command= self.lista_curvas.yview)
        self.label_curvas = tk.Label()
        x_label = max(0, pos_x - 25)
        y_label = max(0, pos_y - 25)
        self.label_curvas.place(x = x_label, y = y_label)
        self.label_curvas['text'] = "Curvas de calibracion canal rojo"
    
    def onselect(self, event):
         widget = event.widget
         selection = widget.curselection()
         self.selection = selection[0]
         self.lista_curvas.place_forget()
         self.label_curvas.place_forget()
         Set_curva_dosimetria(self.objeto_ppal, self.selection, self.ids)

#
# Clase que obtiene los valores de DoN a partir de los valores de pixel.
# matriz_pixel es un array unidimensional de valores de pixeles
# pixel_velo tiene dimension 1
#
class Pixel_to_DoN_1:
    def __init__(self, pixel_velo, matriz_pixel):
        buffer = np.zeros((len(matriz_pixel)))
        for i in range(len(matriz_pixel)):
                if matriz_pixel[i] <= 0: # Para salvar la discontinuidad del logaritmo
                    matriz_pixel[i] = 1
                buffer[i] = -math.log10(matriz_pixel[i]) + math.log10(pixel_velo)
                #print(buffer[i][j])
                if buffer[i] < 0: # Para valores de pixeles menores al del velo
                    buffer[i] = 0
        
        self.DoN = buffer
  
#
# Clase que obtiene los valores de DoN a partir de los valores de pixel.
# matriz_pixel es un array bidimensional de valores de pixeles
# pixel_velo tiene dimension 1
#
class Pixel_to_DoN_2:
    def __init__(self, pixel_velo, matriz_pixel):
        buffer = np.zeros((len(matriz_pixel), len(matriz_pixel[0])))
        for i in range(len(matriz_pixel)):
            for j in range(len(matriz_pixel[0])):
                if matriz_pixel[i][j] <= 0: # Para salvar la discontinuidad del logaritmo
                    matriz_pixel[i][j] = 1
                buffer[i][j] = -math.log10(matriz_pixel[i][j]) + math.log10(pixel_velo)
                if buffer[i][j] < 0: # Para valores de pixeles menores al del velo
                    buffer[i][j] = 0
        self.DoN = buffer

#
# Clase que devuelve el valor medio por canal, de una submatriz3D definida apartir de (x0, y0)
# como esquina superior izquierda y (x1, y1) como esquina inferior derecha de una matriz de datos 3D.
# Los puntos de las rectas correspondientes a (x1, y1) entran en el estudio.
class Get_media_recorte:
    def __init__(self, matriz3D, x0, y0, x1, y1):
        self.num_puntos = 0
        self.media = np.zeros((3))
        
        for k in range(3):
            suma_k = 0
            contador = 0
            for i in range(y0, y1+1):
                for j in range(x0, x1+1):
                    suma_k += matriz3D[k][i][j]
                    contador += 1
            self.media[k] = suma_k/contador
            self.num_puntos = contador

#
# Clase paar dibujar un cuadrado discontinuo sobre el canvas indicado centrado en el aposición del ratón
#
# ventana   -> Tamaño del rectángulo
class Dibuja_rectangulo_canvas:
    def __init__(self, canvas, ancho_canvas, alto_canvas, x, y, ventana, flag):
        self.punto_a_x = max(x - ventana, 0)
        self.punto_b_x = min(x + ventana, ancho_canvas-1)
        self.punto_a_y = max(y - ventana, 0)
        self.punto_b_y = min(y + ventana, alto_canvas-1)
        
        if flag == True:
            self.rectangulo = canvas.create_rectangle(self.punto_a_x, self.punto_a_y, self.punto_b_x, self.punto_b_y,
                                                      dash=(4, 2), width=1, stipple="gray50", tags='rectangle'
                                                      )

#
# Clase para pasar coordenadas de widget a matriz 2D de datos
#
class  Coordenadas_widget_to_imagen:
    def __init__(self, x_widget, y_widget, width_widget, height_widget, width_imagen, height_imagen):
        self.x = int(x_widget*width_imagen/width_widget)
        self.y = int(y_widget*height_imagen/height_widget)


#
# Clase paar obtener el factor pixel/mm de un fichero tif
#
class Get_dpi_tif:
    def __init__(self, fichero_imagen):
        im = Image.open(fichero_imagen)
        info = im.info
        dpiy, dpix = info['dpi'][:2]
        if dpix != dpiy:
            messagebox.showinfo(message="¡¡¡dpi diferente en \'x\' (" + str(dpix) + ") y en \'y\' (" + str(dpiy) + ")!!!" , title="Aviso!!!")
            sys.exit()
        self.factor_pixel_mm = dpix/25.4

#
# Clase que transforma una matriz3D_datos (3, n, m) a una matriz3D_imagen (n, m, 3)
# Principalmente para representacion visual
#
class matriz_datos_2_matriz_im:
    def __init__(self, matriz_datos):
        self.matriz_im = np.zeros((len(matriz_datos[0]), len(matriz_datos[0][0]), 3))
        self.matriz_im[:, :, 0] = matriz_datos[2]
        self.matriz_im[:, :, 1] = matriz_datos[1]
        self.matriz_im[:, :, 2] = matriz_datos[0]

#
# Clase que lee un archivo tipo "pnc" y devuelve:
#   factor_pixel_mm (self.pixel_size_mm)
#   anchura y altura en pixeles (self.width_imagen_pixel y self.height_imagen_pixel)
#   matriz de dosis3D (self.dosis_onc_3D, copiando dosis2D en los tres canales)
#
class Lee_PNC:
    def __init__(self, fichero_pnc):
        f = open(fichero_pnc, "r")
              
        while True:
            linea = f.readline()
            if linea == "":
                break
            dummy = linea[0:1]
            if dummy == ",":
                datax = linea.split(",")
                X0 = float(datax[1])
                X1 = float(datax[2])
                self.width_imagen_pixel = len(datax) - 2
                self.pixel_size_mm = ""
                valor = (X1 - X0)*10.0
                valor = 1 / valor
                self.pixel_size_mm = float(round(valor, 3))

                contador = 0
                while True:
                    linea = f.readline()
                    datay = linea.split(",")
                    if linea != "":
                        contador += 1
                    else:
                        break
                break
        self.height_imagen_pixel = contador
        f.close()
       
        self.dosis_pnc_3D = np.zeros((3, self.height_imagen_pixel, self.width_imagen_pixel))
        
        f = open(fichero_pnc, "r")
              
        while True:
            linea = f.readline()
            dummy = linea[0:1]
            if dummy == ",":
               for i in range(0, self.height_imagen_pixel):
                    linea = f.readline()
                    datay = linea.split(",")
                    for j in range(0, self.width_imagen_pixel):
                        self.dosis_pnc_3D[0][i][j] = float(datay[j+1])
               break
        f.close()

        for i in range(1, 3):
            self.dosis_pnc_3D[i] = self.dosis_pnc_3D[0]

#
# Clase para obtener la extensión del path de un fichero
#
class Get_extension:
    def __init__(self, cadena):
        self.path =""
        long_i = len(cadena)
        long_f = 0
        
        buffer = cadena
        while True:
            pos = buffer.find(".")
            if pos != -1:
                buffer = buffer[pos+1:]
            else:
                long_f = len(buffer)
                break
        buffer = cadena[long_i - long_f:]
        self.extension = buffer

#
# Clase para definir la estructura de los dos objetos principales, uno para la izquierda y otro para la derecha.
# Contiene las funciones de actualizacion y query para cada uno de los elementos.
#
# tipo                  -> pnc ó tif
# factor_pixel_mm       -> Relación pixel por mm para la imagen
# cal                   -> Array de tres enteros que relacion las curvas de calibración con cada canal de la imagen (0 R, 1 G y 2 B).
#                          -1 significa "canal no calibrado". Solo en tif
# cal_id                -> Array de tres strings para ids curvas de calibracion. Solo en tif
# pixel_velo            -> Array de tres enteros para los valores de pixeles del velo en cada canal. Solo en tif
# dosis_linealizacion   -> Dosis conocida a efectos de linealización. Sólo en tif
# pixel_dosis_conocida  -> Array de tres enteros para los valores de pixeles de dosis conocida en cada canal. Solo en tif
# modelo_ajuste         -> Modelo de ajuste de la curva (don, dosis). Utiliza el mismo para los tres canales de la imagen
#               0 -> Polinomio de grado tres. Uno por canal
#               1 -> Percolacion. Uno por canal
# modelo_ajuste_id      -> Etiqueta de identificacion del modelo
# width                 -> Ancho en pixeles de la imagen actual
# height                -> Alto en pixeles de la imagen actual
# pixel                 -> Matriz en pixeles de la imagen actual. Tridimensional. Canal primera dimension
# don                   -> Matriz en densidad optica neta de la imagen actual. Tridimensional. Canal primera dimension
# dosis                 -> Matriz en dosis de la imagen actual. Tridimensional. Canal primera dimension
#
class Objeto_Ppal:
    def __init__(self):
        self.tipo = None
        self.factor_pixel_mm = None
        self.cal = np.zeros((3)) # Solo en tif
        self.cal_id = ["", "", ""] # Solo en tif
        self.pixel_velo = np.zeros((3)) # Solo en tif
        self.pixel_dosis_conocida = np.zeros((3)) # Solo en tif
        self.dosis_linealizacion = "" # Solo en tif
        self.modelo_ajuste = -1
        self.modelo_ajuste_id = ""
        self.width = -1
        self.height = -1
        self.pixel = None # Solo en tif
        self.don = None # Solo en tif
        self.dosis = None
        self.fdo = None  # Film Dose Object. Solo en tif
        
        for i in range(3):
            self.cal[i] = -1
            self.pixel_velo[i] = -1
            self.pixel_dosis_conocida[i] = -1
        
    def set_tipo(self, tipo):
        self.tipo = tipo
        
    def get_tipo(self):
        return self.tipo
        
    def set_factor_pixel_mm(self, factor_pixel_mm):
        self.factor_pixel_mm = factor_pixel_mm
        
    def get_factor_pixel_mm(self):
        return self.factor_pixel_mm
    
    def set_width(self, width):
        self.width = width
        
    def get_width(self):
        return self.width
    
    def set_height(self, height):
        self.height = height
        
    def get_height(self):
        return self.height
    
    def set_cal(self, canal, curva):
        self.cal[canal] = curva
        
    def set_cal_id(self, canal, idd):
        self.cal_id[canal] = idd
    
    def get_cal_id(self, canal):
        return self.cal_id[canal]
    
    def get_cal(self, canal):
        return self.cal[canal]
    
    def set_modelo_ajuste(self, modelo_ajuste):
        self.modelo_ajuste = modelo_ajuste
        
    def get_modelo_ajuste(self):
        return self.modelo_ajuste
    
    def set_modelo_ajuste_id(self, modelo_ajuste_id):
        self.modelo_ajuste_id = modelo_ajuste_id
        
    def get_modelo_ajuste_id(self):
        return self.modelo_ajuste_id
    
    def set_pixel_velo(self, canal, velo):
        self.pixel_velo[canal] = velo
        
    def get_pixel_velo(self, canal):
        if canal != 3:
            return self.pixel_velo[canal]
        else:
            return self.pixel_velo
        
    def set_pixel_dosis_conocida(self, canal, pixel_dosis_conocida):
        self.pixel_dosis_conocida[canal] = pixel_dosis_conocida
        
    def get_pixel_dosis_conocida(self, canal):
        if canal != 3:
            return self.pixel_dosis_conocida[canal]
        else:
            return self.pixel_dosis_conocida
    
    def set_dosis_linealizacion(self, dosis_linealizacion):
        self.dosis_linealizacion = dosis_linealizacion
        
    def get_dosis_linealizacion(self):
        return self.dosis_linealizacion
    
    def set_dosis(self, dosis):
        self.dosis = dosis
        
    def get_dosis(self, query):
        if query == 0:
            return self.dosis[0]
        if query == 1:
            return self.dosis[1]
        if query == 2:
            return self.dosis[2]
        if query == 3:
            return self.dosis
        
    def set_pixel(self, pixel):
        self.pixel = pixel
        
    def get_pixel(self, query):
        if query == 0:
            return self.pixel[0]
        if query == 1:
            return self.pixel[1]
        if query == 2:
            return self.pixel[2]
        if query == 3:
            return self.pixel
        
    def set_don(self, don):
        self.don = don
        
    def get_don(self, query):
        if query == 0:
            return self.don[0]
        if query == 1:
            return self.don[1]
        if query == 2:
            return self.don[2]
        if query == 3:
            return self.don

class MainFrame(tk.Tk):
    
    def __init__(self):
        #Variables acumuladas para ventana ecualizacion
        self.desplazamiento_lateral = 0
        self.desplazamiento_vertical = 0
        self.angulo_deg_giro = 0
        
        # Define escritorio principal
        self.Title = "JJGamma v.1.0"
        ancho_ventana = 1100
        alto_ventana = 600
        padx_ppal = 5
        pady_ppal = 5
        geometria_ventana = str(ancho_ventana + 3*padx_ppal) + "x" + str(alto_ventana + 2*pady_ppal)
        self.ancho_canvas_ppal = int(ancho_ventana/2)
        self.alto_canvas_ppal = int(0.9*alto_ventana)
        
        self.root=tk.Tk()
        self.root.geometry(geometria_ventana)
        self.root.title(self.Title)
        
        self.canvas1_ppal = tk.Canvas(self.root, bg ="blue", height=self.alto_canvas_ppal, width=self.ancho_canvas_ppal)
        self.canvas1_ppal.grid()
        self.canvas1_ppal.place(x=padx_ppal, y=pady_ppal)
        self.canvas1_ppal.bind("<Motion>", self.mouse_move_1)
        self.canvas1_ppal.bind("<Button-1>", self.mouse_clicked_1)
        self.canvas1_ppal.bind("<B1-Motion>", self._on_drag_1)
        self.canvas1_ppal.bind("<ButtonRelease-1>", self._on_drop_1)
        self.canvas2_ppal = tk.Canvas(self.root, bg = "red", height=self.alto_canvas_ppal, width=self.ancho_canvas_ppal)
        self.canvas2_ppal.grid()
        self.canvas2_ppal.place(x=self.ancho_canvas_ppal + padx_ppal, y=pady_ppal)
        self.canvas2_ppal.bind("<Motion>", self.mouse_move_2)
        self.canvas2_ppal.bind("<Button-1>", self.mouse_clicked_2)
        self.canvas2_ppal.bind("<B1-Motion>", self._on_drag_2)
        self.canvas2_ppal.bind("<ButtonRelease-1>", self._on_drop_2)
        
        self.popup_menu = tk.Menu(self.root, tearoff=0)
        self.popup_menu.add_command(label="Gestión Curvas Calibración", command=lambda:self.gestion())
        self.popup_menu.add_command(label="Abrir Imagen", command=lambda:self.abrir_selected(foco))
        self.popup_menu.add_command(label="Flat Scanner", command=lambda:self.flat_scanner(foco))
        self.popup_menu.add_command(label="Establecer Pixel Velo", command=lambda:self.establecer_velo(foco))
        self.popup_menu.add_command(label="Establecer Pixel Dosis conocida (200 cGy)", command  = lambda: self.establecer_dosis_conocida(foco))
        self.popup_menu.add_command(label="Seleccionar Curva de Calibración", command  = lambda: self.establecer_curva_calibracion(foco))
        self.popup_menu.add_command(label="Seleccionar Modelo de Calibración", command  = lambda: self.establecer_modelo_calibracion(foco))
        self.popup_menu.add_command(label="Linealización Dosimétrica", command  = lambda: self.linealizacion(foco))
        self.popup_menu.add_command(label="Multicanal_Micke_2011", command  = lambda: self.multicanal(foco))
        self.popup_menu.add_command(label="Flip Horizontal", command  = lambda: self.flip(foco, 0))
        self.popup_menu.add_command(label="Flip Vertical", command  = lambda: self.flip(foco, 1))
        self.popup_menu.add_command(label="Rotar 90º respecto del centro", command  = lambda: self.rota90(foco, -1))
        self.popup_menu.add_command(label="Rotar 270º respecto del centro", command  = lambda: self.rota90(foco, 1))
        self.popup_menu.add_command(label="Desplazar Horizontal", command  = lambda: self.traslacion(foco, 1))
        self.popup_menu.add_command(label="Desplazar Vertical", command  = lambda: self.traslacion(foco, -1))
        self.popup_menu.add_command(label="Rotar cualquier angulo respecto del centro", command  = lambda: self.rotar_centro(foco))
        self.popup_menu.add_command(label="Recortar", command  = lambda: self.recortar(foco))
        self.popup_menu.add_command(label="Guardar dosis en formato pnc", command  = lambda: self.guardar(foco))
        self.popup_menu.add_command(label="Ver perfiles", command  = lambda: self.perfiles(foco))
        self.canvas1_ppal.bind("<Button-3>", self.popup)
        self.canvas2_ppal.bind("<Button-3>", self.popup)
        
        self.btn1 = tk.Button(self.root, text="Ecualización", command = lambda: self.new_windowe(self.root))
        self.btn1.place(x=230, y=560)
        self.btn2 = tk.Button(self.root, text="Recortar ambas imágenes", command = lambda: self.recortar_ambas())
        self.btn2.place(x=320, y=560)
        self.btn3 = tk.Button(self.root, text="Diferencia ambas imágenes", command = lambda: self.new_windowd(self.root))
        self.btn3.place(x=480, y=560)
        self.btn4 = tk.Button(self.root, text="Gamma ambas imágenes", command = lambda: self.new_windowg(self.root))
        self.btn4.place(x=650, y=560)
        
        self.o1 = None
        self.o2 = None
        
        self.root.mainloop()
        
    def new_windowe(self, master):
        # Comprueba que ambas imagenes estén cargadas y en dosis
        flag_comprobacion_e = True
        if self.o1 != None and self.o2 != None:
            if self.o1.get_tipo() != "pnc" and (self.o1.get_cal(0) == -1 or self.o1.get_modelo_ajuste() == -1):
                flag_comprobacion_e = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 1!!!" , title="Aviso!!!")
                self.close_windowse_salir()
            if self.o2.get_tipo() != "pnc" and (self.o2.get_cal(0) == -1 or self.o2.get_modelo_ajuste() == -1):
                flag_comprobacion_e = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 2!!!" , title="Aviso!!!")
                self.close_windowse_salir()
                
            if flag_comprobacion_e == True:
                n_cols_e = 15
                n_rows_e = 8
                ancho_canvas_e_p = 75
                alto_canvas_e_p = 75
                ancho_canvas_e_g = int(ancho_canvas_e_p * 5.5)
                alto_canvas_e_g = int(alto_canvas_e_p * 5.5)
                ancho_canvas_e_m = ancho_canvas_e_p * 8
                alto_canvas_e_m = int(alto_canvas_e_p * 3.25)
                padx_e = 3
                pady_e = 3
                ancho_ventana_e = n_cols_e*ancho_canvas_e_p + (n_cols_e+2)*padx_e
                alto_ventana_e = n_rows_e*alto_canvas_e_p + (n_rows_e+2)*pady_e
                self.ancho_canvas_e_g  = ancho_canvas_e_g
                self.alto_canvas_e_g  = alto_canvas_e_g
                self.ancho_canvas_e_m = ancho_canvas_e_m
                self.alto_canvas_e_m = alto_canvas_e_m
        
                self.newWindowe = tk.Toplevel(master)
                geom = str(ancho_ventana_e) + "x" + str(alto_ventana_e)
                self.newWindowe.geometry(geom)
                self.newWindowe.grid()
        
                self.Canvas1 = tk.Canvas(self.newWindowe, bg = "blue", height=alto_canvas_e_p, width=ancho_canvas_e_p)
                self.Canvas2 = tk.Canvas(self.newWindowe, bg = 'red', height=alto_canvas_e_p, width=ancho_canvas_e_p)
                self.Canvas3 = tk.Canvas(self.newWindowe, bg = 'green', height=alto_canvas_e_g, width=ancho_canvas_e_g)
                self.Canvas4 = tk.Canvas(self.newWindowe, bg = 'green', height=alto_canvas_e_m, width=ancho_canvas_e_m)
                self.Canvas5 = tk.Canvas(self.newWindowe, bg = 'green', height=alto_canvas_e_m, width=ancho_canvas_e_m)
                self.Canvas1.place(x=padx_e, y=pady_e)
                self.Canvas2.place(x=ancho_canvas_e_g-ancho_canvas_e_p+padx_e, y = pady_e)
                self.Canvas3.place(x=padx_e, y=alto_canvas_e_p + 2*pady_e)
                self.Canvas4.place(x=ancho_canvas_e_g + ancho_canvas_e_p + 2*padx_e, y=pady_e)
                self.Canvas5.place(x=ancho_canvas_e_g + ancho_canvas_e_p + 2*padx_e, y=2*   pady_e + alto_canvas_e_m)
                self.Canvas3.bind("<Button-1>", self.mouse_clicked_e)
                self.Canvas3.bind("<Motion>", self.mouse_move_e)
        
                self.L1 = tk.Label(self.newWindowe)
                self.L1.place(x = int(ancho_canvas_e_g*0.4),y = int(alto_canvas_e_p*0.2))
                self.L1['text'] = "Transparencia"
                self.valor_alpha = tk.StringVar()
                self.E_alpha = tk.Entry(self.newWindowe, width=3, textvariable=self.valor_alpha)
                self.E_alpha.place(x = int(ancho_canvas_e_g*0.4)+30, y = int(alto_canvas_e_p*0.45)+5)
                self.valor_alpha.set("0.5")
                self.B0L1 = tk.Button(self.newWindowe, text = '-', width = 2, command = lambda: self.blend("menos", ancho_canvas_e_g, alto_canvas_e_g))
                self.B0L1.place(x = int(ancho_canvas_e_g*0.35), y = int(alto_canvas_e_p*0.45))
                self.B1L1 = tk.Button(self.newWindowe, text = '+', width = 2, command = lambda: self.blend("mas", ancho_canvas_e_g, alto_canvas_e_g))
                self.B1L1.place(x = int(ancho_canvas_e_g*0.58), y = int(alto_canvas_e_p*0.45))
                self.L2 = tk.Label(self.newWindowe)
                self.L2.place(x = int(ancho_canvas_e_g*0.2)-15,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)-5)
                self.L2['text'] = "Desplazamientos (pixeles)"
                self.valor_despl = tk.StringVar()
                self.E_despl = tk.Entry(self.newWindowe, width=3, textvariable=self.valor_despl)
                self.E_despl.place(x = int(ancho_canvas_e_g*0.2)+30, y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+57)
                self.valor_despl.set("10")
                self.B0L2 = tk.Button(self.newWindowe, text = 'Left', width = 4, command = lambda: self.move(int(self.valor_despl.get())*-1, 0, 'dust1.png', 'dust2.png'))
                self.B0L2.place(x = int(ancho_canvas_e_g*0.2)-20,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+55)
                self.B1L2 = tk.Button(self.newWindowe, text = 'Up', width = 4, command = lambda: self.move(0, int(self.valor_despl.get())*-1,'dust1.png', 'dust2.png'))
                self.B1L2.place(x = int(ancho_canvas_e_g*0.2)+20,y =int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+20)
                self.B2L2 = tk.Button(self.newWindowe, text = 'Right', width = 4, command = lambda: self.move(int(self.valor_despl.get()), 0, 'dust1.png', 'dust2.png'))
                self.B2L2.place(x = int(ancho_canvas_e_g*0.2)+60,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+55)
                self.B3L2 = tk.Button(self.newWindowe, text = 'Down', width = 4, command = lambda: self.move(0, int(self.valor_despl.get()), 'dust1.png', 'dust2.png'))
                self.B3L2.place(x = int(ancho_canvas_e_g*0.2)+20,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+85)
                self.L3 = tk.Label(self.newWindowe)
                self.L3.place(x = int(ancho_canvas_e_g*0.8)-15,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)-5)
                self.L3['text'] = "Giros (º)"
                self.B0L3 = tk.Button(self.newWindowe, text = '-', width = 4, command = lambda: self.giro("menos", float(self.valor_giro.get()), 'dust1.png', 'dust2.png'))
                self.B0L3.place(x = int(ancho_canvas_e_g*0.2)+195,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+25)
                self.B1L3 = tk.Button(self.newWindowe, text = '+', width = 4, command = lambda: self.giro("mas", float(self.valor_giro.get()), 'dust1.png', 'dust2.png'))
                self.B1L3.place(x = int(ancho_canvas_e_g*0.2)+270,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+25)
        
                self.valor_giro = tk.StringVar()
                self.E_giro = tk.Entry(self.newWindowe, width=3, textvariable=self.valor_giro)
                self.E_giro.place(x = int(ancho_canvas_e_g*0.2)+240, y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+30)
                self.valor_giro.set("0.5")
        
                self.L_indice_canal = tk.Label(self.newWindowe)
                self.L_indice_canal.place(x = int(ancho_canvas_e_g*0.8) + 100 ,y =  int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+65)
                self.L_indice_canal['text'] = "Canal dosis Imagen_2"
                self.indice_canal = tk.StringVar()
                self.E_indice_canal = tk.Entry(self.newWindowe, width=3, textvariable=self.indice_canal)
                self.E_indice_canal.place(x = int(ancho_canvas_e_g*0.8)+140, y =  int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+90 )
                self.indice_canal.set("0")
        
                self.L_factor_dosis = tk.Label(self.newWindowe)
                self.L_factor_dosis.place(x = int(ancho_canvas_e_g*0.8) - 55 ,y =  int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+65)
                self.L_factor_dosis['text'] = "Factor dosis Imagen_2"
                self.factor_dosis_o2 = tk.StringVar()
                self.E_factor_dosis_o2 = tk.Entry(self.newWindowe, width=6, textvariable=self.factor_dosis_o2)
                self.E_factor_dosis_o2.place(x = int(ancho_canvas_e_g*0.8)-10, y =  int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)+90)
        
                self.factor_dosis_o2.set("1.0")
        
                self.aceptarButton = tk.Button(self.newWindowe, text = 'Aceptar', width = 6, command = self.close_windowse_aceptar)
                self.aceptarButton.place(x=ancho_ventana_e - 125, y=alto_ventana_e - 50)
        
                self.salirButton = tk.Button(self.newWindowe, text = 'Salir', width = 5, command = self.close_windowse_salir)
                self.salirButton.place(x=ancho_ventana_e - 65, y=alto_ventana_e - 50)
                self.newWindowe.title(self.Title)
        
                self.valor_iso = tk.StringVar()
                self.E_iso = tk.Entry(self.newWindowe, width=3, textvariable=self.valor_iso)
                self.E_iso.place(x = int(ancho_canvas_e_g*0.8)+430, y =  int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02) + 20)
                self.valor_iso.set("180")
                self.Li = tk.Label(self.newWindowe)
                self.Li.place(x = int(ancho_canvas_e_g*0.8)+390,y = int((alto_canvas_e_g + alto_canvas_e_p + 2*pady_e)*1.02)-5)
                self.Li['text'] = "Isodosis a representar"
        
                master.withdraw()
        
                # Lienzos. Prepara ambas imagenes para que tengan el mismo numero de pixeles y el mismo factor pixel/mm
                il = Lienzo(self.o1, self.o2)
                if self.o1.get_tipo() == "pnc":
                    self.o1.set_factor_pixel_mm(float(il.f))
                    self.o1.set_width(int(il.lw_pixel_new))
                    self.o1.set_height(int(il.lh_pixel_new))
                    self.o1.set_dosis(il.dosis_1_new)
                else:
                    self.o1.set_factor_pixel_mm(float(il.f))
                    self.o1.set_width(int(il.lw_pixel_new))
                    self.o1.set_height(int(il.lh_pixel_new))
                    self.o1.set_pixel(il.pixel_1_new)
                    self.o1.set_don(il.don_1_new)
                    self.o1.set_dosis(il.dosis_1_new)
                if self.o2.get_tipo() == "pnc":
                    self.o2.set_factor_pixel_mm(float(il.f))
                    self.o2.set_width(int(il.lw_pixel_new))
                    self.o2.set_height(int(il.lh_pixel_new))
                    self.o2.set_dosis(il.dosis_2_new)
                else:
                    self.o2.set_factor_pixel_mm(float(il.f))
                    self.o2.set_width(int(il.lw_pixel_new))
                    self.o2.set_height(int(il.lh_pixel_new))
                    self.o2.set_pixel(il.pixel_2_new)
                    self.o2.set_don(il.don_2_new)
                    self.o2.set_dosis(il.dosis_2_new)
                    
                self.x_clicked = -1
                self.y_clicked = -1    
                   
                # Muestra imagenes de objetos ppales pequeños de windows_e
                if self.o1.get_tipo() == 'pnc':
                    dst_im = matriz_datos_2_matriz_im(self.o1.get_dosis(3)).matriz_im
                    self.Mostrar(dst_im, 'dust1.png', self.Canvas1, 1, ancho_canvas_e_p, alto_canvas_e_p, "pnc")
                else:
                    dst_im = matriz_datos_2_matriz_im(self.o1.get_pixel(3)).matriz_im
                    self.Mostrar(dst_im, 'dust1.png', self.Canvas1, 1, ancho_canvas_e_p, alto_canvas_e_p, "tif")
            
                if self.o2.get_tipo() == 'pnc':
                    dst_im = matriz_datos_2_matriz_im(self.o2.get_dosis(3)).matriz_im
                    self.Mostrar(dst_im, 'dust2.png', self.Canvas2, 2, ancho_canvas_e_p, alto_canvas_e_p, "pnc")
                else:
                    dst_im = matriz_datos_2_matriz_im(self.o2.get_pixel(3)).matriz_im
                    self.Mostrar(dst_im, 'dust2.png', self.Canvas2, 2, ancho_canvas_e_p, alto_canvas_e_p, "tif")
                    
                #Muestra el alpha en canvas3
                gi = Get_isodosis_2D(self.o1.get_dosis(0), float(self.valor_iso.get()), 2, (255, 0, 0), True)
                cv2.imwrite('dust1.png', gi.imagen_isodosis)
                gi = Get_isodosis_2D(self.o2.get_dosis(0), float(self.valor_iso.get())*float(self.factor_dosis_o2.get()), 2, (0, 0, 255), True)
                cv2.imwrite('dust2.png', gi.imagen_isodosis)
                im1 = cv2.imread('dust1.png')  
                im2 = cv2.imread('dust2.png')
                dst = Alfa(im1, im2, float(self.valor_alpha.get())).dst
                cv2.imwrite('dust3.png', dst)
                self.Mostrare('dust3.png', self.Canvas3, 3, ancho_canvas_e_g, alto_canvas_e_g)
        else:
            messagebox.showinfo(message="¡¡¡Deben estar cargadas las dos imagenes!!!" , title="Aviso!!!")
    
    def Mostrare(self, fichero, canvas, indice, ancho_canvas, alto_canvas):
        img = Image.open(fichero)
        img = img.resize((ancho_canvas, alto_canvas), Image.ANTIALIAS)
        if indice == 1:
            self.img1e = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img1e)
            canvas.itemconfig(canvas_img, image=self.img1e)
        if indice == 2:
            self.img2e = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img2e)
            canvas.itemconfig(canvas_img, image=self.img2e)
        if indice == 3:
            self.img3e = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img3e)
            canvas.itemconfig(canvas_img, image=self.img3e)
        if indice == 4:
            self.img4e = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img4e)
            canvas.itemconfig(canvas_img, image=self.img4e)
        if indice == 5:
            self.img5e = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img5e)
            canvas.itemconfig(canvas_img, image=self.img5e)
    
    def blend(self, tipo, ancho_canvas, alto_canvas):
        valor = float(self.valor_alpha.get())
        if tipo == "mas":
            valor += 0.1
            valor = min(1.0, valor)
        else:
            valor -= 0.1
            valor = max(0.0, valor)
        self.valor_alpha.set(str(valor))
        canal = int(self.indice_canal.get())
        gi = Get_isodosis(self.o1.get_dosis(3), canal, float(self.valor_iso.get()), 2, (255, 0, 0), True)
        cv2.imwrite('dust1.png', gi.imagen_isodosis)
        gi = Get_isodosis(self.o2.get_dosis(3), canal, float(self.valor_iso.get())*float(self.factor_dosis_o2.get()), 2, (255, 0, 0), True)
        cv2.imwrite('dust2.png', gi.imagen_isodosis)
        im1 = cv2.imread('dust1.png')  
        im2 = cv2.imread('dust2.png')
        dstg = Rotar_respecto_any_point(im2, int(self.o2.get_width()/2), int(self.o2.get_height()/2), float(self.angulo_deg_giro), 1.0).dst
        dstt = Traslacion(dstg, int(self.desplazamiento_lateral), int(self.desplazamiento_vertical)).dst
        dstt_im = matriz_datos_2_matriz_im(dstt).matriz_im
        
        dst = Alfa(im1, dstt_im, float(self.valor_alpha.get())).dst
        cv2.imwrite('dust3.png', dst)
        self.Mostrare('dust3.png', self.Canvas3, 3, ancho_canvas, alto_canvas)
    
    def move(self, valor_x, valor_y, fichero1, fichero2):
        self.desplazamiento_lateral += valor_x
        self.desplazamiento_vertical += valor_y
        self.Canvas3.delete("all")
        canal = int(self.indice_canal.get())
        gi = Get_isodosis(self.o1.get_dosis(3), canal, float(self.valor_iso.get()), 2, (255, 0, 0), True)
        cv2.imwrite(fichero1, gi.imagen_isodosis)
        im1 = cv2.imread(fichero1)
        gi = Get_isodosis(self.o2.get_dosis(3), canal, float(self.valor_iso.get())*float(self.factor_dosis_o2.get()), 2, (255, 0, 0), True)
        cv2.imwrite(fichero2, gi.imagen_isodosis)
        im2 = cv2.imread(fichero2)
      
        dstg = Rotar_respecto_any_point(im2, int(self.o2.get_width()/2), int(self.o2.get_height()/2), float(self.angulo_deg_giro), 1.0).dst
        dstt = Traslacion(dstg, int(self.desplazamiento_lateral), int(self.desplazamiento_vertical)).dst
        dstt_im = matriz_datos_2_matriz_im(dstt).matriz_im
        
        dst = Alfa(im1, dstt_im, float(self.valor_alpha.get())).dst
        cv2.imwrite('dust3.png', dst)
        self.Mostrare('dust3.png', self.Canvas3, 3, self.ancho_canvas_e_g, self.alto_canvas_e_g)
        
        if self.x_clicked == -1:
            self.x_clicked = self.ancho_canvas_e_g/2
            
        if self.y_clicked == -1:
            self.y_clicked = self.alto_canvas_e_g/2
            
        cc = Coordenadas_widget_to_imagen(int(self.x_clicked), int(self.y_clicked), self.ancho_canvas_e_g, self.alto_canvas_e_g, self.o1.get_width(),self.o1.get_height())
    
        pe = Perfil_ecualizacion(int(cc.x), int(cc.y), self.desplazamiento_lateral, self.desplazamiento_vertical, float(self.angulo_deg_giro), 
                            self.o1, self.o2, canal, float(self.factor_dosis_o2.get()))
        buffer_x1 = pe.px_1
        buffer_y1 = pe.py_1
        buffer_x2 = pe.px_2
        buffer_y2 = pe.py_2
        
        if canal == 0:
            cadena = "Rojo"
        else:
            if canal == 1:
                cadena = "Verde"
            else:
                cadena = "azul"
        
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis (" + str(cadena) + ")")
        plt.title("Perfil Vertical " + str(cc.x))
        plt.plot(buffer_x1, color = 'blue')
        plt.plot(buffer_x2, color = 'red')
        plt.savefig('salidav.png', bbox_inches='tight')
        self.Mostrare('salidav.png', self.Canvas5, 5, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis")
        plt.title("Perfil Horizontal " + str(cc.y))
        plt.plot(buffer_y1, color = 'blue')
        plt.plot(buffer_y2, color = 'red')
        plt.savefig('salidah.png', bbox_inches='tight')
        self.Mostrare('salidah.png', self.Canvas4, 4, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
        
    def giro(self, tipo, angulo_deg, fichero1, fichero2):
        if tipo == "mas":
            self.angulo_deg_giro -= angulo_deg
        else:
            self.angulo_deg_giro += angulo_deg
        self.Canvas3.delete("all")
        canal = int(self.indice_canal.get())
        gi = Get_isodosis(self.o1.get_dosis(3), canal, float(self.valor_iso.get()), 2, (255, 0, 0), True)
        cv2.imwrite(fichero1, gi.imagen_isodosis)
        gi = Get_isodosis(self.o2.get_dosis(3), canal, float(self.valor_iso.get())*float(self.factor_dosis_o2.get()), 2, (255, 0, 0), True)
        cv2.imwrite(fichero2, gi.imagen_isodosis)
        im1 = cv2.imread(fichero1)  
        im2 = cv2.imread(fichero2)
      
        dstg = Rotar_respecto_any_point(im2, int(self.o2.get_width()/2), int(self.o2.get_height()/2), float(self.angulo_deg_giro), 1.0).dst
        dstt = Traslacion(dstg, int(self.desplazamiento_lateral), int(self.desplazamiento_vertical)).dst
        dstt_im = matriz_datos_2_matriz_im(dstt).matriz_im
        
        dst = Alfa(im1, dstt_im, float(self.valor_alpha.get())).dst
        cv2.imwrite('dust3.png', dst)
        self.Mostrare('dust3.png', self.Canvas3, 3, self.ancho_canvas_e_g, self.alto_canvas_e_g)
        
        if self.x_clicked == -1:
            self.x_clicked = self.ancho_canvas_e_g/2
            
        if self.y_clicked == -1:
            self.y_clicked = self.alto_canvas_e_g/2
            
        cc = Coordenadas_widget_to_imagen(int(self.x_clicked), int(self.x_clicked), self.ancho_canvas_e_g, self.alto_canvas_e_g, self.o1.get_width(),self.o1.get_height())
        
        pe = Perfil_ecualizacion(int(cc.x), int(cc.y), self.desplazamiento_lateral, self.desplazamiento_vertical, float(self.angulo_deg_giro), 
                            self.o1, self.o2, canal, float(self.factor_dosis_o2.get()))
        buffer_x1 = pe.px_1
        buffer_y1 = pe.py_1
        buffer_x2 = pe.px_2
        buffer_y2 = pe.py_2
        
        if canal == 0:
            cadena = "Rojo"
        else:
            if canal == 1:
                cadena = "Verde"
            else:
                cadena = "azul"
        
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis (" + str(cadena) + ")")
        plt.title("Perfil Vertical " + str(cc.x))
        plt.plot(buffer_x1, color = 'blue')
        plt.plot(buffer_x2, color = 'red')
        plt.savefig('salidav.png', bbox_inches='tight')
        self.Mostrare('salidav.png', self.Canvas5, 5, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis")
        plt.title("Perfil Horizontal " + str(cc.y))
        plt.plot(buffer_y1, color = 'blue')
        plt.plot(buffer_y2, color = 'red')
        plt.savefig('salidah.png', bbox_inches='tight')
        self.Mostrare('salidah.png', self.Canvas4, 4, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
    
    def close_windowse_aceptar(self):
        #Aplicar transformacion a imagen2 y actualizar datos de pixel, don y dosis si procede 
        buffer_3D = self.o2.get_dosis(3)
        M = cv2.getRotationMatrix2D((int(self.o2.get_width()/2), int(self.o2.get_height()/2)), float(self.angulo_deg_giro), 1.0)
        rows = len(buffer_3D[0])
        cols = len(buffer_3D[0][0])
        dstg = np.zeros((3, rows, cols))
        for i in range(3):
            dstg[i] = cv2.warpAffine(buffer_3D[i],M,(cols,rows))
        M = np.float32([[1,0,int(self.desplazamiento_lateral)],[0,1,int(self.desplazamiento_vertical)]])
        dstt = np.zeros((3, rows, cols))
        for i in range(3):
            dstt[i] = cv2.warpAffine(dstg[i], M, (cols,rows))
        self.o2.set_dosis(dstt)
            
        if self.o2.get_tipo() != 'pnc':
            buffer_3D = self.o2.get_pixel(3)
            M = cv2.getRotationMatrix2D((int(self.o2.get_width()/2), int(self.o2.get_height()/2)), float(self.angulo_deg_giro), 1.0)
            rows = len(buffer_3D[0])
            cols = len(buffer_3D[0][0])
            dstg = np.zeros((3, rows, cols))
            for i in range(3):
                dstg[i] = cv2.warpAffine(buffer_3D[i],M,(cols,rows))
            M = np.float32([[1,0,int(self.desplazamiento_lateral)],[0,1,int(self.desplazamiento_vertical)]])
            dstt = np.zeros((3, rows, cols))
            for i in range(3):
                dstt[i] = cv2.warpAffine(dstg[i], M, (cols,rows))
            self.o2.set_pixel(dstt)
            
            buffer_3D = self.o2.get_don(3)
            M = cv2.getRotationMatrix2D((int(self.o2.get_width()/2), int(self.o2.get_height()/2)), float(self.angulo_deg_giro), 1.0)
            rows = len(buffer_3D[0])
            cols = len(buffer_3D[0][0])
            dstg = np.zeros((3, rows, cols))
            for i in range(3):
                dstg[i] = cv2.warpAffine(buffer_3D[i],M,(cols,rows))
            M = np.float32([[1,0,int(self.desplazamiento_lateral)],[0,1,int(self.desplazamiento_vertical)]])
            dstt = np.zeros((3, rows, cols))
            for i in range(3):
                dstt[i] = cv2.warpAffine(dstg[i], M, (cols,rows))
            self.o2.set_don(dstt)
        
        #Presentar en escritorio principal
        if self.o1.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o1.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o1.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
            
        if self.o2.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o2.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o2.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
            
        self.desplazamiento_lateral = 0
        self.desplazamiento_vertical = 0
        self.angulo_deg_giro = 0
        self.newWindowe.destroy()
        self.root.deiconify()
    
    def mouse_clicked_e(self, event):
        cc = Coordenadas_widget_to_imagen(event.x, event.y, self.ancho_canvas_e_g, self.alto_canvas_e_g, self.o1.get_width(),self.o1.get_height())
        self.x_clicked = cc.x
        self.y_clicked = cc.x
        canal = int(self.indice_canal.get())
        if canal == 0:
            cadena = "Rojo"
        else:
            if canal == 1:
                cadena = "Verde"
            else:
                cadena = "Azul"
                
        self.newWindowe.title("(" + str(event.x) + ", " + str(event.y) + "): Tamaño. (" + str(self.o1.get_width()) + "," + str(self.o1.get_height()) +
                        "): (" + str(cc.x) + ", " + str(cc.y) + "): Canal " + str(cadena))
        
        self.Canvas3.delete("all")
        gi = Get_isodosis(self.o1.get_dosis(3), canal, float(self.valor_iso.get()), 2, (255, 0, 0), True)
        cv2.imwrite('dust1.png', gi.imagen_isodosis)
        im1 = cv2.imread('dust1.png')
        gi = Get_isodosis(self.o2.get_dosis(3), canal, float(self.valor_iso.get())*float(self.factor_dosis_o2.get()), 2, (255, 0, 0), True)
        cv2.imwrite('dust2.png', gi.imagen_isodosis)
        im2 = cv2.imread('dust2.png')
      
        dstg = Rotar_respecto_any_point(im2, int(self.o2.get_width()/2), int(self.o2.get_height()/2), float(self.angulo_deg_giro), 1.0).dst
        dstt = Traslacion(dstg, int(self.desplazamiento_lateral), int(self.desplazamiento_vertical)).dst
        dstt_im = matriz_datos_2_matriz_im(dstt).matriz_im
        
        dst = Alfa(im1, dstt_im, float(self.valor_alpha.get())).dst
        cv2.imwrite('dust3.png', dst)
        self.Mostrare('dust3.png', self.Canvas3, 3, self.ancho_canvas_e_g, self.alto_canvas_e_g)
        
        self.Canvas3.create_line(event.x, 0, event.x, self.alto_canvas_e_g, dash=(4, 2))
        self.Canvas3.create_line(0, event.y, self.ancho_canvas_e_g, event.y, dash=(4, 2))
    
        pe = Perfil_ecualizacion(int(cc.x), int(cc.y), self.desplazamiento_lateral, self.desplazamiento_vertical, float(self.angulo_deg_giro), 
                            self.o1, self.o2, canal, float(self.factor_dosis_o2.get()))
        buffer_x1 = pe.px_1
        buffer_y1 = pe.py_1
        buffer_x2 = pe.px_2
        buffer_y2 = pe.py_2
        
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis (" + str(cadena) + ")")
        plt.title("Perfil Vertical " + str(cc.x))
        plt.plot(buffer_x1, color = 'blue')
        plt.plot(buffer_x2, color = 'red')
        plt.savefig('salidav.png', bbox_inches='tight')
        self.Mostrare('salidav.png', self.Canvas5, 5, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
        plt.xlabel("pixel)")
        plt.ylabel("dosis")
        plt.title("Perfil Horizontal " + str(cc.y))
        plt.plot(buffer_y1, color = 'blue')
        plt.plot(buffer_y2, color = 'red')
        plt.savefig('salidah.png', bbox_inches='tight')
        self.Mostrare('salidah.png', self.Canvas4, 4, self.ancho_canvas_e_m, self.alto_canvas_e_m)
        plt.clf()
    
    def mouse_move_e(self, event):
        cc = Coordenadas_widget_to_imagen(event.x, event.y, self.ancho_canvas_e_g, self.alto_canvas_e_g, self.o1.get_width(),self.o1.get_height())
        canal = int(self.indice_canal.get())
        if canal == 0:
            cadena = "Rojo"
        else:
            if canal == 1:
                cadena = "Verde"
            else:
                cadena = "Azul"
                
        self.newWindowe.title("(" + str(event.x) + ", " + str(event.y) + "): Tamaño. (" + str(self.o1.get_width()) + "," + str(self.o1.get_height()) +
                        "): (" + str(cc.x) + ", " + str(cc.y) + "): Canal " + str(cadena))
        
        self.Canvas3.delete("all")
        self.Mostrare('dust3.png', self.Canvas3, 3, self.ancho_canvas_e_g, self.alto_canvas_e_g)
        self.Canvas3.create_line(event.x, 0, event.x, self.alto_canvas_e_g, dash=(4, 2))
        self.Canvas3.create_line(0, event.y, self.ancho_canvas_e_g, event.y, dash=(4, 2))
        
    def new_windowd(self, master):
        n_cols_d = 2
        n_rows_d = 1
        ancho_canvas_d_g = 700
        alto_canvas_d_g = 500
        padx_d = 3
        pady_d = 3
        ancho_ventana_d = 800
        alto_ventana_d = n_rows_d*alto_canvas_d_g + (n_rows_d+2)*pady_d
        self.ancho_canvas_d_g  = ancho_canvas_d_g
        self.alto_canvas_d_g  = alto_canvas_d_g
        
        # Comprueba que ambas imagenes estén cargadas, en dosis y están ecualizadas
        flag_comprobacion_d = True
        if self.o1 != None and self.o2 != None:
            if self.o1.get_tipo() != "pnc" and (self.o1.get_cal(0) == -1 or self.o1.get_modelo_ajuste() == -1):
                flag_comprobacion_d = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 1!!!" , title="Aviso!!!")
                self.close_windowsd_salir()
            if self.o2.get_tipo() != "pnc" and (self.o2.get_cal(0) == -1 or self.o2.get_modelo_ajuste() == -1):
                flag_comprobacion_d = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 2!!!" , title="Aviso!!!")
                self.close_windowsd_salir()
            if self.o1.get_factor_pixel_mm() != self.o2.get_factor_pixel_mm():
                flag_comprobacion_d = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsd_salir()
            if self.o1.get_height() != self.o2.get_height():
                flag_comprobacion_d = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsd_salir()
            if self.o1.get_width() != self.o2.get_width():
                flag_comprobacion_d = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsd_salir()
                
            if flag_comprobacion_d == True:
                messagebox.showinfo(message="Se va a calcular la diferencia!!!" , title="Aviso!!!")
                self.newWindowd = tk.Toplevel(master)
                geom = str(ancho_ventana_d) + "x" + str(alto_ventana_d)
                self.newWindowd.geometry(geom)
                self.newWindowd.grid()
                self.Canvas1d = tk.Canvas(self.newWindowd, bg = 'green', height=alto_canvas_d_g, width=ancho_canvas_d_g)
                self.Canvas1d.place(x=padx_d, y=pady_d)
                self.Canvas1d.bind("<Motion>", self.mouse_move_d)
                
                self.salirButtond = tk.Button(self.newWindowd, text = 'Salir', width = 5, command = self.close_windowsd_salir)
                self.salirButtond.place(x=ancho_ventana_d - 65, y=alto_ventana_d - 50)
                
                self.Lcanald = tk.Label(self.newWindowd)
                self.Lcanald.place(x = 740, y = 10)
                self.Lcanald['text'] = "Canal"
                self.canald = tk.StringVar()
                self.Ecanald = tk.Entry(self.newWindowd, width=3, textvariable=self.canald)
                self.Ecanald.place(x = 750, y = 30)
                self.canald.set("0")
            
                master.withdraw()

                #Diferencia
                self.dd = Get_diferencia(self.o1.get_dosis(int(self.canald.get())), self.o2.get_dosis(int(self.canald.get())))
                cv2.imwrite('diferencia.png', self.dd.matriz_diferencia)
                self.Mostrare('diferencia.png', self.Canvas1d, 1, self.ancho_canvas_d_g, self.alto_canvas_d_g)

        else:
            messagebox.showinfo(message="¡¡¡Imagenes no cargadas!!!" , title="Aviso!!!")
    
    def mouse_move_d(self, event):
        cc = Coordenadas_widget_to_imagen(event.x, event.y, self.ancho_canvas_d_g, self.alto_canvas_d_g, self.o1.get_width(),self.o1.get_height())
        canal = int(self.canald.get())
        if canal == 0:
            cadena = "Rojo"
        else:
            if canal == 1:
                cadena = "Verde"
            else:
                cadena = "Azul"
                
        #Diferencia
        self.dd = Get_diferencia(self.o1.get_dosis(int(self.canald.get())), self.o2.get_dosis(int(self.canald.get())))
        cv2.imwrite('diferencia.png', self.dd.matriz_diferencia)
        self.Mostrare('diferencia.png', self.Canvas1d, 1, self.ancho_canvas_d_g, self.alto_canvas_d_g)
        
        cc.x = min(cc.x, self.o1.get_width()-1)
        cc.y = min(cc.y, self.o1.get_height()-1)
        self.newWindowd.title("(" + str(event.x) + ", " + str(event.y) + "): Tamaño. (" + str(self.o1.get_width()) + ", " + str(self.o1.get_height()) +
                        "): (" + str(cc.x) + ", " + str(cc.y) + "): Canal " + str(cadena) + ": Diferencia. " + str(self.dd.matriz_diferencia[int(cc.y)][int(cc.x)]))
    
    def new_windowg(self, master):
        n_cols_g = 2
        n_rows_g = 1
        ancho_canvas_g_g = 500
        alto_canvas_g_g = 500
        padx_g = 3
        pady_g = 3
        ancho_ventana_g = n_cols_g*ancho_canvas_g_g + (n_cols_g+2)*padx_g
        alto_ventana_g = n_rows_g*alto_canvas_g_g + (n_rows_g+2)*pady_g
        self.ancho_canvas_g_g  = ancho_canvas_g_g
        self.alto_canvas_g_g  = alto_canvas_g_g
        
        rango = [0.5, 0.75, 1.0, 1.5]
        colores = [(255, 0, 0), (0, 255, 0), (125, 125, 125), (0, 255, 255), (0, 0, 255)]
        nombres = ["Azul", "Verde", "Gris", "Amarillo", "Rojo"]
       
        # Comprueba que ambas imagenes estén cargadas, en dosis y están ecualizadas
        flag_comprobacion_g = True
        if self.o1 != None and self.o2 != None:
            if self.o1.get_tipo() != "pnc" and (self.o1.get_cal(0) == -1 or self.o1.get_modelo_ajuste() == -1):
                flag_comprobacion_g = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 1!!!" , title="Aviso!!!")
                self.close_windowsg_salir()
            if self.o2.get_tipo() != "pnc" and (self.o2.get_cal(0) == -1 or self.o2.get_modelo_ajuste() == -1):
                flag_comprobacion_g = False
                messagebox.showinfo(message="¡¡¡las imagenes tif deben estar calibradas. Imagen 2!!!" , title="Aviso!!!")
                self.close_windowsg_salir()
            if self.o1.get_factor_pixel_mm() != self.o2.get_factor_pixel_mm():
                flag_comprobacion_g = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsg_salir()
            if self.o1.get_height() != self.o2.get_height():
                flag_comprobacion_g = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsg_salir()
            if self.o1.get_width() != self.o2.get_width():
                flag_comprobacion_g = False
                messagebox.showinfo(message="¡¡¡Ecualiza las imagenes!!!" , title="Aviso!!!")
                self.close_windowsg_salir()
        
            if flag_comprobacion_g == True:
                messagebox.showinfo(message="Se va a calcular la gamma!!!" , title="Aviso!!!")
                self.newWindowg = tk.Toplevel(master)
                geom = str(ancho_ventana_g) + "x" + str(alto_ventana_g)
                self.newWindowg.geometry(geom)
                self.newWindowg.grid()
                self.Canvas1g = tk.Canvas(self.newWindowg, bg = 'green', height=alto_canvas_g_g, width=ancho_canvas_g_g)
                self.Canvas1g.place(x=padx_g, y=pady_g)
                self.Canvas1g.bind("<Motion>", self.mouse_move_g)
        
            self.Ltol_dosis = tk.Label(self.newWindowg)
            self.Ltol_dosis.place(x = ancho_canvas_g_g + 30, y = 10)
            self.Ltol_dosis['text'] = "Tolerancia dosis (cGy)"
            self.Ltol_distancia = tk.Label(self.newWindowg)
            self.Ltol_distancia.place(x = ancho_canvas_g_g + 180,y = 10)
            self.Ltol_distancia['text'] = "Tolerancia distancia (mm)"
            self.tol_dosis = tk.StringVar()
            self.Etol_dosis = tk.Entry(self.newWindowg, width=3, textvariable=self.tol_dosis)
            self.Etol_dosis.place(x = ancho_canvas_g_g + 70, y = 30)
            self.tol_dosis.set("2")
            self.tol_distancia = tk.StringVar()
            self.Etol_distancia = tk.Entry(self.newWindowg, width=3, textvariable=self.tol_distancia)
            self.Etol_distancia.place(x = ancho_canvas_g_g + 220, y = 30)
            self.tol_distancia.set("2")
            self.Lcanal = tk.Label(self.newWindowg)
            self.Lcanal.place(x = ancho_canvas_g_g + 410, y = 10)
            self.Lcanal['text'] = "Canal"
            self.canalg = tk.StringVar()
            self.Ecanalg = tk.Entry(self.newWindowg, width=3, textvariable=self.canalg)
            self.Ecanalg.place(x = ancho_canvas_g_g + 420, y = 30)
            self.canalg.set("0")
            self.Lcodigo0 = tk.Label(self.newWindowg)
            self.Lcodigo0.place(x = ancho_canvas_g_g + 50, y = 70)
            self.Lcodigo0['text'] = "Interpretacion codigo colores:"
            
            
            self.Lcodigo1 = tk.Label(self.newWindowg)
            self.Lcodigo1.place(x = ancho_canvas_g_g + 70, y = 90)
            cadena = "<----> "
            for i in rango:
                cadena += str(i) + " <---> "
            self.Lcodigo1['text'] = cadena
            
            self.Lcodigo2 = tk.Label(self.newWindowg)
            self.Lcodigo2.place(x = ancho_canvas_g_g + 80, y = 110)
            cadena = ""
            for i in nombres:
                cadena += str(i) + "        "
            self.Lcodigo2['text'] = cadena
            
            self.Sipasa = tk.Label(self.newWindowg)
            self.Sipasa.place(x = ancho_canvas_g_g + 50, y = 130)
            self.Sipasa['text'] = ""
            
            self.Nopasa = tk.Label(self.newWindowg)
            self.Nopasa.place(x = ancho_canvas_g_g + 50, y = 150)
            self.Nopasa['text'] = ""
        
            self.salirButtong = tk.Button(self.newWindowg, text = 'Salir', width = 5, command = self.close_windowsg_salir)
            self.salirButtong.place(x=ancho_ventana_g - 65, y=alto_ventana_g - 50)
            self.actualizarButtong = tk.Button(self.newWindowg, text = 'Actualizar gamma', width = 20, command = lambda: self.regamma(ancho_canvas_g_g, alto_canvas_g_g,
                                                                                                                             colores, rango, nombres))
            self.actualizarButtong.place(x=ancho_ventana_g - 350, y=alto_ventana_g - 335)
        
            master.withdraw()
        
            self.gg = Get_gamma(self.o1.get_dosis(int(self.canalg.get())), self.o2.get_dosis(int(self.canalg.get())),
                           int(self.tol_dosis.get()), int(int(self.tol_distancia.get())*float(self.o1.get_factor_pixel_mm())), 
                           colores, rango, nombres)
            self.Sipasa['text'] = "% Puntos cumplen criterio: " + str(round(self.gg.pasa_tantos_por_cien, 2))
            self.Nopasa['text'] = "% Puntos no cumplen criterio: " + str(round(self.gg.no_pasa_tantos_por_cien, 2))
            cv2.imwrite("gamma.png", self.gg.imagen_gamma)
            self.Mostrare('gamma.png', self.Canvas1g, 1, self.ancho_canvas_g_g, self.alto_canvas_g_g)
            
        else:
            messagebox.showinfo(message="¡¡¡Imagenes no cargadas!!!" , title="Aviso!!!")
            
    def close_windowse_salir(self):#Presentar en escritorio principal
        if self.o1.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o1.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o1.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
            
        if self.o2.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o2.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o2.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
        
        self.desplazamiento_lateral = 0
        self.desplazamiento_vertical = 0
        self.angulo_deg_giro = 0
        self.newWindowe.destroy()
        self.root.deiconify()
        
    def regamma(self, ancho_canvas_g_g, alto_canvas_g_g, colores, rango, nombres):
        self.gg = Get_gamma(self.o1.get_dosis(int(self.canalg.get())), self.o2.get_dosis(int(self.canalg.get())),
                           int(self.tol_dosis.get()), int(int(self.tol_distancia.get())*float(self.o1.get_factor_pixel_mm())), 
                           colores, rango, nombres)
        self.Sipasa['text'] = "% Puntos cumplen criterio: " + str(round(self.gg.pasa_tantos_por_cien, 2))
        self.Nopasa['text'] = "% Puntos no cumplen criterio: " + str(round(self.gg.no_pasa_tantos_por_cien, 2))
        cv2.imwrite("gamma.png", self.gg.imagen_gamma)
        self.Canvas1g.delete("all")
        self.Mostrare('gamma.png', self.Canvas1g, 1, ancho_canvas_g_g, alto_canvas_g_g)
        
    def close_windowsd_salir(self):
        if self.o1.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o1.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o1.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
            
        if self.o2.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o2.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o2.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
        
        self.newWindowd.destroy()
        self.root.deiconify()
        
    def close_windowsg_salir(self):
        if self.o1.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o1.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o1.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
            
        if self.o2.get_tipo() == 'pnc':
            dst_im = matriz_datos_2_matriz_im(self.o2.get_dosis(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
        else:
            dst_im = matriz_datos_2_matriz_im(self.o2.get_pixel(3)).matriz_im
            self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
        
        self.newWindowg.destroy()
        self.root.deiconify()
        
    def mouse_move_g(self, event):
        cc = Coordenadas_widget_to_imagen(event.x, event.y, self.ancho_canvas_g_g, self.alto_canvas_g_g, self.o1.get_width(),self.o1.get_height())
                
        cc.x = min(cc.x, self.o1.get_width()-1)
        cc.y = min(cc.y, self.o1.get_height()-1)
        self.newWindowg.title("(" + str(event.x) + " , " + str(event.y) + "): Tamaño. (" + str(self.o1.get_width()) + ", " + str(self.o1.get_height()) +
                        "): (" + str(cc.x) + ", " + str(cc.y) + "): Indice_gamma. " + str(self.gg.valores_gamma[int(cc.y)][int(cc.x)]))
    
    def recortar_ambas(self):
        global flag_recorte
        global flag_recorte_ambas
        
        flag_check = True
        if self.o1 == None or self.o2 == None:
            flag_check = False
            messagebox.showinfo(message="¡¡¡Imagenes no cargadas!!!" , title="Aviso!!!")
        else:
            if self.o1.get_width() != self.o2.get_width() or self.o1.get_height() != self.o2.get_height():
                flag_check = False
                messagebox.showinfo(message="¡¡¡Imagenes no ecualizadas!!!" , title="Aviso!!!")
            
            if self.o1.get_factor_pixel_mm() != self.o2.get_factor_pixel_mm() :
                flag_check = False
                messagebox.showinfo(message="¡¡¡Imagenes no ecualizadas!!!" , title="Aviso!!!")    
            
            if self.o1.get_tipo() != "pnc" and (self.o1.get_cal(0) == -1 or self.o1.get_modelo_ajuste() == -1):
                flag_check = False
                messagebox.showinfo(message="¡¡¡Imagen Izda no calibrada!!!" , title="Aviso!!!")    
            
            if self.o2.get_tipo() != "pnc" and (self.o2.get_cal(0) == -1 or self.o2.get_modelo_ajuste() == -1):
                flag_check = False
                messagebox.showinfo(message="¡¡¡Imagen Dcha no calibrada!!!" , title="Aviso!!!")
        
        if flag_check == True:
            flag_recorte_ambas = True
            flag_recorte = True
        
    def popup(self, event):
        global foco
        
        try:
            self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
            buffer = str(event.widget)
            buffer = buffer[len(buffer)-1]
            if buffer == "2":
                foco = 2
            else:
                foco = 1
        finally:
            self.popup_menu.grab_release()
        
    def gestion(self):
        pass
        
    def abrir_selected(self, foco):
        
        if foco == 1:
            answer = askopenfilename(filetypes=[("Image files", ".tif .pnc")])
            if answer != "":
                self.o1 = Objeto_Ppal()
                extension = Get_extension(answer).extension
                if extension == 'pnc':
                    lp = Lee_PNC(answer)
                    im_datos = lp.dosis_pnc_3D
                    im = matriz_datos_2_matriz_im(im_datos).matriz_im
                    self.Mostrar(im, 'dust1.png', self.canvas1_ppal, foco, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
                    self.o1.set_tipo('pnc')
                    self.o1.set_factor_pixel_mm(lp.pixel_size_mm)
                    self.o1.set_height(lp.height_imagen_pixel)
                    self.o1.set_width(lp.width_imagen_pixel)
                    self.o1.set_dosis(lp.dosis_pnc_3D)
                    self.o1.set_don(None)
                    self.o1.set_pixel(None)
                else:                    
                    im = cv2.imread(answer, -1)
                    self.o1.fdo = FilmDoseClass.FilmDose(answer)

                    median_image = cv2.medianBlur(im, 3)
                    self.Mostrar(median_image, 'dust1.png', self.canvas1_ppal, foco, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                    self.o1.set_tipo('tif')
                    height, width = im.shape[:2]
                    self.o1.set_height(height)
                    self.o1.set_width(width)
                    self.o1.set_factor_pixel_mm(Get_dpi_tif(answer).factor_pixel_mm)
                    self.o1.set_cal(0, -1)
                    self.o1.set_cal(1, -1)
                    self.o1.set_cal(2, -1)
                    pixel = np.zeros((3, self.o1.get_height(), self.o1.get_width()))
                    pixel[0] = median_image[:, :, 2]
                    pixel[1] = median_image[:, :, 1]
                    pixel[2] = median_image[:, :, 0]
                    self.o1.set_dosis(None)
                    self.o1.set_don(None)
                    self.o1.set_pixel(pixel)
                    
        if foco == 2:
            answer = askopenfilename(filetypes=[("Image files", ".tif .pnc")])
            if answer != "":
                self.o2 = Objeto_Ppal()
                extension = Get_extension(answer).extension
                if extension == 'pnc':
                    lp = Lee_PNC(answer)
                    im_datos = lp.dosis_pnc_3D
                    im = matriz_datos_2_matriz_im(im_datos).matriz_im
                    self.Mostrar(im, 'dust2.png', self.canvas2_ppal, foco, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
                    self.o2.set_tipo('pnc')
                    self.o2.set_factor_pixel_mm(lp.pixel_size_mm)
                    self.o2.set_height(lp.height_imagen_pixel)
                    self.o2.set_width(lp.width_imagen_pixel)
                    self.o2.set_dosis(lp.dosis_pnc_3D)
                    self.o2.set_don(None)
                    self.o2.set_pixel(None)
                else:
                    im = cv2.imread(answer, -1)
                    self.o2.fdo = FilmDoseClass.FilmDose(answer)
                    self.o2.fdo.AddCalibrationFromFile('_')
                    median_image = cv2.medianBlur(im, 3)
                    self.Mostrar2(self.o2.fdo.get_nod(), self.canvas2_ppal, foco, self.ancho_canvas_ppal,
                                  self.alto_canvas_ppal, "tif")
                    #self.Mostrar(median_image, 'dust2.png', self.canvas2_ppal, foco, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                    self.o2.set_tipo('tif')
                    height, width = im.shape[:2]
                    self.o2.set_height(height)
                    self.o2.set_width(width)
                    self.o2.set_factor_pixel_mm(Get_dpi_tif(answer).factor_pixel_mm)
                    self.o2.set_cal(0, -1)
                    self.o2.set_cal(1, -1)
                    self.o2.set_cal(2, -1)
                    pixel = np.zeros((3, self.o2.get_height(), self.o2.get_width()))
                    pixel[0] = median_image[:, :, 2]
                    pixel[1] = median_image[:, :, 1]
                    pixel[2] = median_image[:, :, 0]
                    self.o2.set_dosis(None)
                    self.o2.set_don(None)
                    self.o2.set_pixel(pixel)
    
    def flat_scanner(self, foco):
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() != "pnc":
                pass
            else:
                messagebox.showinfo(message="¡¡¡ Sólo para imágenes tiff!!!" , title="Aviso!!!")
                
        if foco == 2 and self.o2 != None:
            if self.o2.get_tipo() != "pnc":
                pass
            else:
                messagebox.showinfo(message="¡¡¡ Sólo para imágenes tiff!!!" , title="Aviso!!!")
                
    def establecer_velo(self, foco):
        global flag_velo_1
        global flag_velo_2
        
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() != "pnc":
                flag_velo_1 = True
                messagebox.showinfo(message="¡¡¡Sitúese sobre el velo y haga click con el botón izquierdo del ratón!!!" , title="Aviso!!!")
            else:
                messagebox.showinfo(message="¡¡¡Sólo para imágenes tiff!!!" , title="Aviso!!!")
                
        if foco == 2 and self.o2 != None:
            if self.o2.get_tipo() != "pnc":
                flag_velo_2 = True
                messagebox.showinfo(message="¡¡¡Sitúese sobre el velo y haga click con el botón izquierdo del ratón!!!" , title="Aviso!!!")
            else:
                messagebox.showinfo(message="¡¡¡Sólo para imágenes tiff!!!" , title="Aviso!!!")
                
    def establecer_dosis_conocida(self, foco):
        global flag_dosis_conocida_1
        global flag_dosis_conocida_2
        
        if foco == 1 and self.o1 != None:
                 if self.o1.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡En las pnc no procede!!!" , title="Aviso!!!")
                 else: # tif
                     flag_dosis_conocida_1 = True
                     valor_dosis_conocida = askstring(title = 'Dosis conocida?', prompt = "(cGy)", initialvalue="200")
                     self.o1.set_dosis_linealizacion(float(valor_dosis_conocida))
                     messagebox.showinfo(message="¡¡¡Sitúese sobre la dosis conocida y haga click con el botón izquierdo del ratón!!!" , title="Aviso!!!")
                     
        if foco == 2 and self.o2 != None:
                 if self.o2.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡En las pnc no procede!!!" , title="Aviso!!!")
                 else: # tif
                     flag_dosis_conocida_2 = True
                     valor_dosis_conocida = askstring(title = 'Dosis conocida?', prompt = "(cGy)", initialvalue="200")
                     self.o2.set_dosis_linealizacion(float(valor_dosis_conocida))
                     messagebox.showinfo(message="¡¡¡Sitúese sobre la dosis conocida y haga click con el botón izquierdo del ratón!!!" , title="Aviso!!!")
        
    def establecer_curva_calibracion(self, foco):
        global fichero_dosimetria
        
        if foco == 1 and self.o1 != None:
                if self.o1.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡Las pnc no necesitan calibracion!!!" , title="Aviso!!!")
                else:
                    if self.o1.get_pixel_velo(0) == -1:                        
                        messagebox.showinfo(message="¡¡¡Primero, establece el velo de la pelicula!!!", title="Aviso!!!")
                    else:
                        messagebox.showinfo(message="Se va a obtener la matriz DoN ...", title="Aviso!!!")
                        don_0 = Pixel_to_DoN_2(self.o1.get_pixel_velo(0), self.o1.get_pixel(0)).DoN
                        don_1 = Pixel_to_DoN_2(self.o1.get_pixel_velo(1), self.o1.get_pixel(1)).DoN
                        don_2 = Pixel_to_DoN_2(self.o1.get_pixel_velo(2), self.o1.get_pixel(2)).DoN
                        don = np.zeros((3, len(don_0), len(don_0[0])))
                        don[0] = don_0
                        don[1] = don_1
                        don[2] = don_2
                        self.o1.set_don(don)
                        messagebox.showinfo(message="Selecciona curva para canal rojo ...(verde y azul consecutivas)", title="Aviso!!!")
                        Selecciona_curva(foco, self.canvas1_ppal, self.o1)
                        
        if foco == 2 and self.o2 != None:
                if self.o2.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡Las pnc no necesitan calibracion!!!" , title="Aviso!!!")
                else:
                    if self.o2.get_pixel_velo(0) == -1:                        
                        messagebox.showinfo(message="¡¡¡Primero, establece el velo de la pelicula!!!", title="Aviso!!!")
                    else:
                        messagebox.showinfo(message="Se va a obtener la matriz DoN ...", title="Aviso!!!")
                        don_0 = Pixel_to_DoN_2(self.o2.get_pixel_velo(0), self.o2.get_pixel(0)).DoN
                        don_1 = Pixel_to_DoN_2(self.o2.get_pixel_velo(1), self.o2.get_pixel(1)).DoN
                        don_2 = Pixel_to_DoN_2(self.o2.get_pixel_velo(2), self.o2.get_pixel(2)).DoN
                        don = np.zeros((3, len(don_0), len(don_0[0])))
                        don[0] = don_0
                        don[1] = don_1
                        don[2] = don_2
                        self.o2.set_don(don)
                        messagebox.showinfo(message="Selecciona curva para canal rojo ...(verde y azul consecutivas)", title="Aviso!!!")
                        Selecciona_curva(foco, self.canvas2_ppal, self.o2)
        
    def establecer_modelo_calibracion(self, foco):
        if foco == 1 and self.o1 != None:
                if self.o1.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡Las pnc no necesitan calibracion!!!" , title="Aviso!!!")
                else:
                    messagebox.showinfo(message="Selecciona modelo de calibración para ajuste Dosis=f(DoNeta)", title="Aviso!!!")
                    Selecciona_modelo(foco, self.canvas1_ppal, self.o1)
                    
        if foco == 2 and self.o2 != None:
                if self.o2.get_tipo() == 'pnc':
                    messagebox.showinfo(message="¡¡¡Las pnc no necesitan calibracion!!!" , title="Aviso!!!")
                else:
                    messagebox.showinfo(message="Selecciona modelo de calibración para ajuste Dosis=f(DoNeta)", title="Aviso!!!")
                    Selecciona_modelo(foco, self.canvas2_ppal, self.o2)
    
    def linealizacion(self, foco):
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() != 'pnc':
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                    if self.o1.get_modelo_ajuste() != -1:
                        if self.o1.get_pixel_dosis_conocida(0) != -1 and self.o1.get_pixel_dosis_conocida(1) != -1 and self.o1.get_pixel_dosis_conocida(2) != -1:
                            gk = Get_K(self.o1)
                            DoN_old = self.o1.get_don(3)
                            DoN_new = np.zeros((3, len(DoN_old[0]), len(DoN_old[0][0])))
                            for k in range(0,3):
                                DoN_new[k] = DoN_old[k] * gk.K[k]
                            self.o1.set_don(DoN_new)
                            messagebox.showinfo(message="K: " + str(gk.K[0]) + " " + str(gk.K[1]) + " " + str(gk.K[2]) , title="Aviso!!!")
                            #new_pixel = DoN_to_Pixel(self.o1.get_don(3), self.o1.get_pixel_velo(3)).pixel
                            #self.o1.set_pixel(new_pixel)
                            lc = Lee_curva_calibracion(self.o1.get_cal(0), self.o1.get_modelo_ajuste(), False, False, False)
                            dosis = DoN3D_to_Dosis(lc, self.o1).Dosis
                            self.o1.set_dosis(dosis)
                            messagebox.showinfo(message="Imagen linealizada!!!", title="Aviso!!!")
                        else:
                            messagebox.showinfo(message="¡¡¡Establece pixel de dosis conocida!!!" , title="Aviso!!!")
                    else:
                        messagebox.showinfo(message="¡¡¡Selecciona modelo de calibración!!!" , title="Aviso!!!")
                else:
                    messagebox.showinfo(message="¡¡¡Imagen no calibrada!!!" , title="Aviso!!!")
            else:
                messagebox.showinfo(message="¡¡¡Las pnc no necesitan linealización!!!" , title="Aviso!!!")
                
        if foco == 2 and self.o2 != None:
            if self.o2.get_tipo() != 'pnc':
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                    if self.o2.get_modelo_ajuste() != -1:
                        if self.o2.get_pixel_dosis_conocida(0) != -1 and self.o2.get_pixel_dosis_conocida(1) != -1 and self.o2.get_pixel_dosis_conocida(2) != -1:
                            gk = Get_K(self.o2)
                            DoN_old = self.o2.get_don(3)
                            DoN_new = np.zeros((3, len(DoN_old[0]), len(DoN_old[0][0])))
                            for k in range(0,3):
                                DoN_new[k] = DoN_old[k] * gk.K[k]
                            self.o2.set_don(DoN_new)
                            messagebox.showinfo(message="K: " + str(gk.K[0]) + " " + str(gk.K[1]) + " " + str(gk.K[2]) , title="Aviso!!!")
                            #new_pixel = DoN_to_Pixel(self.o2.get_don(3), self.o2.get_pixel_velo(3)).pixel
                            #self.o2.set_pixel(new_pixel)
                            lc = Lee_curva_calibracion(self.o2.get_cal(0), self.o2.get_modelo_ajuste(), False, False, False)
                            dosis = DoN3D_to_Dosis(lc, self.o2).Dosis
                            self.o2.set_dosis(dosis)
                            messagebox.showinfo(message="Imagen linealizada!!!", title="Aviso!!!")
                        else:
                          messagebox.showinfo(message="¡¡¡Establece pixel de dosis conocida!!!" , title="Aviso!!!")  
                    else:
                        messagebox.showinfo(message="¡¡¡Selecciona modelo de calibración!!!" , title="Aviso!!!")
                else:
                    messagebox.showinfo(message="¡¡¡Imagen no calibrada!!!" , title="Aviso!!!")
            else:
                messagebox.showinfo(message="¡¡¡Las pnc no necesitan linealización!!!" , title="Aviso!!!")
    
    def multicanal(self, foco):
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() != 'pnc':
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1 and self.o1.get_modelo_ajuste() != -1:
                    #factor = Multicanal_1(self.o1).factor
                    #new_doN = self.o1.get_don(3)*factor
                    #self.o1.set_don(new_doN)
                    lc = Lee_curva_calibracion(self.o1.get_cal(0), self.o1.get_modelo_ajuste(), False, False, False)
                    dosis = DoN3D_to_Dosis(lc, self.o1).Dosis
                    self.o1.set_dosis(dosis)
                    messagebox.showinfo(message="Multicanal aplicado!!!", title="Aviso!!!")
                    #new_pixel = DoN_to_Pixel(new_doN, self.o1.get_pixel_velo(3)).pixel
                    #self.o1.set_pixel(new_pixel)
                else:
                    messagebox.showinfo(message="¡¡¡La imagen tif debe estar calibrada!!!" , title="Aviso!!!")
            else:
                messagebox.showinfo(message="¡¡¡Las pnc no aplican multicanal!!!" , title="Aviso!!!")
                
            if foco == 2 and self.o2 != None:
                if self.o2.get_tipo() != 'pnc':
                    pass
                else:
                    messagebox.showinfo(message="¡¡¡Las pnc no aplican multicanal!!!" , title="Aviso!!!")
    
    def flip(self, foco, args):
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() == "pnc":
                dst = Flip(self.o1.get_dosis(3), args).dst
                self.o1.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Flip(self.o1.get_pixel(3), args).dst
                self.o1.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                   dst = Flip(self.o1.get_don(3), args).dst
                   self.o1.set_don(dst)
                   if self.o1.get_modelo_ajuste() != -1:
                       dst = Flip(self.o1.get_dosis(3), args).dst
                       self.o1.set_dosis(dst)
                
        if foco == 2 and self.o2 != None:
            if self.o2.get_tipo() == "pnc":
                dst = Flip(self.o2.get_dosis(3), args).dst
                self.o2.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Flip(self.o2.get_pixel(3), args).dst
                self.o2.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                   dst = Flip(self.o2.get_don(3), args).dst
                   self.o2.set_don(dst)
                   if self.o2.get_modelo_ajuste() != -1:
                       dst = Flip(self.o2.get_dosis(3), args).dst
                       self.o2.set_dosis(dst)
    
    def rota90(self, foco, args):
        if foco == 1 and self.o1 != None:
            if self.o1.get_tipo() == "pnc":
                dst = Rota_90(self.o1.get_dosis(3), args).dst
                self.o1.set_height(len(dst[0]))
                self.o1.set_width(len(dst[0][0]))
                self.o1.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Rota_90(self.o1.get_pixel(3), args).dst
                self.o1.set_height(len(dst[0]))
                self.o1.set_width(len(dst[0][0]))
                self.o1.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                   dst = Rota_90(self.o1.get_don(3), args).dst
                   self.o1.set_don(dst)
                   if self.o1.get_modelo_ajuste() != -1:
                       dst = Rota_90(self.o1.get_dosis(3), args).dst
                       self.o1.set_dosis(dst)
                
        if foco == 2 and self.o2 != None:
            if self.o2.get_tipo() == "pnc":
                dst = Rota_90(self.o2.get_dosis(3), args).dst
                self.o2.set_height(len(dst[0]))
                self.o2.set_width(len(dst[0][0]))
                self.o2.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Rota_90(self.o2.get_pixel(3), args).dst
                self.o2.set_height(len(dst[0]))
                self.o2.set_width(len(dst[0][0]))
                self.o2.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                   dst = Rota_90(self.o2.get_don(3), args).dst
                   self.o2.set_don(dst)
                   if self.o2.get_modelo_ajuste() != -1:
                       dst = Rota_90(self.o2.get_dosis(3), args).dst
                       self.o2.set_dosis(dst)
    
    def traslacion(self, foco, args):
        if foco == 1 and self.o1 != None:
            if args == 1:
                cadena = "Horizontal"
            else:
                cadena = "Vertical"
            desplazamiento_pixeles = askstring(cadena, 'Desplazamiento (pixeles)?')
            if self.o1.get_tipo() == 'pnc':
                if args == 1:
                    dst = Traslacion(self.o1.get_dosis(3), desplazamiento_pixeles, 0).dst
                else:
                    dst = Traslacion(self.o1.get_dosis(3), 0, desplazamiento_pixeles).dst
                self.o1.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                if args == 1:
                    dst = Traslacion(self.o1.get_pixel(3), desplazamiento_pixeles, 0).dst
                else:
                    dst = Traslacion(self.o1.get_pixel(3), 0, desplazamiento_pixeles).dst
                self.o1.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                    if args == 1:
                        dst = Traslacion(self.o1.get_don(3), desplazamiento_pixeles, 0).dst
                    else:
                        dst = Traslacion(self.o1.get_don(3), 0, desplazamiento_pixeles).dst
                    self.o1.set_don(dst)
                    if self.o1.get_modelo_ajuste() != -1:
                        if args == 1:
                            dst = Traslacion(self.o1.get_dosis(3), desplazamiento_pixeles, 0).dst
                        else:
                            dst = Traslacion(self.o1.get_dosis(3), 0, desplazamiento_pixeles).dst
                        self.o1.set_dosis(dst)
                    
        if foco == 2 and self.o2 != None:
                if args == 1:
                    cadena = "Horizontal"
                else:
                    cadena = "Vertical"
                desplazamiento_pixeles = askstring(cadena, 'Desplazamiento (pixeles)?')
                if self.o2.get_tipo() == 'pnc':
                    if args == 1:
                        dst = Traslacion(self.o2.get_dosis(3), desplazamiento_pixeles, 0).dst
                    else:
                        dst = Traslacion(self.o2.get_dosis(3), 0, desplazamiento_pixeles).dst
                    self.o2.set_dosis(dst)
                    dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                    self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
                else:
                    if args == 1:
                        dst = Traslacion(self.o2.get_pixel(3), desplazamiento_pixeles, 0).dst
                    else:
                        dst = Traslacion(self.o2.get_pixel(3), 0, desplazamiento_pixeles).dst
                self.o2.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                    if args == 1:
                        dst = Traslacion(self.o2.get_don(3), desplazamiento_pixeles, 0).dst
                    else:
                        dst = Traslacion(self.o2.get_don(3), 0, desplazamiento_pixeles).dst
                    self.o2.set_don(dst)
                    if self.o2.get_modelo_ajuste() != -1:
                        if args == 1:
                            dst = Traslacion(self.o2.get_dosis(3), desplazamiento_pixeles, 0).dst
                        else:
                            dst = Traslacion(self.o2.get_dosis(3), 0, desplazamiento_pixeles).dst
                        self.o2.set_dosis(dst)
                    
    def rotar_centro(self, foco):
        if foco == 1 and self.o1 != None:
            angulo_deg = askstring("Rotar respecto centro", 'Angulo (º)?')
            if self.o1.get_tipo() == 'pnc':
                dst = Rotar_centro(self.o1.get_dosis(3), float(angulo_deg)).dst
                self.o1.set_height(len(dst[0]))
                self.o1.set_width(len(dst[0][0]))
                self.o1.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Rotar_centro(self.o1.get_pixel(3), float(angulo_deg)).dst
                self.o1.set_height(len(dst[0]))
                self.o1.set_width(len(dst[0][0]))
                self.o1.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                    dst = Rotar_centro(self.o1.get_don(3), float(angulo_deg)).dst
                    self.o1.set_don(dst)
                    if self.o1.get_modelo_ajuste() != -1:
                        dst = Rotar_centro(self.o1.get_dosis(3), float(angulo_deg)).dst
                        self.o1.set_dosis(dst)
                        
        if foco == 2 and self.o2 != None:
            angulo_deg = askstring("Rotar respecto centro", 'Angulo (º)?')
            if self.o2.get_tipo() == 'pnc':
                dst = Rotar_centro(self.o2.get_dosis(3), float(angulo_deg)).dst
                self.o2.set_height(len(dst[0]))
                self.o2.set_width(len(dst[0][0]))
                self.o2.set_dosis(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Rotar_centro(self.o2.get_pixel(3), float(angulo_deg)).dst
                self.o2.set_height(len(dst[0]))
                self.o2.set_width(len(dst[0][0]))
                self.o2.set_pixel(dst)
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                    dst = Rotar_centro(self.o2.get_don(3), float(angulo_deg)).dst
                    self.o2.set_don(dst)
                    if self.o2.get_modelo_ajuste() != -1:
                        dst = Rotar_centro(self.o2.get_dosis(3), float(angulo_deg)).dst
                        self.o2.set_dosis(dst)
                
    def recortar(self, foco):
        global flag_recorte
        
        if foco == 1 and self.o1 != None: # Canvas de la izquierda
            flag_recorte = True
                     
        if foco == 2 and self.o2 != None: # Canvas de la derecha
            flag_recorte = True
        
    def guardar(self, foco):
        if foco == 1 and self.o1 != None:
            canal = askstring(title = 'Canal de dosis?', prompt = "Rojo, Verde, Azul -> (0, 1, 2)", initialvalue="0")
            if self.o1.get_tipo() == "pnc":
                answer = asksaveasfilename(filetypes=[("Image files", ".pnc")])
                if answer != "":
                    Escribe_PNC(self.o1, canal, answer)
                    messagebox.showinfo(message=answer + " creado!!!", title="Aviso!!!")
            else:
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1 and self.o1.get_modelo_ajuste() != -1:
                    answer = asksaveasfilename(filetypes=[("Image files", ".pnc")])
                    if answer != "":
                        Escribe_PNC(self.o1, canal, answer)
                        messagebox.showinfo(message=answer + " creado!!!", title="Aviso!!!")
                else:
                    messagebox.showinfo(message="Imagen no calibrada!!!", title="Aviso!!!")
        
    def perfiles(self, foco):
        global flag_perfil
        
        flag_perfil = True
        
    def mouse_move_1(self, event):
         global flag_velo_1
         global flag_dosis_conocida_1
         global rectangulo_velo_200
        
         if self.o1 != None:
             if flag_velo_1 == True or flag_dosis_conocida_1 == True:
                 if rectangulo_velo_200 != None:
                     self.canvas1_ppal.delete(rectangulo_velo_200)
                 dr = Dibuja_rectangulo_canvas(self.canvas1_ppal, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(), event.x, event.y, 25, True)
                 rectangulo_velo_200 = dr.rectangulo
             pos_imagen = Coordenadas_widget_to_imagen(event.x, event.y, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                       self.o1.get_width(), self.o1.get_height())
             if self.o1.get_tipo() == "pnc":
                 self.root.title(self.Title + ". Imagen 1. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                               ", " + str(pos_imagen.y) + "): Dosis. " + str(int(self.o1.dosis[0][pos_imagen.y][pos_imagen.x])) + " " +
                               str(int(self.o1.dosis[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o1.dosis[2][pos_imagen.y][pos_imagen.x])) +
                               ": Tamaño. (" + str(self.o1.get_width()) + ", " + str(self.o1.get_height()) + "): Factor Pixel/mm. " + str(round(self.o1.factor_pixel_mm, 3)))
             else:
                if (self.o1.get_cal(0) == -1 and self.o1.get_cal(1) == -1 and self.o1.get_cal(2) == -1) or self.o1.get_modelo_ajuste() == -1:
                     self.root.title(self.Title + ". Imagen 1. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                                   ", " + str(pos_imagen.y) + "): RGB. " + str(int(self.o1.pixel[0][pos_imagen.y][pos_imagen.x])) + " " +
                                   str(int(self.o1.pixel[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o1.pixel[2][pos_imagen.y][pos_imagen.x])) + 
                                   ": Tamaño. (" + str(self.o1.get_width()) + ", " + str(self.o1.get_height()) + "): Factor Pixel/mm. " + str(round(self.o1.factor_pixel_mm, 3)))
                else:
                    self.root.title(self.Title + ". Imagen 1. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                                   ", " + str(pos_imagen.y) + "): Dosis. " + str(int(self.o1.dosis[0][pos_imagen.y][pos_imagen.x])) + " " +
                                   str(int(self.o1.dosis[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o1.dosis[2][pos_imagen.y][pos_imagen.x])) + 
                                   ": Tamaño. (" + str(self.o1.get_width()) + ", " + str(self.o1.get_height()) + "): Factor Pixel/mm. " + str(round(self.o1.factor_pixel_mm, 3)) +
                                   ": Cal. " + self.o1.get_cal_id(0) + ", " + self.o1.get_cal_id(1) + ", " + self.o1.get_cal_id(2) + ": Modelo. " + modelos_ajuste_ids[self.o1.get_modelo_ajuste()])
            
    def mouse_move_2(self, event):
         global flag_velo_2
         global flag_dosis_conocida_2
         global rectangulo_velo_200
         global modelos_ajuste_ids
        
         if self.o2 != None:
             if flag_velo_2 == True or flag_dosis_conocida_2 == True:
                 if rectangulo_velo_200 != None:
                     self.canvas2_ppal.delete(rectangulo_velo_200)
                 dr = Dibuja_rectangulo_canvas(self.canvas2_ppal, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(), event.x, event.y, 25, True)
                 rectangulo_velo_200 = dr.rectangulo
             pos_imagen = Coordenadas_widget_to_imagen(event.x, event.y, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                       self.o2.get_width(), self.o2.get_height())
             if self.o2.get_tipo() == "pnc":
                 self.root.title(self.Title + ". Imagen 2. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                               ", " + str(pos_imagen.y) + "): Dosis. " + str(int(self.o2.dosis[0][pos_imagen.y][pos_imagen.x])) + " " +
                               str(int(self.o2.dosis[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o2.dosis[2][pos_imagen.y][pos_imagen.x])) +
                               ": Tamaño. (" + str(self.o2.get_width()) + ", " + str(self.o2.get_height()) + "): Factor Pixel/mm. " + str(round(self.o2.factor_pixel_mm, 3)))
             else:
                if (self.o2.get_cal(0) == -1 and self.o2.get_cal(1) == -1 and self.o2.get_cal(2) == -1) or self.o2.get_modelo_ajuste() == -1:
                     self.root.title(self.Title + ". Imagen 2. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                                   ", " + str(pos_imagen.y) + "): RGB. " + str(int(self.o2.pixel[0][pos_imagen.y][pos_imagen.x])) + " " +
                                   str(int(self.o2.pixel[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o2.pixel[2][pos_imagen.y][pos_imagen.x])) + 
                                   ": Tamaño. (" + str(self.o2.get_width()) + ", " + str(self.o2.get_height()) + "): Factor Pixel/mm. " + str(round(self.o2.factor_pixel_mm, 3)))
                else:
                    self.root.title(self.Title + ". Imagen 2. (" + str(event.x) + ", " + str(event.y) + "): (" + str(pos_imagen.x) +
                                   ", " + str(pos_imagen.y) + "): Dosis. " + str(int(self.o2.dosis[0][pos_imagen.y][pos_imagen.x])) + " " +
                                   str(int(self.o2.dosis[1][pos_imagen.y][pos_imagen.x])) + " " + str(int(self.o2.dosis[2][pos_imagen.y][pos_imagen.x])) + 
                                   ": Tamaño. (" + str(self.o2.get_width()) + ", " + str(self.o2.get_height()) + "): Factor Pixel/mm. " + str(round(self.o2.factor_pixel_mm, 3)) +
                                   ": Cal. " + self.o2.get_cal_id(0) + ", " + self.o2.get_cal_id(1) + ", " + self.o2.get_cal_id(2) + ": Modelo. " + modelos_ajuste_ids[self.o2.get_modelo_ajuste()])
        
    def mouse_clicked_1(self, event):
        global flag_velo_1
        global flag_dosis_conocida_1
        global rectangulo_velo_200
        global flag_perfil
        global flag_recorte
        global _start
        global _end
        
        if self.o1 != None:
            if flag_velo_1 == True or flag_dosis_conocida_1 == True:
                pos_imagen = Coordenadas_widget_to_imagen(event.x, event.y, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                       self.o1.get_width(), self.o1.get_height())
                dr = Dibuja_rectangulo_canvas(self.canvas1_ppal, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(), event.x, event.y, 25, False)
                x0_canvas = dr.punto_a_x
                y0_canvas = dr.punto_a_y
                x1_canvas = dr.punto_b_x
                y1_canvas = dr.punto_b_y
                pos_imagen_0 = Coordenadas_widget_to_imagen(x0_canvas, y0_canvas, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                       self.o1.get_width(), self.o1.get_height())
                pos_imagen_1 = Coordenadas_widget_to_imagen(x1_canvas, y1_canvas, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                            self.o1.get_width(), self.o1.get_height())
            
                media = Get_media_recorte(self.o1.pixel, pos_imagen_0.x, pos_imagen_0.y, pos_imagen_1.x, pos_imagen_1.y).media
            
                if flag_velo_1 == True:
                    self.o1.set_pixel_velo(0, int(media[0]))
                    self.o1.set_pixel_velo(1, int(media[1]))
                    self.o1.set_pixel_velo(2, int(media[2]))
            
                    self.canvas1_ppal.delete(rectangulo_velo_200)
                    flag_velo_1 = False
                    messagebox.showinfo(message="¡¡¡Velo establecido!!! (" + str(self.o1.get_pixel_velo(0)) + ", " +
                                        str(self.o1.get_pixel_velo(1)) + " " + str(self.o1.get_pixel_velo(2)) + ")", title="Aviso!!!")
                else:
                    self.o1.set_pixel_dosis_conocida(0, int(media[0]))
                    self.o1.set_pixel_dosis_conocida(1, int(media[1]))
                    self.o1.set_pixel_dosis_conocida(2, int(media[2]))
            
                    self.canvas1_ppal.delete(rectangulo_velo_200)
                    flag_dosis_conocida_1 = False
                    messagebox.showinfo(message="¡¡¡Pixel de dosis conocida establecido!!! (" + str(self.o1.get_pixel_dosis_conocida(0)) + ", " +
                                        str(self.o1.get_pixel_dosis_conocida(1)) + " " + str(self.o1.get_pixel_dosis_conocida(2)) + ")", title="Aviso!!!")
                
            if flag_perfil == True:
                x0, y0 = (self.canvas1_ppal.canvasx(event.x), self.canvas1_ppal.canvasy(event.y))
                if self.o1.get_tipo() == 'pnc':
                    pc = Coordenadas_widget_to_imagen(x0, y0, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                self.o1.get_width(), self.o1.get_height())
                    ph = Get_profile(self.o1.dosis[0], pc.y, 0, 0)
                    plt.clf()
                    plt.plot(ph.profile_x, ph.profile_y)
                    plt.title("Horizontal por y =" + str(pc.y) + " pixel")
                    plt.xlabel("pixel")
                    plt.ylabel("dosis")
                    plt.show()
                    
                    pv = Get_profile(self.o1.dosis[0], pc.x, 1, 0)
                    plt.clf()
                    plt.plot(pv.profile_x, pv.profile_y)
                    plt.title("Vertical por x =" + str(pc.x) + " pixel")
                    plt.xlabel("pixel")
                    plt.ylabel("dosis")
                    plt.show()
                else:
                    pc = Coordenadas_widget_to_imagen(x0, y0, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                                self.o1.get_width(), self.o1.get_height())
                    if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1 and self.o1.get_modelo_ajuste() != -1:
                        ph0 = Get_profile(self.o1.dosis[0], pc.y, 0, 0)
                        ph1 = Get_profile(self.o1.dosis[1], pc.y, 0, 0)
                        ph2 = Get_profile(self.o1.dosis[2], pc.y, 0, 0)
                        plt.clf()
                        plt.plot(ph0.profile_x, ph0.profile_y, color = "red")
                        plt.plot(ph1.profile_x, ph1.profile_y, color = "green")
                        plt.plot(ph2.profile_x, ph2.profile_y, color = "blue")
                        plt.title("Horizontal por y =" + str(pc.y) + " pixel")
                        plt.xlabel("pixel")
                        plt.ylabel("dosis")
                        plt.show()
                        
                        pv0 = Get_profile(self.o1.dosis[0], pc.x, 1, 0)
                        pv1 = Get_profile(self.o1.dosis[1], pc.x, 1, 0)
                        pv2 = Get_profile(self.o1.dosis[2], pc.x, 1, 0)
                        plt.clf()
                        plt.plot(pv0.profile_x, pv0.profile_y, color = "red")
                        plt.plot(pv1.profile_x, pv1.profile_y, color = "green")
                        plt.plot(pv2.profile_x, pv2.profile_y, color = "blue")
                        plt.title("Vertical por x =" + str(pc.x) + " pixel")
                        plt.xlabel("pixel")
                        plt.ylabel("dosis")
                        plt.show()
                    else:
                        messagebox.showinfo(message="Imagen no calibrada!!!", title="Aviso!!!")
                        
                flag_perfil = False
            
            if flag_recorte == True:
                _start = (self.canvas1_ppal.canvasx(event.x), self.canvas1_ppal.canvasy(event.y))
                _end = None
            
    def mouse_clicked_2(self, event):
        global flag_velo_2
        global flag_dosis_conocida_2
        global rectangulo_velo_200
        global flag_perfil
        global flag_recorte
        global _start
        global _end
        
        if self.o2 != None:
            if flag_velo_2 == True or flag_dosis_conocida_2 == True:
                pos_imagen = Coordenadas_widget_to_imagen(event.x, event.y, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                       self.o2.get_width(), self.o2.get_height())
                dr = Dibuja_rectangulo_canvas(self.canvas2_ppal, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(), event.x, event.y, 25, False)
                x0_canvas = dr.punto_a_x
                y0_canvas = dr.punto_a_y
                x1_canvas = dr.punto_b_x
                y1_canvas = dr.punto_b_y
                pos_imagen_0 = Coordenadas_widget_to_imagen(x0_canvas, y0_canvas, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                            self.o2.get_width(), self.o2.get_height())
                pos_imagen_1 = Coordenadas_widget_to_imagen(x1_canvas, y1_canvas, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                            self.o2.get_width(), self.o2.get_height())
            
                media = Get_media_recorte(self.o2.pixel, pos_imagen_0.x, pos_imagen_0.y, pos_imagen_1.x, pos_imagen_1.y).media
                
                if flag_velo_2 == True:
                    self.o2.set_pixel_velo(0, int(media[0]))
                    self.o2.set_pixel_velo(1, int(media[1]))
                    self.o2.set_pixel_velo(2, int(media[2]))
            
                    self.canvas2_ppal.delete(rectangulo_velo_200)
                    flag_velo_2 = False
                    messagebox.showinfo(message="¡¡¡Velo establecido!!! (" + str(self.o2.get_pixel_velo(0)) + ", " +
                                        str(self.o2.get_pixel_velo(1)) + " " + str(self.o2.get_pixel_velo(2)) + ")", title="Aviso!!!")
                else:
                    self.o2.set_pixel_dosis_conocida(0, int(media[0]))
                    self.o2.set_pixel_dosis_conocida(1, int(media[1]))
                    self.o2.set_pixel_dosis_conocida(2, int(media[2]))
            
                    self.canvas2_ppal.delete(rectangulo_velo_200)
                    flag_dosis_conocida_2 = False
                    messagebox.showinfo(message="¡¡¡Pixel de dosis conocida establecido!!! (" + str(self.o2.get_pixel_dosis_conocida(0)) + ", " +
                                        str(self.o2.get_pixel_dosis_conocida(1)) + " " + str(self.o2.get_pixel_dosis_conocida(2)) + ")", title="Aviso!!!")
    
            if flag_perfil == True:
                x0, y0 = (self.canvas2_ppal.canvasx(event.x), self.canvas2_ppal.canvasy(event.y))
                if self.o2.get_tipo() == 'pnc':
                    pc = Coordenadas_widget_to_imagen(x0, y0, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                self.o2.get_width(), self.o2.get_height())
                    ph = Get_profile(self.o2.dosis[0], pc.y, 0, 0)
                    plt.clf()
                    plt.plot(ph.profile_x, ph.profile_y)
                    plt.title("Horizontal por y =" + str(pc.y) + " pixel")
                    plt.xlabel("pixel")
                    plt.ylabel("dosis")
                    plt.show()
                    
                    pv = Get_profile(self.o2.dosis[0], pc.x, 1, 0)
                    plt.clf()
                    plt.plot(pv.profile_x, pv.profile_y)
                    plt.title("Vertical por x =" + str(pc.x) + " pixel")
                    plt.xlabel("pixel")
                    plt.ylabel("dosis")
                    plt.show()
                else:
                    pc = Coordenadas_widget_to_imagen(x0, y0, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                                self.o2.get_width(), self.o2.get_height())
                    if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1 and self.o2.get_modelo_ajuste() != -1:
                        ph0 = Get_profile(self.o2.dosis[0], pc.y, 0, 0)
                        ph1 = Get_profile(self.o2.dosis[1], pc.y, 0, 0)
                        ph2 = Get_profile(self.o2.dosis[2], pc.y, 0, 0)
                        plt.clf()
                        plt.plot(ph0.profile_x, ph0.profile_y, color = "red")
                        plt.plot(ph1.profile_x, ph1.profile_y, color = "green")
                        plt.plot(ph2.profile_x, ph2.profile_y, color = "blue")
                        plt.title("Horizontal por y =" + str(pc.y) + " pixel")
                        plt.xlabel("pixel")
                        plt.ylabel("dosis")
                        plt.show()
                        
                        pv0 = Get_profile(self.o2.dosis[0], pc.x, 1, 0)
                        pv1 = Get_profile(self.o2.dosis[1], pc.x, 1, 0)
                        pv2 = Get_profile(self.o2.dosis[2], pc.x, 1, 0)
                        plt.clf()
                        plt.plot(pv0.profile_x, pv0.profile_y, color = "red")
                        plt.plot(pv1.profile_x, pv1.profile_y, color = "green")
                        plt.plot(pv2.profile_x, pv2.profile_y, color = "blue")
                        plt.title("Vertical por x =" + str(pc.x) + " pixel")
                        plt.xlabel("pixel")
                        plt.ylabel("dosis")
                        plt.show()
                    else:
                        messagebox.showinfo(message="Imagen no calibrada!!!", title="Aviso!!!")
                        
                flag_perfil = False
                
            if flag_recorte == True:
                _start = (self.canvas2_ppal.canvasx(event.x), self.canvas2_ppal.canvasy(event.y))
                _end = None
    
    def _on_drag_1(self, event):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        if flag_recorte == True:
            x0, y0 = _start
            ex, ey = self.canvas1_ppal.canvasx(event.x), self.canvas1_ppal.canvasy(event.y)
            _end = (ex, ey)
            self._draw_rectangle_1()
            if flag_recorte_ambas == True:
                self._draw_rectangle_2()
            
    def _on_drag_2(self, event):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        if flag_recorte == True:
            x0, y0 = _start
            ex, ey = self.canvas2_ppal.canvasx(event.x), self.canvas2_ppal.canvasy(event.y)
            _end = (ex, ey)
            self._draw_rectangle_2()
            if flag_recorte_ambas == True:
                self._draw_rectangle_1()
            
    def _draw_rectangle_1(self):
        global _start
        global _end
        
        self.canvas1_ppal.delete("rectangle")
        
        if _end is None or _start is None:    
            return None

        x0, y0 = _start
        x1, y1 = _end

        self.canvas1_ppal.create_rectangle(x0, y0, x1, y1, fill="#18c194",
                                width=1, stipple="gray50", tags='rectangle'
                                )
        
    def _draw_rectangle_2(self):
        global _start
        global _end
        
        self.canvas2_ppal.delete("rectangle")
        
        if _end is None or _start is None:    
            return None

        x0, y0 = _start
        x1, y1 = _end

        self.canvas2_ppal.create_rectangle(x0, y0, x1, y1, fill="#18c194",
                                width=1, stipple="gray50", tags='rectangle'
                                )
    
    def _on_drop_1(self, event):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        if flag_recorte == True:
        
            # Hay que pasar de las coordenads del canvas a las de la imagen
            x0b, y0b = _start
            cs = Coordenadas_widget_to_imagen(x0b, y0b, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                              self.o1.get_width(), self.o1.get_height())
            x1b, y1b = _end
            ce = Coordenadas_widget_to_imagen(x1b, y1b, self.canvas1_ppal.winfo_width(), self.canvas1_ppal.winfo_height(),
                                              self.o1.get_width(), self.o1.get_height())
            if _end is None:
                pass
            
            else:
            
                # Acotar límites de selección a la imagen
                img_x = self.o1.get_width()
                img_y = self.o1.get_height()

                #x0, y0 = _start
                x0 = cs.x
                y0 = cs.y
                x0 = img_x if x0 > img_x else 0 if x0 < 0 else x0
                y0 = img_y if y0 > img_y else 0 if y0 < 0 else y0 
                _start = (x0, y0)

                #x1, y1 = _end
                x1 = ce.x
                y1 = ce.y
                x1 = img_x if x1 > img_x else 0 if x1 < 0 else x1
                y1 = img_y if y1 > img_y else 0 if y1 < 0 else y1       
                _end = (x1, y1)

                # Normalizado para obtener vertice superior izquierdo e inferior derecho
                if x0 > x1:
                    if y0 < y1: # _start es el vértice superior derecho
                        _start = (x1, y0)
                        _end = (x0, y1)
                    else:       # _start es el vértice inferior derecho
                        _start, _end = _end, _start
                else:
                    if y0 > y1:  # _start es el vértice inferior izquierdo
                        _start = (x0, y1)
                        _end = (x1, y0)
                    
            # Redibujar rectágulo
            self._draw_rectangle_1()
    
            self._crop_1()
            self.canvas1_ppal.delete("rectangle")
            if flag_recorte_ambas == True:
                self._draw_rectangle_2()
                self._crop_2()
                self.canvas2_ppal.delete("rectangle")
            
    def _on_drop_2(self, event):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        if flag_recorte == True:
        
            # Hay que pasar de las coordenads del canvas a las de la imagen
            x0b, y0b = _start
            cs = Coordenadas_widget_to_imagen(x0b, y0b, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                              self.o2.get_width(), self.o2.get_height())
            x1b, y1b = _end
            ce = Coordenadas_widget_to_imagen(x1b, y1b, self.canvas2_ppal.winfo_width(), self.canvas2_ppal.winfo_height(),
                                              self.o2.get_width(), self.o2.get_height())
            if _end is None:
                pass
            
            else:
            
                # Acotar límites de selección a la imagen
                img_x = self.o2.get_width()
                img_y = self.o2.get_height()

                #x0, y0 = _start
                x0 = cs.x
                y0 = cs.y
                x0 = img_x if x0 > img_x else 0 if x0 < 0 else x0
                y0 = img_y if y0 > img_y else 0 if y0 < 0 else y0 
                _start = (x0, y0)

                #x1, y1 = _end
                x1 = ce.x
                y1 = ce.y
                x1 = img_x if x1 > img_x else 0 if x1 < 0 else x1
                y1 = img_y if y1 > img_y else 0 if y1 < 0 else y1       
                _end = (x1, y1)

                # Normalizado para obtener vertice superior izquierdo e inferior derecho
                if x0 > x1:
                    if y0 < y1: # _start es el vértice superior derecho
                        _start = (x1, y0)
                        _end = (x0, y1)
                    else:       # _start es el vértice inferior derecho
                        _start, _end = _end, _start
                else:
                    if y0 > y1:  # _start es el vértice inferior izquierdo
                        _start = (x0, y1)
                        _end = (x1, y0)
                    
            # Redibujar rectágulo
            self._draw_rectangle_2()
    
            self._crop_2()
            self.canvas2_ppal.delete("rectangle")
            if flag_recorte_ambas == True:
                self._draw_rectangle_1()
                self._crop_1()
                self.canvas1_ppal.delete("rectangle")
            
    def _crop_1(self):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        flag_recorte = False
        X0, Y0 = _start
        X1, Y1 = _end
       
        if  self.o1 != None:
            if self.o1.get_tipo() == 'pnc':
                dst = Cut(self.o1.get_dosis(3), X0, Y0, X1, Y1).dst
                self.o1.set_dosis(dst)
                self.o1.set_width(len(dst[0][0]))
                self.o1.set_height(len(dst[0]))
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Cut(self.o1.get_pixel(3), X0, Y0, X1, Y1).dst
                self.o1.set_pixel(dst)
                self.o1.set_height(len(dst[0]))
                self.o1.set_width(len(dst[0][0]))
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust1.png', self.canvas1_ppal, 1, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o1.get_cal(0) != -1 and self.o1.get_cal(1) != -1 and self.o1.get_cal(2) != -1:
                    dst = Cut(self.o1.get_don(3), X0, Y0, X1, Y1).dst
                    self.o1.set_don(dst)
                    if self.o1.get_modelo_ajuste() != -1:
                        dst = Cut(self.o1.get_dosis(3), X0, Y0, X1, Y1).dst
                        self.o1.set_dosis(dst)
            if flag_recorte_ambas == True:
                flag_recorte_ambas = False
                self._crop_2()
    
    def _crop_2(self):
        global _start
        global _end
        global flag_recorte
        global flag_recorte_ambas
        
        flag_recorte = False
        X0, Y0 = _start
        X1, Y1 = _end
       
        if  self.o2 != None:
            if self.o2.get_tipo() == 'pnc':
                dst = Cut(self.o2.get_dosis(3), X0, Y0, X1, Y1).dst
                self.o2.set_dosis(dst)
                self.o2.set_width(len(dst[0][0]))
                self.o2.set_height(len(dst[0]))
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "pnc")
            else:
                dst = Cut(self.o2.get_pixel(3), X0, Y0, X1, Y1).dst
                self.o2.set_pixel(dst)
                self.o2.set_height(len(dst[0]))
                self.o2.set_width(len(dst[0][0]))
                dst_im = matriz_datos_2_matriz_im(dst).matriz_im
                self.Mostrar(dst_im, 'dust2.png', self.canvas2_ppal, 2, self.ancho_canvas_ppal, self.alto_canvas_ppal, "tif")
                if self.o2.get_cal(0) != -1 and self.o2.get_cal(1) != -1 and self.o2.get_cal(2) != -1:
                    dst = Cut(self.o2.get_don(3), X0, Y0, X1, Y1).dst
                    self.o2.set_don(dst)
                    if self.o2.get_modelo_ajuste() != -1:
                        dst = Cut(self.o2.get_dosis(3), X0, Y0, X1, Y1).dst
                        self.o2.set_dosis(dst)
            if flag_recorte_ambas == True:
                flag_recorte_ambas = False
                self._crop_1()
            
    # Funcion para representar en el canvas directamente una matriz3D formato imagen, con el canal en la
    # tercera dimension. 
    def Mostrar2(self, matriz_imagen, canvas, foco, ancho_canvas, alto_canvas, tipo):
        max = np.max(matriz_imagen)
        #if max > 0:
        matriz_imagen = matriz_imagen*255/max
        matriz_imagen = matriz_imagen.astype(np.uint8)
        img = Image.fromarray(matriz_imagen)
        img = img.resize((ancho_canvas, alto_canvas), Image.ANTIALIAS)
        if foco == 1:
            self.img1 = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor="nw", image=self.img1)
        if foco == 2:
            self.img2 = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor="nw", image=self.img2)


    def Mostrar(self, matriz3D_im, fichero_imagen, canvas, foco, ancho_canvas, alto_canvas, tipo):
        if foco ==1:
            maximo = np.amax(matriz3D_im)
            if maximo > 255 and tipo != "pnc":
                cv2.imwrite(fichero_imagen, matriz3D_im.astype(np.uint16))
            else:
                cv2.imwrite(fichero_imagen, matriz3D_im)
            img = Image.open(fichero_imagen)
            img = img.resize((ancho_canvas, alto_canvas), Image.ANTIALIAS)
            self.img1 = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img1)
        if foco == 2:
            maximo = np.amax(matriz3D_im)
            if maximo > 255 and tipo != "pnc":
                cv2.imwrite(fichero_imagen, matriz3D_im.astype(np.uint16))
            else:
                cv2.imwrite(fichero_imagen, matriz3D_im)
            img = Image.open(fichero_imagen)
            img = img.resize((ancho_canvas, alto_canvas), Image.ANTIALIAS)
            self.img2 = ImageTk.PhotoImage(img)
            canvas_img = canvas.create_image(0, 0, anchor = "nw", image=self.img2)
    
if __name__ == "__main__":
    main_frame = MainFrame()