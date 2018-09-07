#####################################################################################################################################################################
# BusquedaFlybyMultiple.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                               #
# La función funciona de la siguiente manera: Dados una serie de tramos entre planetas (Ej: Tierra-Venus, Venus-Júpiter), unas restricciones de tiempo y DeltaV, así# 
# como un rango de fechas de salida a calcular, el programa calculará la mejor fecha posible para la misión de asistencia gravitatoria múltiple a realizar en 	    #
# cuestión del DeltaV, suponiendo asistencias gravitatorias con impulso en el periapsis. El método es el utilizado en el libro que se menciona al final de este	    #
# párrafo introductorio. Para realizar este método se ha escrito una función recursiva que calcule tramo a tramo los DeltaV, de forma que compare todas y cada una  #
# de las posibilidades entre las fechas escogidas para un paso de tiempo deltat.										                                            #
# Algoritmo obtenido de:																	                                                                        #
# A.V. Labunsky, O.V. Papkov, K.G. Sukhanov, Multiple Gravity Assist Interplanetary Trajectories, First Edition, ESI Book Series, Amsterdan (The Netherlands), 1998.#
# Funciones necesarias: Clasif.py, Lambert.py, Kepler3.py, Julianday.py, ModBas.py, Gregoriandate.py								    							#
#####################################################################################################################################################################

import numpy as np
from math import *

#Cálculo del flyby múltiple.
tiempo_inicio=time()

#Se debe especificar los tramos como un array según la codificación explicada en Clasif.py del libro Multiple Gravity Assist Interplanetary Trajectories. En este 
#caso el ejemplo es Tierra-Venus-Júpiter.
tramos=[[10003,10002],[10002,10003],[10003,10005]]
N=len(tramos) #Número de segmentos a calcular (trayectorias de un planeta a otro).
lamb="prograde"

#Rango de días de salida a estudiar.
dia_inic=Julianday(1,1,2020)
dia_fin=Julianday(31,12,2035)
#Step utilizado de fechas.
Deltat=2 #días

#Restricción de tiempo de vuelo máximo por segmento en días. El tiempo de vuelo se variará desde 100 días hasta la restricción de tiempo de vuelo de ese tramo.
tiempomax=[700,700,1500]
tofmax=np.sum(tiempomax)
#Restricción de DeltaV por tramo. Se debe especificar en cada segmento ambos planetas. En este caso (DeltaVmaxi=[[4.1,0.4],[0.4,0.4],[0.4,7.7]]) se tiene que el 
#DeltaV de salida de la Tierra debe ser máximo 4.1 km/s, 0.4 km/s en el flyby en Venus y de 7.7 km/s en la puesta en órbita en Júpiter.
DeltaVmaxi=[[4.1,0.4],[0.4,0.4],[0.4,7.7]]
#Vector en el que se van añadiendo los diferentes Delta-V de cada una de las trayectorias. Tiene todos los vectores DeltaV calculados.
ImpulsoTotal=[]
#Vector en el que se van añadiendo los diferentes tiempos de vuelo de cada una de las trayectorias. Tiene todos los vectores tof calculados.
TOF=[]
#Vector en el que se van añadiendo los días de salida de cada una de las trayectorias válidas.
vector_dia_salida_total=[]
terminado=0
#Vector de Delta-V para tener todos los Delta-V de una trayectoria completa.
deltav=np.zeros(N+1)
#Vector de tof para tener todos los tof de una trayectoria completa.
tof=np.zeros(N)
#Vector de días de salida para una trayectoria válida.
vector_dia_salida=np.zeros(N)
#Dias de salida/llegada a cada planeta.
tiempos=np.zeros(N+1)

solucion=[0,[0,0,0],0,0,0,0]
velinf=[0,0,0]
solucionado=False
contador_soluciones=0

def buscarSolucion(indice = 0): #La variable indice es la que indica el tramo en el que se encuentra calculando el programa. 0 es el primero, 1 el segundo... etc.
    
    global solucion, velinf, solucionado, contador_soluciones,deltav,tof,vector_dia_salida,ImpulsoTotal,vector_dia_salida_total,TOF
    solucionado = False
	
	#A partir de aquí empezarán las iteraciones de fechas de salida para el tramo en cuestión. Como la fecha de salida solo puede variar en el primer planeta, esto
	#solo se realiza para el primer tramo.
    for j in range(ceil((dia_fin-dia_inic)/Deltat)):
        if indice==0:
            tiempos[0]=dia_inic+Deltat*j
            print("t",indice+1,"=", Gregoriandate(tiempos[0]))
        elif indice>0 and j>0:
            return False
            
		#A continuación iterará con diferentes fechas de llegada, calculando utilizando la función de ModBas.py el DeltaV dependiendo de si es flyby o Hohmann/no-Hohmann.
        for k in range(ceil((tiempomax[indice])/Deltat)):

            if indice==0:
                tipo1="directo"
                tipo2="flyby"
                Vinf1x=0
                Vinf1y=0
                Vinf1z=0
				#Para el caso de salida desde la Tierra, se comienza con un tiempo de vuelo de 100 días. Esto se puede variar.
                tiempos[indice+1]=dia_inic+Deltat*j+(Deltat*k)+100
            elif indice==N-1:
                tipo2="directo"
                tipo1="flyby"
                Vinf1x=velinf[0]
                Vinf1y=velinf[1]
                Vinf1z=velinf[2]
				#Para el caso del tramo de llegada a Júpiter, se comienza con 700 días de tiempo de vuelo. Esto se puede variar.
                tiempos[indice+1]=tiempos[indice]+(Deltat*k)+700
            else:
                tipo1="flyby"
                tipo2="flyby"
                Vinf1x=velinf[0]
                Vinf1y=velinf[1]
                Vinf1z=velinf[2]
				#Para el caso de la llegada a Venus, se comienza con 50 días de tiempo de vuelo. Esto se puede variar.
                tiempos[indice+1]=tiempos[indice]+Deltat*k+50
            
			#En el caso de que se llegue a la fecha de salida límite y se supere la restricción de tiempo de vuelo se terminará el programa.
            if tiempos[0]>=dia_fin-Deltat and ((tiempos[indice+1]-tiempos[0])>=tofmax or tiempos[indice+1]-tiempos[indice]>=tiempomax[indice]):
                return False
			#En el caso de que se supere el tiempo de vuelo máximo se pasará a la fecha de salida.
            elif (tiempos[indice+1]-tiempos[0])>tofmax or tiempos[indice+1]-tiempos[indice]>=tiempomax[indice]:
                break
			#En caso contrario, calculará el DeltaV. Se han agregado unas excepciones que puede llegar a dar el problema de Lambert por su naturaleza para que las ignore.
            else:
                try:
                    solucion=ModBas(tramos[indice][0],tramos[indice][1],tiempos[indice],tiempos[indice+1],tipo1,tipo2,lamb,Vinf1x,Vinf1y,Vinf1z)
                except ValueError:
                    continue
                except OverflowError:
                    continue
                except ZeroDivisionError:
                    continue
				
				#Altura mínima de flyby y máxima (Esfera de influencia).
                altura1min=Clasif(tramos[indice][0],tiempos[indice])[0]+200
                altura1max=Clasif(tramos[indice][0],tiempos[indice])[2]
				
				#Se comprueban las restricciones.
                if solucion[4]<=DeltaVmaxi[indice][0] and tiempos[indice+1]-tiempos[indice]<=tiempomax[indice]:
                    validacaso=False
					#Solución válida para tramo intermedio.
                    if indice!=0 and solucion[2]>=altura1min and solucion[2]<=altura1max and indice!=len(tramos)-1:
                        deltav[indice]=solucion[4]
                        velinf=solucion[1]
                        validacaso=True
					#Solución válida para el primer tramo.
                    elif indice==0:
                        validacaso=True
                        deltav[0]=solucion[4]
                        velinf=solucion[1]

					#Solución válida para el último tramo.
                    elif indice==N-1 and solucion[2]>=altura1min and solucion[2]<=altura1max:
                        validacaso=True
                        ###########################################################################################################################################
                        #Esto no está totalmente automatizado, así que se debe cambiar manualmente el append en TOF y vector_dia_salida_total en caso de que se   #
						#aumente o disminuya el número de tramos para que contemple todos los tramos.															  #
                        ###########################################################################################################################################
                        deltav[N-1]=solucion[4]
                        deltav[N]=solucion[5]
                        ImpulsoTotal.append(np.sum(deltav))
                        TOF.append([tiempos[1]-tiempos[0],tiempos[2]-tiempos[1],tiempos[3]-tiempos[2]])
                        vector_dia_salida_total.append([tiempos[0],tiempos[1],tiempos[2],tiempos[3]])
                        contador_soluciones=contador_soluciones+1
                    if validacaso:
                        if indice<len(tramos)-1:
							#Si encuentra solución y no es el último tramo pasará al siguiente.
                            solucionado=buscarSolucion(indice+1)
                        else:
							#Si encuentra solución para el último tramo imprimirá que ha encontrado una solución.
                            print("Solución encontrada.")
                            solucionado=True
                        
    return solucionado

#Se llama a la función para el índice 0 (no es necesario especificarlo en la función).
buscarSolucion ()
#Una vez terminada la búsqueda de soluciones, buscará la más óptima en términos de DeltaV.
if contador_soluciones>0:
    Mejor=-1
    Mejorindice=-1
    if len(ImpulsoTotal)>0:
        for i in range(0,ceil((len(ImpulsoTotal)-1)/2)):
            inegativo=(i+1)*(-1)
            if ImpulsoTotal[i]<Mejor or Mejor==-1:
                Mejor=ImpulsoTotal[i]
                Mejorindice=i
            if ImpulsoTotal[inegativo]<Mejor or Mejor==-1:
                Mejor=ImpulsoTotal[inegativo]
                Mejorindice=inegativo
        print("Mejor DeltaV total=",Mejor)
        print("Mejor TOF=", TOF[Mejorindice])
        print("días=", vector_dia_salida_total[Mejorindice])
        print("Terminado.")
        print("Mejorindice=",Mejorindice)
        print("Número de soluciones encontradas:",contador_soluciones)
else:
    print("No hay solución.") 
    
tiempo_fin=time()
tiempo_ejec = (tiempo_fin - tiempo_inicio)/3600 #en horas
print("Tiempo de ejecución:", tiempo_ejec, "horas.")
