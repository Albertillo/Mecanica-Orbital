#####################################################################################################################################################################
# ModBas.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                        		#
# La función funciona de la siguiente manera: Dados los códigos de dos planetas, las fechas de salida del planeta 1 y llegada al planeta 2, si las órbitas son		#
# prograde o retrograde, si se realiza sobre ese planeta en cuestión una asistencia gravitatoria o una transferencia tipo Hohmann/no-Hohmann y la velocidad   		#
# hiperbólica de exceso en la llegada al planeta 1 en caso de que se realice un flyby, el módulo calculará automáticamente el DeltaV necesario para la misión.		#
# Las asistencias gravitatorias las calculará suponiendo un impulso en el periapsis, de forma que el valor de DeltaV que dará el programa será el impulso     		#
# necesario en el periapsis de la órbita hiperbólica de flyby. Además, en el caso del flyby, sólo se calculará el DeltaV si se realiza en el primer planeta y no en #
# el segundo. Esto es debido a que es necesario conocer la velocidad hiperbólica de exceso en la salida, y si el planeta de flyby es el segundo, este dato no se    #
# conoce.																																							#
# Algoritmo obtenido de:																																			#
# A.V. Labunsky, O.V. Papkov, K.G. Sukhanov, Multiple Gravity Assist Interplanetary Trajectories, First Edition, ESI Book Series, Amsterdan (The Netherlands), 1998.#
# Funciones necesarias: Clasif.py, Lambert.py, Kepler3.py   														#												#
#####################################################################################################################################################################

import numpy as np
from math import *
from scipy import optimize

#Módulo básico para flyby múltiple para un segmento (es decir, para el cálculo desde un planeta hasta el siguiente). Libro Multiple Gravity Assist Interplanetary Trajectories.
Vinf1=[0,0,0]
Vinf2=[0,0,0]
a1=0
a2=0
A=0
def ModBas(N1,N2,t1,t2,tipo1,tipo2,lamb,Vinf1x,Vinf1y,Vinf1z): #N1 y N2 cuerpos de inicio y destino según la clasificación del libro (se explica en Clasif.py), t1 y t2 
#fechas de salida y llegada en días julianos, Vinf la Vinf de llegada al planeta 1, si es 0 es que es el planeta de inicio. Tipo="flyby" o "directo". Lamb="prograde" 
#o "retrograde"
    
    global a1, a2, Vinf1, Vinf2,A
    #Tiempo de vuelo
	vart=t2-t1
	
    mu=132712439935.5 #km^3/s^2
    
    cuerpo1=[]
    cuerpo2=[]
    
    #Identificación de cuerpos de salida y llegada (radio del cuerpo, mu y rsoi)
    cuerpo1=Clasif(N1,t1)
    cuerpo2=Clasif(N2,t2)
    R01=cuerpo1[0]
    mu1=cuerpo1[1]
    rsoi1=cuerpo1[2]
    a1=cuerpo1[3]
    e1=cuerpo1[4]
    i1=cuerpo1[5]
    omegar1=cuerpo1[6]
    Omega1=cuerpo1[7]
    L1=cuerpo1[8]
    
    R02=cuerpo2[0]
    mu2=cuerpo2[1]
    rsoi2=cuerpo2[2]
    a2=cuerpo2[3]
    e2=cuerpo2[4]
    i2=cuerpo2[5]
    omegar2=cuerpo2[6]
    Omega2=cuerpo2[7]
    L2=cuerpo2[8]
    
    cuerpo1=[]
    cuerpo2=[]
	
    #Cálculo de las posiciones y velocidades de los cuerpos 1 y 2
    cuerpo1=Kepler3("Sol",a1,e1,i1,omegar1,Omega1,L1)
    cuerpo2=Kepler3("Sol",a2,e2,i2,omegar2,Omega2,L2)

    r1=np.array([cuerpo1[0],cuerpo1[1],cuerpo1[2]])
    vp1=np.array([cuerpo1[3],cuerpo1[4],cuerpo1[5]])
    r2=np.array([cuerpo2[0],cuerpo2[1],cuerpo2[2]])
    vp2=np.array([cuerpo2[3],cuerpo2[4],cuerpo2[5]])

    #Órbita de transferencia:
    transf=Lambert("Sol",r1[0],r1[1],r1[2],r2[0],r2[1],r2[2],vart,lamb)
    v1=transf[0]
    v2=transf[1]
	
	#Si se alcanza el máximo de iteraciones dará esta solución debido a que diverge:
    if transf[2]==10000:
        DeltaV=transf[2]
        vinf_approach2=[10000,10000,10000]
        rp1=0
        rp2=0
        DeltaV1=10000
        DeltaV2=10000
		
	#En caso contrario, se calculará el DeltaV.
    else:
	
        #Cálculo de Delta V en ambos cuerpos

        #Punto 1
        if tipo1=="directo": #"directo" se refiere a no-Hohmann. "flyby" a una asistencia gravitatoria en ese planeta.

            #radio de órbita en el cuerpo 1. Se toma 200 km sobre la superficie del planeta para los planetas interiores y 800000 km para Júpiter.
            if (N1-10000)//1==3 or (N1-10000)//1==2 or (N1-10000)//1==4:
                altura1=R01+200
            elif (N1-10000)//1==5:
                altura1=800000
			#Velocidad hiperbólica de escape
            vinf1=v1-vp1
            modvinf1=np.linalg.norm(vinf1)
			#semieje mayor y excentricidad de la órbita de escape.
            a_esc=-mu1/modvinf1**2
            e_esc=1-altura1/a_esc
			#Velocidad en periapsis de la órbita de escape.
            v_peri_esc=sqrt(mu1/altura1)*sqrt(2-altura1/a_esc)
			#Velocidad de aparcamiento.
            v_aparc=sqrt(mu1/altura1)
			#DeltaV de salida.
            DeltaV1=v_peri_esc-v_aparc

            rp1=0

        elif tipo1=="flyby":
            #Según el libro Multiple Gravity Assist Interplanetary Trajectories.

            #radio de órbita en el cuerpo 1 
            if (N1-10000)//1==3 or (N1-10000)//1==2 or (N1-10000)//1==4:
                altura1=R01+200
            elif (N1-10000)//1==5:
                altura1=800000
				
            vinf1=0
			#Velocidad hiperbólica de exceso al entrar en la órbita hiperbólica alrededor de planeta 1.
            Vinf1=np.array([Vinf1x,Vinf1y,Vinf1z])
			#Velocidad hiperbólica de salida.
            Vinf2=v1-vp1
            #Cálculo de las excentricidades de la órbita hiperbólica e1 y e2 (antes y después del impulso en periapsis).
            A=np.dot(Vinf1,Vinf2)/(np.linalg.norm(Vinf1)*np.linalg.norm(Vinf2))
            a1=-mu1/np.linalg.norm(Vinf1)**2
            a2=-mu1/np.linalg.norm(Vinf2)**2
            
			#Resolución de sistema de ecuaciones para obtener las excentricidades de la órbita de llegada y salida de flyby.
            def f(variables) :
                (e1,e2) = variables

                first_eq = a1*(1-e1)-a2*(1-e2)
                second_eq = A+cos(acos(-1/e1)+acos(-1/e1))
                return [first_eq, second_eq]

            valorese = optimize.fsolve(f, (1.1,1.4) )
            e1=valorese[0]
            e2=valorese[1]
            rp1=a1*(1-e1) #rp2=rp1 en cada flyby.
            
			#DeltaV necesario para el flyby (impulso en periapsis).
            DeltaV1=abs(sqrt(2*mu1/rp1+np.linalg.norm(Vinf2)**2)-sqrt(2*mu1/rp1+np.linalg.norm(Vinf1)**2))
        else:
            print("ERROR: valor tipo1 incorrecto.")

        #Punto 2

        if tipo2=="directo":
			
			#Velocidad hiperbólica de exceso en la llegada al planeta.
            vinf2=vp2-v2

            #Altura de órbita en 2. Mismo criterio que en el caso anterior.
            if (N2-10000)//1==3 or (N2-10000)//1==2 or (N2-10000)//1==4:
                altura2=R02+200
            elif (N2-10000)//1==5:
                altura2=800000
			
			#Módulo de la velocidad hiperbólica de exceso.
            modvinf2=np.linalg.norm(vinf2)
			
			#Semieje mayor y excentricidad de la órbita hiperbólica de llegada.
            a_esc2=-mu2/modvinf2**2
            e_esc2=1-altura2/a_esc2
			#Velocidad en periapsis de la órbita de llegada.
            v_peri_esc2=sqrt(mu2/altura2)*sqrt(2-altura2/a_esc2)
			#Velocidad orbital a la altura deseada.
            v_aparc2=sqrt(mu2/altura2)
			#DeltaV de llegada.
            DeltaV2=v_peri_esc2-v_aparc2

            rp2=0
            vinf_approach2=[0,0,0]
        elif tipo2=="flyby":
            #Se calcula la velocidad hiperbólica de exceso en la llegada para su uso en el siguiente tramo.

            vinf_approach2=vp2-v2
            rp2=0
        else:
            print("ERROR: valor tipo2 incorrecto.")
        if tipo2=="flyby":
            DeltaV2=0
		#DeltaV total.
        DeltaV=DeltaV1+DeltaV2
		
    return DeltaV,vinf_approach2,rp1,rp2,DeltaV1,DeltaV2,vinf1,r1,v1
