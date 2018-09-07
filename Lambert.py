###############################################################################################################################################################
# Lambert.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                       #
# La función funciona de la siguiente manera: Se debe especificar un foco, dos vectores de posición, un tiempo de vuelo y el tipo de órbita (prograde o 	  #
# retrograde). Con esto resolverá el problema de Lambert y devolverá los vectorse de velocidad en los puntos de los vectores de posición. Los focos 		  #
# disponibles en este momento son "Sol", "Tierra" y "Jupiter", pero pueden añadirse facilmente introduciendo las diferentes "mu" en el primer if.			  #
# Algoritmo de Howard D. Curtis, Orbital Mechanics for Engineering Students, First Edition, Elsevier Butterworth-Heinemann, Oxford (UK), 2005.                #
# Función necesaria: StumpffS(z).py y StumpffC(z).py                                                                                                          #
###############################################################################################################################################################

import numpy as np
from math import *

#Problema de Lambert según Orbital Mechanics for Engineering Students, Algoritmo 5.2.
def Lambert(focpoint,r1x,r1y,r1z,r2x,r2y,r2z,tof,tipo): #r en km. tiempo de vuelo "tof" en días. tipo="retrograde" o "prograde". focpoint debe ser o "Sol" o "Tierra" o "Jupiter".

    if focpoint=="Sol":
        mu=132712439935.5 #km^3/s^2
    elif focpoint=="Tierra":
        mu=398600.4 #km^3/s^2
    elif focpoint=="Jupiter":
        mu=126711995.4 #km^3/s^2
    else:
        print("ERROR, FOCO DE LA ÓRBITA NO VÁLIDO.")
		
    #Cálculo del módulo de los vectores y generación de vectores.
    r1=np.array([r1x,r1y,r1z])
    r2=np.array([r2x,r2y,r2z])
    modr1=np.linalg.norm(r1)
    modr2=np.linalg.norm(r2)
    
    #Cálculo del Deltatheta para trayectorias prograde y retrograde:
    if tipo=="prograde":
        if np.cross(r1,r2)[2]>=0:
            vartheta=acos(np.dot(r1,r2)/(modr1*modr2))
        else:
            vartheta=2*pi-acos(np.dot(r1,r2)/(modr1*modr2))
    elif tipo=="retrograde":
        if np.cross(r1,r2)[2]>=0:
            vartheta=2*pi-acos(np.dot(r1,r2)/(modr1*modr2))
        else:
            vartheta=acos(np.dot(r1,r2)/(modr1*modr2))
    else:
        print("Error, tipo de órbita no válida. Escribe prograde o retrograde.")
    #Cálculo A:
    A=sin(vartheta)*sqrt(modr1*modr2/(1-cos(vartheta)))
    
	#Se obtiene z por iteración. Si z<0 la órbita es una hipérbola, si z=0 parábola, si z>0 elipse. 
    #se utiliza ahora ratio=10 para que entre en el bucle
    ratio=10
    #z inicial
    z=0
    iteraciones=0
	
    while abs(ratio)>1e-8:
        y=modr1+modr2+A*(z*S(z)-1)/sqrt(C(z))
        y0=modr1+modr2-A/sqrt(C(0))
        iteraciones=iteraciones+1
        
		#Ecuación 5.40
        Fun=sqrt(y*y*y/(C(z)*C(z)*C(z)))*S(z)+A*sqrt(y)-sqrt(mu)*tof*86400
        
        #Ecuación 5.43
        if z==0:
            Funp=sqrt(2)*sqrt(y0*y0*y0)/40+A*(sqrt(y0)+A*sqrt(1/(2*y0)))/8
        else:
            Funp=sqrt(y*y*y/(C(z)*C(z)*C(z)))*((C(z)-3*S(z)/(2*C(z)))/(2*z)+3*S(z)*S(z)/(4*C(z)))+A*(3*S(z)*sqrt(y)/C(z)+A*sqrt(C(z)/y))/8
        ratio=Fun/Funp
        
        if abs(ratio)>1e-8:
            #Ecuación 5.45
            z=z-ratio    
        if iteraciones>1000:
            SaltarIteracion=10000
            break
    #Ecuación 5.38
    if iteraciones<=1000:
        SaltarIteracion=-1
    if z==0:
        y=modr1+modr2-A/sqrt(C(0))
    else:
        y=modr1+modr2+A*(z*S(z)-1)/sqrt(C(z))
    #Cálculo de f,g y gdot con ecuación 5.46
    f=1-y/modr1
    g=A*sqrt(y/mu)
    gp=1-y/modr2
    #Cálculo de v1 y v2 con las ecuaciones 5.28 y 5.29
    v1=(r2-f*r1)/g
    v2=(gp*r2-r1)/g
    return v1,v2,SaltarIteracion
