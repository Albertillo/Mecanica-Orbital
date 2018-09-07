###############################################################################################################################################################
# Kepler2.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                       #
# La función funciona de la siguiente manera: Se debe especificar un foco de la órbita y sus elementos orbitales (semieje mayor "a" (km), excentricidad "e",  #
# inclinación "i" (rad), longitud del perihelio "omegar" (rad), longitud del nodo ascendente "Omega" (rad), anomalía verdadera "theta" (rad).                 #
# Con esto calculará los vectores de posición y velocidad del cuerpo en coordenadas heliocéntricas. Los focos disponibles en este momento son "Sol", "Tierra" #
# y "Jupiter", pero pueden añadirse facilmente introduciendo las diferentes "mu" en el primer if.											                  #
# Algoritmo de Howard D. Curtis, Orbital Mechanics for Engineering Students, First Edition, Elsevier Butterworth-Heinemann, Oxford (UK), 2005.                #
# Función necesaria: angulomenor2pi.py                                                                                                                        #
###############################################################################################################################################################

import numpy as np
from math import *

def Kepler2(focpoint,a,e,i,omegar,Omega,theta): 
#Semieje mayor "a" en km, inclinación "i" en rad, longitud del perihelio "omegar" en rad, longitud del nodo ascendente "Omega" en rad, anomalía verdadera "theta" en rad.
#focpoint puede tener los valores "Sol", "Tierra" o "Jupiter"
#La función devolverá 6 valores: las tres primeras serán las componentes de las coordenadas heliocéntricas del vector posición, y las tres últimas las componentes del vector velocidad.

    if focpoint=="Sol":
        mu=132712439935.5 #km^3/s^2
    elif focpoint=="Tierra":
        mu=398600.4 #km^3/s^2
    elif focpoint=="Jupiter":
        mu=126711995.4 #km^3/s^2
    else:
        print("ERROR, FOCO DE LA ÓRBITA NO VÁLIDO.")
    #Momento angular específico.
    h=sqrt(mu*a*(1-pow(e,2)))
	#Argumento del perihelio.
    omega=omegar-Omega
	#En caso de que omega>2pi
    omega=angulomenor2pi(omega)
    #Cálculo de la posición
    rp=np.array([(h**2/mu)*(1/(1 + e*cos(theta)))*cos(theta),(h**2/mu)*(1/(1 + e*cos(theta)))*sin(theta),0])
    
    #Cálculo de la velocidad
    vp=np.array([(mu/h)*(-sin(theta)),(mu/h)*(e+cos(theta)),0]);
    
    #Matriz de rotación de perifocal a heliocéntrico
    rot1=np.array([[cos(Omega),sin(Omega),0],[-sin(Omega),cos(Omega),0],[0,0,1]])
    rot2=np.array([[1,0,0],[0,cos(i),sin(i)],[0,-sin(i),cos(i)]])
    rot3=np.array([[cos(omega),sin(omega),0],[-sin(omega),cos(omega),0],[0,0,1]])
    matrizrot=np.dot(np.dot(rot3,rot2),rot1)
    matrizrot=np.transpose(matrizrot)
	
    #Vectores de velocidad y posición en coordenadas heliocéntricas
    vhelio=np.dot(matrizrot,vp)
    V_X=vhelio[0]
    V_Y=vhelio[1]
    V_Z=vhelio[2]
    rhelio=np.dot(matrizrot,rp)
    X=rhelio[0]
    Y=rhelio[1]
    Z=rhelio[2]
    
    return X,Y,Z,V_X,V_Y,V_Z
