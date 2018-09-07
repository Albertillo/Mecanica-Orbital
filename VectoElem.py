###############################################################################################################################################################
# VectoElem.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                     #
# La función funciona de la siguiente manera: Se debe especificar un foco de la órbita y los vectores de posición y velocidad. La función calculará los       #
# elementos orbitales: Semieje mayor "a" (km), inclinación "i" (rad), longitud del nodo ascendente "Omega" (rad), excentricidad "mode", argumento del         #
# periastro "omega" (rad) y anomalía verdadera "theta" (rad). Los focos disponibles en este momento son "Sol", "Tierra" y "Jupiter", pero pueden añadirse     #
# facilmente introduciendo las diferentes "mu" en el primer if.											  													  #
# Algoritmo de Howard D. Curtis, Orbital Mechanics for Engineering Students, First Edition, Elsevier Butterworth-Heinemann, Oxford (UK), 2005.                #
###############################################################################################################################################################

import numpy as np
from math import *

#Elementos orbitales desde vectores de posición y velocidad.
def VectoElem(focpoint,rx,ry,rz,vx,vy,vz): #Posición en km y velocidad en km/s. focpoint solo admite "Sol", "Tierra" o "Jupiter".
    if focpoint=="Sol":
        mu=132712439935.5 #km^3/s^2
    elif focpoint=="Tierra":
        mu=398600.4 #km^3/s^2
    elif focpoint=="Jupiter":
        mu=126711995.4 #km^3/s^2
    else:
        print("ERROR, FOCO DE LA ÓRBITA NO VÁLIDO.")
    r=np.array([rx,ry,rz])
    v=np.array([vx,vy,vz])
    modr=np.linalg.norm(r)
    modv=np.linalg.norm(v)
    a=mu/(2*mu/modr-modv*modv)
    h=np.cross(r,v)
    i=np.arccos(h[2]/np.linalg.norm(h))
    N=np.cross([0,0,1],h)
    
    if N[1]>=0:
        Omega=np.arccos(N[0]/np.linalg.norm(N))
    else:
        Omega=2*pi-np.arccos(N[0]/np.linalg.norm(N))
    vr=np.dot(r,v)/modr
    e=1/mu*((modv*modv-mu/modr)*r-modr*vr*v)
    mode=np.linalg.norm(e)
    if mode>=0:
        omega=np.arccos(np.dot(N,e)/(np.linalg.norm(N)*mode))
    else:
        omega=2*pi-np.arccos(np.dot(N,e)/(np.linalg.norm(N)*mode))
    if vr>=0:
        theta=np.arccos(np.dot(e,r)/(modr*mode))
    else:
        theta=2*pi-np.arccos(np.dot(e,r)/(modr*mode))
    return a,i,Omega,mode,omega,theta
