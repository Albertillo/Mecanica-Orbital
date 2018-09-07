#####################################################################################################################################################################
# Kepler.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                              #
# La función funciona de la siguiente manera: Se debe especificar un código del planeta para obtener los elementos orbitales del mismo dado un día juliano "t"      #
# Algoritmo obtenido de:																	    #
# A.V. Labunsky, O.V. Papkov, K.G. Sukhanov, Multiple Gravity Assist Interplanetary Trajectories, First Edition, ESI Book Series, Amsterdan (The Netherlands), 1998.#
# Función necesaria: angulomenor2pi.py                                                                                                                        	    #
#####################################################################################################################################################################

from math import *

def Clasif(N,t): #Con esto se obtienen los elementos orbitales de los planetas según la clasificación del libro Multiple Gravity Assist Interplanteary Trajectories
    #El código N del planeta es el explicado en el libro Multiple Gravity Assist Interplanteary Trajectories. Consiste de un número de 5 cifras. La primera muestra el 
	#tipo de objeto: 1 planeta, 2 asteroide, 3 cometa, 4 satélite natural y 5 satélite artificial. Si hubiera un subtipo sería el segundo número. Los restantes denotan el cuerpo en sí. Por ejemplo, los planetas se numeran
	#de más cercano a más lejano al Sol, siendo 1 Mercurio, 2 Venus, 3 Tierra... etc. Por tanto, la Tierra sería el 10003, Neptuno el 10008, la Luna 43001 (4= satélite artificial, 3=del planeta Tierra, 001=el primero y único en este caso)... etc.
	#Los que se han agregado a la función son: Venus=10002, Tierra=10003, Marte=10004 y Júpiter=10005.
	
	AU2km=149597871 #factor de conversión AU->km.
	s2deg=1/3600 #segundos a grados
    G=6.674e-11/pow(1000,3) #km^3/kgs^2
	
	#Masa del sol y siglos que han pasado desde J2000 para el día especificado.
    ms=1.9891e30 #kg
    T0=(t-Julianday(1,1,2000))/36525
	
    if N//10000==1:
        if (N-10000)//1==2:
            #print("Venus")
            radio=6051.8 #km
            mu=324858.8 #km^3/s^2
            rsoi=0.616e6 #km
            a=(0.72333199+0.00000092*T0)*AU2km #km
            e=0.00677323-0.00004938*T0
            i=(3.39471-2.86*s2deg*T0)*pi/180 #rad
            i=angulomenor360(i)
            omegar=(131.53298-108.8*s2deg*T0)*pi/180#rad
            omegar=angulomenor360(omegar)
            Omega=(76.68069-996.89*s2deg*T0)*pi/180 #rad
            Omega=angulomenor360(Omega)
            L=(181.97973+210664136.06*s2deg*T0)*pi/180 #rad
            L=angulomenor360(L)

        elif (N-10000)//1==3:
            #print("Tierra")
            radio=6378.14 #km
            mu=398600.4 #km^3/s^2
            rsoi=0.924e6 #km
            a=(1.00000011-0.00000005*T0)*AU2km #km
            e=0.01671022-0.00003804*T0
            i=(0.00005-46.94*s2deg*T0)*pi/180 #rad
            i=angulomenor360(i)
            omegar=(102.94719+1198.28*s2deg*T0)*pi/180#rad
            omegar=angulomenor360(omegar)
            Omega=(-11.26064-18228.25*s2deg*T0)*pi/180 #rad
            Omega=angulomenor360(Omega)
            L=(100.46435+129597740.63*s2deg*T0)*pi/180 #rad
            L=angulomenor360(L)

        elif (N-10000)//1==4:
            #print("Marte")
            radio=3397.0 #km
            mu=42828.3 #km^3/s^2
            rsoi=0.577e6 #km
            a=(1.52366231-0.00007221*T0)*AU2km #km
            e=0.09341233+0.00011902*T0
            i=(1.85061-25.47*s2deg*T0)*pi/180 #rad
            i=angulomenor360(i)
            omegar=(336.04084+1560.78*s2deg*T0)*pi/180 #rad
            omegar=angulomenor360(omegar)
            Omega=(49.57854-1020.19*s2deg*T0)*pi/180 #rad
            Omega=angulomenor360(Omega)
            L=(355.45332+68905103.78*s2deg*T0)*pi/180
            L=angulomenor360(L)

        elif (N-10000)//1==5:
            #print("Júpiter")
            radio=71492 #km
            mu=126711995.4 #km^3/s^2
            rsoi=48.157e6 #km
            a=(5.20336301+0.00060737*T0)*AU2km #km
            e=0.04839266-0.00012880*T0
            i=(1.30530-4.15*s2deg*T0)*pi/180 #rad
            i=angulomenor360(i)
            omegar=(14.75385+839.93*s2deg*T0)*pi/180 #rad
            omegar=angulomenor360(omegar)
            Omega=(100.55615+1217.17*s2deg*T0)*pi/180 #rad
            Omega=angulomenor360(Omega)
            L=(34.40438+10925078.35*s2deg*T0)*pi/180
            L=angulomenor360(L)

        else:
            print("ERROR: código del cuerpo celeste no válido.")
    elif N//10000==4:
        if (N-45000)//1==2:
            print("Europa")
        else:
            print("ERROR: código del cuerpo celeste no válido.")
    else:
        print("ERROR: código del cuerpo celeste no válido.")
    return radio,mu,rsoi,a,e,i,omegar,Omega,L
