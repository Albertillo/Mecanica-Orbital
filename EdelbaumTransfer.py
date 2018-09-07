#########################################################################################################################################################################
# EdelbaumTransfer.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                        #
# La función funciona de la siguiente manera: Se debe especificar una aceleración constante, un foco de la órbita mediante un mu, un semieje mayor inicial y final y    #
# una inclinación inicial y final de órbita. A partir de esto, se calculará mediante la teoría de Edelbaum el DeltaV total, el tiempo de transferencia y diferentes     #
# parámetros relevantes como el beta o el semieje mayor con respecto al tiempo. También se obtendrá el vector posición en cada momento. 								#
# Algoritmo obtenido de "Vladimir A. Chobotov, Orbital Mechanics, Third Edition, AIAA Education Series, Virginia (USA), 2002."
#########################################################################################################################################################################

from math import *

def EdelbaumTransfer(f,mu,a0,afi,i0,ifi): #f es la aceleración en km/s^2, mu la constante del foco en km^3/s^2, a0 y afi los semiejes mayores iniciales y finales (km)
#e i0 e ifi las inclinaciones iniciales y finales (rad), respectivamente.
	
	#variación total de inclinación y velocidades orbitales de inicio y fin.
	deltaif=ifi-i0
	v0=(mu/a0)**0.5
	vfi=(mu/afi)**0.5
	
	#beta0
	beta0=atan2(sin(pi*deltaif/2),v0/vfi-cos(pi*deltaif/2))
	#DeltaV total.
	deltaV=v0*cos(beta0)-v0*sin(beta0)/tan(0.5*deltaif*pi+beta0)
	#Tiempo de vuelo.
	tf=deltaV2/f
	
	DeltaV=[]
	beta=[]
	V=[]
	lambdav=[]
	deltai=[]
	tiempo=[]

	t=0
	contador=0
	rx=[]
	ry=[]
	rz=[]
	theta=[]
	dt=tf/10000 #Paso de tiempo. Si es demasiado alto para la misión a estudiar se recomienda disminuir este valor.
	a=[]

	while t<tf:
		
		#Se calculan los diferentes parámetros de interés en la transferencia de Edelbaum y se adjuntan a un array.
		tiempo.append(t/86400)
		DeltaV.append(f*t)
		beta.append(atan2(v0*sin(beta0),v0*cos(beta0)-f*t))
		ft=f*cos(beta[contador])
		V.append(sqrt(v0**2-2*v0*f*t*cos(beta0)+(f*t)**2))
		lambdav.append(cos(beta[contador])/f)
		deltai.append(2/pi*(atan2(f*t-v0*cos(beta0),v0*sin(beta0))+pi/2-beta0))
		a.append(mu/V[contador]**2)
		t=t+dt
		#P=periodo
		P=2*pi*sqrt(a[contador]**3/mu)
		if contador==0:
			theta.append(0)
		else:
			#Se aumenta el theta en función del periodo debido a que se supone un cambio en la órbita muy bajo.
			theta.append(angulomenor360(2*pi*dt/P+theta[contador-1]))
		#Se calcula el vector de posición en cada momento t.
		if ifi>i0:
			vecposvel=Kepler2("Sol",a[contador],0,i0+deltai[contador],0,0,theta[contador])
		else:
			vecposvel=Kepler2("Sol",a[contador],0,i0-deltai[contador],0,0,theta[contador])
		rx.append(vecposvel[0])
		ry.append(vecposvel[1])
		rz.append(vecposvel[2])
		contador=contador+1
	#Se imprime el tiempo de transferencia.
	if tf<60:
		print("Tiempo total de transferencia=", tf,"s")
	elif tf/60<60:
		print("Tiempo total de transferencia=", tf/60,"min")
	elif tf/3600<24:
		print("Tiempo total de transferencia=", tf/3600,"horas")
	elif tf/86400<365.25:
		print("Tiempo total de transferencia=", tf/86400,"días")
	else:
		print("Tiempo total de transferencia=", tf/86400/365.25,"años")
	#Se imprime el DeltaV total, beta0 y beta final.
	print("DeltaV=", deltaV,"km/s")
	print("beta0=",beta0*180/pi,"º")
	print("betaf=",beta[len(beta)-1]*180/pi,"º")
	
return tf, deltaV, beta0, (beta[len(beta)-1]*180/pi),rx,ry,rz #Devuelve datos de interés. Modificable.