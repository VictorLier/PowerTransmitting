import math
import sympy as sp


# Given parameters
power = 5500 #[W]

rpm = 2900 #[rpm]

M = 1.2 #Mass of belt [kg/m^2]

T = 0.0013 #Thickness [m]

Strength = 170000 #Strength of belt [N/m] - Width of belt

Elong = 3900 #Newton required for 1% elongation [N/m]

W = 0.042 #Width of belt [m]

mu = 0.33

D1 = 0.224 #[m]
D2 = 0.355 #[m]

a = 0.4 #[m]

#A - Rotation speed of the output

OutputRotationSpeed = D1/D2 * rpm

def RPMtoRadians(RPM: float):
    '''
    Output the rpm as [rad/s]
    '''
    radians = RPM * 2*math.pi/60
    return radians


print("The rotational speed of the putpu is", OutputRotationSpeed, "rpm")


#B - The belt speed
v = math.pi*D1 * rpm/60

print("The Belt speed is", v, "[m/s]")


#C - Centrifugal force

q = M*W #belt weight pr. meter [kg/m]

CentrifugalForce = q*v**2

print("The centrifugal force is", CentrifugalForce, "[N]")


#D - Total belt forces


omega1 = RPMtoRadians(rpm)
omega2 = RPMtoRadians(OutputRotationSpeed)

def WrapAngle(D1: float, D2: float, Distance: float):
    '''
    Calculates the wrap angle with the diameter and the distance between two pullies.
    Outputs Alpha1, Alpha2 [rad]
    '''
    a = math.asin((D2-D1)/(2*Distance))
    alpha1 = math.pi-2*a
    alpha2 = math.pi+2*a
    return alpha1, alpha2

#Finding F1


alpha1, alpha2 =WrapAngle(D1,D2,a)


T1 = power/omega1

T2 = T1*(D2/D1)

F2 = (T1/(D1/2))/(math.exp(mu*alpha1)-1)
F1 = F2 * math.exp(mu*alpha1)

F1 = F1 + CentrifugalForce
F2 = F2 + CentrifugalForce


print("The total belt forces F1 =",F1,"[N] and F2 =",F2, "[N]")

#E - Pretension

Fp = (F1+F2)/2

print("The preload is",Fp,"[N]")


#F - Belt forces at 20 Nm output

M2 = 20 #[Nm]

M1 = M2*D1/D2

F1n = M2/D2 + Fp

F2n = Fp - M2/D2

print("The belt forces at M2 = 20 Nm is: F1 =",F1n,"[N] andF2 =", F2n, "[N]")


#G - Unloaded length of the belt

def BeltLength(Dl: float, Ds: float, L: float):
    '''
    Returns the belt length from the large and small diameter [Dl and Ds] and length between centers
    '''

    Length = ((Dl+Ds)*math.pi/2)+(Dl-Ds)*math.asin((Dl-Ds)/(2*L))+2*math.sqrt((L**2-0.25*(Dl-Ds)**2))
    return Length


BeltL = BeltLength(D2,D1,a)

LoadElong = Elong * W
Elongation = (Fp/LoadElong)/100

StartBeltLenght = BeltL/(Elongation+1)

print("The starting belt length should be ",StartBeltLenght,"[m]")


#H