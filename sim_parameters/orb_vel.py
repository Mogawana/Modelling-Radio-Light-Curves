import astropy.units as unt
import astropy.constants as const
from numpy import sqrt, pi



M = 25*unt.Msun
P = 4000*unt.yr
e = 0.5

a = (P**2 * M)**(1/3)
r = a.value*unt.AU

R10 = (15*r)/M.value
R15 = (10*r)/M.value

V10 = (2*pi*R10)/P
V15 = (2*pi*R15)/P

#print(const.G)
#print(1*unt.AU, '=', 1*unt.AU.to(unt.m))
#print(M.to(unt.kg))

v = sqrt( (const.G * M.to(unt.kg)) / r.to(unt.m))

print("for circular orbit:")

print("\nthe total separation distance/semi-major axis is a = {}".format(r))
print("distance from COM to the B1 companion is r10 = {}".format(R10))
print("distance from COM to the red supergiant is r15 = {}".format(R15))
print("velocity of the B1 companion is v10 = {}".format(V10.to(unt.km/unt.s)))
print("velocity of the red supergiant is v15 = {}".format(V15.to(unt.km/unt.s)))
print("relative speed of the two bodies is v = {}".format(v.to(unt.km/unt.s)))

print(" \n\n\nfor elliptic obits:")

M = 25*unt.Msun
P = 4000*unt.yr
e = 0.5

a = (P**2 * M)**(1/3)
r = a.value*unt.AU

R10 = r*(1+e)
R15 = r*(1-e)
R = R10 + R15

V10 = sqrt((const.G * M.to(unt.kg))*((2/R10.to(unt.m)) - (1/r.to(unt.m))))
V15 = sqrt((const.G * M.to(unt.kg))*((2/R15.to(unt.m)) - (1/r.to(unt.m))))


print("\nthe semi-major axis is a = {}".format(r))
print("the total separation distance is r = {}".format(R))
print("distance from COM to the B1 companion is r10 = {}".format(R10))
print("distance from COM to the red supergiant is r15 = {}".format(R15))
print("velocity of the B1 companion is v10 = {}".format(V10.to(unt.km/unt.s)))
print("velocity of the red supergiant is v15 = {}".format(V15.to(unt.km/unt.s)))

