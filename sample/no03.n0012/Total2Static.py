#!/usr/bin/python
import math
Rgas=287.
gamma=1.4

Pt=25.0081
Tt=232.6485
Mach=.4033

r_ts=1.0+(gamma-1.0)/2*Mach**2
Ts = Tt/r_ts
Pt*=6894.757
Ps = Pt/r_ts**(gamma/(gamma-1.0))

print "total  pressure(bar)=", Pt/1e5
print "static temperature(K)=",Ts
print "static pressure(bar)=", Ps/1e5

a = math.sqrt(gamma*Rgas*Ts)
Uin = a*Mach

print "inflow velocity(m/s)=",Uin
