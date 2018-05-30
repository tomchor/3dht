#------
# General
const Ndim = 3
Nprocs = nprocs()
#------

#------
# Domain setup and and IC
U_dim=.5 # m/s
L_dim=0.05 # m
N=2^5
ν_dim=1.5e-5 # m^2/s
Prod_dim = 1. # m^2/s^3
#------

#------
L_scale=L_dim # m
U_scale=0.5 # m/s
#------

#------
L=L_dim/L_scale
ν=ν_dim/(L_scale*U_scale)
Prod=Prod_dim/Prod_dim
k_peak = (2*π)/(0.8*L) # From Pope
#------

#------
# Time setup
dt = .001
dt_out = 0.05
T_f = 10.
out_T = 0:dt_out:T_f
Nt = trunc(Int,T_f/dt)
out_n = trunc.(Int,out_T./dt)
rms=0.001
#------

IC="iso"

