#------
# General
const Ndim = 3
Nprocs = nprocs()
#------

#------
# Domain setup and and IC
L_dim=0.05 # m
N=2^7
ν_dim=1.5e-5 # m^2/s
#------

#------
L_scale=L_dim # m
U_scale=2.0 # m/s
#------

#------
L=L_dim/L_scale
ν=ν_dim/(L_scale*U_scale)
#------

#------
# Time setup
dt = .001
dt_out = 0.05
T_f = 5.
out_T = 0:dt_out:T_f
Nt = trunc(Int,T_f/dt)
out_n = trunc.(Int,out_T./dt)
rms=0.01
#------

IC="jet"

