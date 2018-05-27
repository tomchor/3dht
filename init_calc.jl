#------
# Calculate derivatives
dUh0dx = np.stack([Kx.*im.*Uh0, Ky.*im.*Uh0], axis=1)
#------

#------
# Calculate nonliner term
NL = -get_nonlinear(Uh0, dUh0dx)
#------

#------
E = exp.(ν*K.^2*dt)
AB3_coeffs = np.stack([E.^(-3)*5/12, -E.^(-2)*4/3, E.^(-1)*23/12], axis=0)
#AB3_coeffs = reshape(AB3_coeffs, (3,1,1,1))
#------

#------
Uh = np.stack([ Uh0 for i in 1:3 ], axis=0)
NL = np.stack([ NL for i in 1:3 ], axis=0)
#------

Q = sum(AB3_coeffs.*NL, 1)[1,:,:,:]

#----
# Shift the positions of previous time steps on last moment
Uh = np.roll(Uh, -1, axis=0)
NL = np.roll(NL, -1, axis=0)
#----

#----
# Calculate current time based on last
Uh[end,:,:,:] = E.^(-1).*Uh[end-1,:,:,:] + dt*Q
#----

