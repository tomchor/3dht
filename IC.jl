include("aux_mod.jl")

U0 = zeros(Ndim, Nx, Ny, Nz)
if IC=="jet"

    u0 = [ zonal_jet(z; wid=.25, a=2*π) for x in x_center, y in y_center, z in z_center ]
    v0 = [ 0.0 for x in x_center, y in y_center, z in z_center ]
    w0 = [ 0.0 for x in x_center, y in y_center, z in z_center ]

elseif IC=="blob"
    u0 = [ gaussian(x, y, z) for x in x_center, y in y_center, z in z_center ]
    v0 = [ gaussian(x, y, z) for x in x_center, y in y_center, z in z_center ]
    w0 = [ gaussian(x, y, z) for x in x_center, y in y_center, z in z_center ]
elseif IC=="sine"
    u0 = [ sin(2*π*x/Lx) for x in x_center, y in y_center ]
    v0 = [ 0.0 for x in x_center, y in y_center ]
elseif IC=="iso"
    a = (2*π)/(0.8*L)
    arg = -2*K.^2/a^2
    Uh0[:,:,:,:] = (U_dim/(2*π)).*(K.^4).*exp.(arg)/U_scale
    U0 = irfft(Uh0, Nx, (2,3,4))
else
    println("IC not coded")
end

if IC=="jet" || IC=="blob" || IC=="sine"
    U0[1,:,:,:] = u0
    U0[2,:,:,:] = v0
    U0[3,:,:,:] = w0
end
U0*=U_dim/((maximum(U0) - minimum(U0))*U_scale)

#----
# Add noise
U0+=randn(size(U0))*rms
DS["U"][:values]=U0[1,:,:,:]; DS["V"][:values]=U0[2,:,:,:]; DS["W"][:values]=U0[3,:,:,:];
#----

#------
# Get 
Uh0 = rplan*U0
#------

Uh0 = rm_div(Uh0)
println("MAX(∇⋅U) = ", maximum(∇(Uh0)))
