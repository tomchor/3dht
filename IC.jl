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
else
    println("IC not coded")
end

U0[1,:,:,:] = u0
U0[2,:,:,:] = v0
U0[3,:,:,:] = w0
U0*=U_dim/U_scale

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
