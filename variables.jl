#-----
# Start definitions and allocations
include("aux_mod.jl")
#-----

#-----
# Define wavenumbers
kx = 2*π*np.fft[:rfftfreq](Nx, d=dx)
ky = 2*π*np.fft[:fftfreq](Ny, d=dx)
kz = 2*π*np.fft[:fftfreq](Nz, d=dx)
dkx=diff(kx)[1]
#-----

#-----
Kx = [ ki for ki in kx, kj in ky, kk in kz ]
Ky = [ kj for ki in kx, kj in ky, kk in kz ]
Kz = [ kk for ki in kx, kj in ky, kk in kz ]

Kx = np.stack([ Kx for i in 1:Ndim ], axis=0)
Ky = np.stack([ Ky for i in 1:Ndim ], axis=0)
Kz = np.stack([ Kz for i in 1:Ndim ], axis=0)

K = sqrt.(Kx.^2 + Ky.^2 + Kz.^2)
#-----

#-----
kxkx_k2 = (Kx.*Kx./K.^2)[1,:,:,:]
kyky_k2 = (Ky.*Ky./K.^2)[1,:,:,:]
kzkz_k2 = (Kz.*Kz./K.^2)[1,:,:,:]
kxky_k2 = (Kx.*Ky./K.^2)[1,:,:,:]
kxkz_k2 = (Kx.*Kz./K.^2)[1,:,:,:]
kykz_k2 = (Ky.*Kz./K.^2)[1,:,:,:]

kxkx_k2[1,1] = 0
kyky_k2[1,1] = 0
kzkz_k2[1,1] = 0
kxky_k2[1,1] = 0
kxkz_k2[1,1] = 0
kykz_k2[1,1] = 0
#-----

#-----
U = zeros(Ndim, Nx, Ny, Nz)
dUhdx = Array{Complex{Float64}}(3, 3, Int(Nx//2+1), Ny, Nz)
coords=("x", "y", "z")
DS = xr.Dataset(Dict("U"=>(coords, U[1,:,:,:]), "V"=>(coords, U[2,:,:,:]), "W"=>(coords, U[3,:,:,:])),
    coords=Dict("x" => x_center, "y" => y_center, "z" => z_center))
#-----


#-----
# Aux variable to enforce non-compressibility in NL term
NL_aux = Array{Complex{Float64}}(Ndim, length(kx), length(ky), length(kz))
NL_aux2 = Array{Complex{Float64}}(Ndim, length(kx), length(ky), length(kz))
#-----


#----
# Calculate production term
f_k = np.stack([ (1- heaviside.(K-k_peak)) for i in 1:3 ], axis=0)
#f_k*= Prod/(2*k_peak)
#-----


#-----
# Define auxiliary padded arrays for de-alising
apad = Array{Complex{Float64}}(Int(ceil(size(kx)[1]*3/2)), Int(ceil(size(ky)[1]*3/2)), Int(ceil(size(kz)[1]*3/2)))
bpad = similar(apad)
phys = Array{Float64}(Int(Nx*3//2+2), Int(Nx*3//2), Int(Nx*3//2))
aphys = similar(phys)
bphys = similar(phys)
#-----
