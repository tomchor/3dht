#-----
# Start definitions and allocations
#-----

#-----
# Define wavenumbers
kx = 2*π*np.fft[:rfftfreq](Nx, d=dx)
ky = 2*π*np.fft[:fftfreq](Ny, d=dx)
#-----

if Ndim==2
    Kx = [ ki for ki in kx, kj in ky ]
    Ky = [ kj for ki in kx, kj in ky ]

    Kx = np.stack([ Kx for i in 1:Ndim ], axis=0)
    Ky = np.stack([ Ky for i in 1:Ndim ], axis=0)

    K = sqrt.(Kx.^2 + Ky.^2)
else
    println("Not coded up")
end

#-----
kxky_k2 = (Kx.*Ky./K.^2)[1,:,:]
kxkx_k2 = (Kx.*Kx./K.^2)[1,:,:]
kyky_k2 = (Ky.*Ky./K.^2)[1,:,:]

kxky_k2[1,1] = 0
kxkx_k2[1,1] = 0
kyky_k2[1,1] = 0
#-----

#-----
# Aux variable to enforce non-compressibility in NL term
NL_aux = Array{Complex{Float64}}(2, length(kx), length(ky))
#-----


#-----
# Define auxiliary padded arrays for de-alising
apad = Array{Complex{Float64}}(Int(ceil(size(kx)[1]*3/2)), Int(ceil(size(ky)[1]*3/2)))
bpad = Array{Complex{Float64}}(Int(ceil(size(kx)[1]*3/2)), Int(ceil(size(ky)[1]*3/2)))
phys = Array{Float64}(Int(Nx*3//2+2), Int(Nx*3//2))
aphys = Array{Float64}(Int(Nx*3//2+2), Int(Nx*3//2))
bphys = Array{Float64}(Int(Nx*3//2+2), Int(Nx*3//2))
#-----
