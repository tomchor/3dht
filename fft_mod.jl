
FFTW.set_num_threads(Nprocs)

#-----
# Set forward fft
U0 = Array{Float64}(Ndim, Nx, Ny, Nz)
rplan = plan_rfft(U0, (2,3,4), flags=FFTW.MEASURE);
rplan1 = plan_rfft(U0[1,:,:,:], (1,2,3), flags=FFTW.MEASURE);
#-----

#-----
# Set inverse FFT
Uh0 = rplan*U0
irplan = plan_irfft(Uh0, Nx, (2,3,4), flags=FFTW.MEASURE);
irplan1 = plan_irfft(Uh0[1,:,:,:], Nx, (1,2,3), flags=FFTW.MEASURE);
#-----

#-----
# Create plan for padded arrays
Uh_p = Array{Complex{Float64}}(Int(ceil(size(Uh0)[2]*3/2)), Int(ceil(size(Uh0)[3]*3/2)), Int(ceil(size(Uh0)[4]*3/2)))
ipadplan = plan_irfft(Uh_p, Int(Nx*3//2+2), (1,2,3), flags=FFTW.MEASURE);
U_p = ipadplan*Uh_p
padplan = plan_rfft(U_p, (1,2,3), flags=FFTW.MEASURE);
#-----


function adv(a::Array, b::Array, apad::Array, bpad::Array, a_phys::Array, b_phys::Array, phys::Array)
    """
    a and b should be in Fourier space
    """
    #----
    # Create padded arrays with 3/2 times the original size and all zeros
    apad[:,:,:] = 0+0*im
    bpad[:,:,:] = 0+0*im
    #----

    #----
    # Fill in the amplitudes that we care about
    halfy = Int(Ny//2)-1
    half2 = round(Int, Ny*3//4)
    a = fftshift(a, (2,3))
    b = fftshift(b, (2,3))
    apad[1:Int(Nx//2)+1, half2-halfy:half2+1+halfy, half2-halfy:half2+1+halfy] = a[:,:,:]
    bpad[1:Int(Nx//2)+1, half2-halfy:half2+1+halfy, half2-halfy:half2+1+halfy] = b[:,:,:]
    apad = fftshift(apad, (2,3))
    bpad = fftshift(bpad, (2,3))
    #----

    #----
    # Transform to physical space, renormalize, multiply, and transform back
    A_mul_B!(a_phys, ipadplan, apad)
    A_mul_B!(b_phys, ipadplan, bpad)
    a_phys .*= length(a_phys)/(Nx*Ny*Nz)
    b_phys .*= length(b_phys)/(Nx*Ny*Nz)

    phys = a_phys.*b_phys
    A_mul_B!(apad, padplan, phys)
    g=fftshift(apad, (2,3))[1:Int(Nx//2)+1, half2-halfy:half2+1+halfy, half2-halfy:half2+1+halfy]
    g=fftshift(g, (2,3))
    #----

    #----
    # Return only the wavenumbers that interest us
    return g
    #----
end


function diff_x(arr::Array)
    return irplan*(im*Kx*arr)
end
function diff_y(arr::Array)
    return irplan*(im*Ky*arr)
end


function adv_aliased(a::Array, b::Array)
    """
    a and b should be in Fourier space
    """

    #----
    # Transform to physical space, renormalize, multiply, and transform back
    a_phys = irplan1*a
    b_phys = irplan1*b

    g = rplan1*(a_phys.*b_phys)
    #----

    #----
    # Return only the wavenumbers that interest us
    return g
    #----
end


