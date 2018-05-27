
function zonal_jet(x::Float64; wid::Float64=.25, center::Float64=.5, a::Float64=2*π, b::Int=10)
    x = (x-center)
    u = sqrt.((1+b^2)/(1+b^2*cos.(a*x)^2)).*cos.(a*x)
    return u
end


function gaussian(x, y; x0=.5, y0=.5, a=.25, H=1)
    x=x-x0
    y=y-y0
    arg = (x^2 + y^2)/a^2
    return H*exp(-arg)
end



function check_CFL(uh::Array; dx=1, dy=1, dt=1, norm=[1,1], reset=false)
    if reset==true
        rm("ke.csv")
    end
    u, v = abs.(irplan*uh)
    CFLx = dt*maximum(u)/dx
    CFLy = dt*maximum(v)/dy
    kea = sum(abs.(uh.^2))*dx*dy
    println("KE, CFL = ", kea/norm[1], ", ", (CFLx+CFLy)/norm[2])
    open("ke.csv", "a") do f
        write(f, kea/norm[1])
    end
    return kea, CFLx + CFLy
end


function rm_div(uh::Array)
    NL_aux[1,:,:] = (1 - kxkx_k2).*uh[1,:,:] - kxky_k2.*uh[2,:,:]
    NL_aux[2,:,:] = - kxky_k2.*uh[1,:,:] + (1 - kyky_k2).*uh[2,:,:]
    return NL_aux
end


function vort(uh::Array)
    return irfft(im*(Kx[1,:,:].*uh[2,:,:] - Ky[1,:,:].*uh[1,:,:]), Nx)
end

function ∇(uh::Array)
    return abs.(im*(Kx[1,:,:].*uh[1,:,:] + Ky[1,:,:].*uh[2,:,:]))
end

