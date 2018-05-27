

function sines(f1, f2, x_center, y_center)
    S1 = [ sin.(2*pi*x*f1) for x in x_center, y in y_center ]
    S2 = [ sin.(2*pi*x*f2) for x in x_center, y in y_center ]
    S12 = [ (1/2)*cos.(2*pi*(f2-f1)*x) for x in x_center, y in y_center ]
    return S1, S2, S12
end


function test_dealias(f1, f2)
    S1, S2, S12 = sines(f1, f2, x_center, y_center)

    deal = adv(rfft(S1), rfft(S2))

    fig, axes = plt.subplots(ncols=3, figsize=(12,4))
    a=axes[1][:pcolormesh](irfft(deal, Nx, (1,2)), vmin=-1,vmax=1); plt.colorbar(a, ax=axes[1])
    a=axes[2][:pcolormesh](S1.*S2, vmin=-1,vmax=1); plt.colorbar(a, ax=axes[2])
    a=axes[3][:pcolormesh](S12, vmin=-1,vmax=1); plt.colorbar(a)

end


function profile_test(n)
    for i = 1:n
        A = randn(100,100,20)
        m = maximum(A)
        Afft = fft(A)
        Am = mapslices(sum, A, 2)
        B = A[:,:,5]
        Bsort = mapslices(sort, B, 1)
        b = rand(100)
        C = B.*b
    end
end

profile_test(1)  # run once to trigger compilation
Profile.clear()  # in case we have any previous profiling data
@profile profile_test(10)

