# module ch_fft
using FFTW
using Images

include("fft.jl")
include("type.jl")
include("array_operation.jl")
include("thermocalc.jl")
include("ch_fft_channel.jl")

# @time begin
  

# nx = 64
# ny = 64
# dx = 1.0
# dy = 1.0
# nstep = 2
# nprint = 100
# dtime = 0.01
# coefA = 1.0
# c0 = 0.6
# mobility = 1.0
# grad_coef = 0.5
# rand_size = 0.02
# seed = 123
# η = 2.0

# ttime = 0.0

# c = thermocalc.make_initial_concentration(c0, nx, ny, seed, rand_size)
# c = array_operation.reject_element(c, (0.001, 0.999))
# (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, η)

# for istep in 1:nstep
#   @time begin
#   ttime = ttime + dtime
#   dfdc = thermocalc.get_free_energy(c)
#   dfdck = fft(dfdc)

#   ck = fft(c)

#   number = dtime * mobility * k2.* dfdck
#   denom = 1.0 .+ dtime * coefA * mobility * grad_coef .* k2 .* k2_anyso
#   ck = (ck - number) ./ denom
#   c = real(ifft(ck))
#   c = array_operation.reject_element(c, (0.001, 0.999))

#   if (istep % nprint == 0) || istep == 1
#     println("step: $istep done")
#     img = Gray.(c')
#     save("result/png/res_$istep.png", img)
#   end
# end
# end


# # # end
