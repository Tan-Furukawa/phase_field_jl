# module ch_fft
using FFTW
using Images
include("fft.jl")
include("array_operation.jl")
include("thermocalc.jl")

Base.@kwdef struct InitialParameter
  nx:: Integer = 64
  ny:: Integer = 64
  dx:: AbstractFloat = 1.0
  dy:: AbstractFloat = 1.0
  nstep:: Integer = 20000
  nprint:: Integer = 100
  dtime:: AbstractFloat = 0.01
  coefA:: AbstractFloat = 1.0
  c0:: AbstractFloat = 0.4
  mobility:: AbstractFloat = 1.0
  grad_coef:: AbstractFloat = 0.5
  rand_size:: AbstractFloat = 0.02
  seed:: Integer = 123
end

struct TmpResult
  c:: AbstractArray
  ttime:: AbstractFloat
end

struct Iterator
  state:: Tuple{Integer, TmpResult}
end



function ch_simple_binary_alloy(iparam:: InitialParameter)
  nx = iparam.nx
  ny = iparam.ny
  nxny = nx * ny
  dx = iparam.dx
  dy = iparam.dy
  nstep = iparam.nstep
  nprint = iparam.nprint
  dtime = iparam.dtime
  coefA = iparam.coefA
  c0 = iparam.c0
  mobility = iparam.mobility
  grad_coef = iparam.grad_coef
  rand_size = iparam.rand_size
  seed = iparam.seed

  ttime = 0.0

  c = thermocalc.make_initial_concentration(c0, nx, ny, seed, rand_size)
  c = array_operation.reject_element(c, (0.001, 0.999))
  (kx, ky, k2, k4) = Myfft.prepare_fft(nx, ny, dx, dy)

  function Base.iterate(iter::Iterator, state = iter.state)
    index = state[1]
    tmpResult = state[2]

    c = tmpResult.c
    ttime = tmpResult.ttime

    ttime = ttime + dtime
    dfdc = thermocalc.get_free_energy(c)
    dfdck = fft(dfdc)
    ck = fft(c)

    number = dtime * mobility * k2.* dfdck
    denom = 1.0 .+ dtime * coefA * mobility * grad_coef .* k4
    ck = (ck - number) ./ denom
    c = real(ifft(ck))
    c = array_operation.reject_element(c, (0.001, 0.999))

    res = TmpResult(c, ttime)

    nstep == index && return nothing
    return ((index, res), (index + 1, res))
  end

  itr = Iterator((1, TmpResult(c, ttime)))
  return itr

end

iparam = InitialParameter()
itr = ch_simple_binary_alloy(iparam)

# for res in itr
#   istep = res[1]
#   if (istep % nprint == 0) || istep == 1
#       println("step: $istep done")
#       img = Gray.(res[2].c)
#       save("result/png/res_$istep.png", img)
#   end
# end

# for istep in 1:nstep
#   ttime = ttime + dtime
#   dfdc = thermocalc.get_free_energy(c)
#   dfdck = fft(dfdc)

#   ck = fft(c)

#   number = dtime * mobility * k2.* dfdck
#   denom = 1.0 .+ dtime * coefA * mobility * grad_coef .* k4
#   ck = (ck - number) ./ denom
#   c = real(ifft(ck))
#   c = array_operation.reject_element(c, (0.001, 0.999))

#   if (istep % nprint == 0) || istep == 1
#     println("step: $istep done")
#     img = Gray.(c)
#     save("result/png/res_$istep.png", img)
#   end
# end


# end
