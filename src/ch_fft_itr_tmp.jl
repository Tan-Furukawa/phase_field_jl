# module ch_fft

module SimpleBinaryAlloy

using FFTW
using ..Param
using ..thermocalc
using ..array_operation

export ChFftSimpleBinaryAlloy
@kwdef mutable struct ChFftSimpleBinaryAlloy
  iparam:: Param.InitialParameter
  ttime:: AbstractFloat = 0
  index:: Integer = 1
  fft_param:: Param.FFTParameter = thermocalc.make_initial_fft_param(iparam)
  c:: AbstractArray = thermocalc.make_initial_concentration_from_iparam(iparam)
end

function Base.iterate(iter::ChFftSimpleBinaryAlloy, state = iter)
  index = state.index

  index < 1 && return nothing
  if index == 1
    state.index = 2
    return state, state
  end

  iparam = state.iparam

  dt = iparam.dtime
  tmpResult = state

  c = tmpResult.c
  ttime = tmpResult.ttime

  ttime = ttime + dt
  dfdc = thermocalc.get_free_energy(c)
  dfdck = fft(dfdc)
  ck = fft(c)

  k2 = state.fft_param.k2
  k4 = state.fft_param.k4

  number = dt * iparam.mobility * k2.* dfdck
  denom = 1.0 .+ dt * iparam.coefA * iparam.mobility * iparam.grad_coef .* k4
  ck = (ck - number) ./ denom
  c = real(ifft(ck))
  c = array_operation.reject_element(c, (0.001, 0.999))

  # res = ChFftSimpleBinaryAlloy(index + 1, c, ttime)
  state.c = c
  state.index = index + 1
  state.ttime = ttime

  iparam.nstep == index && return nothing
  return state, state
end

using Images
iparam = Param.InitialParameter()
itr = ChFftSimpleBinaryAlloy(iparam = iparam, index = 1)

for res in itr
  istep = res.index
  if (istep % iparam.nprint == 0) || istep == 1
      println("step: $istep done")
      img = Gray.(res.c)
      save("result/png/res_$istep.png", img)
  end
end

end
