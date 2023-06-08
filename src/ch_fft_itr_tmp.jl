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

# arr = 1:100  # 1から10,000までの範囲の配列
# show(arr, all=true, limit=true)

function Base.iterate(iter::ChFftSimpleBinaryAlloy, state = iter)
  index = state.index
  
  index < 1 && return nothing
  if index == 1
    state.index = 2
    return state, state
  end

  iparam = state.iparam

  dt = iparam.dtime

  ttime = state.ttime

  ttime = ttime + dt
  dfdc = thermocalc.get_free_energy(state.c)
  dfdck = fft(dfdc)
  ck = fft(state.c)

  k2 = state.fft_param.k2
  k4 = state.fft_param.k4

  # イテレータにするとここで毎回アロケーションが発生しちゃう！！
  # number = dt * iparam.mobility * k2 .* dfdck
  # denom = 1.0 .+ dt * iparam.coefA * iparam.mobility * iparam.grad_coef .* k4

  ck = (ck - dt * iparam.mobility * k2 .* dfdck) ./ (1.0 .+ dt * iparam.coefA * iparam.mobility * iparam.grad_coef .* k4)
  state.c = real(ifft(ck))
  state.c = array_operation.reject_element(state.c, (0.001, 0.999))

  state.index = index + 1
  state.ttime = ttime

  iparam.nstep == index && return nothing

  return state, state
end

function my_function(x, y)
  z = x + y
  return z
end

function zundokochannel()
  return Channel{Tuple{Integer, String}}(32) do channel
      i = 1
      while true
          i = i + 1
          put!(channel, (i, "ズン"))
      end
      put!(channel, (i, "キ・ヨ・シ！"))
  end
end

d = zundokochannel()
println(d)
for (i, z) in d
  println(i)
  if i == 100 || i == 102
    break
  end
end



# iparam = Param.InitialParameter()
# itr = ChFftSimpleBinaryAlloy(iparam = iparam)
# iterate(itr)
# iterate(itr)

end

using .SimpleBinaryAlloy
iparam = Param.InitialParameter()
itr = ChFftSimpleBinaryAlloy(iparam = iparam)
iterate(itr)
@time begin
iterate(itr)
end



# @kwdef struct dd
#   c:: AbstractArray = zeros(100,100) .+ 1
# end

# d = dd()

# @time begin
# d1 = d.c
# end #389 allocations

# @time begin
#   # d1 = d.c
#   e = d1 .* 1
# end #81.86 k allocations


# d2 = zeros(100,100) .+ 1
# @time begin
#   h = d2 .* d2
# end

# using Images
# @time begin
# iparam = Param.InitialParameter()
# itr = ChFftSimpleBinaryAlloy(iparam = iparam)

# for res in itr
#   istep = res.index
#   if (istep % iparam.nprint == 0) || istep == 1
#       println("step: $istep done")
#       img = Gray.(res.c)
#       save("result/png/res_$istep.png", img)
#   end
# end

# end

# end

