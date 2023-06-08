include("fft.jl")
include("type.jl")
include("array_operation.jl")
include("thermocalc.jl")

# module SimpleBinaryAlloy

using FFTW
using ..Param
using ..thermocalc
using ..array_operation
using ..Myfft

export PhaseFiledResult
mutable struct PhaseFiledResult
  index:: Integer
  c:: AbstractMatrix
  ttime:: AbstractFloat
  is_last_object:: Bool
end

export ch_simple_binary_alloy
function ch_simple_binary_alloy(iparam:: Param.InitialParameter, pre_result::Union{PhaseFiledResult, Nothing} = nothing)
  return Channel{Tuple{Integer, PhaseFiledResult}}(0) do channel
    nx = iparam.nxny[1]; ny = iparam.nxny[2] 
    dx = iparam.dxdy[1]; dy = iparam.dxdy[2] 

    res = if (pre_result !== nothing) 
      println("hey")
      put!(channel, (1, pre_result))
      pre_result
    else
      c = thermocalc.make_initial_concentration(iparam.c0, iparam.nxny[1], iparam.nxny[2], iparam.seed, iparam.initial_noise_size)
      c = array_operation.reject_element(c, (0.001, 0.999))
      PhaseFiledResult(1, c, 0.0, false)
    end

    c = res.c

    (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, iparam.η)
    ttime = res.ttime
    nstep = iparam.nstep
    dt= iparam.dtime
    mobility = iparam.mobility
    coefA = iparam.coefA; grad_coef = iparam.grad_coef

    for istep in 1:nstep

      if nstep == istep 
        res.is_last_object = true
      end

      res.index == 1 && put!(channel, (istep, res))
      # istep == 1 && put!(channel, res)

      ttime = ttime + dt
      dfdc = thermocalc.get_free_energy(c)
      dfdck = fft(dfdc)
    
      ck = fft(c)
    
      number = dt * mobility * k2.* dfdck
      denom = 1.0 .+ dt * coefA * mobility * grad_coef .* k2 .* k2_anyso
      ck = (ck - number) ./ denom
      c = real(ifft(ck))
      c = array_operation.reject_element(c, (0.001, 0.999))
      res.c = c
      res.index = res.index + 1
      res.ttime = res.ttime
      put!(channel, (istep + 1, res))
    end
  end
end

# end


using Images


iparam = Param.InitialParameter(nstep = 10000, η = 2)
itr = ch_simple_binary_alloy(iparam)

res = nothing;
@time begin
for (i, k) in itr
    if (i % iparam.nprint == 0) || i == 1
        j = k.index
        println("step: $(j) done")
        img = Gray.(k.c')
        save("result/png/res_$(j).png", img)
    end
    if (k.is_last_object) 
      break
    end
end
end

iparam = Param.InitialParameter(nstep = 10000, η = 1)
itr2 = ch_simple_binary_alloy(iparam, take!(itr)[2])
for (i, k) in itr2
  if (i % iparam.nprint == 0) || i == 1
      j = k.index
      println("step: $(j) done")
      img = Gray.(k.c')
      save("result/png/res_$(j).png", img)
  end
end
