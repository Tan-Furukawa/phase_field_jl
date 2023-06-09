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

export ch_simple_binary_alloy
function ch_simple_binary_alloy(iparam:: Param.InitialParameter, pre_result::Union{Tuple{Integer, Param.PhaseFiledResult}, Nothing} = nothing)
  return Channel{Tuple{Integer, Param.PhaseFiledResult}}(0) do channel
    nx, ny = iparam.nxny
    dx, dy = iparam.dxdy
    nstep = iparam.nstep

    res = if (pre_result !== nothing) 
      put!(channel, (1, pre_result))
      _, phase_filed_res = pre_result
      phase_filed_res
    else
      c = thermocalc.make_initial_concentration(iparam.c0, nx, ny, iparam.seed, iparam.initial_noise_size)
      c = array_operation.reject_element!(c, (0.001, 0.999))
      Param.PhaseFiledResult(1, c, 0.0, false, Array{Float64}(undef, nstep))
    end

    c = res.c

    (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, iparam.η)
    ttime = res.ttime
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
      thermocalc.make_random_noise!(c, nx, ny, istep, iparam.noise_per_step)
      array_operation.reject_element!(c, (0.001, 0.999))
      res.c = c
      res.index = res.index + 1
      res.ttime = res.ttime
      res.free_energy[istep] = thermocalc.calculate_free_energy(nx, ny, c, grad_coef)
      put!(channel, (istep + 1, res))
    end
  end
end
# end


function normalize_c(c)
  return 0.3 .* c + 0.7 .* (1.0 .- c)
end

iparam = Param.InitialParameter(nstep = 20000, η = 1.0, noise_per_step = 0.00 , mobility = 1.0)
itr = ch_simple_binary_alloy(iparam)
# itr = SimpleBinaryAlloy.ch_simple_binary_alloy(iparam)


using Colors, Plots
@time begin
for (i, k) in itr
    if (i % iparam.nprint == 0) || i == 1
        j = k.index
        c = k.c
        println("step: $(j) done")
        img = Gray.(c')
        # save("result/png/res_$(j).png", img)
    end
    if k.is_last_object
      break
    end
end
end


using PyPlot
(_, res) = take!(itr)

(nx, ny) = iparam.nxny

# x = 1:nx
# y = 1:ny
# pcolormesh(res.c, cmap="prism")
# colorbar()

# using PyPlot

# prs = [1 2 3;4 5 6;7 8 9]

# prsview = zeros(size(prs))
# for k = 1:size(prs,1)
#     prsview[k,:] = prs[size(prs,1)-k+1,:]
# end

# pcolormesh(prsview,cmap="Blues")
# colorbar()
# savefig("Matrix.png")

# plot((res.c))

# using Plots
# plot(1:iparam.nstep, res.free_energy)

# iparam = Param.InitialParameter(nstep = 10000, η = 2, noise_per_step = 0.02, mobility = 1)
# sum(res.c)/ length(res.c)
# res.c = normalize_c(res.c)
# itr2 = ch_simple_binary_alloy(iparam, res)
# for (i, k) in itr2
#   if (i % iparam.nprint == 0) || i == 1
#       j = k.index
#       println("step: $(j) done")
#       img = Gray.(k.c')
#       save("result/png/res_$(j).png", img)
#   end
# end
