module SimpleBinaryAlloy

using FFTW
using ..Param
using ..Thermocalc
using ..ArrayOperation
using ..Myfft


export calculate_next_step_consentration!
function calculate_next_step_consentration!(c::Matrix{Float64}, fft_param::Param.FFTParameter, iparam::Param.InitialParameter)
  k2 = fft_param.k2
  k2_anyso = fft_param.k_anysotropy
  dfdc = Thermocalc.get_free_energy(c)
  dfdck = fft(dfdc)
  ck = fft(c)
  dt = iparam.dtime
  mobility = iparam.mobility
  coefA = iparam.coefA
  grad_coef = iparam.grad_coef
  (nx, ny) = iparam.nxny

  number = dt * mobility * k2.* dfdck
  denom = 1.0 .+ dt * coefA * mobility * grad_coef .* k2 .* k2_anyso
  ck = (ck - number) ./ denom
  c = real(ifft(ck))
  c = Thermocalc.normalize_to_bulk!(c, iparam.c0)
  Thermocalc.make_random_noise!(c, nx, ny, 123, iparam.noise_per_step)
  ArrayOperation.reject_element!(c, (0.001, 0.999))

  return c
end

export ch_simple_binary_alloy
function ch_simple_binary_alloy(iparam:: Param.InitialParameter, pre_result::Union{Tuple{Integer, Param.PhaseFiledResult}, Nothing} = nothing)
  return Channel{Tuple{Integer, Param.PhaseFiledResult}}(0) do channel
    nx, ny = iparam.nxny
    dx, dy = iparam.dxdy
    nstep = iparam.nstep
    η = iparam.η

    res = if (pre_result !== nothing) 
      put!(channel, pre_result)
      _, phase_filed_res = pre_result
      phase_filed_res
    else
      c = Thermocalc.make_initial_concentration(iparam.c0, nx, ny, iparam.seed, iparam.initial_noise_size)
      c = ArrayOperation.reject_element!(c, (0.001, 0.999))
      Param.PhaseFiledResult(1, c, 0.0, false, Array{Float64}(undef, nstep))
    end

    c = res.c

    fft_params = Myfft.prepare_fft(nx, ny, dx, dy, η)
    # (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, η)
    ttime = res.ttime
    dt= iparam.dtime
    grad_coef = iparam.grad_coef

    for istep in 1:nstep
      c = res.c

      if nstep == istep 
        res.is_last_object = true
      end

      res.index == 1 && put!(channel, (istep, res))
      # istep == 1 && put!(channel, res)

      ttime = ttime + dt

      c = Thermocalc.make_random_noise!(c, nx, ny, 1234, iparam.noise_per_step)
      c = calculate_next_step_consentration!(c, fft_params, iparam)

      res.c = c
      res.index = res.index + 1
      res.ttime = res.ttime
      res.free_energy[istep] = Thermocalc.calculate_free_energy(nx, ny, c, grad_coef)
      put!(channel, (istep + 1, res))
    end
  end
end
end


# function normalize_c(c)
#   return 0.3 .* c + 0.7 .* (1.0 .- c)
# end

# iparam = Param.InitialParameter(nstep = 20000, η = 1.0, noise_per_step = 0.00 , mobility = 10)
# itr = ch_simple_binary_alloy(iparam)
# # itr = SimpleBinaryAlloy.ch_simple_binary_alloy(iparam)
# using Colors, Plots
# using Dates


# anim = Animation()
# @time begin
# for (i, k) in itr
#     if (i % iparam.nprint == 0) || i == 1
#         j = k.index
#         c = k.c
#         println("step: $(j) done")
#         plt = plot(Gray.(c'))
#         frame(anim, plt)
#         # save("result/png/res_$(j).png", img)
#     end
#     if k.is_last_object
#       break
#     end
# end
# end

# date_str = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS-SSS")
# gif(anim, "result/gif/$(date_str).gif", fps = 30)


# (nx, ny) = iparam.nxny


# using Plots

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
