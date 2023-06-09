module NonLinearBinary
using Plots, Colors
using Dates
using FFTW
using Images
using ..Param
using ..ArrayOperation
using ..Thermocalc
using ..Myfft
# include("fft.jl")
# include("array_operation.jl")
# include("Thermocalc.jl")
# include("type.jl")
# iparam = Param.InitialParameter()

export calculate_next_step_consentration!
function calculate_next_step_consentration!(c::Matrix{Float64}, fft_param::Param.FFTParameter, iparam::Param.InitialParameter)
  kx = fft_param.kx
  ky = fft_param.ky
  k2 = fft_param.k2
  k2_anyso = fft_param.k_anysotropy
  nx,ny = iparam.nxny
  c0 = iparam.c0
  w = iparam.w
  dt = iparam.dtime
  noise_per_step = iparam.noise_per_step
  g = c .* (1.0 .- c)

  ck = fft(c)
  gk = fft(g)

  ∂ₓc_k = ck .* kx .* im
  ∂ₓc = ifft(∂ₓc_k)
  ∂ₓg_k = gk .* kx * im
  ∂ₓg = ifft(∂ₓg_k)

  ∂yc_k = ck .* ky * im
  ∂yc = ifft(∂yc_k)
  ∂yg_k = gk .* ky * im
  ∂yg = ifft(∂yg_k)

  Δc_k = -k2 .* ck
  Δc = ifft(Δc_k)

  ∇_g∇c_k = fft(∂ₓg .* ∂ₓc + ∂yg .* ∂yc + g .* Δc)
  term1 = - 2.0w .* ∇_g∇c_k

  Δηc_k = - k2_anyso .* ck
  # Δηc = ifft(Δηc_k)

  ∂ₓΔηc_k = Δηc_k .* kx * im
  ∂ₓΔηc = ifft(∂ₓΔηc_k)
  ∂yΔηc_k = Δηc_k .* ky * im
  ∂yΔηc = ifft(∂yΔηc_k)

  ΔΔηc_k = k2_anyso .* k2 .* ck
  ΔΔηc = ifft(ΔΔηc_k)

  term2 = fft(∂ₓg .* ∂ₓΔηc + ∂yg .* ∂yΔηc + g .* ΔΔηc)

  # semi implicit (not work)
  # ck = ck + dt * (term1 - term2) ./ (1.0 .+ dt * k2) 

  # implict (work)
  ck = ck + dt * (Δc_k + term1 - term2)

  c = real(ifft(ck))
  c = ArrayOperation.reject_element!(c, (0.001, 0.999))
  c = Thermocalc.normalize_to_bulk!(c, c0)
  Thermocalc.make_random_noise!(c, nx, ny, 1234, noise_per_step)
  ArrayOperation.reject_element!(c, (0.001, 0.999))
  return c
end

export ch_nonlinear_binary_alloy
function ch_nonlinear_binary_alloy(iparam:: Param.InitialParameter, pre_result::Union{Tuple{Integer, Param.PhaseFiledResult}, Nothing} = nothing)
  return Channel{Tuple{Integer, Param.PhaseFiledResult}}(0) do channel
    nx, ny = iparam.nxny
    dx, dy = iparam.dxdy
    nstep = iparam.nstep
    dt = iparam.dtime
    c0 = iparam.c0
    mobility = iparam.mobility
    initial_noise_size = iparam.initial_noise_size
    seed = iparam.seed
    η = iparam.η
    w =  iparam.w
    noise_per_step = iparam.noise_per_step

    res = if (pre_result !== nothing) 
      put!(channel, (1, pre_result))
      _, phase_filed_res = pre_result
      phase_filed_res
    else
      c = Thermocalc.make_initial_concentration(c0, nx, ny, seed, initial_noise_size)
      c = ArrayOperation.reject_element!(c, (0.001, 0.999))
      Param.PhaseFiledResult(1, c, 0.0, false, Array{Float64}(undef, nstep))
    end

    c = res.c

    fft_params = Myfft.prepare_fft(nx, ny, dx, dy, η)
    ttime = res.ttime

    for istep in 1:nstep
      if nstep == istep 
        res.is_last_object = true
      end

      res.index == 1 && put!(channel, (istep, res))

      ttime = ttime + dt
      c = calculate_next_step_consentration!(c, fft_params, iparam)

      res.c = c
      res.index = res.index + 1
      res.ttime = res.ttime
      put!(channel, (istep + 1, res))

      # if (istep % nprint == 0) || istep == 1
      #   println("step: $istep done")
      #   img = Gray.(c')
      #   save("result/png/res_$istep.png", img)
      # end
  end
end

end
end

# if abspath(PROGRAM_FILE) == @__FILE__
#   using .NonLinearBinary
#   print("hello")
# end

# print(abspath(PROGRAM_FILE))


iparam = Main.Param.InitialParameter(nstep = 10000, η = 1.0, noise_per_step = 0.00, dtime = 0.01)
itr = Main.NonLinearBinary.ch_nonlinear_binary_alloy(iparam)
# take!(itr)

using Plots, Colors
using Dates
anim = Animation()
@time begin
for (i, k) in itr
    if (i % iparam.nprint == 0) || i == 1
        j = k.index
        c = k.c
        println("step: $(j) done")
        img = Gray.(c')
        plt = plot(img)
        frame(anim, plt)
        save("result/png/res_$(j).png", img)
    end
    if k.is_last_object
      break
    end
end
end

date_str = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS-SSS")
gif(anim, "result/gif/$(date_str).gif", fps = 1000)


# c = Thermocalc.make_initial_concentration(c0, nx, ny, seed, rand_size)
# c = ArrayOperation.reject_element(c, (0.001, 0.999))
# (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, η)

# function ts()
#   return Channel{String}(0) do channel
#     i = 0
#     while true
#       xxxx()
#       i = i + 1
#       put!(channel, "this is $(i)")
#     end
#   end
# end

# itr = ts()
# take!(itr)