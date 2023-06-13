module NonLinearBinary

using Plots, Colors
using Dates
using FFTW
using Images
using ..Param
using ..ArrayOperation
using ..Thermocalc
using ..Myfft

export calculate_free_energy
function calculate_free_energy(c, η, w)
  # ∑∑f(cᵢⱼ) + 1/2 * ∑∑(η f(cᵢ₊₁ⱼ) - f(cᵢⱼ))² +  (f(cᵢⱼ₊₁) - f(cᵢⱼ))²)
  d = 0.001
  return sum(
    c .* log.(max.(c, d)) + (1.0 .- c) .* log.(max.(1.0 .- c, d)) 
   + w * c .* (1.0 .- c)
  ) + 1 / 2 * (
    sum((c[:, 2:end] - c[:, 1:end-1]).^2) * η + 
    sum((c[2:end, :] - c[1:end-1, :]).^2)
  )
end


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

  # explict (work)
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
      res.free_energy[istep] = log(calculate_free_energy(c, η, w))
      put!(channel, (istep + 1, res))
    end

  end
end

end
