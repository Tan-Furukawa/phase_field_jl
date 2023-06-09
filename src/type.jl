module Param

export f
function f()
end

export InitialParameter
Base.@kwdef struct InitialParameter
  nxny:: Tuple{Integer, Integer} = (64, 64)
  dxdy:: Tuple{Float64, Float64} = (1.0, 1.0)
  nstep:: Integer = 20000
  nprint:: Integer = 100
  c0:: Float64 = 0.4
  η:: Float64 = 2.0
  noise_per_step:: Float64 = 0.0
  grad_coef:: Float64 = 0.3
  dtime:: Float64 = 0.01
  coefA:: Float64 = 1.0
  mobility:: Float64 = 1.0
  initial_noise_size:: Float64 = 0.02
  seed:: Integer = 123
  w:: Float64 = 3.0
end

export FFTParameter
struct FFTParameter
  kx::Matrix{Float64}
  ky::Matrix{Float64}
  k2::Matrix{Float64}
  k4::Matrix{Float64}
  k_anysotropy:: Matrix{Float64}
end

export PhaseFiledResult
mutable struct PhaseFiledResult
  index:: Integer
  c:: Matrix{Float64}
  ttime:: Float64
  is_last_object:: Bool
  free_energy:: Array{Float64}
end


end