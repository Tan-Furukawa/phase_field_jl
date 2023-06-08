module Param

export f
function f()
end

export InitialParameter
Base.@kwdef struct InitialParameter
  nxny:: Tuple{Integer, Integer} = (64, 64)
  dxdy:: Tuple{AbstractFloat, AbstractFloat} = (1.0, 1.0)
  nstep:: Integer = 20000
  nprint:: Integer = 100
  c0:: AbstractFloat = 0.4
  η:: AbstractFloat = 2.0
  grad_coef:: AbstractFloat = 0.5
  dtime:: AbstractFloat = 0.01
  coefA:: AbstractFloat = 1.0
  mobility:: AbstractFloat = 1.0
  initial_noise_size:: AbstractFloat = 0.02
  seed:: Integer = 123
end

export FFTParameter
struct FFTParameter
  kx::Matrix{AbstractFloat}
  ky::Matrix{AbstractFloat}
  k2::Matrix{AbstractFloat}
  k4::Matrix{AbstractFloat}
  k_anysotropy:: Matrix{AbstractFloat}
end

end