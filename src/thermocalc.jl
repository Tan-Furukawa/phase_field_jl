# include("array_operation.jl")
# include("type.jl")
# include("fft.jl")


module thermocalc
using ..Myfft
using ..Param
using ..array_operation

using Random
export make_initial_concentration
function make_initial_concentration(c0::Float64, nx::Int, ny::Int, seed::Int, rand_size::Float64)
  Random.seed!(seed)
  c = rand(nx, ny)
  c = c0 .+ (0.5 .- c) * rand_size
  c = array_operation.reject_element(c, (0.001, 0.999))
  return c
end

export get_free_energy
function get_free_energy(c:: AbstractArray)
  return  2.0 .* c .* (1.0 .- c).^2 .- 2.0 .* c.^2 .* (1.0 .- c)
end

export make_initial_concentration_from_iparam
function make_initial_concentration_from_iparam(iparam:: Param.InitialParameter)::AbstractArray
  nx = iparam.nxny[1]
  ny = iparam.nxny[2]
  # set initial consentration
  c = thermocalc.make_initial_concentration(iparam.c0, nx, ny, iparam.seed, iparam.initial_noise_size)
  return c
end

export make_initial_fft_param
function make_initial_fft_param(iparam:: Param.InitialParameter)::Param.FFTParameter
  nx = iparam.nxny[1]
  ny = iparam.nxny[2]
  dx = iparam.dxdy[1]
  dy = iparam.dxdy[2]
  (kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, iparam.η)
  return Param.FFTParameter(kx, ky, k2, k4, k2_anyso)
end

end