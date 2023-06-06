module thermocalc
using Random

export make_initial_concentration
function make_initial_concentration(c0::Float64, nx::Int, ny::Int, seed::Int, rand_size::Float64)
  Random.seed!(seed)
  c = rand(nx, ny)
  c = c0 .+ (0.5 .- c) * rand_size
  return c
end

export get_free_energy
function get_free_energy(c:: AbstractArray)
  return  2.0 .* c .* (1.0 .- c).^2 .- 2.0 .* c.^2 .* (1.0 .- c)
end

end