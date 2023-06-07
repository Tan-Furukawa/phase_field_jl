# module ch_fft
using FFTW
using Images
include("fft.jl")
include("array_operation.jl")
include("thermocalc.jl")

nx = 64
ny = 64
dx = 1.0
dy = 1.0
nstep = 4000
nprint = 100
dtime = 0.01
coefA = 1.0
c0 = 0.4
mobility = 1.0
grad_coef = 0.5
rand_size = 0.2
seed = 123

η = 2.0
ttime = 0.0
w =  2.5

c = thermocalc.make_initial_concentration(c0, nx, ny, seed, rand_size)
c = array_operation.reject_element(c, (0.001, 0.999))
(kx, ky, k2, k4, k2_anyso) = Myfft.prepare_fft(nx, ny, dx, dy, η)

istep = 1
for istep in 1:nstep
  ttime = ttime + dtime
  c[:,1] = c[:,end]
  c[1,:] = c[end,:]

  gg = c .* (1.0 .- c)
  # sg = vcat(gg[2:nx,:], gg[1,:])
  # dg = sg - gg

  gk = fft(gg)
  ck = fft(c)
  g = ifft(gk)

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
  Δηc = ifft(Δηc_k)

  ∂ₓΔηc_k = Δηc_k .* kx * im
  ∂ₓΔηc = ifft(∂ₓΔηc_k)
  ∂yΔηc_k = Δηc_k .* ky * im
  ∂yΔηc = ifft(∂yΔηc_k)

  ΔΔηc_k = k2_anyso .* k2 .* ck
  ΔΔηc = ifft(ΔΔηc_k)

  term2 = fft(∂ₓg .* ∂ₓΔηc + ∂yg .* ∂yΔηc + g .* ΔΔηc)

  # ck = ck + dtime * (term1 - term2) ./ (1.0 .+ dtime * k2) 
  ck = ck + dtime * (term1 - term2 + Δc_k) 

  ck = ck + dtime * (term1 + Δc_k) 

  # OK; ck = ck + dtime * Δc_k

  # istep == 10 && break

  c = real(ifft(ck))

  c = array_operation.reject_element(c, (0.001, 0.999))

  if (istep % nprint == 0) || istep == 1
    println("step: $istep done")
    img = Gray.(c')
    save("result/png/res_$istep.png", img)
  end
end

