module Myfft

export prepare_fft
function prepare_fft(nx, ny, dx, dy, η)
  nx21 = Int(floor(nx / 2) + 1)
  ny21 = Int(floor(ny / 2) + 1)
  nx2 = Int(nx + 2)
  ny2 = Int(ny + 2)

  # L = nx * dx = nx; dkx = 2π/L
  dkx = 2.0π / (nx * dx)
  dky = 2.0π / (ny * dy)

  kx = zeros(nx + 1)
  ky = zeros(ny + 1)
  k2 = zeros(nx, ny)
  k2_anyso = zeros(nx, ny)
  kx_mat = zeros(nx, ny)
  ky_mat = zeros(nx, ny)

  # kx = 2πm/L (m ∈ {z|z ∈ Z and -nx/2 <= z <= nx/2})
  # e.g. when nx = 4 then 
  #   L = 4
  #   kx = (0, π/2, -π, -π/2, 0)
  for i in 1:nx21
    fk1 = (i - 1.0) * dkx
    kx[i] = fk1
    kx[nx2-i] = -fk1
  end

  for i in 1:ny21
    fk1 = (i - 1.0) * dky
    ky[i] = fk1
    ky[ny2-i] = -fk1
  end

  for i in 1:nx
    for j in 1:ny
      k2[i, j] = kx[i]^2 + ky[j]^2
      kx_mat[i, j] = kx[i]
      ky_mat[i, j] = ky[j]
      k2_anyso[i, j] = η * kx[i]^2 + ky[j]^2
    end
  end

  k4 = k2 .^ 2

  return (kx_mat, ky_mat, k2, k4, k2_anyso)
end

# (kx, ky, k2, k4) = prepare_fft(4, 4, 1, 1)
end