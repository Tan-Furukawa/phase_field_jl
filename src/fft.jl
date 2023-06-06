module Myfft

export prepare_fft
function prepare_fft(nx, ny, dx, dy, η)
  nx21 = Int(floor(nx / 2) + 1)
  ny21 = Int(floor(ny / 2) + 1)
  nx2 = Int(nx + 2)
  ny2 = Int(ny + 2)

  # L = nx * dx; delkx = 2π/L
  delkx = 2.0π / (nx * dx)
  delky = 2.0π / (ny * dy)

  kx = zeros(nx + 1)
  ky = zeros(ny + 1)
  k2 = zeros(nx, ny)
  k2_anyso = zeros(nx, ny)

  # kx = (0, π/2, -π, -π/2, 0)
  for i in 1:nx21
    fk1 = (i - 1.0) * delkx
    kx[i] = fk1
    kx[nx2-i] = -fk1
  end

  for i in 1:ny21
    fk1 = (i - 1.0) * delky
    ky[i] = fk1
    ky[ny2-i] = -fk1
  end

  for i in 1:nx
    for j in 1:ny
      k2[i, j] = kx[i]^2 + ky[j]^2
    end
  end

  for i in 1:nx
    for j in 1:ny
      k2_anyso[i, j] = η * kx[i]^2 + ky[j]^2
    end
  end

  k4 = k2 .^ 2

  return (kx, ky, k2, k4, k2_anyso)
end

# (kx, ky, k2, k4) = prepare_fft(4, 4, 1, 1)
end