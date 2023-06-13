iparam = Main.Param.InitialParameter(nstep = 20000, η = 2.0, noise_per_step = 0.00, dtime = 0.01, w = 3, nxny = (128, 128))
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
gif(anim, "result/gif/$(date_str).gif", fps = 100)
_, res = take!(itr)
plot(1:iparam.nstep, res.free_energy)


# (_, res) = take!(itr)
# plot(1:iparam.nstep, res.free_energy)