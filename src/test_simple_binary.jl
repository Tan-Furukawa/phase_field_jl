
iparam = Main.Param.InitialParameter(nstep = 20000, η = 4, noise_per_step = 0.01 , mobility = 20, nxny=(64, 64))
itr = Main.SimpleBinaryAlloy.ch_simple_binary_alloy(iparam)
# itr = SimpleBinaryAlloy.ch_simple_binary_alloy(iparam)
using Colors, Plots
using Dates

function normalize_c(c)
  return 0.2 .* c + 0.8 .* (1.0 .- c)
end

anim = Animation()
@time begin
for (i, k) in itr
    if (i % iparam.nprint == 0) || i == 1
        j = k.index
        c = normalize_c(k.c)
        img = Gray.(c')
        println("step: $(j) done")
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
gif(anim, "result/gif/$(date_str).gif", fps = 10)

iparam = Main.Param.InitialParameter(nstep = 2000, η = 3, noise_per_step = 0.01 , mobility=0.1, nxny=(64, 64))

res = take!(itr)
res[2].c = normalize_c(res[2].c)
itr2 = Main.SimpleBinaryAlloy.ch_simple_binary_alloy(iparam, res)

for (i, k) in itr2
    if (i % iparam.nprint == 0) || i == 1
        j = k.index
        c = k.c
        img = Gray.(c')
        println("step: $(j) done")
        plt = plot(img)
        frame(anim, plt)
        save("result/png/res_$(j).png", img)
    end
    # if k.is_last_object
    #   break
    # end
end

date_str = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS-SSS")
gif(anim, "result/gif/$(date_str).gif", fps = 10)