using ch_fft
using Test

@testset "ch_fft" begin
    using Plots, Colors, Dates

    anim = Animation()
    iparam = ch_fft.Param.InitialParameter(nstep = 1000, η = 1.0, noise_per_step = 0.00 , mobility = 1)
    itr = ch_fft.NonLinearBinary.ch_nonlinear_binary_alloy(iparam)

    @time begin
    for (i, k) in itr
        if (i % iparam.nprint == 0) || i == 1
            j = k.index
            c = k.c
            println("step: $(j) done")
            plt = plot(Gray.(c'))
            frame(anim, plt)
            # save("result/png/res_$(j).png", img)
        end
        if k.is_last_object
          break
        end
    end
    end

    date_str = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS-SSS")
    gif(anim, "result/gif/$(date_str).gif", fps = 30)
end
