# computes $\frac{J_1(2\pi r)}{2\pi r}$ where $J_1$ - Bessel function of first kind
jinc(r) = real(2*besselj(1,2π*r) / (2π*r))

using Roots

jinc_zeros = fzeros(jinc, 0, 5)
println("jinc_zeros : $jinc_zeros")

using Gadfly

my_plot = plot(jinc, 0, 5
        , Guide.XLabel("X"), Guide.YLabel("jinc")
        , Guide.xticks(ticks=[0, 0.61, 1.12, 1.62, 2.12, 2.62, 3.12, 3.62, 4.12, 4.62] ))

draw(PDF("images/jinc.pdf", 5inch, 3inch), my_plot)
draw(PNG("images/jinc.png", 5inch, 3inch), my_plot)
