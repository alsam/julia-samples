using Calculus
using PGFPlots
using QuadGK

# WARNING: integrate(f,a,b) is deprecated, use (quadgk(f,a,b))[1] instead
# f(x) = integrate(z -> besselj(1, z), 0.0, x)

f(x) = QuadGK.quadgk(z -> besselj(1, z), 0.0, x)[1]

function my_draw(f)
    my_plot_pgf = Axis([
                        Plots.Linear(f, (0, 100), style="smooth", legendentry=L"\int_{z=0}^{z=\textbf{r}}J_1(z)dz"),
                    ], legendPos="north west")
    
    save("images/Int_J1_tex3.tex", my_plot_pgf)
end

@time my_draw(f)
@time my_draw(f)

y = z -> besselj(1, z)
g(x) = QuadGK.quadgk(y, 0.0, x)[1]

@time my_draw(g)
@time my_draw(g)

