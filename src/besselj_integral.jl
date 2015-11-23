using Calculus
using Gadfly
using LaTeXStrings # for L"" strings

# WARNING: integrate(f,a,b) is deprecated, use (quadgk(f,a,b))[1] instead
f(x) = integrate(z -> besselj(1, z), 0.0, x) 

my_plot_pgf = plot(f, 0, 100
            , Guide.XLabel("r"
                    , orientation=:horizontal)
#            , Guide.YLabel(L"\int_{z=0}^{z=\textbf{r}}J_1(z)dz"
            , Guide.YLabel("\$\\int_{z=0}^{z=\\textbf{r}}J_1(z)dz\$"
                    , orientation=:vertical))

draw(PGF("images/Int_J1_tex2.tex", 5inch, 3inch), my_plot_pgf)

# see [PGF LaTeX output seems to be broken #753](https://github.com/dcjones/Gadfly.jl/issues/753)
