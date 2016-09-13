using Calculus
using Gadfly
using LaTeXStrings # for L"" strings

# WARNING: integrate(f,a,b) is deprecated, use (quadgk(f,a,b))[1] instead
# f(x) = integrate(z -> besselj(1, z), 0.0, x)

f(x) = quadgk(z -> besselj(1, z), 0.0, x)[1]

function my_draw(f)
    my_plot_pgf = plot(f, 0, 100
                , Guide.XLabel("r"
                        , orientation=:horizontal)
                , Guide.YLabel(String(L"\int_{z=0}^{z=\textbf{r}}J_1(z)dz") # seems Gadfly doesn't support LaTeX L"" anymore
                                                                    # workaround: edit manually PGF .tex file
                #, Guide.YLabel("∫₀ʳJ₁(z)dz" # it is very ugly
                        , orientation=:vertical))
    
    draw(PGF("images/Int_J1_tex2.tex", 5inch, 3inch), my_plot_pgf)
end

# Gadfly used to support LaTeX labels but doesn't support anymore
# [Does Gadfly support LaTeX in plot titles and labels?](http://stackoverflow.com/questions/36880531/does-gadfly-support-latex-in-plot-titles-and-labels)
# see also [PGF LaTeX output seems to be broken #753](https://github.com/dcjones/Gadfly.jl/issues/753)

@time my_draw(f)
@time my_draw(f)

#using FastAnonymous
#y = @anon z -> besselj(1, z)
# FastAnonymous doesn't work for 0.5, 0.6

y = z -> besselj(1, z)
g(x) = quadgk(y, 0.0, x)[1]

@time my_draw(g)
@time my_draw(g)

