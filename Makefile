# julia samples

jinc:
	julia src/jinc.jl

besselj_integral:
	julia src/besselj_integral.jl
	xelatex -output-directory=images images/Int_J1_tex2.tex

besselj_integral_pgf:
	julia src/besselj_integral_pgfplots.jl
	xelatex -output-directory=images images/Int_J1_tex3.tex
qp:
	julia src/qp.jl benchmarks/3QP/toy1 toy1.out

