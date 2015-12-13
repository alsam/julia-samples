# julia samples

jinc:
	julia src/jinc.jl

besselj_integral:
	julia src/besselj_integral.jl
	xelatex -output-directory=images images/Int_J1_tex2.tex

qp:
	julia src/qp.jl benchmarks/3QP/toy1 toy1.out

