__precompile__()
function exec()
    resize!(ARGS, 2)
    ARGS[1] = "../benchmarks/3QP/toy1"
    ARGS[2] = "toy1.out"
    include("qp.jl")
end

@time exec()
@time exec()
