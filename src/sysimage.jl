using PackageCompiler
create_sysimage(:IterativeSolvers, sysimage_path="sys_iterative_solvers.so", precompile_execution_file="precompile_qp.jl")
#create_sysimage(:IterativeSolvers, sysimage_path="sys_iterative_solvers.so", precompile_execution_file="custom_sysimage.jl")

