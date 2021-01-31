__precompile__()
function exec()
    resize!(ARGS, 2)
    ARGS[1] = "../benchmarks/3QP/toy1"
    ARGS[2] = "toy1.out"
    include("qp.jl")
end

@time exec()
@time exec()

Base.init_depot_path()
Base.init_load_path()

@eval Module() begin
    Base.include(@__MODULE__, "qp.jl")
    for (pkgid, mod) in Base.loaded_modules
        if !(pkgid.name in ("Main", "Core", "Base"))
            eval(@__MODULE__, :(const $(Symbol(mod)) = $mod))
        end
    end
    for statement in readlines("app_precompile.jl")
        try
            Base.include_string(@__MODULE__, statement)
        catch
            # See julia issue #28808
            Core.println("failed to compile statement: ", statement)
        end
    end
end # module

empty!(LOAD_PATH)
empty!(DEPOT_PATH)
