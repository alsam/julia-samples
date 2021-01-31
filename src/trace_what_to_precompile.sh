#!/bin/sh
julia --startup-file=no --trace-compile=app_precompile.jl precompile_qp.jl
