# The MIT License (MIT)
# 
# Copyright (c) 2014 Alexander Samoilov
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

using IterativeSolvers
using Base.Test

type Pos
    x::Float64
    y::Float64

    Pos(x::Float64 = -1.0, y::Float64 = -1.0) = new(x, y)
end

type Gate
    id::UInt
    pos::Pos
    nets::Vector{UInt}

    Gate(id::UInt, nets) = new(id, Pos(), nets)
end

type Pad
    pos::Pos
    net_num::UInt

    Pad(pos::Pos, net::UInt) = new(pos, net)
end

type Net
    weight::Float64

    # convention : positive number - gate number, negative one - inverted pad number
    gates::Vector{UInt}
    pads::Vector{UInt}

    gates_set::IntSet
    pads_set::IntSet

    Net(weight::Float64 = 1.0) = new(weight, UInt[], UInt[], IntSet(), IntSet())
end

function net_is_empty(n::Net)
    length(gates) == 0 && length(pads) == 0
end

function read_input(fname::String)

    is = open(fname, "r")
    line = readline(is)
    nums = split(line)
    (num_gates::UInt, num_nets::UInt) = (parse(nums[1]), parse(nums[2]))
    println("-I- read num_gates: $num_gates num_nets: $num_nets")

    # read gates
    gates = Array(Gate, num_gates)
    for i = 1:num_gates
        line = readline(is)
        nums = split(line)
        (gate_num, num_nets_connected) = (parse(nums[1]), parse(nums[2]))
        @assert(gate_num == i, "gate_num: $gate_num i: $i")
        nets = [parse(nums[j + 2]) for j = 1:num_nets_connected]
        gates[i] = Gate(i, nets)
    end

    # read pads
    line = readline(is)
    num_pads::UInt = parse(line)
    println("-I- read num_pads: $num_pads")
    pads = Array(Pad, num_pads)
    for i = 1:num_pads
        line = readline(is)
        nums = split(line)
        (pad_num::UInt, net_num::UInt, x::Float64, y::Float64) = (parse(nums[1]), parse(nums[2]), parse(nums[3]), parse(nums[4]))
        @assert(pad_num == i, "pad_num: $pad_num i: $i")
        pads[i] = Pad(Pos(x,y), net_num)

    end

    close(is)
    return (num_gates, num_nets, num_pads, gates, pads)
end

function fill_net_sets(n::Net)
    n.gates_set = IntSet(n.gates)
    n.pads_set  = IntSet(n.pads)
end

function create_nets(num_nets::UInt, gates::Vector{Gate}, pads::Vector{Pad})
    num_gates::UInt = length(gates)
    num_pads::UInt  = length(pads)
    nets = [Net() for i = 1:num_nets]
    for i = 1:num_gates
        gate = gates[i]
        for j = 1:length(gate.nets)
            net_no = gate.nets[j]
            push!(nets[net_no].gates, i)
        end
    end

    for i = 1:num_pads
        push!(nets[pads[i].net_num].pads, i)
    end

    for i = 1:num_nets
        fill_net_sets(nets[i])
    end

    return nets
end

function qp_step(nets::Vector{Net}, gates::Vector{Gate}, pads::Vector{Pad})
    num_nets   ::UInt = length(nets)
    num_gates  ::UInt = length(gates)
    num_pads   ::UInt = length(pads)

    # connectivity matrix
    C = spzeros(num_gates, num_gates)
    for i = 1:num_nets
        net = nets[i]
        gs = net.gates
        ps = net.pads
        clique_size = length(gs)
        @printf("-D- NET : %d\n",i)
        for j = 1:clique_size
            for k = j+1:clique_size
                jc = gs[j]
                kc = gs[k]
                @printf("-D-      Conn %d <-> %d\n", jc, kc)
                C[jc, kc] = C[kc, jc] = net.weight
            end
        end
    end

    println("-D- C : \n$C")

    A = spzeros(num_gates, num_gates)
    # use dense vectors for rhs 
    bx = zeros(num_gates)
    by = zeros(num_gates)

    A = -C
    for i = 1:num_gates
        # see https://github.com/JuliaLang/julia/issues/6480
        diag_val::Float64 = sum(C[i, :])[1]
        pads_contrib::Float64 = 0.0
        b_accum = (0.0, 0.0)
        adjacent_nets = gates[i].nets
        for j = 1:length(adjacent_nets)
            net_num = adjacent_nets[j]
            the_net = nets[net_num]
            ps = the_net.pads
            net_weight = the_net.weight
            for pad_idx in ps
                pad = pads[pad_idx]
                pads_contrib += net_weight
                (x, y) = (pad.pos.x, pad.pos.y)
                b_accum = (b_accum[1] + net_weight * x, b_accum[2] + net_weight * y)
            end
        end
        @printf("-D- i: %i diag_val: %f pads_contrib: %f\n",i,diag_val,pads_contrib)
        A[i, i] = diag_val + pads_contrib
        (bx[i], by[i]) = b_accum
    end

    println("-D- A : \n$A")
    println("-D- bx : \n$bx")
    println("-D- by : \n$by")

    maxiter = 200
    tol = 1e-12
    x,ch = cg(A,bx,tol=tol,maxiter=maxiter)
    #see https://github.com/JuliaLang/julia/issues/6485#issuecomment-40063871
    @test_approx_eq_eps A*x bx cond(full(A))*sqrt(tol)
    @test ch.isconverged
    println("ch_x : $ch")
    
    y,ch = cg(A,by,tol=tol,maxiter=maxiter)
    @test_approx_eq_eps A*y by cond(full(A))*sqrt(tol)
    @test ch.isconverged
    println("ch_y : $ch")

    return [(i,Pos(x[i], y[i])) for i = 1:num_gates]
end

my_isless(a::Tuple{UInt,Pos}, b::Tuple{UInt,Pos}) = a[2].x * 1e6 + a[2].y < b[2].x * 1e6 + b[2].y

typealias IndexedPosVec Vector{Tuple{UInt,Pos}}

function assignment_cut(gate_pos::IndexedPosVec)
    num_gates = length(gate_pos)

    # avoid using anonymous functions
    sorted_gate_pos = sort!(gate_pos, alg = QuickSort, lt = my_isless)
    println("-D- sorted 1st plc : \n$sorted_gate_pos")

    median::Int = div(num_gates, 2)
    return ( sorted_gate_pos[1:median], sorted_gate_pos[median+1:num_gates] )
end

function left_side_containment(left_pos  ::IndexedPosVec,
                               right_pos ::IndexedPosVec,
                               nets      ::Vector{Net},
                               gates     ::Vector{Gate},
                               pads      ::Vector{Pad})
    num_gates = length(gates)
    left_gates = [deepcopy(gates[left_pos[i][1]]) for i = 1:length(left_pos)]
    left_gates_renumbering = Array(Int, num_gates)
    fill!(left_gates_renumbering, -1)
    for i = 1:length(left_gates)
        left_gates_renumbering[left_gates[i].id] = i
    end
    left_pads = deepcopy(pads)
    println("-D- left_gates : $left_gates")
    println("-D- left_gates_renumbering : $left_gates_renumbering")
    left_gates_indices  = IntSet([left_pos[i][1] for i = 1:length(left_pos)])
    println("-D- left gates indices : $left_gates_indices")
    right_gates_indices = IntSet([right_pos[i][1] for i = 1:length(right_pos)])
    println("-D- right gates indices : $right_gates_indices")

    const median::Float64 = 50.0

    num_nets = length(nets)
    left_nets = [Net() for i = 1:num_nets]
    for i = 1:num_nets
        orig_net = nets[i]
        if isempty(intersect(orig_net.gates_set,left_gates_indices))
            println("-D- the orig_net $orig_net DOESN'T belongs to left partition, EXCLUDE it")
        else
            net_no::UInt = i
            println("-D- the orig_net $orig_net belongs to left partition, fix it and include")
            orig_pads_indices = orig_net.pads
            new_pads_indices = UInt[]
            # process orig pads
            for orig_pad_idx in orig_pads_indices
                orig_pad = pads[orig_pad_idx]
                # fix abscissa if needed
                if orig_pad.pos.x <= median # ok, add this pad_idx AS IS
                    push!(new_pads_indices, orig_pad_idx)
                else
                    new_pad = deepcopy(orig_pad)
                    new_pad.pos.x = median
                    push!(left_pads, new_pad)
                    new_pad_idx::UInt = length(left_pads)
                    push!(new_pads_indices, new_pad_idx)
                end
            end
            orig_gates_indices = orig_net.gates
            new_gates_indices = UInt[]
            # process orig gates
            for orig_gate_idx in orig_gates_indices
                if orig_gate_idx in left_gates_indices # Bingo!
                    @assert(left_gates_renumbering[orig_gate_idx] > 0)
                    remapping::UInt = left_gates_renumbering[orig_gate_idx]
                    push!(new_gates_indices, remapping)
                else # gate belongs to right partition, convert it to pad
                    orig_gate = gates[orig_gate_idx]
                    new_pad = Pad(Pos(median, orig_gate.pos.y), net_no)
                    push!(left_pads, new_pad)
                    new_pad_idx::UInt = length(left_pads)
                    push!(new_pads_indices, new_pad_idx)
                end
            end
            the_net = left_nets[net_no]
            the_net.pads = new_pads_indices
            the_net.gates = new_gates_indices
            fill_net_sets(the_net)
        end
    end

    (left_nets, left_gates, left_pads)
end

function right_side_containment(left_pos  ::IndexedPosVec,
                                right_pos ::IndexedPosVec,
                                nets      ::Vector{Net},
                                gates     ::Vector{Gate},
                                pads      ::Vector{Pad})
    num_gates = length(gates)
    right_gates = [deepcopy(gates[right_pos[i][1]]) for i = 1:length(right_pos)]
    right_pads = deepcopy(pads)
    println("-D- right_gates : $right_gates")
    left_gates_indices  = IntSet([left_pos[i][1] for i = 1:length(left_pos)])
    println("-D- left gates indices : $left_gates_indices")
    right_gates_indices = IntSet([right_pos[i][1] for i = 1:length(right_pos)])
    println("-D- right gates indices : $right_gates_indices")
    right_gates_renumbering = Array(Int, num_gates)
    fill!(right_gates_renumbering, -1)
    for i = 1:length(right_gates)
        right_gates_renumbering[right_gates[i].id] = i
    end
    println("-D- right_gates_renumbering : $right_gates_renumbering")

    const median::Float64 = 50.0

    num_nets = length(nets)
    right_nets = [Net() for i = 1:num_nets]
    for i = 1:num_nets
        orig_net = nets[i]
        if isempty(intersect(orig_net.gates_set,right_gates_indices))
            println("-D- the orig_net $orig_net DOESN'T belongs to right partition, EXCLUDE it")
        else
            net_no::UInt = i
            println("-D- the orig_net $orig_net belongs to right partition, fix it and include")
            orig_pads_indices = orig_net.pads
            new_pads_indices = UInt[]
            # process orig pads
            for orig_pad_idx in orig_pads_indices
                orig_pad = pads[orig_pad_idx]
                # fix abscissa if needed
                if orig_pad.pos.x >= median # ok, add this pad_idx AS IS
                    push!(new_pads_indices, orig_pad_idx)
                else
                    new_pad = deepcopy(orig_pad)
                    new_pad.pos.x = median
                    push!(right_pads, new_pad)
                    new_pad_idx::UInt = length(right_pads)
                    push!(new_pads_indices, new_pad_idx)
                end
            end
            orig_gates_indices = orig_net.gates
            new_gates_indices = UInt[]
            # process orig gates
            for orig_gate_idx in orig_gates_indices
                if orig_gate_idx in right_gates_indices # Bingo!
                    @assert(right_gates_renumbering[orig_gate_idx] > 0)
                    remapping::UInt = right_gates_renumbering[orig_gate_idx]
                    push!(new_gates_indices, remapping)
                else # gate belongs to right partition, convert it to pad
                    orig_gate = gates[orig_gate_idx]
                    new_pad = Pad(Pos(median, orig_gate.pos.y), net_no)
                    push!(right_pads, new_pad)
                    new_pad_idx::UInt = length(right_pads)
                    push!(new_pads_indices, new_pad_idx)
                end
            end
            the_net = right_nets[net_no]
            the_net.pads = new_pads_indices
            the_net.gates = new_gates_indices
            fill_net_sets(the_net)
        end
    end

    (right_nets, right_gates, right_pads)

end

function main()

    # for x in ARGS
    #   println(x)
    # end

    @assert(length(ARGS) >= 1)
    inp_name = ARGS[1]
    println("-I- reading from $inp_name")

    (num_gates::UInt, num_nets::UInt, num_pads::UInt, gates, pads) = read_input(inp_name)
    println("-I- num_gates: $num_gates num_nets: $num_nets num_pads: $num_pads")
    println("-D- gates: $gates")
    println("-D- pads: $pads")

    nets = create_nets(num_nets, gates, pads)
    println("-D- nets: $nets")
    gate_pos = qp_step(nets, gates, pads)
    println("-D- 1st plc : \n$gate_pos")

    # update gates positions after 1st placement
    for i = 1:num_gates
        gates[i].pos = gate_pos[i][2]
    end
    println("-D- all gates after 1st placement: \n$gates")

    (left,right) = assignment_cut(gate_pos)
    @assert(length(left) + length(right) == num_gates)
    @printf("-D- nam_gates %d length of left part %d length of right part %d\n",
            num_gates, length(left), length(right))
    println("-D- left part $left")
    println("-D- right part $right")

    (left_nets::Vector{Net}, left_gates::Vector{Gate}, left_pads::Vector{Pad}) = left_side_containment(left, right, nets, gates, pads)
    println("-D- left_gates: $left_gates")
    println("-D- left_pads: $left_pads")
    #left_nets = create_nets(new_nets_no, left_gates, left_pads)
    println("-D- left_nets: $left_nets")
    gate_pos_2 = qp_step(left_nets, left_gates, left_pads)
    println("-D- 2nd plc : \n$gate_pos_2")

    # update gates positions after 2nd placement
    for i = 1:length(left_gates)
        left_gate = left_gates[i]
        global_gates_idx = left_gate.id
        gates[global_gates_idx].pos = gate_pos_2[i][2]
    end

    println("-D- all gates after 2nd placement: \n$gates")

    (right_nets::Vector{Net}, right_gates::Vector{Gate}, right_pads::Vector{Pad}) = right_side_containment(left, right, nets, gates, pads)
    println("-D- right_gates: $right_gates")
    println("-D- right_pads: $right_pads")
    println("-D- right_nets: $right_nets")
    gate_pos_3 = qp_step(right_nets, right_gates, right_pads)
    println("-D- 3d plc : \n$gate_pos_3")

    # update gates positions after 3d placement
    for i = 1:length(right_gates)
        right_gate = right_gates[i]
        global_gates_idx = right_gate.id
        gates[global_gates_idx].pos = gate_pos_3[i][2]
    end

    println("-D- all gates after 3d placement: \n$gates")

    if length(ARGS) >= 2 # write to file with the name given by ARGS[2]
        out_fname::String = ARGS[2]
        os = open(out_fname, "w")
        for i = 1:num_gates
            gate = gates[i]
            pos = gate.pos
            @printf(os, "%d %f %f\n", i, pos.x, pos.y)
        end
        close(os)
    end
end

main()
