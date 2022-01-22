struct Leaf
    clause::Vector{Int}
    can_split::Bool
    captured::BitArray
    num_predicted::Float64
    num_captured::Int
    function Leaf(clause::Vector{Int64}, can_split::Bool, captured::BitVector, y::BitVector)
        s1 = sum(captured .* y)
        s2 = sum(captured)
        new(clause, can_split, captured, max(s1, s2 - s1), s2)
    end
    function Leaf(l::Leaf, b::Bool)
        new(l.clause, b, l.captured, l.num_predicted, l.num_captured)
    end
end
struct Tree
    leaves::Vector{Leaf}
    node::Int
    miss::Float64
    penalty::Float64
end


function sorttree(clauses::Vector{Vector{Int}})
    if isempty(clauses[1])
        return (clauses)
    end
    most = map(clauses) do x
        abs.(x)
    end |> t -> intersect(t...)
    i = most[1]
    function removeele(leaves, i)
        map(t -> setdiff(t, [i]), filter(x -> i in x, leaves))
    end
    union([[i, x...] for x in sorttree(removeele(clauses, i))], [[-i, x...] for x in sorttree(removeele(clauses, -i))])
end

function Base.print(t::Tree)
    print("Node:", node(t.leaves), " Miss:", t.miss, "\n\r")
    if (isempty(t.leaves))
        print("  R\n\r")
    else
        t = map(x -> x.clause, t.leaves)
        t = sorttree(t)
        sort!(t, rev = true)
        h = max(map(x -> length(x), t)...)
        width = 2
        print(repeat(" ", width * (2^(h - 0) - 1)) * " R\n\r")
        for i in 1:h
            completenode = map(x -> length(x) >= i ? [string(x[i])] : repeat(["  "], 2^(i - length(x))), t)
            tt = reduce((x, y) -> begin
                    (" " in y[1]) || y[1] != y[2] ? append!(x, [y[1]]) : x
                end, zip(completenode[2:length(completenode)], completenode[1:length(completenode)-1]), init = [completenode[1]])
            strings = map(x -> length(x) < width ? repeat(" ", width - length(x)) * x : x, Iterators.Flatten(tt))
            print(repeat(" ", width * (2^(h - i) - 1)))
            for s in strings
                print(s)
                print(repeat(" ", width * (2^(h - i + 1) - 1)))
            end
            print("\n\r")
        end
    end
    nothing
end

function node(leaves::Vector{Leaf})
    if (isempty(leaves))
        return (0)
    end
    t = map(x -> x.clause, leaves)
    sort!(t, rev = true)
    h = max(map(x -> length(x), t)...)
    n = 0
    for i in 1:h
        completenode = map(t -> t[i], filter(x -> length(x) >= i, t))
        tt = reduce((x, y) -> begin
                y[1] != y[2] ? append!(x, [y[1]]) : x
            end, zip(completenode[2:length(completenode)], completenode[1:length(completenode)-1]), init = [completenode[1]])
        n += length(tt)
    end
    n
end