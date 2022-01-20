module Optimizer
    using ..OSDT: Tree,Leaf,node
    using Test
    mutable struct Fitter
        X::Matrix{Int64}
        y::Array{Int64}
        lamb::Float64
        queue::Vector{Tree}
        B_c::Float64
        numdata::Int
        BestTree::Tree
        function Fitter(X::Matrix{Int64},y::Array{Int64},lamb::Float64)
            temp=new(X,y,lamb,Tree[],1,size(X)[1])
            temp.BestTree=Tree([Leaf(Int[],true,repeat([0],temp.numdata),y)],temp)
            push!(temp.queue,temp.BestTree)
            temp
        end
    end

    function get_capture(clause::Vector{Int},fitter::Fitter)
        res=repeat([1],fitter.numdata)
        for i in clause
            if i<0
                res.*=-view(fitter.X, :, -i).+1
            else
                res.*=view(fitter.X, :, i)
            end
        end
        res
    end
    function Leaf(clause::Vector{Int64},can_split::Bool,fitter::Fitter)
        Leaf(clause,can_split,get_capture(clause,fitter),fitter.y)
    end
    function get_splitable_leaves!(fitter::Fitter,t::Tree,nrule::Int)
        if !isempty(filter(x->x.can_split,t.leaves))
            function leave2list(x::Leaf)
                tt=(x.can_split ? [[Leaf([x.clause...,i],true,fitter),Leaf([x.clause...,-i],true,fitter)] for i in setdiff(1:nrule,abs.(x.clause))] : [[x]])
                tt
            end
            leave2list(t.leaves[1])
            combo=map(leave2list,t.leaves)
            trees=map(x->Tree(collect(Iterators.flatten(x)),fitter),Iterators.product(combo...))
            for tree in trees
                treecombo=map(x::Leaf->x.can_split ? [Leaf(x,true),Leaf(x,false)] : [x],tree.leaves)
                inserttrees=map(x->Tree([x...],fitter),collect(Iterators.product(treecombo...)))
                push!(fitter.queue,inserttrees...)
            end
        else
        end
    end


    function bbound!(fitter::Fitter)
        X_train=fitter.X
        y_train=fitter.y
        nrule=size(X_train)[2]
        queue=fitter.queue
        while size(queue)[1]>0
            tree=pop!(queue)
            if tree.penalty<fitter.B_c
                fitter.B_c=tree.penalty
                fitter.BestTree=tree
            end
            get_splitable_leaves!(fitter,tree,nrule)
        end
        print(fitter.BestTree)
    end

    function Tree(leaves::Vector{Leaf},fitter::Fitter)
        ac=map(leaves) do x
            x.num_predicted
        end|>sum|>t->t/fitter.numdata
        Tree(leaves,1-ac,1-ac+fitter.lamb*node(leaves))
    end

    export Fitter,bbound!

    @test get_capture([-1,2],Fitter([1 0],[1],0.2))==[0]
    @test get_capture([-1,2],Fitter([1 0;0 1;],[1,0],0.2))==[0,1]
end

