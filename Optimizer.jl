module Optimizer
    using ..OSDT: Tree,Leaf,node
    using Test
    mutable struct Fitter
        X::BitMatrix
        y::BitArray
        lamb::Float64
        queue::Vector{Tree}
        B_c::Float64
        numdata::Int
        BestTree::Tree
        LeafCache::Dict{Tuple{Vector{Int},Bool},Leaf}
        TreeCache::Dict{Vector{Tuple{Vector{Int},Bool}},Bool}
        function Fitter(X::BitMatrix,y::BitArray,lamb::Float64)
            temp=new(X,y,lamb,Tree[],1,size(X)[1])
            temp.BestTree=Tree([Leaf(Int[],true,temp)],temp,0)
            temp.LeafCache=Dict()
            temp.TreeCache=Dict()
            push!(temp.queue,temp.BestTree)
            temp
        end
    end

    function get_capture(clause::Vector{Int},fitter::Fitter)
        res=trues(fitter.numdata)
        for i in clause
            if i<0
                res.*=view(fitter.X, :, -i).==0
            else
                res.*=view(fitter.X, :, i).==1
            end
        end
        res
    end
    function Leaf(clause::Vector{Int64},can_split::Bool,fitter::Fitter)
        Leaf(clause,can_split,get_capture(clause,fitter),fitter.y)
    end
    function Leaf(clausepair::Tuple{Vector{Int64}, Bool},fitter::Fitter)
        (clause,can_split)=clausepair
        clausepair=(sort!(clause),can_split)
        if haskey(fitter.LeafCache,clausepair)
            fitter.LeafCache[clausepair]
        else
            l=Leaf(clause,can_split,fitter)
            fitter.LeafCache[clausepair]=l
            l
        end
    end
    function pushtofitter(fitter::Fitter,tree::Tree)
        sorttree=map(x->(x.clause,x.can_split),tree.leaves)
        sort!(sorttree)
        if !haskey(fitter.TreeCache,sorttree)
            fitter.TreeCache[sorttree]=true
            push!(fitter.queue,tree)
        end
    end
    function get_splitable_leaves!(fitter::Fitter,t::Tree,nrule::Int)
        num_tosplit=length(filter(x->x.can_split,t.leaves))
        if num_tosplit!=0
            function leave2list(x::Leaf)
                tt=(x.can_split ? [[([x.clause...,i],true),([x.clause...,-i],true,)] for i in setdiff(1:nrule,abs.(x.clause))] : [[(x.clause,false)]])
                tt
            end
            combo=map(leave2list,t.leaves)
            trees=map(x->collect(Iterators.flatten(x)),Iterators.product(combo...))
            for tree in trees
                treecombos=map(tree)do (x,y)
                    if y
                        length(x)==nrule ? [(x,false)] : [(x,true),(x,false)]
                    else
                        [(x,y)]
                    end 
                end
                for treecombo in Iterators.product(treecombos...)
                    x=map(treecombo)do clpair
                        Leaf(clpair,fitter)                       
                    end
                    inserttree=Tree([x...],fitter,t.node+num_tosplit)
                    pushtofitter(fitter,inserttree)
                end
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
        n=node(leaves)
        ac=map(leaves) do x
            x.num_predicted
        end|>sum|>t->t/fitter.numdata
        Tree(leaves,n,1-ac,1-ac+fitter.lamb*n)
    end

    function Tree(leaves::Vector{Leaf},fitter::Fitter,node::Int)
        ac=map(leaves) do x
            x.num_predicted
        end|>sum|>t->t/fitter.numdata
        Tree(leaves,node,1-ac,1-ac+fitter.lamb*node)
    end

    export Fitter,bbound!

    # @test get_capture([-1,2],Fitter([1 0],[1],0.2))==[0]
    # @test get_capture([-1,2],Fitter([1 0;0 1;],[1,0],0.2))==[0,1]
end

