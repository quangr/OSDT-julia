include("OSDT.jl")
using .OSDT:bbound!,Fitter,Leaf
X=[1 1 0;1 0 1;0 1 0]
y=[0;1;1]
f=Fitter(X,y,0.05)
bbound!(f)