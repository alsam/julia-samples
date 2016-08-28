using ApproxFun

x = Fun(identity, [0,100])

g = sqrt(x)

g2 = g(2)

x = Fun(identity, [0.1,1.])

w = Fun(identity,[0.1,1.])

T = Fun(identity,[0.1,1.])

ν = Fun(identity,[0.1,1.])

#kern = ν^4 * besselj(1, 2π*√(x*x+w*w/4) )^2

kern = Fun( x -> (ν^4 * besselj(1, 2π*√(x*x+w*w/4) )^2 / (x^2 + w^2/4) - T)^2 )
