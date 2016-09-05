using Calculus

f(k) = quadgk(η -> 1 / √(1 - k^2 * cos(η)^2), 0, π / 2 )[1]
g(k) = quadgk(η -> 1 / √(1 - k^2 * sin(η)^2), 0, π / 2 )[1]
δ(x) = f(x) - g(x)


for x = 0.01 : 0.02 : 0.99
    del = δ(x)
    println("δ($x) = $del")
end

