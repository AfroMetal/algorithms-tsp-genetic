type Tour
    route::Array{Int16}
    distance::Float32

    function Tour(r::Array{Int16})
        if length(r) < 1
            error("route has no cities")
        end
        new(r, distance([Int16(1) ; r ; Int16(1)]))
    end
end

import Base.*
function Base.:*(tour1::Tour, tour2::Tour)
    i = rand(1:n-2)
    j = rand(i+1:n-1)

    route = zeros(Int16, n-1)
    # copy middle part
    route[i:j] = tour1.route[i:j]

    # fill other places
    spare = 1
    for city in tour2.route
        if !(city in route)
            for idx in spare:n-1
                if route[idx] == 0
                    route[idx] = city
                    spare = idx
                    break
                end
            end
        end
    end

    if 0.5 > rand()
        Tour(route)
    else
        Tour(mutate(route, 0.01))
    end
end

function distance(solution::Array)
    cost = Float32(0)
    for i in 1:length(solution)-1
        cost += d[solution[i],solution[i+1]]
    end
    return Float32(cost)
end

# function aproximateMST()
#     T = Array{Tuple{Int16, Int16}}(0)
#     X = [1]
#
#     while length(X) != n
#         crossing = Array{Tuple{Int16, Int16}}(0)
#         for x in X
#             for y in 1:n
#                 if !(y in X) && d[x,y] != 0
#                     crossing = [crossing ; (Int16(x), Int16(y))]
#                 end
#             end
#         end
#
#         edge = sort(crossing, by=e->d[e[1],e[2]])[1]
#         T = [T ; edge]
#         X = [X ; edge[2]]
#     end
#
#     visited = Int16[]
#
#     function dfs(tree::Array{Tuple{Int16, Int16}}, v::Int16)
#         visited = [visited ; v]
#         for edge in [e for e in tree if v in e]
#             w = edge[1]
#             if edge[1] == v
#                 w = edge[2]
#             end
#             if !(w in visited)
#                 dfs(tree, w)
#             end
#         end
#     end
#
#     dfs(T,Int16(1))
#     return visited[2:n]
# end

function aproximategreedy(from::Int16)
    city = from
    visited = Dict(zip(collect(Int16, 1:n), falses(n)))
    visited[1] = true
    visited[city] = true
    route = unique(Int16[1 ; city])
    while length(route) != n
        city = sort([c for c in 2:n if !visited[c]], by=x->d[city,x])[1]
        visited[city] = true
        route = Int16[route ; city]
    end
    return route[2:n]
end

function hillclimbing(route::Array{Int16}, rounds::Int)
    route = [Int16(1) ; route ; Int16(1)]
    for _ in 1:rounds
        for idx1 in 2:n-3
            idx2 = rand(idx1+2:n-1)
            if d[route[idx1-1], route[idx1]] -
            d[route[idx1], route[idx1+1]] +
            d[route[idx1-1], route[idx2]] +
            d[route[idx2], route[idx1+1]] -
            d[route[idx2-1], route[idx2]] -
            d[route[idx2], route[idx2+1]] +
            d[route[idx2-1], route[idx1]] +
            d[route[idx1], route[idx2+1]] < 0
                route[idx1], route[idx2] = route[idx2], route[idx1]
            end
        end
    end
    return route[2:n]
end

function mutate(route::Array{Int16}, prob::Float64)
    for i in eachindex(route)
        if prob > rand()
            idx = rand(1:length(route))
            route[i], route[idx] = route[idx], route[i]
        end
    end
    return route
end

tic()
tic()

const start = time()
const DEBUG = length(ARGS) >= 1 && ARGS[1] == "1"

# read data size N
const n = parse(Int32, readline(STDIN))

# 2xN node's coordinates matrix
data = zeros(Float32, n, 2)

# read file data to data matrix
for i in 1:n
    v = readline(STDIN)
    (i_v, x_v, y_v) = split(v)
    data[parse(Int16, i_v),1] = parse(Float32, x_v)
    data[parse(Int16, i_v),2] = parse(Float32, y_v)
end

const maxtime = parse(Int64, readline(STDIN))

# calculate nodes distances as NxN matrix
d = zeros(Float32, n, n)
for i = 1:n
    for j = i+1:n
        d[i,j] = norm(data[i,1:2] - data[j,1:2])
        d[j,i] = d[i,j]
    end
end
# free data memory
data = 0
gc()

if n >= 20
    const k = 2*Int16(ceil(sqrt(n)))
    const maxiters = Int(floor(3.0e7*round(sqrt(n))/n^2))
else
    const k = 30
    const maxiters = 100
end

aproxsol = aproximategreedy(Int16(1))

population = [Tour(mutate(hillclimbing(aproxsol, Int(ceil(maxiters/2))), 0.01)) for _ in 1:k]
population[1] = Tour(aproxsol)
aproxsol = 0
sort!(population, by=x->x.distance)
gc()

newpopulation = Array{Tour}(k)
bestpopulation = copy(population)

noimprov = 0
totalnoimprov = 0
iters = 0
bestsol = population[1]

const initialcost = bestsol.distance

if DEBUG
    toc()
end

function endprogram()
    if DEBUG
        @printf("initial cost:\t%.3f\n", initialcost)
    	@printf("best cost:\t%.3f\n", bestsol.distance)
        @printf("%.2f%% improvement\n", 100.0-100*bestsol.distance/initialcost)
    	@printf("maxiters:\t%d\n", maxiters)
    	@printf("iterations:\t%d\n", iters)
    	toc()
    else
        @printf("%.5f\n", bestsol.distance)
        for city in [Int16(1) ; bestsol.route ; Int16(1)]
            print(STDERR, city)
            print(STDERR, " ")
        end
        println(STDERR)
    end
end
atexit(endprogram)

function checktime()
    if time() - start > maxtime - 1.0
        exit()
    end
end

while iters < maxiters
    checktime()
    iters += 1

    for i in 1:k
        pop1 = zero(Int16)
        pop2 = zero(Int16)

        while pop1 == 0
            pop1 = findfirst(x -> max(0.05, 0.3-x/k) > rand(), eachindex(population))
        end
        while pop2 == 0 || pop2 == pop1
            pop2 = findfirst(x -> max(0.05, 0.3-x/k) > rand(), eachindex(population))
        end
        # println((pop1,pop2))

        ofs1 = population[pop1] * population[pop2]
        checktime()
        ofs2 = population[pop2] * population[pop1]

        if ofs1.distance <= ofs2.distance
            newpopulation[i] = ofs1
        else
            newpopulation[i] = ofs2
        end
        checktime()
    end

    population = sort(newpopulation, by=x->x.distance)

    if bestsol.distance > population[1].distance
        bestsol = population[1]
        checktime()
        # if DEBUG
        #     @printf("%.3f\t%d\n", bestsol.distance, iters)
        # end
        bestpopulation = copy(population)
        noimprov = 0
        totalnoimprov = 0
    else
        noimprov += 1
        totalnoimprov += 1
        if noimprov > max(20, maxiters / 20)
            checktime()
            # if DEBUG
            #     @printf("no improvement for %d, mutating...\n", noimprov)
            # end
            population = sort([Tour(hillclimbing(mutate(p.route, 0.01), totalnoimprov)) for p in bestpopulation], by=x->x.distance)
            noimprov = 0
            checktime()
        end
        if totalnoimprov > max(50, maxiters / 5)
            exit()
        end
    end
end
