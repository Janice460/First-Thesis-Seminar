# Packages
import Pkg
Pkg.add("GLPK")
Pkg.add("JuMP")
using GLPK 
using JuMP
Pkg.add("Polyhedra")
Pkg.add("CDDLib")
using Polyhedra, CDDLib
import Printf

# ====================== MODEL =======================
#    min  cx + dy
#    s.a. Ax + By >= e
#              Cx <= f
#             x,y >= 0

# Data
c = [1, 3]
d = [-2, -3]
nvarx = length(c)
nvary = length(d)
e = [-2; -3]
A = [1 -3;
    -1 -3]
B = [1 -2; 
    -1 -1]
C = [-1 0;
    0 -1]
f = [-0.5; -0.5]
M = -100
maxiter = 100
epsi = 1e-6

function MaestroRestricto(M,c, C,f)
    MR = Model(GLPK.Optimizer)
    @variable(MR, x[1:nvarx] >= 0)
    @variable(MR, alpha >= M)
    @objective(MR, Min, c' * x + alpha)
    @constraint(MR, C * x .<= f)
    return MR
end

function SubDual(x)
    SPD = Model(GLPK.Optimizer)
    lb = length(e)
    @variable(SPD, lambda[1:lb] >= 0)
    c2 = e - A * x
    @objective(SPD, Max, c2' * lambda)
    @constraint(SPD, cons, B' * lambda .<= d )
    return SPD
end

function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

function Benders(epsi, maxiter, M, C,f,d)
    println("Iteration  Lower Bound  Upper Bound          Gap")
    MR = MaestroRestricto(M,c, C,f)
    x = MR[:x]
    alpha = MR[:alpha]
    xgen = Matrix(undef, 0, nvarx)
    for k in 1:maxiter
        # Solve restrict master
        optimize!(MR)
        lower_bound = objective_value(MR)
        xk = value.(x)
        xgen = vcat(xgen,reshape(xk,(1,nvarx)))

        # Solve subproblem
        SubD = SubDual(xk)
        lambda = SubD[:lambda]
        optimize!(SubD)
        statusSubp = string(MOI.get(SubD, MOI.TerminationStatus()))

        if statusSubp == "OPTIMAL"
            obj = objective_value(SubD)
            λ = value.(lambda)
            upper_bound = c' * xk + obj
            gap = (upper_bound - lower_bound)
            print_iteration(k, lower_bound, upper_bound, gap)
            if gap < epsi
                println("Terminating with the optimal solution with gap ", gap)
                if gap == 0
                    println("The objective value is ", upper_bound)
                else
                    println("The objective value is near to our numerical results, is more than ", lower_bound, " and less than", upper_bound)
                end
                break
            end
            Alamb = A' * λ
            cut = @constraint(MR, alpha >=  e'*λ - Alamb' * x)
            @info "Adding the optimality cut $(cut)"
        elseif statusSubp == "INFEASIBLE"
            return print("El problema principal es infactible")
            break
        else # unbounded case
            poly = polyhedron(SubD, CDDLib.Library(:exact))
            hrep(poly)
            B2 = rays(vrep(poly))

            # Extreme directions
            extdirec = Matrix(undef, 0, 9)
            for j in B2
                m = []
                for i in 1:9
                    push!(m, vec(j)[i])
                end
                m = reshape(m,(1,54))
                extdirec = vcat(extdirec, m) # almacena las direcciones extremas en una matriz
            end 
            # already have extreme rays 
            nrays = size(extdirec)[1] # number of rays
            MR = MaestroRestricto(M,c, C,f)
            x = MR[:x]
            alpha = MR[:alpha]
            for g in 1:nrays
                ray = extdirec[g,:]
                Aray = A'*vec(ray)
                bray = b' * ray
                cut2 = @constraint(MR, 0 >= Aray' * x - bray )
                @info "Adding the feasibility cut $(cut2)"
            end
        end
    end

    print("x values throught iterations: \n",xgen, "\n")
    print("Optimal value x= ", vec(xgen[end,:]))
end

Benders(epsi, maxiter, M, C,f,d)