
# lanczos example from :
# https://journals.aps.org/prx/abstract/10.1103/PhysRevX.9.041017
# figure 2, "X in XX"

include("pauli_strings.jl")
import .pauli_strings as ps
using PyPlot
using ProgressBars


# XX hamiltonian: 1d chain with XX+YY interraction
function XX(N)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += ('X',j,'X',j+1)
        H += ('Z',j,'Z',j+1)
    end
    return H
end

# X local operator: X operator on each site
function X(N)
    H = ps.Operator(N)
    for j in 1:N
        H += ('X',j)
    end
    return H
end


N = 32 # system size
H = XX(N) #hamiltonian
O = X(N) #operator


ioff()


# heisenberg evolution of the operator using rk4
# plot tr(O*O(t))
function evolve(H, O, M, eps, dt)
    ts = range(0, stop=2, step=dt)
    echo = []
    O0 = deepcopy(O)
    for t in ProgressBar(ts)
        push!(echo, ps.trace(O*ps.dagger(O0))/ps.trace(O0*O0))
        O = ps.rk4(H, O, dt; heisenberg=true)
        O = ps.truncate(O, M)
        O = ps.cutoff(O, eps)
    end
    plot(ts, echo, label="cutoff=$eps")
end

plt.cla()
# time evolve O for different cutoff values
for eps in (0.1,0.06,0.02,0.01)
    evolve(H, O, 100, eps, 0.1)
end

legend()
title("N=$N")
xlabel("t")
ylabel("tr(O(0)*O(t))")
savefig("time_evolve_example.png")
show()
