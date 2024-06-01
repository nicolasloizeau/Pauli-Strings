# "implementation of ref https://arxiv.org/pdf/1812.08657.pdf using pauli strings"


module pauli_lanczos

include("pauli_strings.jl")
import ..pauli_strings as ps
using ProgressBars


"in ref https://arxiv.org/pdf/1812.08657.pdf the norm is defined with a /Tr(1) (above eq 3)"
function norm2(H, N)
    return ps.opnorm(H)/sqrt(2^N)
end


"
H : hamiltonian MPO
O : operator MPO
nmax : numer of lanczos steps
maxl : maximum pauli string length. Used by truncate at every step
epsilon : cutoff value. At every step a cutoff is applied to the pauli strings coeficients
https://journals.aps.org/prx/pdf/10.1103/PhysRevX.9.041017 equation 4
"
function lanczos(H, O, nmax, maxl, epsilon)
    N = H.N
    O0 = deepcopy(O)
    b = norm2(ps.com(H, O0), N)
    O1 = ps.com(H, O0)/b
    bs = [b]
    for n in ProgressBar(0:nmax)
        LHO = ps.com(H, O1)
        A = LHO-b*O0
        b = norm2(A, N)
        O = A/b
        O = ps.truncate(O, maxl)
        O = ps.cutoff(O, epsilon)
        O0 = deepcopy(O1)
        O1 = deepcopy(O)
        push!(bs, b)
    end
    return bs
end

end
