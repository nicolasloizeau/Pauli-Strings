

module models
include("pauli_strings.jl")
import ..pauli_strings as ps
using Random

"""XX hamiltonian https://arxiv.org/abs/1812.08657"""
function XX(N)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += ('X',j,'X',j+1)
        H += ('Z',j,'Z',j+1)
    end
    return H
end

"""X operator https://arxiv.org/abs/1812.08657"""
function X(N)
    H = ps.Operator(N)
    for j in 1:N
        H += ('X',j)
    end
    return H
end

"""return an operator corresponding to a majorana encoded as a pauli string"""
function majorana(N::Int, k::Int)
    s = fill(0, N)
    if k%2 == 0
        s[div(k, 2)+1] = 1
    else
        s[div(k, 2)+1] = 2
    end
    s[1:div(k, 2)] .+= 3
    H = ps.Operator(N)
    return H+s
end

"""return the 2N majoranas supported on N qubits"""
function majoranas(N::Int)
    ms = []
    for k in 0:N*2-1
        push!(ms, majorana(N,k))
    end
    return ms
end

"""indices for the syk summation i<j<k<l"""
function syk_indexes(N)
    ind = []
    for i in 1:N*2
        for j in 1:i-1
            for k in 1:j-1
                for l in 1:k-1
                    push!(ind, (i,j,k,l))
                end
            end
        end
    end
    return ind
end

"""generate a random SYK hamiltonian supported on N qubits with J random normal(0,1)"""
function sykH(N::Int)
    majs = majoranas(N)
    H = ps.Operator(N)
    for (i, j, k, l) in syk_indexes(N)
        H += randn()*majs[i]*majs[j]*majs[k]*majs[l]
    end
    return H
end


"""hamiltonian from fig 4b of https://arxiv.org/pdf/1812.08657"""
function chaotic_chain(N::Int)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += ('X',j,'X',j+1)
    end
    for j in 1:N
        H += (-1.05,'Z',j)
        H += (0.5,'X',j)
    end
    return H
end

"""operator from fig 4b of https://arxiv.org/pdf/1812.08657"""
function chaotic_chain_op(N::Int)
    H = ps.Operator(N)
    for j in 1:(N - 1)
        H += (-1.05,'X',j,'X',j+1)
    end
    for j in 1:N
        H += ('Z',j)
    end
    return H
end



"""hamiltonian from eq 20 of https://arxiv.org/pdf/1812.08657"""
function XZ2D(N::Int)
    if sqrt(N) != floor(sqrt(N))
        throw(ArgumentError("N is not the square of an integer"))
    end
    L::Int = sqrt(N)
    H = ps.Operator(N)
    for x in 1:L-1
        for y in 1:L-1
            # convert x,y in qubit index
            i = L*(y-1)+x
            j = L*(y-1)+(x+1)
            H += ('X',i,'Z',j)
            i = L*(y-1)+x
            j = L*y+x
            H += ('Z',i,'X',j)
        end
    end
    return H
end

"""operator for use with XZ2D"""
function X00(N)
    O = ps.Operator(N)
    O += ('X',1)
    return O
end

end
