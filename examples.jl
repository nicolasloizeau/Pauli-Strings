
include("pauli_strings.jl")
import .pauli_strings as ps

# initiatialize an empty operator of N=3 qubits
o = ps.Operator(3)

# add a pauli strings to the operator
o += "XYZ"

# add a pauli string with a coeficient (o += 2XYZ)
o += 2, "XYZ"
o += 2, [1,0,3]


# add a 2-local term
# a term is of the form (J, a, i, b, j), a and b are sites operator types, i, and j are sites
# the following line means o+= 2 Z_2 Z_3
o += 2, "Z", 2, "Z", 3
o += 4, "S+", 1, "S-", 2

# no J is interpreted as J=1
o += "Y", 2, "Y", 3


# add a 1-local term
# add a Z operator at site 2
o += 1, "Z", 2

# multiply two operators
o1 = ps.Operator(3)
o2 = ps.Operator(3)
o1 += "XYZ"
o2 += "XXX"
o3 = o1*o2

# add and subtract
o3 = o1+o2
o3 = o1-o2

# adding a scalar is equivalent to adding the unit times the scalar
o3 = o1+2

# multiply operator by a scalar
o = 5*o

# trace:
println(ps.trace(o))

# frobenius norm:
println(ps.opnorm(o))

# conjugate transpose:
println(ps.dagger(o))

# number of terms in o:
println(length(o))

#convert an operator to a list of coeficients and pauli strings :
coefs, strings = ps.op_to_strings(o)

#cutoff: remove terms with coef less that cutoff. Return a copy
println(o)
println(ps.cutoff(o, 6))

#truncate: remove terms with pauli lengh less that truncate. Return a copy
println(ps.truncate(o, 2))
