using SparseArrays
using Convex
using SCS

# CHSH

GId = sparse([1,2,3,4,5], [1,2,3,4,5], [1,1,1,1,1])
GA1 = sparse([1], [2,1], [1,1], 5, 5)
GA2 = sparse([1], [3,1], [1,1], 5, 5)
GB1 = sparse([1], [4,1], [1,1], 5, 5)
GB2 = sparse([1], [5,1], [1,1], 5, 5)
GA1A2 = sparse([2], [3,2], [1,1], 5, 5)
GA1B1 = sparse([2], [4,2], [1,1], 5, 5)
GA1B2 = sparse([2], [5,2], [1,1], 5, 5)
GA2B1 = sparse([3], [4,3], [1,1], 5, 5)
GA2B2 = sparse([3], [5,3], [1,1], 5, 5)
GB1B2 = sparse([4], [5,4], [1,1], 5, 5)

A1 = Variable()
A2 = Variable()
B1 = Variable()
B2 = Variable()
A1A2 = Variable()
A1B1 = Variable()
A1B2 = Variable()
A2B1 = Variable()
A2B2 = Variable()
B1B2 = Variable()

objective = A1B1 + A1B2 + A2B1 - A2B2

Gamma = (GId + GA1*A1 + GA2*A2 + GB1*B1 + GB2*B2
         + GA1A2*A1A2 + GA1B1*A1B1 + GA1B2*A1B2
         + GA2B1*A2B1 + GA2B2*A2B2 + GB1B2*B1B2)

constraints = [(Gamma in :SDP)]

problem = maximize(objective, constraints)

solve!(problem, SCS.Optimizer)

println(evaluate(objective))
