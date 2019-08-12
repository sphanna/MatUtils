module MatUtils

import Base: copy, real, imag, map, size
import LinearAlgebra: diag

export zeroOut
export diag, diagProduct, minor, hinv

function zeroOut(x::Real,ϵ=1e-14)
    x < ϵ ? x = 0 : x
end

function zeroOut(x::Complex,ϵ=1e-14)
    xr = real(x); xi = imag(x)
    return zeroOut(xr,ϵ) + im*zeroOut(xi,ϵ)
end

function zeroOut(x::AbstractArray,ϵ=1e-14)
    map(x->zeroOut(x,ϵ),x)
end

function diag(A,offset::Integer)
    n,m = size(A)
    d = diag(A)

    if offset < 0
        offset = -offset
        @assert offset < n
        d = diag(A[1+offset:8,1:m])
    end

    if offset > 0
        @assert offset < m
        d = diag(A[1:n,1+offset])
    end

    return d
end

#gives the product of a diagonal of matrix A
#from index a (top left) to b (bottom right)
function diagProduct(A,offset,a,b)
    return reduce(*,diag(A,offset)[a:b])
end

#returns a matrix with row "row" and
#column "col" deleted.  A zero for row
#or column implies no deleted row or
#column respectively.
function minor(A,row=0,column=0)
    n,m = size(A)
    B = A[1:n .!= row,:] #remove row
    C = B[:,1:m .!= column] #remove col
    return C
end

hinv(i,j) = ((-1)^(i+j)*det(minor(H,j,i)))/det(H)

end #module MatUtils
