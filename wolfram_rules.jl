using GLMakie

# ----- CA Logic -----

function ruletable(rule::Int)
    @assert 0 ≤ rule ≤ 255 "Rule number must be between 0 and 255!"
    return ntuple(i -> Bool((rule >> (i-1)) & 1), 8)
end

function step(row::AbstractVector{Bool}, ruletable::NTuple{8, Bool}; wrap::Bool=true)
    n = length(row)
    next = similar(row)
    @inbounds for i in 1:n
        l = wrap ? row[mod1(i-1, n)] : (i == 1 ? false : row[i-1])
        c = row[i]
        r = wrap ? row[mod1(i+1, n)] : (i == n ? false : row[i+1])
        idx = (l << 2) | (c << 1) | r
        next[i] = ruletable[idx+1]
    end
    return next
end


