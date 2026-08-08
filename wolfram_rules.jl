using GLMakie

# ----- CA Logic -----

function ruletable(rule::Int)
    @assert 0 ≤ rule ≤ 255 "Rule number must be between 0 and 255!"
    return ntuple(i -> Bool((rule >> (i-1)) & 1), 8)
end

function step!(next::AbstractVector{Bool}, curr::AbstractVector{Bool}, table::NTuple{8, Bool}; 
                wrap::Bool=true)
    n = length(curr)
    @inbounds for i in 1:n
        l = wrap ? curr[mod1(i-1, n)] : (i == 1 ? false : curr[i-1])
        c = curr[i]
        r = wrap ? curr[mod1(i+1, n)] : (i == n ? false : curr[i+1])
        idx = (l << 2) | (c << 1) | r
        next[i] = table[idx+1]
    end
    return next
end