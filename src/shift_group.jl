##
## OBSOLETE
##

export generators, addgen, two_products, make_generators,
    verify_generated

function invel(g)
    return (6 - g[1], 6 - g[2])
end

# Need to include identity (0, 0) because
# This includes the generators themselves
# in the product of pairs of them.
function make_generators()
    _generators = (
        (0,0),
        (1,1),
        (3,5),
        (1,3),
        (3,4),
        (2,3)
    )
    return tuple(_generators..., invel.(_generators)...,)
end

const generators = make_generators()

function addgen(g1, g2)
    (a1, b1) = g1
    (a2, b2) = g2
    return ((a1 + a2) % 6, (b1 + b2) % 6)
end

function two_products()
    twoels = [addgen(x, y) for x in generators for y in generators];
    return unique(sort(twoels))
end

function verify_generated()
    @assert length(two_products()) == 35
end
