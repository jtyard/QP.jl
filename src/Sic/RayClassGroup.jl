# It might be nice to have an explicit type for ray class groups that e.g. ideals can be coerced into.  I didn't get far with this though.

using Oscar

mutable struct RayClassGroup
    C  # The group
    ph # C -> ideals 
end

