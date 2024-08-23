


@attributes mutable struct MatPolyMap{
    D <: MatPolyRing,
    C <: NCRing,
    U,
    V} <: Map{D, C, Map, MatPolyMap}

  domain::D
  codomain::C
  coeff_map::U
  img_gens::Vector{V}
  temp_ring           # temporary ring used when evaluating maps

  function MatPolyMap{D, C, U, V}(domain::D,
                                codomain::C,
                                coeff_map::U,
                                img_gens::Vector{V}) where {D, C, U, V}
      @assert V === elem_type(C)
      for g in img_gens
        @assert parent(g) === codomain "elements does not have the correct parent"
      end
    return new{D, C, U, V}(domain, codomain, coeff_map, img_gens)
  end
end

function MatPolyMap(d::D, c::C, cm::U, ig::Vector{V}) where {D, C, U, V}
  return MatPolyMap{D, C, U, V}(d, c, cm, ig)
end

################################################################################
#
#  Field access
#
################################################################################

domain(f::MatPolyMap) = f.domain
codomain(f::MatPolyMap) = f.codomain
