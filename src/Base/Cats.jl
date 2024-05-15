# Experimental attempts at working with additive categories 
#
# Preadditive category = Ab category.  Then (finite) product = coproduct = biproduct ⊕.
# Additive category = preadditive category with ⊕ defined everywhere.
# Preabelian category = additive category with all kernels. 
# Abelian category = preabelian category with all monics and epis normal.

# Here are basic types that can be instantiated to compute in a preadditive category 

#
# Can I use Catlab?
# i.e. 
# using Catlab
# @instance Catlab.Theories.ThMonoidalCategory
# @instance Catlab.Theories.ThAdditiveCategory

# For now get it working with oplus and otimes
#using Oscar
#import Oscar.⊗
#import Oscar.⊕

#export QuantumSystem

#abstract type QuantumSystem end

#using Catlab

mutable struct Qudit 
    d
    basis
    operator_basis
end

mutable struct Map{T}
    dom::T
    codom::T
end

#@instance Catlab.Theories.ThAdditiveCategory{Qudit,QuditMap}

mutable struct Morphism end

mutable struct Object end

mutable struct Sum{T}
    o::Tuple{T}
end



#mutable struct Tensor{T} <: T
#    o::Vector{T}
#end    


