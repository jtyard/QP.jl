using Oscar

export fiducial, FIDUCIALS_DIR

FIDUCIALS_DIR = (@__DIR__)*"/fiducials/"

# Returns a fiducial in dimension d - todo: return more like SICdata, rcf, complex conjugation?
function fiducial(d::Int)
    if d == 2
        F = cyclotomic_field(12)[1]
        return (identity_matrix(F,2) + map(F,heis(0,1,2) + heis(1,0,2) + heis(1,1,2))/sqrt(F(3)))/2
    elseif d == 3
        F = cyclotomic_field(3)[1]
        return F[0;1;-1]*F[0 1 -1]/2
    elseif d == 4
        return fiducial("4a")
    elseif d == 5
        return fiducial("5a")
    elseif d == 7
        return fiducial("7b")
    else
        error("Not implemented")
    end
end

# e.g. fiducial("7b") returns the Scott-Grassl 7b fiducial
function fiducial(label::String)
    if label == "4a"
        include(FIDUCIALS_DIR*"4a.jl")
        return Phi
    elseif label == "5a"
        include(FIDUCIALS_DIR*"5a.jl")
        return Phi
    elseif    label == "7a"
        include(FIDUCIALS_DIR*"7ab.jl")
        return Psia
    elseif label == "7b"
        include(FIDUCIALS_DIR*"7ab.jl")
        return Psib
    else
        error("Not implemented")
    end
end

