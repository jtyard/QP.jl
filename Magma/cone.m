// assumes d is odd

f0 := f(0,0);
light := [f(i,i) : i in [1..(d-1)/2]] cat [f(i,-i) : i in [1..(d-1)/2]];
forwp := [fp(i,j) : j in [-i+1..i-1], i in [0..(d-1)/2]];
forwm := [fm(i,j) : j in [-i+1..i-1], i in [0..(d-1)/2]];

Ilight := ideal<R|light>;
Iforwp := ideal<R|forwp>;
Iforwm := ideal<R|forwm>;
Icone := Ilight + Iforwp;


monomials := MonomialsOfDegree(R,2);
N := #monomials;

R2 := VectorSpace(Q,N);
function Vf(f)
    return R2 ! [MonomialCoefficient(f,m) : m in monomials];
end function;


function Vfspan(fs)
    return sub<R2|[Vf(f) : f in fs]>;
end function;

function Vdim(fs)
    return Dimension(Vfspan(fs));
end function;

RP0 := ideal<R|f(1,1), f(2,2), f(3,3), f(1,6),f(2,5),f(1,4),f(2,1) + f(3,-2) + f(3,-1), f(2,-1)+f(3,2)+f(3,1)>;

RHT := ideal<R|f(1,1) + f(2,2) + f(3,3), f(1,6) + f(2,5) + f(3,4), f(0,1) + f(0,2) + f(0,3), f(2,1)+f(3,-2)+f(3,-3),f(2,-1)+f(3,1)+f(3,2)>;

t2 := tt(2);
tV := ZeroMatrix(Q,N,N);

