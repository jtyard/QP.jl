// Computing number of stabilizer codes of a given rank
nmax := 10;
size := func<n|2^(n^2)*&*[4^m-1 : m in [1..n]]>;

function numcodes(n,k)
    return size(n)/(size(k)*size(n-k));
end function;

for n in [2..nmax] do
    nums :=  [numcodes(n,k) : k in [1..Floor(n/2)]];
    print n, nums;
end for;