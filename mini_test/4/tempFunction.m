function output = tempFunction(q)
    Tmm = 2000;
    [Tmax,Err,ci] = findTmax(q)
    output = Tmax - Tmm;
end