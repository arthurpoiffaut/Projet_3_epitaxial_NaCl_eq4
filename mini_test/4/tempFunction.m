function output = tempFunction(q)
    Tmm=2250;
    [Tmax,Err,ci,T] = findTmax(q);
    output = Tmax - Tmm;
end