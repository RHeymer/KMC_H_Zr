%% function vec2dist converts vector x=[da,db,dc] to a distance in metres

function dist = vec2dist(x)
    a0 = 3.232e-10; %lattice parameter in metres
    
    da = x(1);
    db = x(2);
    dc = x(3);
    
    dx = a0*(da+(0.5*db));
    dy = a0*(sqrt(3)/2)*db;
    dz = a0*(2/3)*sqrt(6)*dc;
    
    dist = sqrt(dx.^2+dy.^2+dz.^2);
end