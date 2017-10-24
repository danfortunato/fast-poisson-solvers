function pass = test_poisson( )

pass = [];
pass = [ pass test_poisson_rectangle()    ];
pass = [ pass test_poisson_cylinder()     ];
pass = [ pass test_poisson_solid_sphere() ];

if ( all(pass) )
    pass = 1;
end

end
