function pass = test_all( )

pass = [];
pass = [ pass test_poisson()    ];
pass = [ pass test_transforms() ];

if ( all(pass) )
    pass = 1;
end

end
