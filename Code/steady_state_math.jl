#=
Solves Q2.Part a)
=#

function fem1d_heat_steady(n,a,b,ua,ub,k,f,x)

    quad_num = 2;
    abscissa = (-0.577350269189625764509148780502, +0.577350269189625764509148780502)
    weight = (1.0, 1.0)

    Amat = zeros( n, n );
    bvec = zeros( n, 1 );

    Amat[1,1] = 1.0;
    bvec[1] = ua;

    for i in 2:n-1
    
        xl = x[i-1]
        xm = x[i]
        xr = x[i+1]

        al,am,ar,bm = 0.0,0.0,0.0,0.0

        for q in 1:quad_num
            # Integrate over the LEFT interval, between XL and XM
            xq = ( ( 1.0 - abscissa[q] ) * xl
                 + ( 1.0 + abscissa[q] ) * xm ) / 2.0
            wq = weight[q] * ( xm - xl ) / 2.0

            vl = ( xm - xq ) / ( xm - xl ); 
            vlp =     - 1.0  / ( xm - xl );

            vm = ( xq - xl ) / ( xm - xl );
            vmp = +1.0       / ( xm - xl );

            vr =  0.0;
            vrp = 0.0;
      
            kxq = k( xq );
            fxq = f( xq );

            al = al + wq * ( kxq * vlp * vmp );
            am = am + wq * ( kxq * vmp * vmp );
            ar = ar + wq * ( kxq * vrp * vmp );
            bm = bm + wq * ( fxq * vm );
        

            #Integrate over the RIGHT interval, between XM and XR
            xq = ( ( 1.0 - abscissa[q] ) * xm
                 + ( 1.0 + abscissa[q] ) * xr ) / 2.0
            wq = weight[q] * ( xr - xm ) / 2.0

            vl = 0.0;
            vlp = 0.0;

            vm = ( xr - xq ) / ( xr - xm );
            vmp = -1.0       / ( xr - xm );

            vr = ( xq - xm ) / ( xr - xm );
            vrp = +1.0  / ( xr - xm );

            kxq = k( xq );
            fxq = f( xq );

            al = al + wq * ( kxq * vlp * vmp );
            am = am + wq * ( kxq * vmp * vmp );
            ar = ar + wq * ( kxq * vrp * vmp );
            bm = bm + wq * ( fxq * vm );

        end

        Amat[i,i-1] = al;
        Amat[i,i] = am;
        Amat[i,i+1] = ar;
        bvec[i] = bm;

    end

    Amat[n,n] = 1.0;
    bvec[n] = ub;

    u = Amat \ bvec;

    return u
end 