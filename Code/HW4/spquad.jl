import Statistics as ST
using FFTW: ifft
using LazyGrids:ndgrid

function matrix_intersection(A,B)

    indices = Int[]
    (m,n) = size(A)
    for i in 1:m
        row = A[i,:]
        for row_B in eachrow(B)
            if row == row_B
                push!(indices,i)
                break
            end
        end
    end

    # return_matrix = Matrix{Float64}(undef,length(indices),n)
    # for i in 1:length(indices)
    #     return_matrix[i,:] = A[indices[i],:]
    # end

    # return return_matrix, indices
    return indices
end



function spquad(dim, ord, bpt = nothing)
    # Compute the sparse grid quadrature abscissae and weights
    # on an orthotope/hyperrectangle using the Clenshaw-Curtis rule.

    # Input Parameters:
    # DIM - number of dimensions
    # ORD - order of the integration rule (Alireza Doostan: level parameter l in the lecture notes)
    # BPT - boundary points (Optional. Defaults to [-1,1]^dim)

    # If used, BPT should be a 2 by DIM matrix containing the
    # endpoints of the hyperrectangle in each column.

    # Example usage:
    # f=@(x) (1+x(:,1)).*exp(x(:,2).*x(:,3))+(1+x(:,2)).*exp(x(:,3).*x(:,1))
    # [x,w]=spquad(3,4,[-1 0 2; 1 1 3]);
    # Q=w'*f(x);

    # written by Greg von Winckel - March 3, 2008
    # Contact: gregory(dot)von-winckel(at)uni-graz(dot)at
    # Check that the number dimensions agree
    if bpt != nothing
        if size(bpt, 2) != dim
            error("dimension mismatch")
        end
    else
        a = [-1,1]'
        bpt = repeat(a,dim)'  #Does what repmat in matlab does
    end
    # Length and midpoint in each dimension
    len = diff(bpt,dims=1)
    mpt = ST.mean(bpt,dims=1)
    if ord == 0
        # zeroth order case is just the midpoint rule
        x = mpt
        w = prod(len)
    elseif ord != 0 && dim == 1
        # Special 1D case
        x = mpt .+ len .* cos.(pi * 2.0 .^ (-ord) .* (0:(2^ord))) ./ 2
        w = clencurt(2^ord + 1) .* len ./ 2
    else
        # node and weight for single point grid
        x0 = mpt
        w0 = prod(len)
        # Construct higher order grid hierarchically
        for j = 1:ord
            x, w = sparsegridnd(dim, j, mpt, len)
            I = matrix_intersection(x, x0)
            w[I] .+= w0
            # for i in I
            #     w[i] += w0
            # end
            x0, w0 = x, w
        end
    end
    return x, w
end


#=
The uniq function is an attempt to replicate the `unique` function in Matlab.
The code was borrowed from here: https://discourse.julialang.org/t/unique-indices-method-similar-to-matlab/34446/7
=#
uniq(A; dims=1) = begin
    @assert ndims(A) ∈ (1, 2)
    slA = ndims(A) > 1 ? eachslice(A; dims) : A
  
    ia = unique(i -> slA[i], axes(A, dims))
    sort!(ia; by=i -> slA[i])
  
    C = stack(slA[ia]; dims)
    slC = ndims(A) > 1 ? eachslice(C; dims) : C
  
    ic = map(r -> findfirst(==(slA[r]), slC), axes(A, dims))
  
    C, ia, ic
  end

function sparsegridnd(n, ord, mpt, len)

    # Get configurations of all possible subgrids
    v = genindex(n, ord);
    vmax = 1 + v[1, end];

    # Compute all orders of one-dimensional quadrature rules
    P = uaf( get_cos_points, (0:vmax-1),1 )
    Q = uaf( diffweight, (0:vmax),1 )

    # Take the union of all possible subgrids
    m = size(v, 1)
    dim_for_cell2mat = 1
    xw = uaf(k -> getpts(P, Q, v[k, :], mpt, len), 1:m, dim_for_cell2mat, 1);
    x = xw[:, 1:n]
    w = xw[:, n+1]

    # Kludge to deal with small errors introduced by ndgrid
    # need a better way to do this
    roundn(a, n) = round(a * 10^n) / 10^n
    x = roundn.(x, 15)

    # Get unique points and weights
    final_x, ii, jj = uniq(x, dims=1)
    rows = length(ii)
    cols = length(jj)

    # Do node condesation for combining weights
    new_w = zeros(rows, cols)
    for index in 1:cols
        i_index = jj[index]
        j_index = index   
        new_w[i_index , j_index] = 1.0
    end
    final_w = new_w*w

    return final_x, final_w
end


function get_cos_points(x)
    r=0:2^x
    r = r * pi*2^(-x*1.0)
    r = cos.(r)
    r = r.*(x!=0)
    if(x==0)
        r = [0.0]
        return r
    else
        return sort(union(r))
    end
end


function getpts(P,Q,v,mpt,len)

    n = length(v)
    grid_points1 = [P[i + 1] for i in v] 
    x = ndgrid(grid_points1...)
    grid_points2 = [Q[i + 1] for i in v] 
    w = ndgrid(grid_points2...)
    d = length(v)

    dim_for_cell2mat = 2
    X = uaf(k -> mpt[k] .+ ( x[k][:] .* (len[k]/2) ), 1:d, dim_for_cell2mat,1)
    W = prod(uaf(k -> w[k][:], 1:d, dim_for_cell2mat, 1), dims=2)
    XW = hcat(X, W)
    
    # println(X)
    # # return X
    # println(W)
    return XW
end


function genindex(n, L1, head=nothing)
    if n == 1
        # v = ones(1,1)*L1
        v = L1
    else
        v = uaf(j -> genindex(n - 1, L1 - j, j), (0:L1)',1,1)
    end
    
    if head !== nothing
        # println("V is ", v)
        new_column = fill(head, size(v, 1))
        v = hcat(new_column, v) # Or v = [new_column v]
        # println("Modified V is ", v)
        # v = [head .+ zeros(size(v, 1)) ; v]
        # v = hcat(fill(head, size(v, 1)), v)
    end
    
    return v
end


function diffweight(ii)
    if ii == 0
        dw = [2.0]
    elseif ii == 1
        dw = [1; -2; 1] / 3
    else
        q(x) = clencurt(2^(x) + 1)
        dw = q(ii)
        dw[1:2:end] .-= q(ii - 1)
    end
    return dw
end

function clencurt(N1)
    if N1 == 1
        w = 2
    else
        N = N1 - 1
        c = zeros(N1,1)
        c[1:2:N1,1] .= (2 ./ [1; 1 .- (2:2:N).^2])
        f = real(ifft([c[1:N1]; c[N:-1:2]],1))
        w = 2 * ([f[1]; 2 .* f[2:N]; f[N1]]) / 2
    end
end


# Shorthand array function with uniform output
# Optional third argument converts output to matrix
function cell2mat(cell_array,dim_for_cell2mat)
    # num_rows = length(cell_array)
    # num_columns = length(cell_array[1])
    # tbr = Matrix{eltype(cell_array[1][1])}(undef, num_rows, num_columns)
    # for i in 1:num_rows
    #     for j in 1:num_columns
    #         tbr[i, j] = cell_array[i][j]
    #     end
    # end
    if(dim_for_cell2mat==2)
        return hcat(cell_array...)
    elseif(dim_for_cell2mat==1)
        return vcat(cell_array...)
    end
end

function uaf(arg, dex, d, varargin=nothing)
    result = map(x->arg(x), dex)
    # println("Before")
    # println(result)
    if varargin != nothing
        result = cell2mat(result,d)
    end
    # println("After")
    # println(result)
    return result
end



#=
dim = 3
ord = 4

a = [-1,1]'
bpt = repeat(a,dim)'  
len = diff(bpt,dims=1)
mpt = ST.mean(bpt,dims=1)


n = dim
v = genindex(n,ord);
vmax = 1+v[1,n]

P = uaf( get_cos_points, (0:vmax-1) )
Q = uaf( diffweight, (0:vmax) )
m = size(v, 1)




=#