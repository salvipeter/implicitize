module MRep

using LinearAlgebra
using Printf

# Bezier-surfaces are represented by a 2D array of control points

const Point = Vector{Float64}

"""
    read_bzr(filename)

Reads a Bezier surface from file.
"""
function read_bzr(filename)
    read_numbers(f, numtype) = map(s -> parse(numtype, s), split(readline(f)))
    open(filename) do f
        n, m = read_numbers(f, Int)
        surface = Array{Point,2}(undef, n + 1, m + 1)
        for i in 1:n+1, j in 1:m+1
            surface[i,j] = read_numbers(f, Float64)
        end
        surface
    end
end

"""
    write_bzr(surface, filename)

Writes a Bezier surface to file.
"""
function write_bzr(surface, filename)
    n, m = size(surface) .- 1
    open(filename, "w") do f
        println(f, "$n $m")
        for i in 1:n+1, j in 1:m+1
            p = surface[i,j]
            println(f, "$(p[1]) $(p[2]) $(p[3])")
        end
    end
end

"""
    bernstein_all(n, u)

Returns all Bernstein polynomials of degree `n` for parameter `u`.
"""
function bernstein_all(n, u)
    coeff = [1.0]
    for j in 1:n
        saved = 0.0
        for k in 1:j
            tmp = coeff[k]
            coeff[k] = saved + tmp * (1.0 - u)
            saved = tmp * u
        end
        push!(coeff, saved)
    end
    coeff
end

"""
    bezier_evaluate(surface, uv)

Evaluates a Bezier surface at `uv`.
"""
function bezier_evaluate(surface, uv)
    n, m = size(surface) .- 1
    coeff_u = bernstein_all(n, uv[1])
    coeff_v = bernstein_all(m, uv[2])
    result = [0., 0, 0]
    for i in 1:n+1, j in 1:m+1
        result += surface[i,j] * coeff_u[i] * coeff_v[j]
    end
    result
end

"""
    bezier_derivative_u(surface)

Computes the control points of the first `u`-derivative of `surface`.
"""
function bezier_derivative_u(surface)
    n, m = size(surface) .- 1
    cnet = Array{Point}(undef, n, m + 1)
    for i in 1:n, j in 1:m+1
        cnet[i,j] = (surface[i+1,j] - surface[i,j]) * n
    end
    cnet
end

"""
    bezier_derivative_u(surface, uv)

Computes the first `u`-derivative of `surface` at `uv`.
"""
bezier_derivative_u(surface, uv) = bezier_evaluate(bezier_derivative_u(surface), uv)

"""
    bezier_derivative_v(surface)

Computes the control points of the first `v`-derivative of `surface`.
"""
function bezier_derivative_v(surface)
    n, m = size(surface) .- 1
    cnet = Array{Point}(undef, n + 1, m)
    for i in 1:n+1, j in 1:m
        cnet[i,j] = (surface[i,j+1] - surface[i,j]) * m
    end
    cnet
end

"""
    bezier_derivative_v(surface, uv)

Computes the first `v`-derivative of `surface` at `uv`.
"""
bezier_derivative_v(surface, uv) = bezier_evaluate(bezier_derivative_v(surface), uv)

"""
    bezier_distance(surface, uv, p)

Computes the unsigned distance between the point `p` and the `surface` evaluated at `uv`.
"""
bezier_distance(surface, uv, p) = norm(p - bezier_evaluate(surface, uv))

"""
    bezier_signed_distance(surface, uv, p)

Computes the signed distance between the point `p` and the `surface` evaluated at `uv`.
"""
function bezier_signed_distance(surface, uv, p)
    q = bezier_evaluate(surface, uv)
    du = bezier_derivative_u(surface, uv)
    dv = bezier_derivative_v(surface, uv)
    n = normalize!(cross(du, dv))
    norm(p - q) * (dot(p - q, n) < 0 ? -1 : 1)
end

"""
    write_obj(verts, tris, filename)

Writes the given vertices and triangles into a Wavefront Obj file.
"""
function write_obj(verts, tris, filename)
    open(filename, "w") do f
        for v in verts
            println(f, "v $(v[1]) $(v[2]) $(v[3])")
        end
        for t in tris
            println(f, "f $(t[1]) $(t[2]) $(t[3])")
        end
    end
end

"""
    write_bezier_mesh(surface, filename, resolution)

Samples `surface` with the given `resolution`, and writes the result
into the specified Wavefront Obj file.
"""
function write_bezier_mesh(surface, filename, resolution)
    samples = [[u, v] for u in range(0.0, stop=1.0, length=resolution)
                      for v in range(0.0, stop=1.0, length=resolution)]
    verts = map(uv -> bezier_evaluate(surface, uv), samples)
    tris = []
    for i in 2:resolution, j in 2:resolution
        index = (j - 1) * resolution + i
        push!(tris, [index - resolution - 1, index - resolution, index])
        push!(tris, [index, index - 1, index - resolution - 1])
    end
    write_obj(verts, tris, filename)
end

"""
    write_cnet(surface, filename)

Write the control network of the Bezier `surface`
into the specified Wavefront .obj file.
"""
function write_cnet(surface, filename)
    d = size(surface, 1) - 1
    open(filename, "w") do f
        for i in 0:d, j in 0:d
            p = surface[i+1,j+1]
            println(f, "v $(p[1]) $(p[2]) $(p[3])")
        end
        for i in 1:d+1, j in 1:d
            index = (i - 1) * (d + 1) + j
            println(f, "l $index $(index+1)")
            index = (j - 1) * (d + 1) + i
            println(f, "l $index $(index+d+1)")
        end
    end
end

"""
    bezier_subdivide_u(surface, u)

Splits a Bezier surface into two parts along the isocurve at parameter `u`.
The return value is `(left, right)`.
"""
function bezier_subdivide_u(surface, u)
    n, m = size(surface) .- 1
    left, right = similar(surface), similar(surface)
    left[1,:], right[end,:] = surface[1,:], surface[end,:]
    tmp = copy(surface)
    for k in 1:n
        for i in 1:n-k+1
            tmp[i,:] = tmp[i,:] * (1 - u) + tmp[i+1,:] * u
        end
        left[k+1,:], right[end-k,:] = tmp[1,:], tmp[n-k+1,:]
    end
    (u < 0 ? reverse(left, dims = 1) : left, u > 1 ? reverse(right, dims = 1) : right)
end

"""
    bezier_subdivide_v(surface, v)

Splits a Bezier surface into two parts along the isocurve at parameter `v`.
The return value is `(left, right)`.
"""
function bezier_subdivide_v(surface, v)
    n, m = size(surface) .- 1
    left, right = similar(surface), similar(surface)
    left[:,1], right[:,end] = surface[:,1], surface[:,end]
    tmp = copy(surface)
    for k in 1:m
        for j in 1:m-k+1
            tmp[:,j] = tmp[:,j] * (1 - v) + tmp[:,j+1] * v
        end
        left[:,k+1], right[:,end-k] = tmp[:,1], tmp[:,m-k+1]
    end
    (v < 0 ? reverse(left, dims = 2) : left, v > 1 ? reverse(right, dims = 2) : right)
end

"""
    bezier_subdivide_uv(surface, uv)

Splits a Bezier surface into four parts along the isocurves at parameter `uv`.
The return value is an array of [top-left, bottom-left, top-right, bottom-right].

TODO: this should be done more efficiently.
"""
function bezier_subdivide_uv(surface, uv)
    left, right = bezier_subdivide_u(surface, uv[1])
    tl, bl = bezier_subdivide_v(left, uv[2])
    tr, br = bezier_subdivide_v(right, uv[2])
    [tl, bl, tr, br]
end

"""
    bezier_bbox(bezier)

Returns the bounding box of a Bezier curve or surface as a pair of points.
"""
function bezier_bbox(bezier)
    (reduce((p,q) -> min.(p, q), bezier),
     reduce((p,q) -> max.(p, q), bezier))
end


# Bezier curve functions

"""
    curve_evaluate(curve, u)

Evaluates a Bezier curve given by its control points `curve` at parameter `u`.
"""
function curve_evaluate(curve, u)
    d = length(curve) - 1
    coeff = bernstein_all(d, u)
    result = [0., 0, 0]
    for i in 1:d+1
        result += curve[i] * coeff[i]
    end
    result
end

"""
    curve_derivative(curve, u)

Computes the first derivative of a Bezier curve
given by its control points `curve` at parameter `u`.
"""
function curve_derivative(curve, u)
    d = length(curve) - 1
    dcp = map(i -> d * (curve[i+1] - curve[i]), 1:d)
    curve_evaluate(dcp, u)
end


# Volume I/O

"""
    sample_volume(f, bbox, resolution)

Samples a function `f` inside the given bounding box `bbox` (a pair of points).
The `resolution` is given as the # of voxels in the shortest edge of the bounding box.
"""
function sample_volume(f, bbox, resolution)
    axis = bbox[2] - bbox[1]
    minlen = minimum(axis)
    size = minlen / resolution
    res = Int.(ceil.(axis / size))
    voxels = Array{Float64}(undef, res[1], res[2], res[3])
    for i in 1:res[1], j in 1:res[2], k in 1:res[3]
        p = bbox[1] + [i, j, k] * size
        voxels[i,j,k] = f(p)
    end
    voxels
end

"""
    write_volume(voxels, filename)

Writes the 3D matrix `voxels` into `filename.raw`
with the metaimage header `filename.mhd`.
"""
function write_volume(voxels, filename)
    n = size(voxels)
    open("$(filename).mhd", "w") do f
        println(f, "NDims = 3")
        println(f, "DimSize = $(n[1]) $(n[2]) $(n[3])")
        println(f, "ElementSize = 1.0 1.0 1.0")
        println(f, "ElementSpacing = 1.0 1.0 1.0")
        println(f, "ElementType = MET_DOUBLE")
        println(f, "ElementByteOrderMSB = False") # LittleEndian
        println(f, "ElementDataFile = $(filename).raw")
    end
    open("$(filename).raw", "w") do f
        for val in voxels
            write(f, val)
        end
    end
end


# Implicit Matrix Representations based on Buse Laurent's paper (CAD 46, pp. 14-24, 2014)

"""
    curve_to_matrix(curve; v, epsilon)

Converts a Bezier `curve` to implicit matrix form.
`epsilon` defines the tolerance for minimal non-zero singular values;
`v` is the degree used by the matrix.

Note that this is a rational curve, all control points have 4 coordinates,
the first being the weight (which the rest are scaled with).
"""
function curve_to_matrix(curve; v = nothing, epsilon = 1.0e-8)
    d = length(curve) - 1
    if v == nothing
        v = d - 1
    end
    S = zeros(d + v + 1, 4(v + 1))
    for c in 0:3, j in 0:v, i in 0:d
        w = binomial(v, j) * binomial(d, i) / binomial(d + v, i + j)
        S[i+j+1,c*(v+1)+j+1] += w * curve[i+1][c+1]
    end
    F = svd(S, full=true)
    from = findfirst(x -> x < epsilon, F.S)
    if from == nothing
        from = d + v + 2
    end
    F.Vt[from:end,:]'
end

"""
    trinomial(d, i, j)

Generalization of binomials; `d` is the degree, the three variables are
`(i, j, d - i - j)`.

The computation involves factorials, which limits its use to low degrees.
"""
trinomial(d, i, j) = factorial(d) / (factorial(i) * factorial(j) * factorial(d - i - j))

"""
    triangle_to_matrix(triangle; v, epsilon)

Converts a Bezier `triangle` to implicit matrix form.
The surface is given as a dictionary, having a `:degree` entry,
as well as `(i,j)` entries for control points with indices `(i,j,d-i-j)`.

`epsilon` defines the tolerance for minimal non-zero singular values;
`v` is the degree used by the matrix.

Note that this is a rational curve, all control points have 4 coordinates,
the first being the weight (which the rest are scaled with).
"""
function triangle_to_matrix(triangle; v = nothing, epsilon = 1.0e-8)
    index(d, i, j) = binomial(d - i + 1, 2) + j + 1
    d = triangle[:degree]
    if v == nothing
        v = 2(d - 1)
    end
    nrow = binomial(d + v + 2, 2)
    ncol = binomial(v + 2, 2)
    S = zeros(nrow, 4 * ncol)
    for c in 0:3, k in 0:v, l in 0:v-k, i in 0:d, j in 0:d-i
        w = trinomial(v, k, l) * trinomial(d, i, j) / trinomial(d + v, i + k, j + l)
        S[index(d+v,i+k,j+l),c*ncol+index(v,k,l)] += w * triangle[i,j][c+1]
    end
    F = svd(S, full=true)
    from = findfirst(x -> x < epsilon, F.S)
    if from == nothing
        from = nrow + 1
    end
    F.Vt[from:end,:]'
end

"""
    surface_to_matrix(surface; v, epsilon)

Converts a Bezier `surface` to implicit matrix form.
`epsilon` defines the tolerance for minimal non-zero singular values;
`v` is a two-element array containing the degrees used by the matrix.

Note that unlike the other two matrix-conversion functions,
this uses non-rational Bezier surfaces of the same format
as the rest of this file.
"""
function surface_to_matrix(surface; v = nothing, epsilon = 1.0e-8)
    index(d, i, j) = i * (d[2] + 1) + j + 1
    d = collect(size(surface) .- 1)
    if v == nothing
        v = [2 * d[1] - 1, d[2] - 1]
    end
    nrow = (v[1] + d[1] + 1) * (v[2] + d[2] + 1)
    ncol = (v[1] + 1) * (v[2] + 1)
    S = zeros(nrow, 4 * ncol)
    for c in 0:3, k in 0:v[1], l in 0:v[2], i in 0:d[1], j in 0:d[2]
        w = binomial(v[1], k) * binomial(v[2], l) * binomial(d[1], i) * binomial(d[2], j) /
            (binomial(v[1] + d[1], i + k) * binomial(v[2] + d[2], j + l))
        p = c == 0 ? 1 : surface[i+1,j+1][c]
        S[index(d+v,i+k,j+l),c*ncol+index(v,k,l)] += w * p
    end
    F = svd(S, full=true)
    from = findfirst(x -> x < epsilon, F.S)
    if from == nothing
        from = nrow + 1
    end
    F.Vt[from:end,:]'
end

"""
    matrix_to_distance(matrix)

Converts an implicit matrix representation into a distance function.
"""
function matrix_to_distance(matrix)
    rows = div(size(matrix, 1), 4)
    M0 = matrix[1:rows,:]
    M1 = matrix[rows+1:2*rows,:]
    M2 = matrix[2*rows+1:3*rows,:]
    M3 = matrix[3*rows+1:end,:]
    function(p)
        M = M0 + p[1] * M1 + p[2] * M2 + p[3] * M3
        prod(svdvals(M))
    end
end

"""
    examples(resolution)

Generates implicit matrix representation examples at a given resolution.
"""
function examples(resolution = 30)
    # curve = [[1,0,0,0], [1,1/3,0,0], [1,2/3,1/3,0], [1,1,1,1]]
    # m = curve_to_matrix(curve)

    # Bezier triangle on the unit sphere
    # triangle = Dict(:degree => 2,
    #                 (0,0) => [1,1,0,0],
    #                 (0,1) => [1,1,0,1],
    #                 (0,2) => [2,0,0,2],
    #                 (1,0) => [1,1,1,0],
    #                 (1,1) => [1,1,1,1],
    #                 (2,0) => [2,0,2,0])
    # m = triangle_to_matrix(triangle)

    surface = read_bzr("test.bzr")
    m = surface_to_matrix(surface)

    f = matrix_to_distance(m)
    bbox = [[0,0,0],[3,3,3]]
    voxels = sample_volume(f, bbox, resolution)
    write_volume(voxels, "/tmp/test")
end

end # module
