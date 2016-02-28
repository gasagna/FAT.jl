module Fields

import Base: call, *, /, +, -, .*, norm, show

using FAT.Meshes
export Field, ScalarField, VectorField, TensorField, inner, grad!, grad, zeroField, curl!


""" Generic tensor field type.

    The parameters R is 1 + the rank of the field, 
    i.e. 1 for scalar field, 2 for vector fields, and 
    3 for tensors. This matches the underlying 
    representation of the data. 

    We use double precision throughout. This could be 
    modified in a later stage.

    The `internalField` and `boundaryField` fields contain
    the solution inside the domain and on its boundary 
    as in output from OpenFoam simulations. 
"""
type Field{R}
    internalField::Array{Float64, R}
    boundaryField::Array{Float64, R}
    mesh::HexaMesh{Float64}
    Field(internalField::Array{Float64, R}, boundaryField::Array{Float64, R}, 
        mesh::HexaMesh{Float64}) = new(internalField, boundaryField, mesh)
end

# an outer constructor is usually required for more complex types
Field{R}(internalField::Array{Float64, R}, 
         boundaryField::Array{Float64, R}, 
         mesh::HexaMesh{Float64}) = Field{R}(internalField, boundaryField, mesh)

# define aliases for scalar, vector and tensor fields
typealias ScalarField Field{1}
typealias VectorField Field{2}
typealias TensorField Field{3}

" Construct a scalar, vector or tensor field of rank `rank` and of dimension `dim` that matches the `mesh` argument "
function zeroField{Tf<:Real}(rank::Integer, mesh::HexaMesh; dim::Integer=2, RealType::Type{Tf}=Float64)
    dim in [2, 3] || error("`dim` must be either 2 or 3")
    internalField = rank == 0 ? Array{Float64}(ncells(mesh)) :
                    rank == 1 ? Array{Float64}(ncells(mesh), dim) :
                    rank == 2 ? Array{Float64}(ncells(mesh), dim, dim) :
                    error("`rank` must be between 0 and 2")
    boundaryField = rank == 0 ? Array{Float64}(nboundaryfaces(mesh)) :
                    rank == 1 ? Array{Float64}(nboundaryfaces(mesh), dim) :
                    rank == 2 ? Array{Float64}(nboundaryfaces(mesh), dim, dim) :
                    error("`rank` must be between 0 and 2")
    return Field{rank+1}(internalField, boundaryField, mesh)
end

function show{R}(io::IO, u::Field{R})
    R == 1 ? print("Scalar field object") :
    R == 2 ? print("Vector field object") :
    R == 3 ? print("Tensor field object") : nothing
    print(" at $(object_id(u)): \n")
    print(" ~ $(size(u.internalField, 1)) cell centre values \n")
    print(" ~ $(size(u.boundaryField, 1)) boundary face values \n")
    print(" ~ on a ")
    show(io, mesh(u); space="  ")
end

# overload call to make aliases working constructors
call(::ScalarField, internalField::Array{Float64, 1},boundaryField::Array{Float64, 1}, 
     mesh::HexaMesh{Float64}) = Field{1}(internalField, boundaryField, mesh)
call(::VectorField, internalField::Array{Float64, 2},boundaryField::Array{Float64, 2}, 
     mesh::HexaMesh{Float64}) = Field{2}(internalField, boundaryField, mesh)
call(::TensorField, internalField::Array{Float64, 3},boundaryField::Array{Float64, 3}, 
     mesh::HexaMesh{Float64}) = Field{3}(internalField, boundaryField, mesh)

# get the mesh
mesh(u::Field) = u.mesh

# Now common operations between fields
# ====================================
*(u::Field, a::Real) = Field(scale(u.internalField, a), scale(u.boundaryField, a), mesh(u))   
*(a::Real, u::Field) = u*a   
/(u::Field, a::Real) = u*inv(a)

# operations on Fields of equal rank. We assume they have 
# same mesh and same amount of data. It would be strange 
# if they did not
for op in [:-, :+]
    @eval $op{R}(u::Field{R}, v::Field{R}) = Field{R}($op(u.internalField, v.internalField), $op(u.boundaryField, v.boundaryField), mesh(u))
end 

# negation
-(u::Field) = Field(-u.internalField, -u.boundaryField, mesh(u))

# multiplication of two ScalarField is product of data
*(u::ScalarField, v::ScalarField) = ScalarField(u.internalField.*v.internalField,
                                                u.boundaryField.*v.boundaryField, mesh(u))

# # dot product of two VectorfieldType is sum over product of components
# function *(u::VectorField, v::VectorField)
#     # we use maximum because
#     internalFieldOut = Array(Float64, maximum(size(u.internalField)))
#     boundaryFieldOut = Array(Float64, maximum(size(u.boundaryField)))
#     nc = ncells(mesh(u))
#     nbf = nboundaryfaces(mesh(u))
#     # note that we sum over all cells and we have a +nc to add the v component
#     @inbounds begin
#         for i in 1:nc
#             internalFieldOut[i] = (u.internalField[i]*v.internalField[i] + 
#                                    u.internalField[i+nc]*v.internalField[i+nc])
#         end
#         for i = 1:nbf
#             boundaryFieldOut[i] = (u.boundaryField[i]*v.boundaryField[i] + 
#                                    u.boundaryField[i+nbf]*v.boundaryField[i+nbf])
#         end
#     end
#     ScalarField(internalFieldOut, boundaryFieldOut, mesh(u))
# end

# inner product between two vector fields is the integral of the 
# dot product but the inner product can be computed without temporaries
function inner(u::VectorField, v::VectorField)
    # We need to sum the u and v components, so we add a skip.
    # Remember that the internalField field is an array of size
    # ncells x 2 or 3, if there are three velocity components
    # This is laid out in memory in column major format
    # and it is a big chunk of memory.
    N = ncells(mesh(u))
    I = 0.0
    @inbounds begin 
        for (i, cell) in enumerate(cells(mesh(u)))
            I += (u.internalField[i]*v.internalField[i] +
                  u.internalField[i+N]*v.internalField[i+N])*volume(cell)
        end
    end
    I
end

# the inner product results in the norm
norm(u::Field) = sqrt(inner(u, u))

# # integral of a scalar field
# function integral(u::ScalarField) 
#     I = 0.0
#     for (i, cell) in enumerate(cells(mesh(u)))
#         @inbounds I += u.internalField[i]*volume(cell)
#     end
#     I
# end
                        
# Dot product of vector field and tensor field is a vector field.
# This is useful to compute the projection onto the convective term.
# This only works in 2D now. Note that we do not fill the boundaryField 
# of the output here, whereas we probably should in some way. The 
# problem is that we do not have the boundaryField for the tensor, and we 
# should probably implement an interpolation scheme to interpolate
# values from the cells to the boundary faces. In reality it should not
# be needed because we would use this product u*V to compute inner(w, u*V)
# which does not involve the solution at the boundary of the domain, but 
# only the internal field.
function *(u::VectorField, V::TensorField) 
    out = zeroField(1, mesh(u))
    @inbounds begin
        for i in 1:ncells(mesh(u))
            out.internalField[i, 1] = (u.internalField[i, 1]*V.internalField[i, 1, 1]+ 
                                       u.internalField[i, 2]*V.internalField[i, 1, 2])
            out.internalField[i, 2] = (u.internalField[i, 1]*V.internalField[i, 2, 1]+
                                       u.internalField[i, 2]*V.internalField[i, 2, 2])
        end
    end
    out
end


" Compute gradient of vector field "
grad(u::VectorField) = grad!(u, zeroField(2, mesh(u)))

" Compute gradient of vector field. In place version, without allocations "
function grad!(u::VectorField, out::TensorField)
    @inbounds begin
        out.internalField[:] = 0.0
        out.boundaryField[:] = 0.0

        # for (faceID, face) in enumerate(boundaryfaces(mesh(u)))
        #     ownID = ownerID(face)
        #     face_value_x = u.boundaryField[faceID - ninternalfaces(mesh(u)), 1]
        #     face_value_y = u.boundaryField[faceID - ninternalfaces(mesh(u)), 2]
        #     out.internalField[ownID, 1, 1] += svec(face, Val{'x'})*face_value_x
        #     out.internalField[ownID, 1, 2] += svec(face, Val{'x'})*face_value_y
        #     out.internalField[ownID, 2, 1] += svec(face, Val{'y'})*face_value_x
        #     out.internalField[ownID, 2, 2] += svec(face, Val{'y'})*face_value_y
        # end
        for (i, face) in enumerate(internalfaces(mesh(u)))
            # get owner and neighbour of the present face
            ownID, neighID = ownerID(face), neighbourID(face)     
            
            # weight for linear interpolation of face centre
            α = mesh(u).αs[i]
            face_value_x = (1.0-α)*u.internalField[ownID, 1] + α*u.internalField[neighID, 1]
            face_value_y = (1.0-α)*u.internalField[ownID, 2] + α*u.internalField[neighID, 2]

            # assign value to owner 
            out.internalField[ownID, 1, 1] += svec(face, Val{:x})*face_value_x
            out.internalField[ownID, 1, 2] += svec(face, Val{:x})*face_value_y
            out.internalField[ownID, 2, 1] += svec(face, Val{:y})*face_value_x
            out.internalField[ownID, 2, 2] += svec(face, Val{:y})*face_value_y
            # and neighbour. The minus sign is because faces are always oriented
            # from the owner to the neighbour, which gives a negative contributions
            # to the gradient for the neighbour cell.
            out.internalField[neighID, 1, 1] -= svec(face, Val{:x})*face_value_x
            out.internalField[neighID, 1, 2] -= svec(face, Val{:x})*face_value_y
            out.internalField[neighID, 2, 1] -= svec(face, Val{:y})*face_value_x
            out.internalField[neighID, 2, 2] -= svec(face, Val{:y})*face_value_y
        end
        # normalise result with volume 
        for (i, cell) in enumerate(cells(mesh(u)))
            # this should become 1:3
            for j = 1:2, k = 1:2
                out.internalField[i, j, k] /= volume(cell)
            end
        end

        # Here we should do interpolation from the cell values to the boundary faces.
        # In some cases the value of the velocity is known at the boundaries, in some
        # others the gradient is known. However, for simplicity of the code, it is probably
        # a better idea to always use interpolation, maybe a good scheme.
    end
    out
end

" Compute scalar vorticity field "
function curl!(u::VectorField, ω::ScalarField, tmp::TensorField)
    tmp = grad!(u, tmp)
    # We now only compute ω_z = ∂v/∂x - ∂u/∂y. In the three dimensional case
    # we should also compute the oher components.
    @inbounds begin
        for i = 1:ninternalfaces(mesh(u))
            ω.internalField[i] = tmp.internalField[i, 2, 1] - tmp.internalField[i, 1, 2]
        end
        for i = 1:nboundaryfaces(mesh(u))
            ω.boundaryField[i] = tmp.boundaryField[i, 2, 1] - tmp.boundaryField[i, 1, 2]
        end
    end
    ω
end




end













