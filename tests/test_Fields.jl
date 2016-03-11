using Base.Test

using FAT.Meshes
using FAT.Fields
using FAT.Simulations

# get the test mesh
sim = SimulationData("ldc_test")

# test constructors
for D in [2, 3]
	zsf = ScalarField(mesh(sim), D, Int32)
	@test ndims(zsf) == D
	@test eltype(zsf) == Int32
	@test length(zsf.internalField) == ncells(mesh(sim))
	@test length(zsf.boundaryField) == nboundaryfaces(mesh(sim))
	@test mesh(zsf) == mesh(sim)
end 
@test_throws ErrorException ScalarField(mesh(sim), 4)
@test_throws ErrorException ScalarField(mesh(sim), 1)

for D in [2, 3]
	zvf = VectorField(mesh(sim), D)
	@test ndims(zvf) == D
	@test eltype(zvf) == Float64
	@test mesh(zvf) == mesh(sim)
	@test length(zvf.scalars) == ndims(zvf)
	for ui in zvf.scalars
		@test ndims(ui) == D
		@test eltype(ui) == Float64
		@test length(ui.internalField) == ncells(mesh(sim))
		@test length(ui.boundaryField) == nboundaryfaces(mesh(sim))
	end
end
@test_throws ErrorException VectorField(mesh(sim), 4)
@test_throws ErrorException VectorField(mesh(sim), 1)

for D in [2, 3]
	ztf = TensorField(mesh(sim), D, Float16)
	@test ndims(ztf) == D
	@test eltype(ztf) == Float16
	@test mesh(ztf) == mesh(sim)
	@test length(ztf.vectors) == ndims(ztf)
	for u in ztf.vectors
		@test ndims(u) == D
		@test eltype(u) == Float16
		@test length(u.scalars) == D
		for ui in u.scalars
			@test ndims(ui) == D
			@test eltype(ui) == Float16
			@test length(ui.internalField) == ncells(mesh(sim))
			@test length(ui.boundaryField) == nboundaryfaces(mesh(sim))
		end
	end
end
@test_throws ErrorException VectorField(mesh(sim), 4)
@test_throws ErrorException VectorField(mesh(sim), 1)

# test similar
for typ in [ScalarField, VectorField, TensorField]
	@eval begin
		a = $(typ)(mesh(sim), 2)
		b = similar(a)
		@test mesh(a) == mesh(b)
		@test mesh(a) === mesh(b)
		@test eltype(a) === eltype(b)
		@test ndims(a) === ndims(b)
	end
end

# test zero
a = ScalarField(mesh(sim), 3)
b = zero(a)
@test sumabs2(b.internalField) == 0
@test sumabs2(b.boundaryField) == 0

a = VectorField(mesh(sim), 3)
b = zero(a)
for d = 1:3
	@test sumabs2(b.scalars[d].internalField) == 0
	@test sumabs2(b.scalars[d].boundaryField) == 0
end

a = TensorField(mesh(sim), 3)
b = zero(a)
for di = 1:3, dj = 1:3
	@test sumabs2(b.vectors[dj].scalars[di].internalField) == 0
	@test sumabs2(b.vectors[dj].scalars[di].boundaryField) == 0
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ScalarField-ScalarField operations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subtraction
u = sim[10.0, :U].scalars[1]
out = u - u
@test sumabs2(out.internalField) == 0
@test sumabs2(out.boundaryField) == 0

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.sub!(u, u, similar(u))
@test sumabs2(out.internalField) == 0
@test sumabs2(out.boundaryField) == 0

# addition
u = sim[10.0, :U].scalars[1]
out = u + u
@test all(out.internalField .== 2*u.internalField)
@test all(out.boundaryField .== 2*u.boundaryField)

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.add!(u, u, similar(u))
@test all(out.internalField .== 2*u.internalField)
@test all(out.boundaryField .== 2*u.boundaryField)

# multiplication
u = sim[10.0, :U].scalars[1]
out = u * u
@test all(out.internalField .== u.internalField.^2.0)
@test all(out.boundaryField .== u.boundaryField.^2.0)

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.mul!(u, u, similar(u))
@test all(out.internalField .== u.internalField.^2.0)
@test all(out.boundaryField .== u.boundaryField.^2.0)

# test muladd, i.e., out = u*v + w
u = sim[10.0, :U].scalars[1]
v = 2.0*u
w = 3.0*u
out = FAT.Fields.muladd!(u, v, w, similar(u))
@test all(out.internalField .== 2*u.internalField.^2 + 3.0*u.internalField)
@test all(out.boundaryField .== 2*u.boundaryField.^2 + 3.0*u.boundaryField)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ScalarField-Real operations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# multiplication
u = sim[10.0, :U].scalars[1]
out = u * 1.5
@test all(out.internalField .== 1.5*u.internalField)
@test all(out.boundaryField .== 1.5*u.boundaryField)

u = sim[10.0, :U].scalars[1]
out = 2.5*u 
@test all(out.internalField .== 2.5*u.internalField)
@test all(out.boundaryField .== 2.5*u.boundaryField)

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.mul!(u, 1.5, similar(u))
@test all(out.internalField .== 1.5*u.internalField)
@test all(out.boundaryField .== 1.5*u.boundaryField)

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.mul!(1.5, u, similar(u))
@test all(out.internalField .== 1.5*u.internalField)
@test all(out.boundaryField .== 1.5*u.boundaryField)

# division
u = sim[10.0, :U].scalars[1]
out = u / 2.0
@test all(out.internalField .== 0.5*u.internalField)
@test all(out.boundaryField .== 0.5*u.boundaryField)

u = sim[10.0, :U].scalars[1]
out = FAT.Fields.div!(u, 3.0, similar(u))
@test all(out.internalField .== u.internalField/3.0)
@test all(out.boundaryField .== u.boundaryField/3.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VectorField-VectorField operations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subtraction 
u = sim[10.0, :U]
out = u - u
for i = 1:2
	@test sumabs2(out.scalars[i].internalField) == 0
	@test sumabs2(out.scalars[i].boundaryField) == 0
end

u = sim[10.0, :U]
out = FAT.Fields.sub!(u, u, similar(u))
for i = 1:2
	@test sumabs2(out.scalars[i].internalField) == 0
	@test sumabs2(out.scalars[i].boundaryField) == 0
end

# addition 
u = sim[10.0, :U]
out = u + u
for i = 1:2
	@test all(out.scalars[i].internalField .== 2*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 2*u.scalars[i].boundaryField)
end

u = sim[10.0, :U]
out = FAT.Fields.add!(u, u, similar(u))
for i = 1:2
	@test all(out.scalars[i].internalField .== 2*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 2*u.scalars[i].boundaryField)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VectorField-Real operations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# multiplication
u = sim[10.0, :U]
out = FAT.Fields.mul!(1.5, u, similar(u))
for i = 1:2
	@test all(out.scalars[i].internalField .== 1.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 1.5*u.scalars[i].boundaryField)
end

u = sim[10.0, :U]
out = FAT.Fields.mul!(u, 1.5, similar(u))
for i = 1:2
	@test all(out.scalars[i].internalField .== 1.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 1.5*u.scalars[i].boundaryField)
end

u = sim[10.0, :U]
out = 1.5*u
for i = 1:2
	@test all(out.scalars[i].internalField .== 1.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 1.5*u.scalars[i].boundaryField)
end

u = sim[10.0, :U]
out = u*1.5
for i = 1:2
	@test all(out.scalars[i].internalField .== 1.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 1.5*u.scalars[i].boundaryField)
end

# division
u = sim[10.0, :U]
out = FAT.Fields.div!(u, 2.0, similar(u))
for i = 1:2
	@test all(out.scalars[i].internalField .== 0.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 0.5*u.scalars[i].boundaryField)
end

u = sim[10.0, :U]
out = u/2.0
for i = 1:2
	@test all(out.scalars[i].internalField .== 0.5*u.scalars[i].internalField)
	@test all(out.scalars[i].boundaryField .== 0.5*u.scalars[i].boundaryField)
end

# ~~~~~~~
# dotgrad
# ~~~~~~~
u = sim[10.0, :U]
∇u = grad(u)
u∇u = dotgrad(u, ∇u)

# this tests the 2D version
out_1 = u.scalars[1]*∇u.vectors[1].scalars[1] + u.scalars[2]*∇u.vectors[1].scalars[2]
@test out_1.internalField == u∇u.scalars[1].internalField
@test out_1.boundaryField == u∇u.scalars[1].boundaryField
out_2 = u.scalars[1]*∇u.vectors[2].scalars[1] + u.scalars[2]*∇u.vectors[2].scalars[2]
@test out_2.internalField == u∇u.scalars[2].internalField
@test out_2.boundaryField == u∇u.scalars[2].boundaryField


# ~~~~~~~~~~~~~~~~~~~~~~
# Inner product and norm
# ~~~~~~~~~~~~~~~~~~~~~~
# VectorField
u = sim[10.0, :U]
u = 0.0*u
@test inner(u, u) == 0

# symmetry
u = sim[10.0, :U]
v = sim[9.0, :U]
@test inner(u, v) == inner(v, u)

# set scalar fields to constants. This means that we do 
# not test the case for which the fields vary in the
# domain, but only the underlying "algorithm".
for d in 1:2
	u.scalars[d].internalField[:] = d
end
@test inner(u, u) ≈ 1^2 + 2^2
@test norm(u) ≈ sqrt(5)

# ScalarField
@test inner(u.scalars[1], u.scalars[1]) ≈ 1.0
@test inner(u.scalars[2], u.scalars[2]) ≈ 4.0
@test inner(u.scalars[1], u.scalars[2]) ≈ 2.0
@test inner(u.scalars[2], u.scalars[1]) ≈ 2.0

@test norm(u) ≈ sqrt(5)


# ~~~~~~~~~~~~~~~~~~~
# Derivative and curl
# ~~~~~~~~~~~~~~~~~~~
# We check against the vorticity field calculated using OpenFoam. This
# should be enough to check the correctness of our algorithm. The rest 
# just tests the interface.
u = sim[1.0, :U]
ω = sim[1.0, :vorticity]
a = curl(u).internalField
b = ω.internalField
@test norm((a-b)./a, Inf) < 1e-8

u = sim[1.0, :U]
ω = sim[1.0, :vorticity]
out = ScalarField(mesh(u), ndims(u), eltype(u))
tmp = ScalarField(mesh(u), ndims(u), eltype(u))
a = curl!(u, out, tmp).internalField
b = ω.internalField
@test norm((a-b)./a, Inf) < 1e-8

# ~~~~~~~~
# Gradient
# ~~~~~~~~
# calculate gradient and then use OpenFoam vorticity field
u = sim[1.0, :U]
∇u = grad(u)
vort = ∇u.vectors[2].scalars[1] - ∇u.vectors[1].scalars[2]
ω = sim[1.0, :vorticity]
a = vort.internalField
b = ω.internalField
@test norm((a-b)./a, Inf) < 1e-8

# use in-place grad!
u = sim[1.0, :U]
∇u = grad!(u, TensorField(mesh(u), ndims(u), eltype(u)))
vort = ∇u.vectors[2].scalars[1] - ∇u.vectors[1].scalars[2]
ω = sim[1.0, :vorticity]
a = vort.internalField
b = ω.internalField
@test norm((a-b)./a, Inf) < 1e-8

# test outputs of gradient and curl
u = sim[1.0, :U]
for el in [mesh, ndims, eltype]
	@test el(grad(u)) == el(u)
	@test el(curl(u)) == el(u)
end