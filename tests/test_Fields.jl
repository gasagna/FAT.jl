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


# # test function, a linear function
# # FIXME: we should test different functions
# f(x, y) = 1 + x - 2y

# # set internal field
# for (i, cell) in enumerate(cells(mesh(sim)))
# 	# the average of the function f for linear function
# 	# is the value at the cell centre.
# 	p.internalField[i] = f(centre(cell).x, centre(cell).y)
# end

# # we should use the faces to compute the gradient, so we 
# # also fill the boundary field
# for (i, face) in enumerate(boundaryfaces(mesh(sim)))
#  	p.boundaryField[i] = f(centre(face).x, centre(face).y)
# end


# # test

# # test derivative
# for (dir, val) in zip([:x, :y], [1.0, -2.0])
# 	dpdi = zeroScalarField(mesh(sim), 2)
# 	FAT.Fields.der!(p, dpdi, dir)
# 	for (i, cell) in enumerate(cells(mesh(sim)))
# 		@test dpdi.internalField[i] ≈ val
# 	end
# end

# # cannot compute derivative with respect to z for 2d field
# dpdi = zeroScalarField(mesh(sim), 2)
# @test_throws ErrorException FAT.Fields.der!(p, dpdi, :z)
# @test_throws ErrorException FAT.Fields.der!(p, dpdi, 3)

# # we have not filled the boundary field by interpolation
# dpdi = zeroScalarField(mesh(sim), 2)
# FAT.Fields.der!(p, dpdi, :x)
# @test all(dpdi.boundaryField .== 0.0)

# # test gradient

# # # test norm and inner product
# # u.internalField[:] = 2.0
# # @test inner(u, u) ≈ 8
# # @test inner(u, u) == norm(u)^2

# # # test zero constructor
# # scalar = zeroField(0, mesh(sim))
# # @test ndims(scalar.internalField) == 1
# # @test ndims(scalar.boundaryField) == 1

# # vector = zeroField(1, mesh(sim))
# # @test ndims(vector.internalField) == 2
# # @test ndims(vector.boundaryField) == 2

# # tensor = zeroField(2, mesh(sim))
# # @test ndims(tensor.internalField) == 3
# # @test ndims(tensor.boundaryField) == 3


# # # test velocity gradient 
# # gradu = zeroField(2, mesh(sim))
# # @time grad!(u, gradu)
# # @time grad!(u, gradu)
# # println(1000*minimum([@elapsed grad!(u, gradu) for i = 1:100]), " ms")


# # @profile minimum([@elapsed grad!(u, gradu) for i = 1:1000])
# # Profile.print()

# # @code_llvm grad!(u, gradu)



# # # # test gradient of a scalar field
# # # uu = ScalarField(u.internalField[1:ncells(FAT.Fields.mesh(u))], 
# # #  				 u.boundaryField[1:nboundaryfaces(FAT.Fields.mesh(u))], 
# # #  				 FAT.Fields.mesh(u))
# # # utmp = zeroField(1, mesh(sim))
# # # println(1000*minimum([@elapsed grad!(uu, utmp) for i = 1:100]), " ms")

# # # # test calculation of vorticity
# # # gradu = zeroField(2, mesh(sim))
# # # ω = zeroField(0, mesh(sim))

# # # import PyPlot: tricontourf, show, pygui, scatter, figure, clf, savefig, plot

# # # x, y, z = cellcentres(mesh(sim))

# # # for (i, v) in enumerate(fields(sim, :U, 100.0:0.5:105.0))
# # # 	ω = curl!(v, ω, gradu)
# # # 	tricontourf(x, y, ω.internalField[:], -15:0.1:15, cmap=PyPlot.cm[:seismic])
# # # 	name = @sprintf "figure-%04d.png" i
# # # 	plot(x, y, "k.", ms=1, alpha=0.2)
# # # 	savefig("vid/$name", dpi=800)
# # # 	println(i)
# # # 	close()
# # # end
# # # # tricontourf(x, y, gradu.internalField[:, 1, 1], -5:0.05:5)
