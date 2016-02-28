using Base.Test

using FAT.Meshes
using FAT.Fields
using FAT.Simulations

# get the test mesh
sim = SimulationData("ldc_test")


# test zero constructors
for D in [2, 3]
	zsf = zeroScalarField(mesh(sim), D; dtype=Int32)
	@test ndims(zsf) == D
	@test eltype(zsf) == Int32
	@test length(zsf.internalField) == ncells(mesh(sim))
	@test length(zsf.boundaryField) == nboundaryfaces(mesh(sim))
	@test mesh(zsf) == mesh(sim)
end 
@test_throws ErrorException zeroScalarField(mesh(sim), 4)
@test_throws ErrorException zeroScalarField(mesh(sim), 1)

for D in [2, 3]
	zvf = zeroVectorField(mesh(sim), D)
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
@test_throws ErrorException zeroVectorField(mesh(sim), 4)
@test_throws ErrorException zeroVectorField(mesh(sim), 1)

for D in [2, 3]
	ztf = zeroTensorField(mesh(sim), D)
	@test ndims(ztf) == D
	@test eltype(ztf) == Float64
	@test mesh(ztf) == mesh(sim)
	@test length(ztf.vectors) == ndims(ztf)
	for u in ztf.vectors
		@test ndims(u) == D
		@test eltype(u) == Float64
		@test length(u.scalars) == D
		for ui in u.scalars
			@test ndims(ui) == D
			@test eltype(ui) == Float64
			@test length(ui.internalField) == ncells(mesh(sim))
			@test length(ui.boundaryField) == nboundaryfaces(mesh(sim))
		end
	end
end
@test_throws ErrorException zeroVectorField(mesh(sim), 4)
@test_throws ErrorException zeroVectorField(mesh(sim), 1)

# get one scalar field
p = sim[10.0, :p]

# test arithmetic
@test (p*1.0).internalField == (1.0*p).internalField
@test (p*1.0).boundaryField == (1.0*p).boundaryField
@test (p/2.0).internalField == (0.5*p).internalField
@test (p/2.0).boundaryField == (0.5*p).boundaryField

# get one scalar field
U = sim[10.0, :U]
for i = 1:ndims(U)
	# test arithmetic
	@test (U*1.0).scalars[i].internalField == (1.0*U).scalars[i].internalField
	@test (U*1.0).scalars[i].boundaryField == (1.0*U).scalars[i].boundaryField
	@test (U/2.0).scalars[i].internalField == (0.5*U).scalars[i].internalField
	@test (U/2.0).scalars[i].boundaryField == (0.5*U).scalars[i].boundaryField
end



# test function, a linear function
# FIXME: we should test different functions
f(x, y) = 1 + x - 2y

# set internal field
for (i, cell) in enumerate(cells(mesh(sim)))
	# the average of the function f for linear function
	# is the value at the cell centre.
	p.internalField[i] = f(centre(cell).x, centre(cell).y)
end

# we should use the faces to compute the gradient, so we 
# also fill the boundary field
for (i, face) in enumerate(boundaryfaces(mesh(sim)))
 	p.boundaryField[i] = f(centre(face).x, centre(face).y)
end


# test

# test derivative
for (dir, val) in zip([:x, :y], [1.0, -2.0])
	dpdi = zeroScalarField(mesh(sim), 2)
	FAT.Fields.der!(p, dpdi, dir)
	for (i, cell) in enumerate(cells(mesh(sim)))
		@test dpdi.internalField[i] ≈ val
	end
end

# cannot compute derivative with respect to z for 2d field
dpdi = zeroScalarField(mesh(sim), 2)
@test_throws ErrorException FAT.Fields.der!(p, dpdi, :z)
@test_throws ErrorException FAT.Fields.der!(p, dpdi, 3)

# we have not filled the boundary field by interpolation
dpdi = zeroScalarField(mesh(sim), 2)
FAT.Fields.der!(p, dpdi, :x)
@test all(dpdi.boundaryField .== 0.0)

# test gradient

# # test norm and inner product
# u.internalField[:] = 2.0
# @test inner(u, u) ≈ 8
# @test inner(u, u) == norm(u)^2

# # test zero constructor
# scalar = zeroField(0, mesh(sim))
# @test ndims(scalar.internalField) == 1
# @test ndims(scalar.boundaryField) == 1

# vector = zeroField(1, mesh(sim))
# @test ndims(vector.internalField) == 2
# @test ndims(vector.boundaryField) == 2

# tensor = zeroField(2, mesh(sim))
# @test ndims(tensor.internalField) == 3
# @test ndims(tensor.boundaryField) == 3


# # test velocity gradient 
# gradu = zeroField(2, mesh(sim))
# @time grad!(u, gradu)
# @time grad!(u, gradu)
# println(1000*minimum([@elapsed grad!(u, gradu) for i = 1:100]), " ms")


# @profile minimum([@elapsed grad!(u, gradu) for i = 1:1000])
# Profile.print()

# @code_llvm grad!(u, gradu)



# # # test gradient of a scalar field
# # uu = ScalarField(u.internalField[1:ncells(FAT.Fields.mesh(u))], 
# #  				 u.boundaryField[1:nboundaryfaces(FAT.Fields.mesh(u))], 
# #  				 FAT.Fields.mesh(u))
# # utmp = zeroField(1, mesh(sim))
# # println(1000*minimum([@elapsed grad!(uu, utmp) for i = 1:100]), " ms")

# # # test calculation of vorticity
# # gradu = zeroField(2, mesh(sim))
# # ω = zeroField(0, mesh(sim))

# # import PyPlot: tricontourf, show, pygui, scatter, figure, clf, savefig, plot

# # x, y, z = cellcentres(mesh(sim))

# # for (i, v) in enumerate(fields(sim, :U, 100.0:0.5:105.0))
# # 	ω = curl!(v, ω, gradu)
# # 	tricontourf(x, y, ω.internalField[:], -15:0.1:15, cmap=PyPlot.cm[:seismic])
# # 	name = @sprintf "figure-%04d.png" i
# # 	plot(x, y, "k.", ms=1, alpha=0.2)
# # 	savefig("vid/$name", dpi=800)
# # 	println(i)
# # 	close()
# # end
# # # tricontourf(x, y, gradu.internalField[:, 1, 1], -5:0.05:5)
