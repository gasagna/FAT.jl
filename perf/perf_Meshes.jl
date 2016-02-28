using FAT.Meshes

casedir = "/Users/davide/Dropbox/Codes/FAT/tests/ldc_large"

@time m = HexaMesh(casedir)
@time m = HexaMesh(casedir)
println(minimum([@elapsed HexaMesh(casedir) for i = 1:10]))

@profile m = HexaMesh(casedir)

Profile.print()