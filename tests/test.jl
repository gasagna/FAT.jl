immutable PolyHedralCell{T<:Real} 
	volume::T                   
end

immutable HexaMesh{T<:Real}
	cells::Vector{PolyHedralCell{T}}
end

const cells = PolyHedralCell{Float64}[PolyHedralCell{Float64}(rand()) for i = 1:1000000]
const mesh = HexaMesh(cells)								   

function test(m)
	a = work(m)
	for i = 2:100
		a += work(m)
	end
	a
end

function work(m)
	I = 0.0
	@inbounds begin 
		for (i, cell) in enumerate(m.cells)
			I += cell.volume
		end
	end
	I
end


function testfast(m)
	a = workfast(m)
	for i = 2:100
		a += workfast(m)
	end
	a
end

function workfast(m)
	I = 0.0
	for i in 1:length(m.cells)
		@inbounds I += m.cells[i].volume
	end
	I
end


println(minimum([@elapsed test(mesh) for i = 1:50]))
println(minimum([@elapsed testfast(mesh) for i = 1:50]))
