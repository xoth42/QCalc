### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d62e8a66-84d9-11f0-13c4-c12995195713
begin 
	using Pkg
	# using QuantumClifford, BPGates
	# using BPGates:toQCcircuit
	# using QuantumOpticsBase
	
	Pkg.activate("../")
	# Pkg.add("ProgressLogging")
end

# ╔═╡ c82d3aa7-2a32-4b34-839d-fbd78464094c
begin 
	using Format,Printf
	function show_pi(n)
		if n == 0
			return "0"
		else
			return "$(format(round(n/pi,digits=2), mixedfraction=true)) π"
		end
	end
	ranges = 0:pi/4:4pi;nothing
	# todo: fraction output
	# format(1//2,mixedfraction=true)
end

# ╔═╡ b11d0525-1bfc-458d-a46a-dd15992b3964
using PlutoUI

# ╔═╡ acfdde1b-e2b1-4f8d-8543-6eda875c1a27
using ProgressLogging

# ╔═╡ dd1dd101-4e38-4b14-8021-2752ec46ac54
begin
	using QuantumClifford, BPGates
	using BPGates:toQCcircuit
	using QuantumOpticsBase
	
	function to_qasm_example(op::CNOTPerm)
	    unknowns = Set() 
	    # here are the easy ops to turn to qasm (simplified but similar to how Quantikz looks at the gates)
	    easy = [tHadamard,tPhase,tId1,tSWAP,tCPHASE,tCNOT,tSWAP,tSWAP*tCNOT*tSWAP]
	
	    # Get cliff ops
	    for cliffordOp in toQCcircuit(op)
	        !isa(cliffordOp, SparseGate) && continue 
	        g = cliffordOp.cliff
	
	        # to illustrate the point, just ignore the easy ones
	        if !(g in easy)
	            push!(unknowns, g)
	            # println("Translation todo:\n$g\n$(g |> Operator)\non $(cliffordOp.indices)")
	        end
	    end
	    return unknowns
	end
	
	unk = Set()
	for i in 1:500
	    union!(unk,to_qasm_example(rand(CNOTPerm,1,2)))
	end
	
	# for g in unk
	#     # println("Translation todo:")
	# 	println(g)
	# 	println(g |> Operator)
	# end
	

end

# ╔═╡ 45b89936-efcf-4f9a-b76a-ec0f9e07dc9b
begin
	using ProgressLogging:@progress
	debug = false
	function find_3params_for_gate(param_ranges,gate,target;err=1e-1)
		tot = 0
		# for 3 params only, gate of 2 by 2
		# param_ranges: [start:step:end, start:step:end, ...]
		# For each parameter to input into the gate, provide a range that that parameter is allowed. 
		# gate will be called gate(param1,param2,...)
		# until the output of the gate call is equivalent to the provided target within an error, ;err = 1e-10
		t_r = real.(target)
		t_i = imag.(target)
		function is_within(matrix,t_r,t_i,err) 
			for i in (1:4)
				# if difference is greater than error, exit
				# real
				r = real(matrix[i]) 
				diff_r = abs(r - t_r[i])
				debug && println(diff_r)
				if diff_r > err
					return false
				end
				# imag
				img = imag(matrix[i]) 
				diff_i = abs(img - t_i[i])
				if diff_i > err
					return false
				end
			
			end
			return true
		end

		found = (guess)-> is_within(guess,t_r,t_i,err)
					
		@progress for i in param_ranges[1]
			for j in param_ranges[2]
				for k in param_ranges[3]
					tot += 1
					debug && println("$i\t$j\t$k")
					g = gate(i,j,k)
					debug && println(g)
					if found(g)
						println("Tries: $tot")
						return ((i,j,k),g)
					end
				end
			end
		end
		println("Tries: $tot")
		return nothing
	end
end

# ╔═╡ 68f5b298-8de1-4a6e-a2b9-fea70ddf9c2a
begin
	# OpenQASM 2.0 U gate: Rz(ϕ)Ry(θ)Rz(λ)
	function U(θ,ϕ,λ)
		return [
			exp(-1/2im * (ϕ + λ))*cos(θ/2)	-exp(-1/2im * (ϕ - λ))*sin(θ/2);
			exp(1/2im * (ϕ - λ))*sin(θ/2)	exp(1/2im * (ϕ + λ))*cos(θ/2)
		]
	end
end

# ╔═╡ 56f7c249-af4b-4ed2-a341-408e87f8bb6b
md"""
* θ: $(@bind θ PlutoUI.Slider(ranges, default=pi,show_value=show_pi))
* ϕ: $(@bind ϕ PlutoUI.Slider(ranges, default=pi,show_value=show_pi)) 
* λ: $(@bind λ PlutoUI.Slider(ranges, default=pi,show_value=show_pi))
"""

# ╔═╡ a5980897-5328-4d11-a22c-f4cb29a672c3
begin
	err = 1e-10
	function norm(z::ComplexF64)
		r = real(z)
		i = imag(z)
		return sqrt(r^2 + i^2)
	end
	function comp_out(z)
		# if number is less than err in magnitude
		if norm(z) < err
			return "0"
		end
		
		str = ""
		r = real(z)
		i = imag(z)
		if r > err || r < -err # real part is nonzero
			str *= @sprintf "%1.3f" r
		end
	 	if i > err || i < -err # i is nonzero
			if str != "" # there was a real number
				str *= " " # add a space
				if i > 0 # if there was a real and the imag is pos, add a +
					str *= "+ " 
				# else
				# 	str *= "- "
				end
			# elseif i < 0 # there was no real number, imag part is negative, need a -
			# 	str *= "- "
			end
			# now add the imag num
			str *= @sprintf "%1.3fim" i
		end
		return str
	end
end

# ╔═╡ 206f6d52-4c8e-499d-8d2f-ca0e7bbe5417
begin
	u = U(θ,ϕ,λ)
	m = comp_out.(u)
end

# ╔═╡ eef4145b-fd7d-46ce-a7bd-a859c895014a
u*sqrt(2)

# ╔═╡ aa03cae3-ec49-42c1-a667-bc0d55a2d6ad
begin
	targets = []
	for n in unk
		push!(targets,n)
	end
	"targets: $(length(targets))"
	
end

# ╔═╡ 3a58c739-f72b-4bd8-a596-07e822281ad7
length(targets) > 1 && md"
* target: $(@bind t PlutoUI.Slider(1:1:length(targets), default=1,show_value=true))
"

# ╔═╡ ab48c4af-80fc-4a49-8a92-9cf7eda20ef2
begin
	stab_op = targets[t]
	op = stab_op |> Operator
	matrix = op.data
	comp_out.(matrix)
end


# ╔═╡ 27b1e0ac-a2c5-4952-8e11-f5ad407e12ff
matrix

# ╔═╡ f5ba370b-82be-4f8a-8d79-e0cce14d1efc
begin
	target_range = 0:pi/32/3:4pi
	target_ranges = (target_range,target_range,target_range)
	# target = matrix
	target = (1/sqrt(2))*[
		1+0im	1+0im;
		1+0im	-1+0im
	]
	out = find_3params_for_gate(target_ranges,U,target)
end

# ╔═╡ 234025cd-9ffa-4958-8a55-8d6cc6e1b144
out == nothing ? "not found" : out

# ╔═╡ Cell order:
# ╠═d62e8a66-84d9-11f0-13c4-c12995195713
# ╠═68f5b298-8de1-4a6e-a2b9-fea70ddf9c2a
# ╟─c82d3aa7-2a32-4b34-839d-fbd78464094c
# ╟─b11d0525-1bfc-458d-a46a-dd15992b3964
# ╠═56f7c249-af4b-4ed2-a341-408e87f8bb6b
# ╟─a5980897-5328-4d11-a22c-f4cb29a672c3
# ╠═206f6d52-4c8e-499d-8d2f-ca0e7bbe5417
# ╠═eef4145b-fd7d-46ce-a7bd-a859c895014a
# ╠═acfdde1b-e2b1-4f8d-8543-6eda875c1a27
# ╟─dd1dd101-4e38-4b14-8021-2752ec46ac54
# ╟─aa03cae3-ec49-42c1-a667-bc0d55a2d6ad
# ╟─3a58c739-f72b-4bd8-a596-07e822281ad7
# ╟─ab48c4af-80fc-4a49-8a92-9cf7eda20ef2
# ╟─27b1e0ac-a2c5-4952-8e11-f5ad407e12ff
# ╟─45b89936-efcf-4f9a-b76a-ec0f9e07dc9b
# ╠═f5ba370b-82be-4f8a-8d79-e0cce14d1efc
# ╠═234025cd-9ffa-4958-8a55-8d6cc6e1b144
