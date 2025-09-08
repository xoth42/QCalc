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
	# function U(θ,ϕ,λ)
	# 	return [
	# 		exp(-1/2im * (ϕ + λ))*cos(θ/2)	-exp(-1/2im * (ϕ - λ))*sin(θ/2);
	# 		exp(1/2im * (ϕ - λ))*sin(θ/2)	exp(1/2im * (ϕ + λ))*cos(θ/2)
	# 	]
	# end

	# This one  deals with the global phase!!
	# OpenQASM 3.0 U gate: exp(i(ϕ+λ)/2) Rz(ϕ)Ry(θ)Rz(λ)
	function U(θ,ϕ,λ)
		return [
			cos(θ/2)	-exp(1im * λ)*sin(θ/2);
			exp(1im * ϕ)*sin(θ/2)	exp(1im * (ϕ + λ))*cos(θ/2)
		]
	end
		
end

# ╔═╡ 59d7a44f-c33b-4051-9a33-39bb4701fcd8
md"set U gate parameters. U(π/2, 0, π) = H"

# ╔═╡ 56f7c249-af4b-4ed2-a341-408e87f8bb6b
md"""
* θ: $(@bind θ PlutoUI.Slider(ranges, default=pi/2,show_value=show_pi))
* ϕ: $(@bind ϕ PlutoUI.Slider(ranges, default=0,show_value=show_pi)) 
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
	function comp_out(z) # Turn a complex number value into a pretty output like 1 + 1im, removing the imaginary part if needed
		# if number is less than err in magnitude
		if norm(z) < err
			return "0"
		end
		
		str = ""
		r = real(z)
		i = imag(z)
		if r > err || r < -err # real part is nonzero
			str *= @sprintf "%1.2f" r
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
			str *= @sprintf "%1.2fim" i
		end
		return str
	end

	function pretty_matrix(m)
		cout = comp_out.(m)
		return "[\t$(cout[1])\t$(cout[3])\n\t$(cout[2])\t$(cout[4])\t]"
	end

end

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
	
	for g in unk
	    # println("Translation todo:")
		# println(g)
		dump(g)
		println(g)
		matrix  = g |> Operator
		println(pretty_matrix(matrix.data), "\n\n")
		
	end
	

end

# ╔═╡ 2cf46b08-0b96-415c-b345-df503ec14af7
md"Resulting matrix:"

# ╔═╡ 206f6d52-4c8e-499d-8d2f-ca0e7bbe5417
begin
	u = U(θ,ϕ,λ)
	println(pretty_matrix(u))
end

# ╔═╡ 3fbff09a-db52-44a9-aa5e-812059d5875a
md"Now trying to find other gates from this basis (brute force)"

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

# ╔═╡ ea29edf4-bc01-487d-be1e-a74321c34686
md"Target to gate to find representation for:"

# ╔═╡ ab48c4af-80fc-4a49-8a92-9cf7eda20ef2
begin
	stab_op = targets[t]
	op = stab_op |> Operator
	matrix = op.data
	comp_out.(matrix)
end


# ╔═╡ e67c9298-5175-4692-93c9-cc5c7b6e33fa
md"QCliff representation:"

# ╔═╡ 300ee489-b091-418e-a13f-28832e33c3e2
dump(stab_op)

# ╔═╡ a514a4bf-3b40-4d1a-b9a5-7ee6374a5f52
stab_op

# ╔═╡ f5ba370b-82be-4f8a-8d79-e0cce14d1efc
begin
	function find_this_gate(gate,grain=pi/2)
		target_range = 0:grain:4pi
		target_ranges = (target_range,target_range,target_range)
		target = gate # set target to the unknown gate
		return find_3params_for_gate(target_ranges,U,target)
	end
	out = find_this_gate(matrix);nothing
end

# ╔═╡ 67477359-634f-49ad-b32d-7190742ad45e
md"Resulting params:  U(θ,ϕ,λ)"

# ╔═╡ cb220082-061e-4475-aad3-7b080a71a201
begin
	
	function pretty_param(param; err=1e-10) # the param will come in (ex: 0, 1.5708, 4.71239...)
		# and we need to return a readable factor of pi
		if param < 0
			return param # should not be negative
		end

		# it is under error cutoff (basically zero)
		if param < err 
			return 0
		end
		
		pis = param / pi # initial factors of pi

		# re-use another func for cleanup
		pis = comp_out(pis + 0im)
				
		return "$(pis)π"
	end

	function pparams(params;init="U")
		s =  "$(init)("
		for i in 1:2
			s *= string(pretty_param(params[i]))
			s *= ", "
		end
		s *=  string(pretty_param(params[3]) )
		s *= ")"
		return s
	end
end

# ╔═╡ 234025cd-9ffa-4958-8a55-8d6cc6e1b144
out == nothing ? "not found" : pparams(out[1])

# ╔═╡ 3473cfca-190f-45f8-968a-0613ee2042ac
md"For this gate:"

# ╔═╡ c4cb1323-42ff-4d38-b76a-e3e8b8928823
begin
	# println(out[2])
	println(pretty_matrix(out[2]))
end

# ╔═╡ Cell order:
# ╠═d62e8a66-84d9-11f0-13c4-c12995195713
# ╟─68f5b298-8de1-4a6e-a2b9-fea70ddf9c2a
# ╟─c82d3aa7-2a32-4b34-839d-fbd78464094c
# ╟─b11d0525-1bfc-458d-a46a-dd15992b3964
# ╟─59d7a44f-c33b-4051-9a33-39bb4701fcd8
# ╟─56f7c249-af4b-4ed2-a341-408e87f8bb6b
# ╟─a5980897-5328-4d11-a22c-f4cb29a672c3
# ╟─2cf46b08-0b96-415c-b345-df503ec14af7
# ╟─206f6d52-4c8e-499d-8d2f-ca0e7bbe5417
# ╟─3fbff09a-db52-44a9-aa5e-812059d5875a
# ╠═dd1dd101-4e38-4b14-8021-2752ec46ac54
# ╟─aa03cae3-ec49-42c1-a667-bc0d55a2d6ad
# ╟─3a58c739-f72b-4bd8-a596-07e822281ad7
# ╟─ea29edf4-bc01-487d-be1e-a74321c34686
# ╟─ab48c4af-80fc-4a49-8a92-9cf7eda20ef2
# ╟─e67c9298-5175-4692-93c9-cc5c7b6e33fa
# ╟─300ee489-b091-418e-a13f-28832e33c3e2
# ╟─a514a4bf-3b40-4d1a-b9a5-7ee6374a5f52
# ╟─45b89936-efcf-4f9a-b76a-ec0f9e07dc9b
# ╟─f5ba370b-82be-4f8a-8d79-e0cce14d1efc
# ╟─67477359-634f-49ad-b32d-7190742ad45e
# ╟─cb220082-061e-4475-aad3-7b080a71a201
# ╠═234025cd-9ffa-4958-8a55-8d6cc6e1b144
# ╟─3473cfca-190f-45f8-968a-0613ee2042ac
# ╟─c4cb1323-42ff-4d38-b76a-e3e8b8928823
