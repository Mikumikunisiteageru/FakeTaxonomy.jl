# src/RandomSpeciesNames.jl

module RandomSpeciesNames

using Unicode

function names2kmers(
		taxonnames::AbstractVector{<:AbstractString}, k::Integer))
	speciesnames = unique(Unicode.normalize.(
		join.(getindex.(split.(taxonnames), [1:2]), ' '), stripmark=true))
	strings = "^"^k .* speciesnames .* "\$"
	kmers = Dict{String, Dict{Char, Int}}()
	for string = strings
		for i = 1 : length(string) - k
			tmer = string[i:i+k-1]
			next = string[i+k]
			get!(get!(kmers, tmer, Dict{Char, Int}()), next, 0)
			kmers[tmer][next] += 1
		end
	end
	return kmers
end

function try_generate(k, kmers)
	str = "^" ^ k
	while ! endswith(str, "\$")
		nextdict = get(kmers, str[end-2:end], nothing)
		nothing === nextdict && return nothing
		chars = collect(keys(nextdict))
		freqs = getindex.([nextdict], chars)
		cufreqs = cumsum(freqs)
		randpos = rand(1:cufreqs[end])
		i = searchsortedfirst(cufreqs, randpos)
		str *= chars[i]
	end
	return str[k+1:end-1]
end

function generate(k, kmers)
	while true
		str = try_generate(k, kmers)
		nothing === str && continue
		count(" ", str) == 1 || continue
		return str
	end
end

end # module
