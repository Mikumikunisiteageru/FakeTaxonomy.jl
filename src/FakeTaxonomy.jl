# src/FakeTaxonomy.jl

module FakeTaxonomy

export FakeTaxonomer, nominate

deaccent(str::AbstractString) = Base.Unicode.normalize(str; stripmark=true)

function specify(name::AbstractString)
	words = split(name, ' ')
	length(words) < 2 && 
		throw(ArgumentError("`$name` is not a legal (infra)species name!"))
	return join(words[1:2], ' ')
end

struct FakeTaxonomer
	k::Int
	kmers::Dict{String,Dict{Char,Int}}
end
function FakeTaxonomer(names::AbstractVector{<:AbstractString}, k::Integer)
	spnames = unique(deaccent.(specify.(names)))
	strings = "^" ^ k .* spnames .* "\$"
	kmers = Dict{String, Dict{Char, Int}}()
	for string = strings
		for i = 1 : length(string) - k
			tmer = string[i:i+k-1]
			next = string[i+k]
			get!(get!(kmers, tmer, Dict{Char, Int}()), next, 0)
			kmers[tmer][next] += 1
		end
	end
	return FakeTaxonomer(k, kmers)
end
function FakeTaxonomer(str::AbstractString)
	ii = isnumeric.(collect(str))
	ii[1] && throw(ArgumentError("Input string malformed!"))
	st = findall(xor.(ii[1:end-1], ii[2:end]))
	kp1 = st[1]
	k = kp1 - 1
	kmers = Dict{String, Dict{Char, Int}}()
	key = "-" ^ kp1
	for (i, (s, t)) = enumerate(zip([1; st.+1], [st; length(str)]))
		if isodd(i)
			key = key[1:k-t+s] * str[s:t]
		else
			get!(kmers, key[1:k], Dict{Char, Int}())
			kmers[key[1:k]][key[kp1]] = parse(Int, str[s:t])
		end
	end
	return FakeTaxonomer(k, kmers)
end

function Base.string(ftn::FakeTaxonomer)
	kmers = ftn.kmers
	code = ""
	current = "!" ^ length(first(keys(kmers)))
	for key = sort!(collect(keys(kmers)), by = key -> replace(key, '^'=>'\0'))
		code *= key[findfirst(map(!=, current, key)):end]
		current = key
		for char = sort!(collect(keys(kmers[key])))
			code *= char * string(kmers[key][char])
		end
	end
	return code
end

function sample(ftn::FakeTaxonomer, start::AbstractString="")
	str = "^" ^ ftn.k * start
	haskey(ftn.kmers, str[end-ftn.k+1:end]) || 
		throw(ArgumentError("No possible species name starts with `$start`!"))
	while ! endswith(str, "\$")
		nextdict = get(ftn.kmers, str[end-ftn.k+1:end], nothing)
		nothing === nextdict && return nothing
		chars = collect(keys(nextdict))
		freqs = getindex.([nextdict], chars)
		cufreqs = cumsum(freqs)
		randpos = rand(1:cufreqs[end])
		i = searchsortedfirst(cufreqs, randpos)
		str *= chars[i]
	end
	return str[ftn.k+1:end-1]
end

function nominate(ftn::FakeTaxonomer, start::AbstractString="")
	while true
		str = sample(ftn, start)
		nothing === str && continue
		count(' ', str) == 1 && return str
	end
end

end # module
