

vectors = Dict{Int, Vector{Vector{Int}}}()

for _=1:1000
  v = Int[]
  for _ = 1:1000
    push!(v, int(rand()*500000))
  end
  fixed_vector = sort(unique(v))
  vectors[length(fixed_vector)] = push!(get(vectors, length(fixed_vector), Vector{Int}[]), fixed_vector)
end
vectors = vectors[1000]


function test_fd(f::Function)
  #
  for i=1:length(vectors)
    for j=1:length(vectors)
      a = f(vectors[i], vectors[j])
    end
  end
end

@time test_fd(fast_delete)
@time test_fd(setdiff)