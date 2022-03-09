@parameters p0 p1 p2 p3
p = [p1 => 1.0, p2 => 2.0, p3 => 1.0]

p_names = first.(p)
p_lookup_table = []

for i in range(1, length(p_names))
    push!(p_lookup_table, string(p_names[i]))
end

println(p_lookup_table)