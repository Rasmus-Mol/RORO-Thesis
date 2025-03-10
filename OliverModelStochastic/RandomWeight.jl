# Code to generate random weight cargo
using Random
using Distributions

# Randomly change n cargos weight
function generate_random_weight_cargo(n_cargo::Int, weight_range::Tuple{Float64, Float64})
    Random.seed!(1234)
    cargo = Vector{Cargo}(undef, n_cargo)
    for i in 1:n_cargo
        cargo[i] = Cargo(
            id=i,
            cargo_type_id=rand(1:5),
            weight=rand(Uniform(weight_range[1], weight_range[2])),
            length=rand(Uniform(2.0, 10.0)),
            width=rand(Uniform(2.0, 10.0)),
            height=rand(Uniform(2.0, 10.0)),
            hazardous=rand() < 0.1
        )
    end
    return cargo
end



