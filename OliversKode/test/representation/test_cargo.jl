using Test

include("../../src/representation/cargo.jl")

@testset "Cargo Tests" begin
    @testset "CargoTypeInfo" begin
        car_info = get_type_info(Car)
        @test car_info.length == 4.5
        @test car_info.width == 1.8
        @test car_info.height == 1.6
        
        container20_info = get_type_info(Container20)
        @test container20_info.length == 20.0
        @test container20_info.width == 2.4
        @test container20_info.height == 2.6
    end

    @testset "Cargo Construction" begin
        cargo = Cargo(
            id="TEST001",
            cargo_type=Car,
            weight=2.0,
            loading_port=1,
            discharge_port=3,
            hazardous=2,
            refers=true
        )
        
        @test cargo.id == "TEST001"
        @test cargo.cargo_type == Car
        @test cargo.weight == 2.0
        @test cargo.priority == 1  # default value
        @test cargo.hazardous == 2
        @test cargo.refers == true
        
        # Test dimension getters
        @test get_length(cargo) == 4.5
        @test get_width(cargo) == 1.8
        @test get_height(cargo) == 1.6
    end

    @testset "CargoCollection" begin
        cargo1 = Cargo(id="T1", cargo_type=Car, weight=2.0, loading_port=1, discharge_port=2, hazardous=1, refers=false)
        cargo2 = Cargo(id="T2", cargo_type=Truck, weight=10.0, loading_port=2, discharge_port=3, hazardous=0, refers=true)
        
        collection = CargoCollection([cargo1, cargo2])
        @test length(collection.items) == 2
        @test collection.total_weight ≈ 12.0
    end

    @testset "StructArray Operations" begin
        cargos = generate_test_collection(5)
        @test length(cargos) == 5
        
        # Test filtering
        cars = get_cargo_by_type(cargos, Car)
        @test all(x -> x == Car, cars.cargo_type)
        
    end

    @testset "Random Generation" begin
        Random.seed!(123)  # For reproducibility
        cargo = generate_random_cargo("TEST001")
        @test cargo.id == "TEST001"
        @test cargo.weight > 0
        @test cargo.loading_port != cargo.discharge_port
        @test cargo.hazardous in 0:3
        @test cargo.refers in [true, false]
        
        # Test dimension validity
        @test get_length(cargo) > 0
        @test get_width(cargo) > 0
        @test get_height(cargo) > 0
        
        collection = generate_test_collection(10)
        @test length(collection) == 10
    end
end