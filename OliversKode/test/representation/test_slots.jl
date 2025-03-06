using Test
using XLSX
using DataFrames  # Add this import
include("../../src/representation/slot.jl")

@testset "Slot Tests" begin
    @testset "Slot Creation" begin
        slot = Slot(
            id=1,
            loadmaster_id=100,
            deck_id=1,
            cargo_type_id=2,
            lcg=133.0,
            tcg=5.0,
            vcg=14.5,
            refrigerated=false,
            length=4.6,
            width=4.5,
            on_deck=true
        )
        
        @test slot.id == 1
        @test slot.loadmaster_id == 100
        @test slot.deck_id == 1
        @test slot.on_deck == true
    end

    @testset "SlotCollection" begin
        slots = generate_test_slots(10, 2)
        @test length(slots) == 10
        @test slots[1].id == 1
        @test slots[end].id == 10
        
        deck_1_slots = get_slots_by_deck(slots, 1)
        @test all(slot -> slot.deck_id == 1, deck_1_slots)
        
        deck_slots = get_deck_slots(slots)
        @test all(slot -> slot.on_deck == true, deck_slots)
    end
end