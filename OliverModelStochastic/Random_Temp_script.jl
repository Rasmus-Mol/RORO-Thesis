@unpack vessel, slots, cargo = problem

deck1 = filter(x -> x.deck_id == 1, slots)
deck2 = filter(x -> x.deck_id == 2, slots)
deck3 = filter(x -> x.deck_id == 3, slots)

types1 = [filter(x -> x.cargo_type_id == i,deck1) for i in 1:4]
length(types1[1])
length(types1[2])
length(types1[3])
length(types1[4])

types2 = [filter(x -> x.cargo_type_id == i,deck2) for i in 1:4]
length(types2[1])
length(types2[2])
length(types2[3])
length(types2[4])

types3 = [filter(x -> x.cargo_type_id == i,deck3) for i in 1:4]
length(types3[1])
length(types3[2])
length(types3[3])
length(types3[4])
