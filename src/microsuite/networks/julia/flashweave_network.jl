using FlashWeave

if length(ARGS) != 4
    error("usage: flashweave_network.jl <sample-feature-table.tsv> <output.edgelist> <sensitive> <heterogeneous>")
end

table_path = ARGS[1]
output = ARGS[2]
sensitive = lowercase(ARGS[3]) in ["1", "true", "t", "yes"]
heterogeneous = lowercase(ARGS[4]) in ["1", "true", "t", "yes"]

network = learn_network(table_path, sensitive=sensitive, heterogeneous=heterogeneous)
save_network(output, network)
