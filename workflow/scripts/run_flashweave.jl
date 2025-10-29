#!/usr/bin/env julia

# Usage: julia flashweave.jl --data singlem_genus_percentage.tsv --meta meta.tsv --alpha 0.001 --n_obs_min -1 --transposed --output output_prefix
# Simple example: julia flashweave.jl --data 967med_MAGs-coverM-relativeabundance.tsv --transposed

using ArgParse
using FlashWeave

# Set command line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--data"
        help = "Path to input data file (e.g., singlem_genus_percentage.tsv)"
        required = true

        "--meta"
        help = "Path to meta data file (optional)"
        required = false

        "--alpha"
        help = "Alpha value for statistical significance"
    arg_type = Float64
    default = 0.001

        "--n_obs_min"
        help = "Minimum number of observations per variable"
    arg_type = Int
    default = -1

        "--transposed"
        help = "Whether data needs to be transposed"
        action = :store_true  # means if present, it will be true

        "--no_heterogeneous"
        help = "Whether to use heterogeneous mode (true/false, default: true)"
        action = :store_true  # means if present, it will be true

        "--outtable"
        help = "output edgelist file"

        "--outgraph"
        help = "output gml graph file"
    end

    return parse_args(s)
end

# Main program
function main()
    args = parse_commandline()
    # heterogeneous parameter：default true. --no_heterogeneous means false
    hetero_val = !get(args, "no_heterogeneous", false)

    println("FlashWeave network inference starting...")

    # Resolve paths and outputs (with sensible defaults)
    data_path = args["data"]
    outgraph = args["outgraph"]
    outtable = args["outtable"]

    # Working/output directory
    out_dir = dirname(outtable)
    if !isdir(out_dir)
        mkdir(out_dir)
    end

    # Define heterogeneous tag
    hetero_tag = hetero_val ? "heteroTRUE" : "heteroFALSE"

    # Read data file
    data = open(data_path, "r")
    data_lines = readlines(data)
    close(data)
    # Remove first header line
    data_lines = data_lines[2:end]
    # Remove lines starting with "Undefined" or "unmapped"
    data_lines = filter(line -> !startswith(line, "Undefined") && !startswith(line, "unmapped"), data_lines)
    # Save new data file
    data_new_path = joinpath(out_dir, "fw_input_abundance.tsv")
    data_new = open(data_new_path, "w")
    write(data_new, join(data_lines, "\n"))
    close(data_new)
    println("Input data file header removed and saved to: $data_new_path")
    println("Using data: $data_path")

    meta_data_path = get(args, "meta", nothing)

    if meta_data_path !== nothing
        println("Using meta data: $meta_data_path")
        netw_results = learn_network(
            data_new_path,
            meta_data_path,
            alpha=args["alpha"],
            n_obs_min=args["n_obs_min"],
            normalize=true,
            transposed=args["transposed"],
            sensitive=true,
            heterogeneous=hetero_val
        )
    else
        println("Running without meta data...")
        netw_results = learn_network(
            data_new_path,
            alpha=args["alpha"],
            n_obs_min=args["n_obs_min"],
            normalize=true,
            transposed=args["transposed"],
            sensitive=true,
            heterogeneous=hetero_val
        )
    end

    # print the learn_network command
    println("learn_network command: learn_network(\"$(data_new_path)\", $(meta_data_path === nothing ? "nothing" : "\"$(meta_data_path)\""), alpha=$(args["alpha"]), n_obs_min=$(args["n_obs_min"]), normalize=true, transposed=$(args["transposed"]), sensitive=true, heterogeneous=$(hetero_val))")

    # Save network outputs
    save_network(outgraph, netw_results)
    edgelist_path = "$(outgraph).edgelist"
    save_network(edgelist_path, netw_results)

    println("Network inference completed.")

    # Generate a more readable tsv files
    edgelist = open(edgelist_path, "r")
    edgelist_lines = readlines(edgelist)
    close(edgelist)
    # Remove first 2 lines
    edgelist_lines = edgelist_lines[3:end]
    # Add header
    header = "Source\tTarget\tFlashweave Weight"
    edgelist_lines = [header; edgelist_lines]
    # Save new edgelist
    edgelist_new_path = outtable
    edgelist_new = open(edgelist_new_path, "w")
    write(edgelist_new, join(edgelist_lines, "\n"))
    close(edgelist_new)

end

main()
