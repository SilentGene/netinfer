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
        default = 0.001

        "--n_obs_min"
        help = "Minimum number of observations per variable"
        default = -1

        "--transposed"
        help = "Whether data needs to be transposed"
        action = :store_true  # means if present, it will be true

        "--no_heterogeneous"
        help = "Whether to use heterogeneous mode (true/false, default: true)"
        action = :store_true  # means if present, it will be true

        "--output"
        help = "Prefix for output files"
    end

    return parse_args(s)
end

# Main program
function main()
    args = parse_commandline()
    # heterogeneous parameter：default true. --no_heterogeneous means false
    hetero_val = !get(args, "no_heterogeneous", false)

    println("FlashWeave network inference starting...")

    # Get current directory
    curr_dir = pwd()
    println("Working directory: $curr_dir")

    data_path = args["data"]

    # Check if input data file exists
    if !isfile(data_path)
        println("Error: Input data file not found: $data_path")
        return
    end

    # Define heterogeneous tag
    hetero_tag = hetero_val ? "heteroTRUE" : "heteroFALSE"
    
    # Create output folder
    base_name_no_ext = splitext(basename(data_path))[1]
    out_dir = joinpath(curr_dir, base_name_no_ext * "_flashweave" * "_" * hetero_tag)
    if !isdir(out_dir)
        mkdir(out_dir)
    end

    # Get alpha value as string
    alpha_str = string(args["alpha"])
    alpha_tag = "-p" * alpha_str
    if args["output"] == nothing
        out_prefix = joinpath(out_dir, "flashweave_" * base_name_no_ext * alpha_tag * "_" * hetero_tag)
    else
        out_prefix = joinpath(out_dir, "flashweave_" * args["output"] * alpha_tag * "_" * hetero_tag)
    end

    # Read data file
    data = open(data_path, "r")
    data_lines = readlines(data)
    close(data)
    # Remove first header line
    data_lines = data_lines[2:end]
    # Remove lines starting with "Undefined" or "unmapped"
    data_lines = filter(line -> !startswith(line, "Undefined") && !startswith(line, "unmapped"), data_lines)
    # Save new data file
    data_new_path = joinpath(out_dir, "$(base_name_no_ext)_remove1n.tsv")
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

    save_network("$(out_prefix).gml", netw_results)
    save_network("$(out_prefix).edgelist", netw_results)

    println("✅ Network inference completed.")
    println("➡ Results saved to: $(out_prefix).gml and $(out_prefix).edgelist")

    # Generate a more readable edgelist.tsv
    edgelist_path = "$(out_prefix).edgelist"
    edgelist = open(edgelist_path, "r")
    edgelist_lines = readlines(edgelist)
    close(edgelist)
    # Remove first 2 lines
    edgelist_lines = edgelist_lines[3:end]
    # Add a new column
    edgelist_lines = map(line -> line * "\t" * split(line, "\t")[1] * "-->" * split(line, "\t")[2], edgelist_lines)
    # Add header
    header = "Source\tTarget\tFlashweave Weight\tSource-->Target"
    edgelist_lines = [header; edgelist_lines]
    # Save new edgelist
    edgelist_new_path = "$(out_prefix)_edges_info.tsv"
    edgelist_new = open(edgelist_new_path, "w")
    write(edgelist_new, join(edgelist_lines, "\n"))
    close(edgelist_new)
    println("➡ Edgelist saved to: $(edgelist_new_path)")

end

main()
