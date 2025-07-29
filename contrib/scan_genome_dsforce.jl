using BioSequences
using FastaIO
using DataFrames
using CSV
using ArgParse
using DSForces

"""
Read the contig with the specified name from the genome FASTA file
"""
function ReadSequence(ifn::String, contig::String)
    c_seq = ""
    FastaReader(ifn) do fr
        for (name, seq) in fr
            S = split(name)
            if S[1] == contig
                c_seq = seq
            end
        end
    end
    @assert length(c_seq) > 0 "Contig $(contig) not found, please try again with a different contig name."
    return c_seq
end


"""
Command line parser
"""
function ParseCommandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "infile"
            help = "Input file with the genome assembly."
            required = true
        "contig"
            help = "Contig to be considered for DS force computation."
            required = true
        "slice_start"
            help = "Start position (1-based) on the contig sequence."
            arg_type = Int
            required = true
        "slice_end"
            help = "End position (1-based) on the contig sequence"
            arg_type = Int
            required = true
        "--outfile"
            help = "Name of the output file."
            default = ""
        "--window_length"
            help = "Length of the sliding window."
            arg_type = Int
            default = 3000
    end
    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
    # argument parsing
    p_args = ParseCommandline()
    infile = p_args["infile"]
    contig = p_args["contig"]
    slice_start = p_args["slice_start"]
    slice_end = p_args["slice_end"]
    window_length = p_args["window_length"]
    outfile = p_args["outfile"]

    if outfile == ""
        outfile = "$(contig)_$(slice_start)-$(slice_end).csv"
    end

    # checks
    @assert slice_start >= 1 "'Slice start' parameter must be positive."

    # prepare seq
    contig_seq = LongDNA{4}(ReadSequence(infile, contig))
    seq = contig_seq[slice_start:min((slice_end + window_length - 1), length(contig_seq))]
    if length(seq) < window_length
        io = open(outfile, "w")
        close(io)
        exit(0)
    end

    # compute DS
    DSFs, range1s, range2s = ComputeDSForce(seq,
        return_LCS_positions=true, sliding_window_length=window_length)
    
    # prepare output as dataframe
    window_starts = [slice_start + i for i in 0:length(seq) - window_length]
    window_ends = window_starts .+ (window_length-1)
    corr_rangeAs = [r1 .+ (slice_start-1) for r1 in range1s]
    corr_rangeBs = [r2 .+ (slice_start-1) for r2 in range2s]
    res_df = DataFrame(["contig" => [contig for _ in 1:length(DSFs)], 
                        "window start:end" => [r1:r2 for (r1, r2) in zip(window_starts, window_ends)], 
                        "DS_force" => round.(DSFs, digits=6), 
                        "seqA start:end" => corr_rangeAs, 
                        "seqB start:end" => corr_rangeBs])
    # write output
    CSV.write(outfile, res_df)
end
