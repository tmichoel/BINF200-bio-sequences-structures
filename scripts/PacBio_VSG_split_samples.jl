using DrWatson
@quickactivate "BINF200"

using BioSequences
using FASTX

# Read all filtered reads in a vector of FASTA records
vsgseq = Vector{FASTARecord}();
FASTAReader( open(datadir("exp_raw", "PacBio_VSG", "PacBio_VSG_filtered_reads.fasta")) ) do reader
    for record in reader
        push!(vsgseq, record)
    end
end

# Parse sample names
parsed_descriptions = split.(description.(vsgseq),"/");
sample_names = map(x -> x[1], parsed_descriptions)

# For each sample name, find the corresponding records and write them to a new file
for nm = String.(unique(sample_names))
    records = vsgseq[contains.(description.(vsgseq), nm)];
    FASTAWriter(open(datadir("exp_pro", "PacBio_VSG", "PacBio_VSG_filtered_reads_" * nm * ".fasta"), "w+")) do writer
        for record in records
            write(writer, record)
        end
    end
end
