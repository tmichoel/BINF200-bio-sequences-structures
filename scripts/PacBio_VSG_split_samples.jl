using DrWatson
@quickactivate "BINF200"

using BioSequences
using FASTX

samplenames = ["3_1", "3_2", "3_3", "3_4", "3_5", "6_1", "6_2", "6_3", "6_4", "6_5", "10_1", "10_2", "10_3", "10_4", "10_5", "12_1", "12_2", "12_3", "12_4", "12_5"]

vsgseq = Vector{FASTARecord}();

FASTAReader( open(datadir("exp_raw", "PacBio_VSG", "PacBio_VSG_filtered_reads.fasta")) ) do reader
    for record in reader
        push!(vsgseq, record)
    end
end
