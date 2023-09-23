### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9ddcb453-89eb-4642-b313-40d9f2548d7c
using DrWatson

# ╔═╡ 4d7d51cc-2751-4a98-bc27-b2e0f21990ec
# ╠═╡ show_logs = false
@quickactivate "BINF200"

# ╔═╡ f1c539ca-6ca9-40f2-9a46-cdccef31c838
begin
	using BioSequences
	using FASTX
	using DataFrames
	using CSV
	using Plots
	using Statistics
	using StatsPlots
	using LaTeXStrings
	using PlutoUI
	using Downloads
end

# ╔═╡ e355a610-4af8-11ee-29f9-d916ed9c2903
md"""
# Application of long read sequencing to determine expressed antigen diversity in *Trypanosoma brucei* infections

## Setup the environment
"""

# ╔═╡ 731d038c-78e0-4740-8f80-174e1c7d70b3
md"""
## Select your sample
"""

# ╔═╡ 0b807b4f-7634-468c-8619-df0348981882
@bind sample MultiSelect(["balbc_3_0", "balbc_3_1", "balbc_3_2", "balbc_3_3", "balbc_3_4", "balbc_6_0", "balbc_6_1", "balbc_6_2", "balbc_6_4", "balbc_6_5", "balbc_10_0", "balbc_10_1", "balbc_10_2", "balbc_10_4", "balbc_10_5","balbc_12_1", "balbc_12_2", "balbc_12_3", "balbc_12_4", "balbc_12_5"])

# ╔═╡ f3637e6c-2984-47e5-a0c0-394f3bdd2cae
begin
	filtered_reads = Vector{FASTARecord}();
   FASTAReader( open(datadir("PacBio_VSG", "filtered_reads", "PacBio_VSG_filtered_reads_" * sample[1] * ".fasta")) ) do reader
	      for record in reader
	         push!(filtered_reads, record);
	      end
	   end
end

# ╔═╡ b9ecaf9c-2c40-494d-83ea-bc9a4bb10056
md"""
## Tasks

### Show that you understand the biological experiment and data

To count the number of sequences in the input file, we read the FASTA file with input sequences into a vector of FASTA records. The length of the vector is the number of sequences: **$(length(filtered_reads))**
"""

# ╔═╡ d187d9f9-e756-484f-acfe-4accf7c59552
md"""
### BLAST the sample sequences against the reference VSG database

The correct blast command is

```
blastn  -query $qf -db $dbf -out $of  -max_target_seqs 1 -outfmt '6 qseqid sseqid score bitscore evalue qlen slen length sstart send nident mismatch gaps positive'
```

where

- `qf` is the input file name PacBio\_VSG\_filtered\_reads\_ $(sample).fasta
- `of` is the output file name PacBio\_VSG\_filtered\_reads\_blastn\_$(sample).txt
- `dbf` is the database file name TREU927-v26\_VSGTranscripts/TREU927-v26\_VSGTranscripts.fasta
"""

# ╔═╡ 6314b615-22a9-4aa4-8b59-614281426e80
begin
	df = DataFrame(CSV.File( datadir("PacBio_VSG", "blastn", "PacBio_VSG_filtered_reads_blastn_" * sample[1] * ".txt"), header=false ));
	   rename!(df, ["qseqid", "sseqid", "score", "bitscore", "evalue", "qlen", "slen", "length", "sstart", "send", "nident", "mismatch", "gaps", "positive"]);
end

# ╔═╡ 47ef717f-fd1f-47c9-a080-894f8eda0149
md"""
### 

### Analyze the BLAST output

Answer the following questions / perform the following tasks:

**1. What is the number of sequences in *your* blast output file? Is it the same as in the input file?**

Import the blast results in a DataFrame. The column names can be constructed from the blastn command. The number of rows of the DataFrame is the number of sequences: **$(nrow(df))**. 
"""

# ╔═╡ d927f106-35c7-4648-8c5b-adc0fd6c6c60
df_qseqid_count = combine(groupby(df[:,1:2], :qseqid), nrow, :sseqid => x -> length(unique(x)))

# ╔═╡ b6b50963-f5ae-4efd-aa42-b156baaefaed
md"""
Despite setting the option `-max_target_seq 1` in the `blastn` command, this number may be *greater* than the number of input sequences. The reason is that sometimes for *one* query sequence, *two* separate alignments are found against the *same* subject sequence. To show this, we compute for each unique query id, its number of alignments (number of rows it appears in in `df`) and the number of unique subject ids in those rows. Grouping the DataFrame by the query sequence id and using the split-apply-combine strategy, we can count the number of subject sequence ids per query in one line. After this computation, we see that now the number of unique query sequences (**$(nrow(df_qseqid_count))**) is indeed equal to the number of input sequences (**$(length(filtered_reads))**).
"""

# ╔═╡ 8db859f2-2794-431c-9a8a-221832296d11
md"""
We can verify that sometimes two alignments were found per query sequence:
"""

# ╔═╡ fccf0588-ca40-46ec-ad64-619eaa4c184a
unique(df_qseqid_count.nrow)

# ╔═╡ 3c3459b2-c3fe-47ef-94af-82d045d858b0
md"""
But the number of unique subject sequences per query sequence was always one:
"""

# ╔═╡ 0789debf-927f-4bb3-8297-a2bcb2ddd43f
unique(df_qseqid_count.sseqid_function)

# ╔═╡ 857f9004-8f56-4e15-91c4-6adcb67d8b6b
md"""
**2. We define the alignment coverage as the percentage of the subject sequence covered by the alignment. Compute the alignment coverage for all sequences from the blastn output.**

We compute the alignment coverage and add it as a column to our main blast DataDrame:
"""

# ╔═╡ ce6bf287-f90f-49ce-a04e-a64d0f693b5f
df.alignment_coverage = df.length ./ df.slen;

# ╔═╡ 5c4e50c4-0204-491b-851f-226e3bb28dd1
md"""
**3. Visualize the distribution of alignment coverages as a histogram.**
"""

# ╔═╡ 0e21763c-fa15-4f19-88d4-2cc68c5bade1
histogram(100*df.alignment_coverage, label="", xlabel="Alignment coverage (%)")

# ╔═╡ 69b35b29-bd8e-4908-9848-b08123992f6f
df_hc = subset(df, :alignment_coverage => x -> x.>= 0.6)

# ╔═╡ 00745bf2-13fe-4ed8-8744-ca2b9b801877
md"""
**4. Remove alignments with coverage less than 60% and verify that each query sequence is now aligned to at most one subject sequence.**

Most alignments have high coverage, but **$(sum(df.alignment_coverage .< 0.6))** have coverage below 60%. We create a new DataFrame keeping only the high-coverage (>= 60%). The number of rows in the filtered DataFrame is **$(nrow(df_hc))**, compared to the total number of input sequences, which was **$(length(filtered_reads))**.
"""

# ╔═╡ 8c27c535-e079-4fda-8bba-2eeaaab4c7ea
plot(df_hc.score[1:100:end], df_hc.bitscore[1:100:end], xlabel="Raw score", ylabel="Bit score", label="")

# ╔═╡ bc8dcde6-ba0a-4a69-a77f-19777724c57c
begin
	i1 = argmin(df.score);
   	i2 = argmax(df.score);
   	x1, y1 = [df.score[i1], df.bitscore[i1]];
   	x2, y2 = [df.score[i2], df.bitscore[i2]];
   	a = (y2 - y1) / (x2 - x1);
   	b = y2 - a * x2;
	λ = a*log(2);
	K = 2^(-b);
end

# ╔═╡ c6a0eade-d189-4de6-b06f-39a23e094708
md"""
**5. The bit score $S'$ is derived from the raw score $S$ using the formula**

```math
S' = \frac{\lambda S - \ln K}{\ln 2}
```

In other words, the bit score is a linear function of the raw score, which we can verify immediately by plotting them against each other. If ``(x_1,y_1)`` and ``(x_2,y_2)`` are two point on a line ``y=a x + b``, the slope ``a`` is found from

```math
a = \frac{y_2 - y_1}{x_2 - x_1}.
```

Once we know the slope, the intercept ``b`` is found from ``b = y_2 - a x_2``. Applied to the alignment scores, we only need two pairs of scores, and might as well take the ones furthest apart. 

Applied to the score formula, we see that the slope parameter is given ``a=\frac{\lambda}{\ln 2}``, which results in the value **``\lambda=``$(λ)**, and the intercept parameter by ``b=-\frac{\ln K}{\ln 2}``, which results in the value **``K=``$(K)**.
"""

# ╔═╡ 44027cb8-0a54-42cb-baa3-0c56cbf661e9
y2-a*x2

# ╔═╡ ce28181d-db1c-4521-a6e3-17ec9386dcd6
 vsg = unique(df_hc.sseqid)

# ╔═╡ 51a53766-c67f-46bf-9724-53b1eb420258
md"""
## Count VSG expression levels

**1. Extract the unique VSG ids in *your* sample**

We obtain the unique VSGs directly from the blast results DataFrame. There are **$(length(vsg))** of them.
"""

# ╔═╡ f9e62367-ca9d-49a4-b576-ddcd3ed3f917
md"""
**2. For each unique VSG: count the number of sequences aliging to that VSG, and the average number of identical matches, mismatches, and gaps for its alignments.**

Grouping the DataFrame by the subject sequence id (that is, by VSG) and using the split-apply-combine strategy, we can count the number of query sequence ids per VSG in one line. 
"""

# ╔═╡ 1f61124e-6fb8-4b60-84ea-f5e95622e186
begin
	df_vsg = combine( groupby(df_hc, :sseqid) , nrow,  :nident => mean, :mismatch => mean, :gaps => mean);
	rename!(df_vsg, :sseqid => "VSG", :nrow => "count");
end

# ╔═╡ d586e716-545a-4e21-8b52-aaea737fcdd4
md"""
**3. Identify the 10 most abundant VSGs in *your* sample and visualize their relative expression levels in a pie chart.**

We sort the VSG DataFrame in reverse order by their abundance count, and then we can take the first 10 elements for plotting.
"""

# ╔═╡ a06a7f65-aedd-40d9-8e03-02143fbf290e
sort!(df_vsg, :count, rev=true)

# ╔═╡ f866826f-c4e2-4111-89b9-22071a53de92
pie(df_vsg.VSG[1:10], df_vsg.count[1:10])

# ╔═╡ b78bb3c7-46c2-4ecd-8c50-f7fd0fa740bb
md"""
### Identify open reading frames

**1. Find the relevant `getorf` command on the [longread-application repository][2] and adapt it to work on *your* input sample. Include the exact `getorf` command in your report.**

The correct `getorf command is

```
/usr/local/emboss/bin/getorf -sequence $qf -outseq $orf -minsize 1200 -find 3 -reverse N
```

where

- `qf` is the input file name PacBio\_VSG\_filtered\_reads\_ $(sample).fasta
- `orf` is the output file name PacBio\_VSG\_filtered\_reads\_ORF\_$(sample).fasta
"""

# ╔═╡ 096919e9-e4a8-40fc-9434-8d1937e31c04
begin
	predicted_orfs = Vector{FASTARecord}();
	   FASTAReader( open(datadir("PacBio_VSG", "orf", "PacBio_VSG_filtered_reads_ORF_" * sample[1] * ".fasta")) ) do reader
	      for record in reader
	         push!(predicted_orfs, record)
	      end
	   end
end

# ╔═╡ 932e630c-47c9-49c0-9139-04eab4effbdf
length(predicted_orfs)

# ╔═╡ aa2901b6-188f-4594-8b7e-b6dea446e0e4
length(filtered_reads)

# ╔═╡ ac37c9f6-b878-4c11-8bf3-5bc03b473443
100 * length(predicted_orfs) / length(filtered_reads)

# ╔═╡ Cell order:
# ╟─e355a610-4af8-11ee-29f9-d916ed9c2903
# ╠═9ddcb453-89eb-4642-b313-40d9f2548d7c
# ╠═4d7d51cc-2751-4a98-bc27-b2e0f21990ec
# ╠═f1c539ca-6ca9-40f2-9a46-cdccef31c838
# ╟─731d038c-78e0-4740-8f80-174e1c7d70b3
# ╟─0b807b4f-7634-468c-8619-df0348981882
# ╟─b9ecaf9c-2c40-494d-83ea-bc9a4bb10056
# ╠═f3637e6c-2984-47e5-a0c0-394f3bdd2cae
# ╟─d187d9f9-e756-484f-acfe-4accf7c59552
# ╟─47ef717f-fd1f-47c9-a080-894f8eda0149
# ╠═6314b615-22a9-4aa4-8b59-614281426e80
# ╟─b6b50963-f5ae-4efd-aa42-b156baaefaed
# ╠═d927f106-35c7-4648-8c5b-adc0fd6c6c60
# ╟─8db859f2-2794-431c-9a8a-221832296d11
# ╠═fccf0588-ca40-46ec-ad64-619eaa4c184a
# ╟─3c3459b2-c3fe-47ef-94af-82d045d858b0
# ╠═0789debf-927f-4bb3-8297-a2bcb2ddd43f
# ╟─857f9004-8f56-4e15-91c4-6adcb67d8b6b
# ╠═ce6bf287-f90f-49ce-a04e-a64d0f693b5f
# ╟─5c4e50c4-0204-491b-851f-226e3bb28dd1
# ╟─0e21763c-fa15-4f19-88d4-2cc68c5bade1
# ╟─00745bf2-13fe-4ed8-8744-ca2b9b801877
# ╠═69b35b29-bd8e-4908-9848-b08123992f6f
# ╟─c6a0eade-d189-4de6-b06f-39a23e094708
# ╠═8c27c535-e079-4fda-8bba-2eeaaab4c7ea
# ╠═bc8dcde6-ba0a-4a69-a77f-19777724c57c
# ╠═44027cb8-0a54-42cb-baa3-0c56cbf661e9
# ╟─51a53766-c67f-46bf-9724-53b1eb420258
# ╠═ce28181d-db1c-4521-a6e3-17ec9386dcd6
# ╟─f9e62367-ca9d-49a4-b576-ddcd3ed3f917
# ╠═1f61124e-6fb8-4b60-84ea-f5e95622e186
# ╟─d586e716-545a-4e21-8b52-aaea737fcdd4
# ╠═a06a7f65-aedd-40d9-8e03-02143fbf290e
# ╠═f866826f-c4e2-4111-89b9-22071a53de92
# ╟─b78bb3c7-46c2-4ecd-8c50-f7fd0fa740bb
# ╠═096919e9-e4a8-40fc-9434-8d1937e31c04
# ╠═932e630c-47c9-49c0-9139-04eab4effbdf
# ╠═aa2901b6-188f-4594-8b7e-b6dea446e0e4
# ╠═ac37c9f6-b878-4c11-8bf3-5bc03b473443
