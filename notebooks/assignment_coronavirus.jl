### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ e04efa7d-af26-486b-99e7-4829b7e5c852
using DrWatson

# ╔═╡ 4da703b9-d9ed-4ccb-8073-09b11fb43b0f
@quickactivate "BINF200"

# ╔═╡ 01724442-fb44-477f-a2a4-e200e400bf2d
begin
	using BioSequences
	using BioAlignments
	using BioSymbols
	using FASTX
	using DataFrames
	using CSV
	using Plots
	using Plots.PlotMeasures
	using Clustering
	using Statistics
	using StatsPlots
	using LaTeXStrings
	using LinearAlgebra
	using Random
	using PlutoUI
end

# ╔═╡ f503ce8e-5d23-11ee-0b33-f73001d28c23
md"""
# Multiple sequence alignment, phylogenetics, and motif analysis of coronavirus genomes

## Setup the environment
"""

# ╔═╡ f58842a5-60f0-4184-9025-053dcd5815d2
md"""
## Data

Download all files in the following OneDrive folder:

[https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/Coronavirus?csf=1&web=1&e=D5umem](https://universityofbergen-my.sharepoint.com/:f:/r/personal/tom_michoel_uib_no/Documents/public/BINF200/Coronavirus?csf=1&web=1&e=D5umem)

The following files are available:

- **protein\_N\_data.fasta** - Sequences of the gene coding for coronavirus nucleocapsid (N) protein in a number of coronaviruses
- **GCA_011537005.1_partial_genomic.fasta** - Part of the BetaCoV/Wuhan/IPBCAMS-WH-02/2019 genome
- **motifCountMatrix.csv** - Count matrix of a sequence motif

They are saved in:
"""

# ╔═╡ cae3ce79-4128-4b44-801e-4c1b677e863c
datadir("Coronavirus")

# ╔═╡ a4516d89-8330-418b-80e1-6318bf8c2c86
md"""
### Multiple sequence alignment and phylogenetic tree construction for the coronavirus nucleocapsid protein



We will focus on different genera of corona viruses namely: alpha, beta, gamma and delta. Their genomes, gene and protein sequences, together with annotations and data reports are available from NCBI:

- Alpha coronavirus: [https://www.ncbi.nlm.nih.gov/datasets/taxonomy/693996/](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/693996/)
- Beta coronavirus: [https://www.ncbi.nlm.nih.gov/datasets/taxonomy/694002/](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/694002/)
- Gamma coronavirus: [https://www.ncbi.nlm.nih.gov/datasets/taxonomy/694013/](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/694013/)
- Delta coronavirus: [https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1159901/](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1159901/)

For simplicity, we will use data from only one important gene, that encodes the coronavirus nucleocapsid (N) protein. This is a structural protein that forms complexes with genomic RNA, interacts with the viral membrane protein during virion assembly and plays a critical role in enhancing the efficiency of virus transcription and assembly. You can read more about it in the paper [*The SARS-CoV-2 Nucleocapsid Protein and Its Role in Viral Structure, Biological Functions, and a Potential Target for Drug or Vaccine Mitigation*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8227405/).

A subset of viruses was used to create **protein\_N\_data.fasta**.

#### Parse the fasta file

How many sequences are contained in the file **protein\_N\_data.fasta**? List the names of the sequences.
"""

# ╔═╡ fbab5158-a837-4c70-ad47-cc526bd7461a
begin
	protein_N_data = Vector{FASTARecord}();
   	FASTAReader( open(datadir("Coronavirus","protein_N_data.fasta")) ) do reader
	      for record in reader
	         push!(protein_N_data, record);
	      end
	   	end
end

# ╔═╡ db865d39-0a20-4103-822e-83d0a5f3d52c
md"""
*There are $(length(protein_N_data)) sequences in the file. Their names are:*
"""

# ╔═╡ 5a1ce827-59d5-4ce8-9f27-1bc4d28c6eb8
for seq = description.(protein_N_data)
	println(seq)
end

# ╔═╡ 73877ecd-a6ff-40f9-a3f1-48dd6ab81bb6
md"""
#### Find protein N in a specific coronavirus genome

From the sequences in **protein\_N\_data.fasta**, find the sequence for which the first letter in its name is closest in the alphabet to the **first letter in your first name**. If there are multiple sequences starting with the same letter, pick on arbitrarily. For the selected sequence:

- Find the assembly accession ID in the table above.
- Go to the NCBI website (cf. links above) and find the corresponding genome assembly.
- What are the genomic coordinates (start and end position) of gene N in this genome? (Hint: follow the RefSeq link)
"""

# ╔═╡ 95d3c40f-c21a-4715-ac75-9ca5c9e9fa44
md"""
#### Multiple sequence alignment

Build a multiple sequence alignment for the **protein\_N\_data** using a multiple sequence alignment tool of your choice. (Hint: check out the services provided by EMBL's European Bioinformatics Institute (EMBL-EBI).)

*A good choice is [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/).*

#### Phylogenetic tree reconstruction

Based on the results from the previous step, build a phylogenetic tree. (Hint: at this stage it is not required to make an “advanced tree”, providing a simple tree is enough). Save the image of the phylogenetic/guide tree.

*After uploading **protein\_N\_data.fasta** and running Clustal Omega, the output page has a tab "Phylogenetic Tree". Here is a screenshot of the tree:*

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_clustalo_phylotree.png"))
"""

# ╔═╡ cbf62395-5f9a-41bf-bade-c32a021fbe60
md"""
#### Interpretation

Based on the results from the previous two steps, what do you see? Elaborate with a small text (3-4 lines): Explain what you observe from the multiple sequence alignment itself (hint: check the number of conserved sites), and give a short interpretation of the phylogenetic tree you have constructed.

### Step-by-step multiple sequence alignment and phylogenetic tree construction using UPGMA

#### Compute pairwise similarities 

Use the Needleman-Wunsch (dynamic programming) pairwise alignment algorithm to build a matrix of global alignment scores for each pair of sequences in **protein\_N\_data.fasta**. You can choose between multiple options:

- Implement the Needleman-Wunsch algorithm yourself. (Hint: You have probably done this already in BINF100) 
- Use an existing implementation of the algorithm. (Hint: Check biopython, biojulia)
- Use the *needleall* command line program from the EMBOSS suite. (Hint: You installed the whole EMBOSS suite for Assignment 1.)
- Use a webserver such as EMBL-EBI's EMBOSS Needle service. (Hint: Manually inputting every pair of sequences will be extremely tedious, though they do provide APIs.)

*We will use the [BioAlignments](https://github.com/BioJulia/BioAlignments.jl) package.*

*First we set the score model. Here we use the popular [EDNAFULL](https://rosalind.info/glossary/dnafull/) substitution matrix:*
"""

# ╔═╡ fe340d20-3ce2-43f9-81aa-185e7b7c7edb
scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);

# ╔═╡ 4eed04b6-39d3-4d9c-9d12-3bb6bf566fec
md"""
*We will store the pairwise similarities in the upper diagonal of a square matrix. Note that in the next task we will need self-alignment scores, and these are therefore included on the diagonal.*
"""

# ╔═╡ af109968-cab5-4603-8a44-a044abf1fae2
function pairalign_vec(seq, scoremodel)
	N = length(seq);
	S = zeros(N,N);
	for i = 1:N
		for j = i:N
			S[i,j] = BioAlignments.score( pairalign(GlobalAlignment(), seq[i], seq[j], scoremodel) )
		end
	end
	return S
end

# ╔═╡ fdf16615-3a4e-4961-bf42-3f78312e6efe
S = pairalign_vec(FASTX.sequence.(protein_N_data), scoremodel)

# ╔═╡ 18bb5a3f-2cb2-4cb3-8e4e-28396f3f6ed3
md"""
#### Generate a pairwise distance matrix

Generate a distance matrix from the score matrix you have created in the previous step. For this task we will use Feng & Doolittle's formulation, and we will compute the distance $D$ using formula:

```math
D = -\log S_{eff} = -\log \frac{S_{obs}-S_{rand}}{S_{max} - S_{rand}}
```

where 

- ``S_{obs}`` is the observed pairwise alignment score
- ``S_{max}`` is the best alignment score for both sequences, obtained by taking the average of the score of aligning either sequence to itself
- ``S_{rand}`` is the expected (average) score for aligning two random sequences of the same length and residue composition, obtained by random shuffling the nucleotide composition of the two sequences.

Compute $S_{rand}$ by taking the average score of **10** pairwise alignments between random sequences with the same sequence compositions as the original sequences.
"""

# ╔═╡ 81b0c557-9965-4ece-b2f4-9e3800bcc387
function seq_freq(seq)
	# we assume our input sequence is a string of DNA letters
	alph = ['A','C','G','T']
	if sort(unique(seq))!=alph
		error("Non-ACGT detected.")
	end
	# create count vector
	c = vec(map(x -> count(x,seq), alph))/length(seq)
end

# ╔═╡ 75c19231-4f64-42c3-831d-5c9d010f652f
function random_score(seq1, seq2, scoremodel)
	nr = 10
	sp1 = SamplerWeighted(dna"ACGT", seq_freq(seq1)[1:3]);
	sp2 = SamplerWeighted(dna"ACGT", seq_freq(seq2)[1:3]);
	s = 0.0
	for k=1:nr
		rs1 = randseq(DNAAlphabet{2}(), sp1, length(seq1))
		rs2 = randseq(DNAAlphabet{2}(), sp2, length(seq2))
		s += BioAlignments.score(pairalign(GlobalAlignment(), rs1, rs2, scoremodel))
	end
	s /= nr
end

# ╔═╡ bac7d343-e83f-4a89-87c1-0cf52e6547e4
function random_pairalign_vec(seq, scoremodel)
	N = length(seq);
	S = zeros(N,N);
	for i = 1:N
		for j = i:N
			S[i,j] = random_score(seq[i], seq[j], scoremodel)
		end
	end
	return S
end

# ╔═╡ 4087cf53-c4c7-49e6-9148-1778f4d33f2c
 Srand = random_pairalign_vec(FASTX.sequence.(protein_N_data), scoremodel)

# ╔═╡ d2a53401-e003-4012-bdbd-5ff7d4712531
md"""
*Now compute ``S_{max}``, the average of the score of aligning either sequence to itself, as a matrix:*
"""

# ╔═╡ 60e20fe5-a39d-4b33-b386-4d8b574293a0
Smax = triu(0.5*(diag(S)' .+ diag(S)))

# ╔═╡ 76448d4e-2894-49a6-8f75-b75361e36927
md"""
*Finally compute the distance matrix:*
"""

# ╔═╡ 6e18f4f5-c51c-464d-b633-aa5a37df3a8a
D = triu( -log.( max.(S .- Srand, 1.) ./ (Smax .- Srand) ) , 1)

# ╔═╡ fd0b5a14-2838-455b-b6af-5af4471e6009
D[2,8]

# ╔═╡ 0d3979ab-8f0e-476f-a623-8d9a15ce9af4
md"""
#### Generate a "guide tree" of phylogenetic relationships

Generate a "guide tree" of phylogenetic relationships from the pairwise distance matrix you have created in the previous step using the UPGMA method. You can choose between multiple options:

- Implement the UPGMA hierarchical clustering algorithm yourself. (Hint: You can represent the tree as a binary tree, either implementing a tree class yourself, or using an existing data structure.)
- Use an existing implementation of the algorithm. (Hint: UPGMA is more commonly known as hierarchical clustering with average linkage. Check SciPy or similar packages for other languages.)

*Here we use the second option:*
"""

# ╔═╡ 1a1df46a-66b0-4733-bebc-46f9896e7766
hc = hclust(D, linkage=:average, branchorder=:optimal, uplo=:U)

# ╔═╡ 053bd29d-475b-44a5-984e-962a2b124f73
md"""
*Set labels:*
"""

# ╔═╡ b22a39b7-a35b-4923-a18c-8fd630bb7fea
xtl = map(x -> String(x), description.(protein_N_data))

# ╔═╡ af76bd93-ad2e-439e-957e-90729c858617
md"""
#### Interpret your results

Visualize your guide tree and compare it to the phylogenetic tree constructed in @sec-phylo-1. Elaborate with a small text (3-4 lines) to explain what you observe.
"""

# ╔═╡ bdb9ea2f-a22d-485c-b838-aec7703836fa
begin
	pl = plot(hc,orientation=:horizontal, xlims=[0, 4])
	yticks!(1:20,xtl[hc.order])
end

# ╔═╡ c4b9078e-bd1a-4ebe-905a-8343fc522a63
md"""
### Sequence motifs

Do simple motif searching on corona virus sequences using the input dataset (**protein\_N\_data.fasta**) we have already analysed.

#### MEME analysis

Connect to the MEME platform at [https://meme-suite.org/](https://meme-suite.org/).

- Find the MEME motif discovery tool.
- Input **protein\_N\_data.fasta** to discover enriched motifs in this set of sequences, allowing for zero or one motif occurrence per sequence and finding upto 5 motifs. Which discovery mode, sequence alphabet, and site distribution options do you select?

Open and download the **MEME HTML output file** and include the sequence logos of the motifs found in your report.

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_meme_motifs_1.png"))

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_meme_motifs_2.png"))

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_meme_motifs_3.png"))

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_meme_motifs_4.png"))

$(PlutoUI.LocalResource("../figures/Coronavirus/protein_N_data_meme_motifs_5.png"))
"""

# ╔═╡ 6ecb20fc-2deb-488c-875a-c4e65616229f
md"""
#### Convert count matrix to PWM

We will work with a 20-nucleotide subset of the first motif found by the MEME software, given by the count matrixin the file **motifCountMatrix.csv**.
"""

# ╔═╡ e3cf83f5-bb34-41f8-bd7c-99e1c077897b
begin
	df = DataFrame(CSV.File(datadir("Coronavirus","motifCountMatrix.csv")))
	rename!(df, :"base\\position" => "base")
	df.base = only.(df.base)
	df
end

# ╔═╡ 50a8d2fd-f1fe-406d-af2f-43dc713efdf4
md"""
1. Compare the count matrix against your sequence logos and mark the 20-nucleotide window corresponding to this count matrix in the right logo.
"""

# ╔═╡ 78ccc4f5-55ab-469d-97a4-f83ed1e69e57
md"""
2. Convert the count matrix to a position-specific probability matrix (PPM) ``P``. To avoid zeros in the PPM, we add *pseudo-counts* and define

   ```math
   P_{k,i} = \frac{\text{Count}_{k,i} + 0.25 * \sqrt{N}}{N + \sqrt{N}},
   ```

   where $\text{Count}_{k,i}$ is the value of the count matrix for nucleotide ``k`` in motif position ``i``, and ``N`` is the number of sequences in **protein_N_data.fasta** (Hint: Count the totals in each column of the count matrix).


"""

# ╔═╡ ebfd05f9-f130-4974-af88-9e2de6220bec
Count = Matrix(df[:,2:end])

# ╔═╡ ab6b63a4-4669-4636-a260-cdb282e76446
N = sum(Count, dims=1)[1]

# ╔═╡ 21196e69-b864-4880-9722-c9362d611d84
N*(N+1)/2

# ╔═╡ 2f193934-df62-4bfe-b03f-1ee583cb6a1c
PPM = (Count .+ 0.25*√N) ./ (N + √N)

# ╔═╡ 64dd2abf-858b-4119-a65e-b14c332ad42f
md"""
3. Convert the PPM matrix to a position-specific weight matrix (PWM) ``W`` using the formula

   ```math
   W_{k,i} = \log_2 \frac{P_{k,i}}{0.25}
   ```

   What would be the value of ``W`` for a random background site with equal counts for all nucleotides and using the pseudo-count formula above to compute the random probabilities?
"""

# ╔═╡ 67a88754-ce20-4b49-9f5f-eb5c87ecafa8
PWM = log2.( PPM ./ 0.25)

# ╔═╡ c828dab8-cbf4-4d42-aa75-39a7503f9d05
md"""
*If all nucleotides are equally frequent, their counts would be ``0.25N``, their probabilities would be ``(0.25N + 0.25\sqrt N)/(N+\sqrt N)=0.25``, and their weights would be ``\log_2 1 = 0``.*
"""

# ╔═╡ a1107731-4b65-4916-9fa1-a152878ec313
md"""
#### Scan a coronavirus genome for motif occurrences

Scan part of the BetaCoV/Wuhan/IPBCAMS-WH-02/2019 genome (the sequence in the file **GCA\_011537005.1\_partial\_genomic.fasta**) and score all possible motif occurrences. Use the sliding window approach presented in the lecture and report (table and figure) both the log-odds score and the odds of each possible motif starting position in the genome sequence.
"""

# ╔═╡ 816a19eb-34ac-493f-80d2-f57055607613
genome = first( FASTAReader( open(datadir("Coronavirus","GCA_011537005.1_partial_genomic.fasta")) ) )

# ╔═╡ 9e1659ed-1716-446c-b68f-d43bd68c2377
genome_seq = FASTX.sequence(genome)

# ╔═╡ 5e0dc69d-24c9-4715-bfc6-136cfb0dd93f
md"""
*Convert the genome sequence into an array of integers where each value corresponds to the row index of the corresponding nucleotide in the PWM:*
"""

# ╔═╡ bf2585bb-7969-4807-ab45-d9dad0a39447
begin
	genome_vec = zeros(Int64,length(genome_seq))
	for k=1:length(df.base)
		genome_vec[findall(df.base[k],genome_seq)] .= k
	end
	genome_vec
end

# ╔═╡ 810edb80-cda4-45fa-8c35-a57c2d5f3c07
md"""
*Function to score a given segment of the sequence*
"""

# ╔═╡ e68228a7-e8e7-49d5-aa77-6322b3de498f
function score_segment(genome_vec, PWM, start)
	ml = size(PWM,2)
	return tr(PWM[genome_vec[start:start+ml-1],1:ml])
end

# ╔═╡ 090aca02-1f39-4a0d-8c2d-c83a31ec92ac
md"""
*Get score for all possible start positions:*
"""

# ╔═╡ 3d07c839-68f3-4534-b270-9108af9201d6
begin
	starts = 1:length(genome_vec)-size(PWM,2)+1
	score = zeros(size(starts))
	for start in starts
		score[start] = score_segment(genome_vec, PWM, start)
	end
end

# ╔═╡ 7f3b38bd-7fb8-417b-819d-701a25b38097
md"""
*Convert score to odss:*
"""

# ╔═╡ 46e549f8-4c99-4294-a4e5-068aafaa1853
odds = 2 .^score

# ╔═╡ 7cfa7582-c62f-4a07-b6bf-9435c5406478
md"""
*Plot the odds ratios:*
"""

# ╔═╡ 17ac9ab5-96c9-4033-9677-5062c4beb700
plot(starts,odds,line=:stem)

# ╔═╡ f348def8-70be-4846-a682-84930dcb2633
md"""
Elaborate with a small text (3-4 lines) to explain what you observe.
"""

# ╔═╡ Cell order:
# ╟─f503ce8e-5d23-11ee-0b33-f73001d28c23
# ╠═e04efa7d-af26-486b-99e7-4829b7e5c852
# ╠═4da703b9-d9ed-4ccb-8073-09b11fb43b0f
# ╠═01724442-fb44-477f-a2a4-e200e400bf2d
# ╟─f58842a5-60f0-4184-9025-053dcd5815d2
# ╠═cae3ce79-4128-4b44-801e-4c1b677e863c
# ╟─a4516d89-8330-418b-80e1-6318bf8c2c86
# ╠═fbab5158-a837-4c70-ad47-cc526bd7461a
# ╟─db865d39-0a20-4103-822e-83d0a5f3d52c
# ╠═5a1ce827-59d5-4ce8-9f27-1bc4d28c6eb8
# ╟─73877ecd-a6ff-40f9-a3f1-48dd6ab81bb6
# ╟─95d3c40f-c21a-4715-ac75-9ca5c9e9fa44
# ╟─cbf62395-5f9a-41bf-bade-c32a021fbe60
# ╠═fe340d20-3ce2-43f9-81aa-185e7b7c7edb
# ╟─4eed04b6-39d3-4d9c-9d12-3bb6bf566fec
# ╠═af109968-cab5-4603-8a44-a044abf1fae2
# ╠═fdf16615-3a4e-4961-bf42-3f78312e6efe
# ╠═21196e69-b864-4880-9722-c9362d611d84
# ╟─18bb5a3f-2cb2-4cb3-8e4e-28396f3f6ed3
# ╠═81b0c557-9965-4ece-b2f4-9e3800bcc387
# ╠═75c19231-4f64-42c3-831d-5c9d010f652f
# ╠═bac7d343-e83f-4a89-87c1-0cf52e6547e4
# ╠═4087cf53-c4c7-49e6-9148-1778f4d33f2c
# ╟─d2a53401-e003-4012-bdbd-5ff7d4712531
# ╠═60e20fe5-a39d-4b33-b386-4d8b574293a0
# ╟─76448d4e-2894-49a6-8f75-b75361e36927
# ╠═6e18f4f5-c51c-464d-b633-aa5a37df3a8a
# ╠═fd0b5a14-2838-455b-b6af-5af4471e6009
# ╟─0d3979ab-8f0e-476f-a623-8d9a15ce9af4
# ╠═1a1df46a-66b0-4733-bebc-46f9896e7766
# ╟─053bd29d-475b-44a5-984e-962a2b124f73
# ╠═b22a39b7-a35b-4923-a18c-8fd630bb7fea
# ╟─af76bd93-ad2e-439e-957e-90729c858617
# ╠═bdb9ea2f-a22d-485c-b838-aec7703836fa
# ╟─c4b9078e-bd1a-4ebe-905a-8343fc522a63
# ╟─6ecb20fc-2deb-488c-875a-c4e65616229f
# ╠═e3cf83f5-bb34-41f8-bd7c-99e1c077897b
# ╟─50a8d2fd-f1fe-406d-af2f-43dc713efdf4
# ╟─78ccc4f5-55ab-469d-97a4-f83ed1e69e57
# ╠═ebfd05f9-f130-4974-af88-9e2de6220bec
# ╠═ab6b63a4-4669-4636-a260-cdb282e76446
# ╠═2f193934-df62-4bfe-b03f-1ee583cb6a1c
# ╟─64dd2abf-858b-4119-a65e-b14c332ad42f
# ╠═67a88754-ce20-4b49-9f5f-eb5c87ecafa8
# ╟─c828dab8-cbf4-4d42-aa75-39a7503f9d05
# ╟─a1107731-4b65-4916-9fa1-a152878ec313
# ╠═816a19eb-34ac-493f-80d2-f57055607613
# ╠═9e1659ed-1716-446c-b68f-d43bd68c2377
# ╟─5e0dc69d-24c9-4715-bfc6-136cfb0dd93f
# ╠═bf2585bb-7969-4807-ab45-d9dad0a39447
# ╟─810edb80-cda4-45fa-8c35-a57c2d5f3c07
# ╠═e68228a7-e8e7-49d5-aa77-6322b3de498f
# ╟─090aca02-1f39-4a0d-8c2d-c83a31ec92ac
# ╠═3d07c839-68f3-4534-b270-9108af9201d6
# ╟─7f3b38bd-7fb8-417b-819d-701a25b38097
# ╠═46e549f8-4c99-4294-a4e5-068aafaa1853
# ╟─7cfa7582-c62f-4a07-b6bf-9435c5406478
# ╠═17ac9ab5-96c9-4033-9677-5062c4beb700
# ╟─f348def8-70be-4846-a682-84930dcb2633
