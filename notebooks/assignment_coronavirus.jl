### A Pluto.jl notebook ###
# v0.19.27

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

- **protein_N_data.fasts** - Sequences of the gene coding for coronavirus nucleocapsid (N) protein in a number of coronaviruses
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

$(PlutoUI.LocalResource("../figures/protein_N_data_clustalo_phylotree.png"))
"""

# ╔═╡ cbf62395-5f9a-41bf-bade-c32a021fbe60
md"""
#### Interpretation

Based on the results from the previous two steps, what do you see? Elaborate with a small text (3-4 lines): Explain what you observe from the multiple sequence alignment itself (hint: check the number of conserved sites), and give a short interpretation of the phylogenetic tree you have constructed.

### Step-by-step multiple sequence alignment and phylogenetic tree construction using UPGMA

#### Compute pairwise similarities 

Use the Needleman-Wunsch (dynamic programming) pairwise alignment algorithm to build a matrix of global alignment scores for each pair of sequences in **protein_N_data.fasta**. You can choose between multiple options:

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
			S[i,j] = score( pairalign(GlobalAlignment(), seq[i], seq[j], scoremodel) )
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
- ``S_{rand}`` is the expected (average) score for aligning two random sequences of the same length and residue composition, obtained by random shuffling the nucleotide composition of the two sequences. (Hint: more info about the Feng & Doolittle can be found at this URL: [https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle](https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle))

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
		s += score(pairalign(GlobalAlignment(), rs1, rs2, scoremodel))
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
D = triu( -log.( (S .- Srand) ./ (Smax .- Srand) ) , 1)

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

# ╔═╡ 404f6d5f-08ac-48c1-9178-73cd0bd3d051
hc.order

# ╔═╡ b22a39b7-a35b-4923-a18c-8fd630bb7fea
description.(protein_N_data)

# ╔═╡ bdb9ea2f-a22d-485c-b838-aec7703836fa
begin
	plot(hc, xrotation=90, bottom_margin = 50mm)
	xticks!(1:20, description.(protein_N_data)[hc.order])
end

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
# ╟─18bb5a3f-2cb2-4cb3-8e4e-28396f3f6ed3
# ╠═81b0c557-9965-4ece-b2f4-9e3800bcc387
# ╠═75c19231-4f64-42c3-831d-5c9d010f652f
# ╠═bac7d343-e83f-4a89-87c1-0cf52e6547e4
# ╠═4087cf53-c4c7-49e6-9148-1778f4d33f2c
# ╟─d2a53401-e003-4012-bdbd-5ff7d4712531
# ╠═60e20fe5-a39d-4b33-b386-4d8b574293a0
# ╟─76448d4e-2894-49a6-8f75-b75361e36927
# ╠═6e18f4f5-c51c-464d-b633-aa5a37df3a8a
# ╟─0d3979ab-8f0e-476f-a623-8d9a15ce9af4
# ╠═1a1df46a-66b0-4733-bebc-46f9896e7766
# ╠═404f6d5f-08ac-48c1-9178-73cd0bd3d051
# ╠═b22a39b7-a35b-4923-a18c-8fd630bb7fea
# ╠═bdb9ea2f-a22d-485c-b838-aec7703836fa
