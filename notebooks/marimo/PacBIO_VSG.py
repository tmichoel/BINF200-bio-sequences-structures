import marimo

__generated_with = "0.14.16"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Application of long read sequencing to determine expressed antigen diversity in *Trypanosoma brucei* infections""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Select your sample""")
    return


@app.cell
def _(os):
    os.getcwd()
    return


@app.cell(hide_code=True)
def _(mo):
    import os
    datadir = os.path.join("data", "PacBio_VSG")
    dropdown = mo.ui.dropdown(["balbc_3_0", "balbc_3_1", "balbc_3_2", "balbc_3_3", "balbc_3_4", "balbc_6_0", "balbc_6_1", "balbc_6_2", "balbc_6_4", "balbc_6_5", "balbc_10_0", "balbc_10_1", "balbc_10_2", "balbc_10_4", "balbc_10_5","balbc_12_1", "balbc_12_2", "balbc_12_3", "balbc_12_4", "balbc_12_5"])
    dropdown
    return datadir, dropdown, os


@app.cell(hide_code=True)
def _(filtered_reads, mo):
    mo.md(
        rf"""
    ## Tasks

    ### Show that you understand the biological experiment and data

    To count the number of sequences in the input file, we read the FASTA file with input sequences into a vector of FASTA records. The length of the vector is the number of sequences: {len(filtered_reads)}
    """
    )
    return


@app.cell
def _(datadir, dropdown, os):
    from Bio import SeqIO
    fastaname = os.path.join(datadir, "filtered_reads", "PacBio_VSG_filtered_reads_" + dropdown.value + ".fasta") 
    filtered_reads = list(SeqIO.parse(fastaname, 'fasta'))
    return SeqIO, filtered_reads


@app.cell(hide_code=True)
def _(mo):
    mo.md(r""" """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
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
    )
    return


@app.cell(hide_code=True)
def _(df, mo):
    mo.md(
        rf"""
    ### Analyze the BLAST output

    Answer the following questions / perform the following tasks:

    **1. What is the number of sequences in *your* blast output file? Is it the same as in the input file?**

    Import the blast results in a DataFrame. The column names can be constructed from the blastn command. The number of rows of the DataFrame is the number of sequences: {len(df)}
    """
    )
    return


@app.cell
def _(datadir, dropdown, os):
    import pandas as pd
    blastname = os.path.join(datadir, "blastn", "PacBio_VSG_filtered_reads_blastn_" + dropdown.value + ".txt")
    df = pd.read_csv(blastname, sep='\t', header=None)
    df.columns = ["qseqid", "sseqid", "score", "bitscore", "evalue", "qlen", "slen", "length", "sstart", "send", "nident", "mismatch", "gaps", "positive"]
    df
    return (df,)


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(df_qseqid_count, filtered_reads, mo):
    mo.md(rf"""Despite setting the option `-max_target_seq 1` in the `blastn` command, this number may be *greater* than the number of input sequences. The reason is that sometimes for *one* query sequence, *two* separate alignments are found against the *same* subject sequence. To show this, we compute for each unique query id, its number of alignments (number of rows it appears in in `df`) and the number of unique subject ids in those rows. Grouping the DataFrame by the query sequence id and using the split-apply-combine strategy, we can count the number of subject sequence ids per query in one line. After this computation, we see that now the number of unique query sequences (**{len(df_qseqid_count)}**) is indeed equal to the number of input sequences (**{len(filtered_reads)}**).""")
    return


@app.cell
def _(df):
    df_qseqid_count = df.groupby("qseqid").agg(num_sseqid=('sseqid', 'nunique'),
        num_align=('sseqid', 'size'))
    df_qseqid_count
    return (df_qseqid_count,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We can verify that sometimes two alignments were found per query sequence:""")
    return


@app.cell
def _(df_qseqid_count):
    df_qseqid_count.num_align.unique()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""But the number of unique subject sequences per query sequence was always one:""")
    return


@app.cell(hide_code=True)
def _(df_qseqid_count):
    df_qseqid_count.num_sseqid.unique()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    **2. We define the alignment coverage as the percentage of the subject sequence covered by the alignment. Compute the alignment coverage for all sequences from the blastn output.**

    We compute the alignment coverage and add it as a column to our main blast DataDrame:
    """
    )
    return


@app.cell
def _(df):
    df['alignment_coverage'] = df['length']/df['slen'] * 100
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""**3. Visualize the distribution of alignment coverages as a histogram.**""")
    return


@app.cell
def _(df):
    import seaborn as sns
    df['alignment_coverage'].hist(bins=50)
    return (sns,)


@app.cell
def _(df):
    sum(df['alignment_coverage'] < 0.6)
    return


@app.cell(hide_code=True)
def _(df, df_hc, filtered_reads, mo):
    mo.md(
        rf"""
    **4. Remove alignments with coverage less than 60% and verify that each query sequence is now aligned to at most one subject sequence.**

    Most alignments have high coverage, but **{sum(df['alignment_coverage'] < 60)}** have coverage below 60%. We create a new DataFrame keeping only the high-coverage (>= 60%). The number of rows in the filtered DataFrame is **{len(df_hc)}**, compared to the total number of input sequences, which was **{len(filtered_reads)}**.
    """
    )
    return


@app.cell
def _(df):
    df_hc = df[df['alignment_coverage'] >= 60]
    df_hc
    return (df_hc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    **5. The bit score $S'$ is derived from the raw score $S$ using the formula**

    $$
    S' = \frac{\lambda S - \ln K}{\ln 2}
    $$

    In other words, the bit score is a linear function of the raw score, which we can verify immediately by plotting them against each other. If $(x_1,y_1)$ and $(x_2,y_2)$ are two point on a line $y=a x + b$, the slope $a$ is found from

    $$
    a = \frac{y_2 - y_1}{x_2 - x_1}.
    $$

    Once we know the slope, the intercept $b$ is found from $b = y_2 - a x_2$. Applied to the alignment scores, we only need two pairs of scores, and might as well take the ones furthest apart. 

    Applied to the score formula, we see that the slope parameter is given $a=\frac{\lambda}{\ln 2}$, which results in the value **$\lambda=$**, and the intercept parameter by $b=-\frac{\ln K}{\ln 2}$.
    """
    )
    return


@app.cell(hide_code=True)
def _(K, lam, mo):
    mo.md(rf"""This results in the values $\lambda=$ {lam:.2f} and $K=$ {K:.2f}.""")
    return


@app.cell
def _(df_hc, sns):
    sns.scatterplot(data=df_hc, x='score', y='bitscore')
    return


@app.cell
def _(df_hc):
    import numpy as np
    i1, i2 = np.argmin(df_hc['score']), np.argmax(df_hc['score'])
    x1, y1 = df_hc['score'].iloc[i1], df_hc['bitscore'].iloc[i1]
    x2, y2 = df_hc['score'].iloc[i2], df_hc['bitscore'].iloc[i2]
    a = (y2 - y1) / (x2 - x1)
    b = y2 - a * x2
    lam = a * np.log(2)
    K = np.exp(-b * np.log(2))
    return K, lam


@app.cell(hide_code=True)
def _(mo, vsg):
    mo.md(
        rf"""
    ## Count VSG expression levels

    **1. Extract the unique VSG ids in *your* sample**

    We obtain the unique VSGs directly from the blast results DataFrame. There are **{len(vsg)}** of them.
    """
    )
    return


@app.cell
def _(df_hc):
    vsg = df_hc.sseqid.unique()
    return (vsg,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    **2. For each unique VSG: count the number of sequences aliging to that VSG, and the average number of identical matches, mismatches, and gaps for its alignments.**

    Grouping the DataFrame by the subject sequence id (that is, by VSG) and using the split-apply-combine strategy, we can count the number of query sequence ids per VSG in one line.
    """
    )
    return


@app.cell
def _(df_hc):
    df_vsg = df_hc.groupby("sseqid").agg(
        count=('qseqid', 'size'),
        avg_nident=('nident', 'mean'),
        avg_mismatch=('mismatch', 'mean'),
        avg_gaps=('gaps', 'mean'))
    df_vsg.sort_values('count', ascending=False, inplace=True)
    df_vsg
    return (df_vsg,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    **3. Identify the 10 most abundant VSGs in *your* sample and visualize their relative expression levels in a pie chart.**

    We take the first 10 rows of the sorted VSG dataframe for plotting.
    """
    )
    return


@app.cell
def _(df_vsg, dropdown):
    import matplotlib.pyplot as plt
    df_vsg_top10 = df_vsg.head(10)
    plt.figure(figsize=(10, 6))
    plt.pie(df_vsg_top10['count'], labels=df_vsg_top10.index, autopct='%1.1f%%', startangle=140)
    plt.title('Top 10 VSGs in Sample ' + dropdown.value)
    plt.axis('equal')  # Equal aspect ratio
    plt.show()
    return


@app.cell(hide_code=True)
def _(dropdown, mo):
    mo.md(
        rf"""
    ### Identify open reading frames

    **1. Find the relevant `getorf` command on the [longread-application repository](https://github.com/siddharthjayaraman/longread-application) and adapt it to work on *your* input sample. Include the exact `getorf` command in your report.**

    The correct `getorf command is

    ```
    /usr/local/emboss/bin/getorf -sequence $qf -outseq $orf -minsize 1200 -find 3 -reverse N
    ```

    where

    - `qf` is the input file name PacBio\_VSG\_filtered\_reads\_{dropdown.value}.fasta
    - `orf` is the output file name PacBio\_VSG\_filtered\_reads\_ORF\_{dropdown.value}.fasta
    """
    )
    return


@app.cell
def _(SeqIO, datadir, dropdown, os):
    orfname = os.path.join(datadir, "orf", "PacBio_VSG_filtered_reads_ORF_" + dropdown.value + ".fasta") 
    predicted_orfs = list(SeqIO.parse(orfname, 'fasta'))
    return (predicted_orfs,)


@app.cell(hide_code=True)
def _(filtered_reads, mo, predicted_orfs):
    mo.md(rf"""We find {len(predicted_orfs)} ORFs in {len(filtered_reads)}, a percentage of {len(predicted_orfs)/len(filtered_reads)*100:.2f}% of the input sequences.""")
    return


if __name__ == "__main__":
    app.run()
