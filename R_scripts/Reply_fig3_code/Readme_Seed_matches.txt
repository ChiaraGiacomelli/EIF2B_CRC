This R scripts were written to search for shRNA seed-matches in Ribo-seq results to account for potential off-target effects of the shRNA sequences. This was in response to a reviewer concerned by the fact that we had used only one shRNA per target subunit.

Each script searches for a different kind of seed match (6mer, 7mer-A1, 7mer-m8, 8mer).
This is defined by what is a seed so that the match is calculated accordingly (ca line 55 of code per script)

Overall structure:
1) Reads in the differential expression tables for 2B1 and 2B4 independently.
NB: These DEseq2 output had been calculated using already the correction for non-targeting control using the design = ~ batch + shRNA + condition + shRNA:condition and extracting the results shRNA2BN.conditionplus (where N was 1 or 4, independently)

2) RPF for 1 and 4, and totals for 1 and 4 are then merged with the TE tables
Repeated columns (ie those present in more than one of the tables) are discarded.
Groups of regulation are also then compared although cumulative groups are not further used in the script - it's just that I was anyway merging data and I thought was good to have an overview.

3) From the FASTA file (gencode.v38.pc_transcripts_filtered) (separately for CDS and 3' UTR) the seed sequence matches are searched and counted for shCTRl, sh2B1, and sh2B4 - the number of seed matches is retained and merged to the merged_TE file.

4) The seed matches are then categorized on whether they are present for all shRNAs, or only subsets. This categorization is briefly explored numerically by printing the summary on screen.

5) From this table, plots are being made with the cumulative fraction of transcripts containing seed matches for CTRL, or 2B1, or 2B4 shRNAs.
The legend enumerates how many transcript have how many matches, as well as the average log2FC of what is indicated in the title (ie Ribosome occupancy displays in the legend the average the log2FC for RPFs).
For the analysis of the matches of shRNA2B1 and 2B4, the respective actual target is filtered out from the merged_TE_seed_matches table prior to plotting given that we are searching for off target effects. I did not test whether keeping the gene was making a difference, I simply removed it because logically for me it made sense.