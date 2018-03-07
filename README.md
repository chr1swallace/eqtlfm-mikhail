# eqtlfm-mikhail
code to use GUESSFM for eqtl fine mapping


The initial preparation of the data is documented in get-data.sh.

Then, each gene was analysed separately.  To see how the code is run, try

```
$ ./run.rb -h
runfalsehelptrueintfalse
Usage: runguess.rb [OPTIONS] COMMAND

  COMMANDS are:
      (run things)
           prep :: prepare GUESS input
           guess :: run guess
           expand :: expand guess output

           clean :: remove slurm files
```

For each gene that was represented multiple times, in `all_tested_affinity_expression_associations.txt`, I checked whether all the snp positions were within 1e+5 of each other or not.  If they were, I merged the region, and if not I conducted multiple analyses per gene. See assoc-processed.tab.

Note for  this region I couldn't find any SNPs in the interval
```{sh}
   snp_chr         gene_id   minpos   maxpos            myid
1:       6 ENSG00000196735 32590893 32592973 ENSG00000196735
```

For each analysis, I took a 2e+5 window either side of the leftmost and rightmost SNP position to construct a region for fine mapping.  Then I tested each SNP for association and if no SNP had a p < 1e-5, I didn't bother to finemap (no information).  For the remainder, I prepared GUESS input, ran GUESS, and expanded the results in R using GUESSFM to generate a posterior over the most likely model space.  The R files prep-files.R and expand-results.R have code to do the first and third steps.

The file `results.tgz` contains these results saved as files "*-snpmod-99.RData".  These can be explored using the GUESSFM library (github.com/chr1swallace/GUESSFM).  But I also generated the most useful summary stats from them.  The files *-nsnp.csv give the posterior for the number of causal snps.  The files *-mppi.csv give the marginal posterior probabilities that each SNP is included in a causal model.  Any value above 0.001 indicates it is plausibly causal, but the larger the stronger.  However, if n SNPs are in tight LD and 1 of them is always in a causal model, then this MPPI is split across them, so that each might be expected to have MPPI =~ 1/n.  Thus lower values can reflect extended LD too.  There are functions in GUESSFM to explore the individual models in more detail for the genes you care most about, but this might be overkill for your aims.

Each results directory also includes a png showing visual qc of the GUESS run.  What we want to see here is that the different coloured lines all overlap, and that there is no trend (in the first 3 graphs) - that is, the chains have converged.  It is worth checking these, particularly if you have an unusually large number of causal variants in the model.
