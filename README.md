# Outbreak in Washington State caused mainly by two seperate introduction

## Results

### Outbreak in Washington State caused by repeated introductions and an early rapid growth

SARS-CoV-2 was repeatedly introduced into Washington State. The outbreak can be separated into two groups that lead to the majority of cases and originated from at least two different introductions into Washington State. The first one was likely introduced at the beginning of February from China. The second one is derived from Europe and was most likely introduced between mid and the end of February. To date, these two local outbreak cluster make up the vast majority of cases in Washington State. Additionally, we see evidence for several additional introductions of lineages into Washington State that are derived from lineages that previously circulated in Europe, as well as from some, were the origin of the lineage is more uncertain.  These lineages were most likely introduced from areas were sampling and sequencing is sparse, which could include other areas of the United States.

When estimating the effective reproduction number through time, we see that the outbreak grew rapidly between mid and end of February. During this time, it was unknown that SARS-CoV-2 was spreading in the USA. The Re then dropped to about 2.5 by the end of February, which is consistent with previously reported values for reproduction numbers in other places.

On the XXth of March, local spread of SARS-CoV-2 was first reported in the Washington State area and the number of confirmed new cases began to flatten soon after. On the XXth of March, several businesses started to institute home office and on the XXth of March, a stay at home order was issued from the Governor.

### Testing of cases accurately reflects trends in new cases.

We next estimated the population dynamics of SARS-CoV-2 in WA from genetic sequence data using a coalescent approach, which conditions on sampling. This means that the information about the population dynamics come from the branching patterns and not how many sequences there are through time.


Additionally, we used a birth death model to estimate the transmission dynamics





## Methods and Materials

### Introductions into Washington State

In order to distinguish between sequences that are connected by local transmission, we cluster all sequences from Washington State together based

We then classify each introduction as either being derived from China or from Europe. To do so, we use the nextstrain pipeline to

### Estimating population dynamics jointly from multiple local outbreak clusters

In order to estimate the population dynamics of the Washington State outbreak, we use a coalescent approach to infer these dynamics jointly from all known local outbreak clusters. To do so, we model the coalescence and migration of lineages within Washington State as structured coalescent process with known migration history. The know migration history here is given by the clustering of sequences into local outbreak clusters. The migration events from anywhere outside WA into WA are always assumed to have happened before the common ancestor of all sequences in each local outbreak cluster. How long before this common ancestor time is inferred during the MCMC.

We then infer the effective population size and rates of introductions through time using a skyline type approach. Effective population sizes and rates of introduction are allowed to change at predefined time points. Between these predefined time points where the rates are estimated, the rates are interpolated. This is equivalent to assuming exponential growth or decline between the effective population sizes at theses time points.

In contrast to backwards in time coalescent approaches, we can consider different local outbreak clusters as independent observations of the same underlying population process.
