# egaDifferential

Compare the interaction graphs of two conditions: gained, lost and shared interactions (optionally tested by Chicdiff), contact-weighted chromatin changes at promoter, direct and indirect (two-hop) level and up to --maxHops graph steps, and their relation to expression changes.

```text
egaDifferential [-h] [--verbose] [--version] --interactionsA
                             INTERACTIONSA [INTERACTIONSA ...] --interactionsB
                             INTERACTIONSB [INTERACTIONSB ...]
                             [--format {ibed,washu,bedpe,mifsud,homer,cool}]
                             [--baits BAITS] [--coolField COOLFIELD]
                             [--maxFdr MAXFDR]
                             [--scoreThreshold SCORETHRESHOLD]
                             [--replicates {all,intersect,union}]
                             [--referenceScore REFERENCESCORE]
                             [--supportThreshold SUPPORTTHRESHOLD]
                             [--supportFloor SUPPORTFLOOR]
                             [--chicdiff CHICDIFF]
                             [--chicdiffPadj CHICDIFFPADJ]
                             [--names NAMES NAMES] --tss TSS
                             [--geneTypes [GENETYPES ...]] --tracksA TRACKSA
                             [TRACKSA ...] --tracksB TRACKSB [TRACKSB ...]
                             [--markSign [MARKSIGN ...]] --expression
                             EXPRESSION [--log2FCColumn LOG2FCCOLUMN]
                             [--padjColumn PADJCOLUMN] [--alpha ALPHA]
                             [--expressionColumns EXPRESSIONCOLUMNS EXPRESSIONCOLUMNS]
                             [--minExpression MINEXPRESSION]
                             [--expressionPseudoCount EXPRESSIONPSEUDOCOUNT]
                             [--geneKey {gene_ids,genes}]
                             [--rnaThreshold RNATHRESHOLD]
                             [--chromatinThreshold CHROMATINTHRESHOLD]
                             [--permutations PERMUTATIONS]
                             [--communityPermutations COMMUNITYPERMUTATIONS]
                             [--genePermutations GENEPERMUTATIONS]
                             [--maxHops MAXHOPS]
                             [--stepFeature {combined,marks,contacts,decomposed}]
                             [--stepWeighting {product,unweighted,normalised}]
                             [--proximitySupport PROXIMITYSUPPORT]
                             [--stateThreshold STATETHRESHOLD]
                             [--stateThresholds [STATETHRESHOLDS ...]]
                             [--stateScale {fold,genomeMean}]
                             [--stateMinFragments STATEMINFRAGMENTS]
                             [--markRoles [MARKROLES ...]]
                             [--roles {declared,learned,both}]
                             [--roleBootstrap ROLEBOOTSTRAP]
                             [--stateTracksA STATETRACKSA [STATETRACKSA ...]]
                             [--stateTracksB STATETRACKSB [STATETRACKSB ...]]
                             [--statePeaksA STATEPEAKSA [STATEPEAKSA ...]]
                             [--statePeaksB STATEPEAKSB [STATEPEAKSB ...]]
                             [--statePeakMinOverlap STATEPEAKMINOVERLAP]
                             [--stateWindow STATEWINDOW]
                             [--stateMinChange STATEMINCHANGE]
                             [--annotation ANNOTATION]
                             [--elements ELEMENTS [ELEMENTS ...]]
                             [--primaryRepressive [PRIMARYREPRESSIVE ...]]
                             [--elementMinGenes ELEMENTMINGENES]
                             [--structuralVariants STRUCTURALVARIANTS [STRUCTURALVARIANTS ...]]
                             [--svContacts SVCONTACTS]
                             [--copyNumber COPYNUMBER [COPYNUMBER ...]]
                             [--svWindow SVWINDOW] [--svMinPairs SVMINPAIRS]
                             [--svMinObsExp SVMINOBSEXP]
                             [--svThreshold SVTHRESHOLD]
                             [--copyNumberThreshold COPYNUMBERTHRESHOLD]
                             [--bootstrap BOOTSTRAP] [--resolution RESOLUTION]
                             [--seed SEED] [--threads THREADS] [--tads TADS]
                             [--tadFormat {auto,hicexplorer,arrowhead,cooltools,bed,bed-boundaries,bed-domains,ega}]
                             [--tadWindow TADWINDOW] [--tadChain TADCHAIN]
                             [--tadRegions TADREGIONS]
                             [--tadMinStrength TADMINSTRENGTH] --outputPrefix
                             OUTPUTPREFIX
```

| Option | Default | Description |
|---|---|---|
| `--verbose` |  | Print progress information |
| `--version` |  | show program's version number and exit |
| `--interactionsA` `INTERACTIONSA ...` (required) |  | Called interactions of condition A, one file per replicate (e.g. CHiCAGO ibed) |
| `--interactionsB` `INTERACTIONSB ...` (required) |  | Called interactions of condition B, one file per replicate |
| `--format`, `-f` `{ibed, washu, bedpe, mifsud, homer, cool}` | `ibed` |  |
| `--baits` `BAITS` |  | Baited fragments of the capture design (BED or CHiCAGO .baitmap); required for --format cool, for other formats it sets the bait status and orientation of each pair |
| `--coolField` `COOLFIELD` | `count` | Pixel column read from cool matrices |
| `--maxFdr` `MAXFDR` |  | Keep pairs with an FDR at most this value (inputs with an FDR column, e.g. HOMER) |
| `--scoreThreshold` `SCORETHRESHOLD` |  | Keep pairs with score >= threshold in both conditions (e.g. CHiCAGO 5) |
| `--replicates` `{all, intersect, union}` | `all` | 'all': pairs present in every replicate and reaching --scoreThreshold in at least one; 'intersect': reaching it in every replicate; 'union': in any replicate |
| `--referenceScore` `REFERENCESCORE` | `5.0` | Score that maps to edge strength 1 (the calling threshold, CHiCAGO 5) |
| `--supportThreshold` `SUPPORTTHRESHOLD` |  | Calling score of the support classes of gained and lost edges (default: --scoreThreshold, else --referenceScore) |
| `--supportFloor` `SUPPORTFLOOR` |  | Gained / lost edges whose pair is listed in the other condition at a score of at least this value (below the calling score) are threshold-sensitive, below it stable; absence from a file whose lowest score is above the floor is not testable (default: 0.6 x the calling score, 3 at CHiCAGO 5). Columns support, support_reason, other_max_score, replicates_floor_<name> of <prefix>_edges.tsv.gz, edge attribute support, <prefix>_support_inputs.tsv (lowest score per input file), <prefix>_support.tsv (edges per status and class) and the support columns of the all-level explanations |
| `--chicdiff` `CHICDIFF` |  | Chicdiff results table (fwrite of <prefix>_results.Rds) testing the change of each interaction |
| `--chicdiffPadj` `CHICDIFFPADJ` | `0.05` |  |
| `--names` `NAMES` | `['A', 'B']` | Condition names |
| `--tss` `TSS` (required) |  | GTF or TSS BED for bait gene annotation |
| `--geneTypes` `GENETYPES ...` | `['protein_coding', 'lincRNA', 'antisense']` |  |
| `--tracksA` `TRACKSA ...` (required) |  | mark=file.bw[,rep.bw] for condition A |
| `--tracksB` `TRACKSB ...` (required) |  | mark=file.bw[,rep.bw] for condition B |
| `--markSign` `MARKSIGN ...` | `['H3K27me3=-1', 'H3K9me3=-1']` | mark=+1\|-1\|0, the expected direction of a mark on expression (default +1); an alias of --markRoles with +1 = active, -1 = repressive, 0 = neutral |
| `--expression` `EXPRESSION` (required) |  | Gene table: either a log2 fold change column (--log2FCColumn) or two expression columns (--expressionColumns) |
| `--log2FCColumn` `LOG2FCCOLUMN` |  | Column holding log2(B/A), e.g. DESeq2 log2FoldChange |
| `--padjColumn` `PADJCOLUMN` |  | Adjusted p-value column; restricts differential genes to padj <= --alpha |
| `--alpha` `ALPHA` | `0.05` |  |
| `--expressionColumns` `EXPRESSIONCOLUMNS` |  | Expression columns of A and B (e.g. RPKM) |
| `--minExpression` `MINEXPRESSION` | `1.0` | With --expressionColumns: keep genes expressed at least this much in A or B |
| `--expressionPseudoCount` `EXPRESSIONPSEUDOCOUNT` | `1.0` |  |
| `--geneKey` `{gene_ids, genes}` | `gene_ids` |  |
| `--rnaThreshold` `RNATHRESHOLD` | `1.0` | \|log2FC\| of a differential gene |
| `--chromatinThreshold` `CHROMATINTHRESHOLD` | `1.0` | \|log2FC\| of a chromatin feature counted as explaining a gene |
| `--permutations` `PERMUTATIONS` | `100` |  |
| `--communityPermutations` `COMMUNITYPERMUTATIONS` | `10000` | Random gene sets per community size for the coherence test of community expression |
| `--genePermutations` `GENEPERMUTATIONS` | `1000` | Permutations of the per-gene test of the indirect level (0 skips it) |
| `--maxHops` `MAXHOPS` | `3` | Graph levels explaining a gene: promoter, 1 hop (direct), 2 hops (indirect), ... up to this many hops (at least 2; level k >= 3 writes columns <mark>_hop<k>, n_hop<k>, ...) |
| `--stepFeature` `{combined, marks, contacts, decomposed}` | `combined` | What the step levels (1 .. --maxHops) measure. 'combined': log2 (S_B + 1) / (S_A + 1) of the contact-weighted mark signal S with the contact strengths r and marks x of each condition (changes of r and of x are confounded); 'marks': the mark term, S with the marks of B versus A under the fixed strengths r_bar = (r_A + r_B) / 2 on the union graph (an edge called in one condition enters with half its strength); 'contacts': the contact term, S with the strengths of B versus A under the fixed marks x_bar = (x_A + x_B) / 2; 'decomposed': the mark term for the cascade, the learned roles, the mirror control and the per-gene permutation test, and in addition the columns <mark>_<level>_combined / _marks / _contacts (and _combined_mirror, _contacts_mirror) and in the all-level table combined_changes, mark_term_changes, contact_term_changes, concordant_combined_R, concordant_mark_term_R, concordant_contact_term_R and driver_R (marks, contacts, both, neither: which term changes concordantly at the level). The promoter levels and contact_<level> are the same under every choice |
| `--stepWeighting` `{product, unweighted, normalised}` | `product` | Weights of the nodes at k steps in S (for every --stepFeature). 'product': the sum over the shortest paths to a node of the product of the strengths r along the path; 'unweighted': every node at k steps enters once (no path multiplicity, no strengths); 'normalised': the product weights divided by their sum over the nodes at k steps, so S is a weighted mean of the mark signal over those nodes (no dependence on the number of nodes and paths). The promoter levels and contact_<level> are the same under every choice |
| `--proximitySupport` `PROXIMITYSUPPORT` |  | Step curve of analysis/hop_distance.py (TSV with steps, concordance, ci_low, ci_high): adds to every explained gene the probability that graph neighbours at its explaining level are closer in the nucleus than distance-matched pairs (proximity_support, proximity_support_ci; 1.0 / NA at the promoter) |
| `--stateThreshold` `STATETHRESHOLD` | `2.0` | Chromatin states: a mark is present on a fragment if its mean signal over the fragment is at least this value (fold enrichment; see --stateScale). Bivalent promoter (bait): an active-role and a repressive-role mark present; poised enhancer (other fragment): an active- or enhancer-role and a repressive-role mark present. Marks with peak calls (--statePeaksA / --statePeaksB) are present where a peak overlaps the fragment instead |
| `--stateThresholds` `STATETHRESHOLDS ...` | `[1.5, 2.0, 3.0]` | Presence thresholds of the state transition report (<prefix>_bivalent_transitions.tsv, <prefix>_state_thresholds.tsv); the cascade uses --stateThreshold |
| `--stateScale` `{fold, genomeMean}` | `fold` | 'fold': the bigWig signal is a fold enrichment over input (e.g. Roadmap .fc.signal); 'genomeMean': coverage tracks without input, divided by their genome-wide mean (sumData / nBasesCovered of the bigWig header) before the threshold |
| `--stateMinFragments` `STATEMINFRAGMENTS` | `50` | Learned states: mark combinations seen on at least this many fragments (both conditions) are kept, the rest pooled |
| `--markRoles` `MARKROLES ...` | `[]` | Declared mark roles, mark=role entries and / or TSV files with the columns mark and role; roles: active, repressive, enhancer, architectural, neutral, unknown (enhancer: active-type at fragments that are not baited promoters; architectural, neutral, unknown: no expected direction). They take precedence over --markSign; marks declared nowhere are active |
| `--roles` `{declared, learned, both}` | `declared` | Roles used by the explanation cascade and the states: declared (--markRoles / --markSign), learned from the data separately at every level (<prefix>_mark_roles.tsv: promoter, step 1 .. --maxHops and the own interactions; at each level the cascade uses the role learned there, a mark without a consistent effect at a level explains no gene there), or both (declared for the main outputs, learned in addition as <prefix>_explanation_cascade_learned.tsv and columns *_learned) |
| `--roleBootstrap` `ROLEBOOTSTRAP` | `1000` | Chromosome block-bootstrap resamples of the learned roles and learned states |
| `--stateTracksA` `STATETRACKSA ...` |  | mark=file.bw for the chromatin states of condition A (e.g. fold enrichment over input), when the --tracksA signal is not on a fold-enrichment scale; marks not listed use --tracksA |
| `--stateTracksB` `STATETRACKSB ...` |  | as --stateTracksA for condition B |
| `--statePeaksA` `STATEPEAKSA ...` |  | mark=peaks.bed[,rep.bed] for the chromatin states of condition A: peak calls of any caller (MACS2 / MACS3 narrowPeak or broadPeak, SEACR, generic BED; columns chrom, start, end are read, several files of a mark are pooled). A mark with peaks is present on a fragment when at least one peak overlaps it by at least --statePeakMinOverlap bp, under both presence statistics and every --stateThresholds value (bivalent promoters, poised enhancers, element classes, learned states); marks without peaks keep the threshold rule. The signal of the mark (--stateTracksA, else --tracksA) still gives the change required by --stateMinChange and the values in the state summaries; the mark changes of the explanation levels always come from the signal tracks. A mark given here needs peaks in --statePeaksB too |
| `--statePeaksB` `STATEPEAKSB ...` |  | as --statePeaksA for condition B |
| `--statePeakMinOverlap` `STATEPEAKMINOVERLAP` | `1` | Overlap in bp of a peak with a fragment that makes the mark present (--statePeaksA / --statePeaksB) |
| `--stateWindow` `STATEWINDOW` | `1000` | Chromatin states are called under two presence statistics of a mark on a fragment, both reported: the mean signal over the fragment (outputs without suffix, columns *_mean) and the maximum of the mean signal over windows of at most this many bp within the fragment (outputs *_max1kb, columns *_max1kb), which keeps narrow peaks on long fragments |
| `--stateMinChange` `STATEMINCHANGE` | `0.5` | A state transition needs every mark crossing the presence threshold to change by at least this log2 fold (pseudo count 1) in the direction of the crossing |
| `--annotation` `ANNOTATION` |  | Element classes: GTF or TSS BED whose transcription start sites make a fragment a promoter, baited or not (gene types --geneTypes; default: --tss) |
| `--elements` `ELEMENTS ...` |  | External element annotation, BED with a name column (e.g. ENCODE cCREs, ChromHMM states): one file for both conditions or one per condition (A B); the agreement of the element classes with it is written to <prefix>_element_agreement.tsv |
| `--primaryRepressive` `PRIMARYREPRESSIVE ...` |  | Element classes: primary repressive marks (e.g. the Polycomb-type mark); a mark that is repressive without being primary makes a heterochromatin-like element (default: the marks declared repressive) |
| `--elementMinGenes` `ELEMENTMINGENES` | `50` | Roles learned per chromatin class of the nodes: classes with at least this many genes at a level |
| `--structuralVariants` `STRUCTURALVARIANTS ...` |  | Breakpoints of structural variants called by an external SV caller (e.g. from whole-genome sequencing): BEDPE with the conditions in which every join is present (column conditions, present_<condition>, pairs_<condition> / obs_exp_<condition>), condition=file.bedpe (every join of the file present in that condition) or condition=table of inter-chromosomal 1 Mb bin pairs with their read support (columns chrom1 start1 chrom2 start2 pairs expected obs_exp). Adds the level 'structural variant' after the promoter levels, flags the graph edges that join the two sides of a breakpoint (sv=True) and writes <prefix>_sv_genes.tsv, _sv_breakpoints.tsv |
| `--svContacts` `SVCONTACTS` |  | Contacts of baited fragments across the junctions (TSV: join, bait_chrom, bait_start, bait_end, partner_chrom, partner_start, partner_end, pairs_<condition>, inter_total_<condition>): the promoter's own contacts with the partner side, written as edges sv=True to <prefix>_sv_edges.tsv and <prefix>_union_sv.graphml. With --format ibed the called inter-chromosomal and long-range pairs across a breakpoint (dropped from the cis graph) are written to <prefix>_sv_called_edges.tsv and added to <prefix>_union_sv.graphml as edges sv=True |
| `--copyNumber` `COPYNUMBER ...` |  | condition=file.bedGraph of the relative copy number per condition (from an external copy number caller, e.g. read depth of whole-genome sequencing): adds the level 'copy number' after the structural variant level |
| `--svWindow` `SVWINDOW` | `1.0` | Mb: promoters within this distance of a breakend, partner nodes within it of the partner breakend |
| `--svMinPairs` `SVMINPAIRS` | `50` | Bin-pair tables of --structuralVariants: read pairs of a breakpoint bin pair |
| `--svMinObsExp` `SVMINOBSEXP` | `20.0` | Bin-pair tables of --structuralVariants: observed / expected of a breakpoint bin pair |
| `--svThreshold` `SVTHRESHOLD` | `1.0` | log2 change of the support (obs/exp) of a join present in both conditions counted as gained or lost |
| `--copyNumberThreshold` `COPYNUMBERTHRESHOLD` | `0.5` | log2 B/A of the relative copy number at the promoter counted as a change |
| `--bootstrap` `BOOTSTRAP` | `1000` | Bootstrap resamples of the mirror comparison |
| `--resolution` `RESOLUTION` | `1.0` | Community resolution parameter |
| `--seed` `SEED` |  |  |
| `--threads` `THREADS` |  | Worker processes for the bigWig reads of the mirrored loci, the per-gene permutation test and the community coherence test (default: the available cores, at most 32); results do not depend on it |
| `--tads` `TADS` |  | TAD / boundary calls: hicFindTADs prefix or file (_domains.bed, _boundaries.bed, _tad_score.bm, _score.bedgraph), Juicer Arrowhead / 4DN 2D domain list, cooltools insulation table, BED of boundaries or domains (optional score column), or egaBoundaries output |
| `--tadFormat` `{auto, hicexplorer, arrowhead, cooltools, bed, bed-boundaries, bed-domains, ega}` | `auto` | Format of --tads |
| `--tadWindow` `TADWINDOW` |  | cooltools: insulation window (default: the largest in the table) |
| `--tadChain` `TADCHAIN` |  | Chain file from the assembly of --tads to that of the graph (UCSC chromosome names) |
| `--tadRegions` `TADREGIONS` |  | BED of named regions (chrom start end name) used to name the partner regions |
| `--tadMinStrength` `TADMINSTRENGTH` |  | Only boundaries at or above this strength percentile (0-1) count as crossed |
| `--outputPrefix`, `-o` `OUTPUTPREFIX` (required) |  | Output prefix; with --tads also <prefix>_tad_crossing_* (egaTadCrossing) |
