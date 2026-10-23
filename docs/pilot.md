# OrthoFinder pilot: findings

> **Note, 2026-10-23.** The repository is now on coordinator 0.14 (pin `cichlid-coordinator>=0.14,<0.15`;
> before that 0.13; the instruction sheet is `docs/orthofinder_repo_spec.md`, Dispatch 5, written against
> 0.13.0). The pilot notes below were written against 0.10.0 and are kept as they were, as the record of the
> pilot. See `docs/SPEC_VERSION.md`.

*Draft, 2026-10-20. Written against cichlid-coordinator v0.10.0 and the agreed settings in
`orthofinder_repo_spec.md` section 4. Sections marked "pending" are filled in from the pilot output.*

## Machine and programs

| | |
|---|---|
| Machine | 32 cores, 31 GB RAM, 32 GB swap; 109 GB free on `/` (the data root `~/Temp/TreeFamData` is on it) |
| Environment | conda `phylo-trees` (coordinator 0.10.0 installed with `pip install -e`) |
| OrthoFinder | 3.1.5 (bioconda; `orthofinder -v` prints `OrthoFinder:v3.1.5`) |
| IQ-TREE | 3.1.3 (`iqtree3`) |
| FAMSA (OrthoFinder's default aligner) | 2.4.1-45c9b2b |
| DIAMOND | 2.2.8 |
| MCL | 22-282 |

The launcher `bin/orthofinder` runs `/usr/bin/env python3`; outside the environment that is the system
Python, without Biopython, and it fails. Run it inside the environment (`conda activate phylo-trees` or
`conda run -n phylo-trees`).

## Pilots

1. **5 species, full proteomes:** Mzeb, Acal, Gacu, Drer, Locu (125,879 proteins), with the agreed settings
   and `-s` given the catalog's tree pruned to the five (`(Locu,(Drer,(Gacu,(Mzeb,Acal))));`). For timing.
2. **5 species, subsampled:** the same species, the genes of 317 orthogroups of 4 to 20 genes from pilot 1
   plus every 500th gene (2,183 proteins). Same settings. For the output structure.

The pilot proteomes were built by the coordinator (`download genomes`, `proteomes`) into a scratch data
root (`CICHLID_DATA_ROOT`), before the cloud copies were reachable; nothing from them was pushed.

## Option names (installed version)

| Sheet | Installed 3.1.5 |
|---|---|
| `-M msa` | `-M msa` (and it is now the default) |
| `-T iqtree` | **`-T iqtree3`**; the choices are `fasttree`, `fasttree_fastest`, `raxml`, `iqtree3` |
| aligner | `-A`, default `famsa` (others: `muscle`, `mafft`) |
| `-s <tree>` | `-s <file>`, "user-specified rooted species tree" |
| `-o <dir>` | `-o <dir>` exists. **OrthoFinder writes a `Results_<Mon><day>` folder inside it** (`-n <name>` appends a name) |
| `-X` | `-X`: "Don't add species names to sequence IDs" |
| threads | `-t` (search, default 32 here), `-a` (analysis) |

## How OrthoFinder runs IQ-TREE (cost)

OrthoFinder's built-in command is `iqtree3 -s INPUT --prefix PATH/IDENTIFIER -quiet` (its `config.json`).
No model is given, so **IQ-TREE runs ModelFinder on every orthogroup**: the logs say "ModelFinder will
test up to 1232 protein models". Orthogroups are processed from the largest down, one IQ-TREE per core.
In pilot 1 the largest orthogroups (100 to 598 sequences, 600 to 1,300 columns) take more than 40 minutes
each. This is the agreed setting (OrthoFinder's defaults); the estimate below is for it as agreed.

**Where IQ-TREE's time goes** (from the logs of the 25 trees finished so far: 5 from pilot 1, 20 from
pilot 2; wall-clock seconds, machine fully loaded):

| Orthogroups | Trees | ModelFinder | Tree search and the rest | ModelFinder's share |
|---|---|---|---|---|
| 100 to 130 sequences | 5 | 579 s | 2,873 s | 17% |
| 10 to 20 sequences | 20 | 518 s | 601 s | 46% |

**Correction (21:30, the largest clusters of pilot 1):** ModelFinder's share is much larger for them, because
it picks FreeRate models (`+R6`, `+R7`) that are expensive to fit:

| Cluster | Sequences × columns | ModelFinder | Model chosen | Status after 7.9 h |
|---|---|---|---|---|
| OG0000000 | 598 × 1,261 | 1 h 57 min | JTT+F+I+G4 | tree search, iteration 213 |
| OG0000002 | 280 × 514 | 4 h 35 min | JTT+I+R7 | tree search, iteration 390 |
| OG0000006 | 182 × 1,053 | 5 h 38 min | VT+F+R6 | tree search, iteration 190 |
| OG0000001 | 576 × 525 | **still running** | | still in ModelFinder |

So for the largest clusters ModelFinder is about 20 to 65% of the time, and a single cluster can hold the
whole run up: the pilot's end is set by OG0000001, not by the thousands of small clusters.

For the first trees finished, ModelFinder was not most of the time: for large orthogroups the tree search dominated (up to 1,292 s for
119 sequences and 574 columns). (Superseded by the whole pilot, below: over all
17,598 runs ModelFinder was about half of IQ-TREE's time.)
ModelFinder chose 14 different models for the 25 orthogroups (JTT, JTTDCMUT, VT, WAG, Dayhoff, LG and
the Q.* matrices, with +G4, +I, +R2, +F); LG+G4 was chosen once.

**If a fixed model is ever agreed** (not the agreed setting; recorded only so the option is known):
OrthoFinder 3.1.5 reads an optional `~/config_orthofinder_user.json` (not present on this machine) and
adds its tools to the built-in ones, so an entry such as `"iqtree3_lg"` with
`iqtree3 -s INPUT --prefix PATH/IDENTIFIER -m LG+G4 -quiet` would be used with `-T iqtree3_lg`, with no
change to OrthoFinder itself. It would need a new run (`-b` can reuse the finished search).

## Pilot 1: results

Finished 2026-10-21 02:04, exit status 0. Kept in `pilot_runs/main_5species/` (outside the repository).

| Step | Result |
|---|---|
| Wall time | **12 h 39 min** on 32 cores (13:25 to 02:04) |
| DIAMOND search and preliminary clusters | about 6 minutes; 22,939 clusters, 17,598 with 3 or more sequences (each gets an IQ-TREE run) |
| Gene trees | 17,598 IQ-TREE runs, **275 hours of IQ-TREE time** (268 CPU-hours); all 32 cores busy for about 9 h, then the largest clusters alone |
| The last trees | OG0000001 (576 sequences): 12 h 33 min, of which 9 h 49 min choosing its model (JTT+R8). OG0000000 (598): 11 h 4 min. They set the end of the run |
| Final steps (splitting, rooting, result folders) | about 1 hour |
| Final orthogroups | 24,275; resolved trees 17,162; alignments 19,941 |
| Peak memory | **3.3 GB** |
| Disk | **2.0 GB**: `WorkingDirectory/` 1.6 GB (IQ-TREE's files in `Alignments_ids/` 1.1 GB, DIAMOND results 146 MB for 25 species pairs, databases 75 MB), alignments 137 MB, orthogroup sequences 133 MB |

### Where IQ-TREE's time went (all 17,598 runs)

| Cluster size | Runs | IQ-TREE time | Share of the time | ModelFinder's share of it |
|---|---|---|---|---|
| 3 to 9 | 16,049 | 43 h | 16% | 67% |
| 10 to 49 | 1,454 | 63 h | 23% | 52% |
| 50 to 199 | 91 | 131 h | 48% | 47% |
| 200 to 700 | 4 | 38 h | 14% | 45% |
| **All** | 17,598 | **275 h** | | **51%** |

Over the whole run **ModelFinder is about half of IQ-TREE's time** (141 of 275 hours). The earlier figures
(17% and 46%, from the first 25 trees) understated it. 91 clusters of 50 to 199 sequences took half the
time; the 16,049 smallest took a sixth.

Time grows with about the square of the number of sequences: n^2.23 over all runs, n^1.95 to n^2.05 when
alignment length is held fixed.

## The 35-species search and clustering (2026-10-21)

The real Run 1 command (`-M msa -T iqtree3 -s <catalog tree>`, 32 threads) on the 35 proteomes from the
coordinator, stopped by a watcher at the start of the tree stage. Kept in `work/diagnostic_35/` (outside
the repository), for reuse with `-b`.

| Step | Result |
|---|---|
| DIAMOND search, 1,225 species pairs | 3 h 14 min (02:34 to 05:48) |
| Clustering (homologues, MCL) | 13 min |
| Preliminary clusters | 33,499 (848,665 proteins, all accounted for); **15,052 with 3+ sequences** get a tree |
| Cluster sizes | median 36 sequences (pilot: 5); largest 5,613 (pilot: 598) |
| Disk | 8.7 GB in `WorkingDirectory/` (mostly the 1,225 search results) |
| Checkpoint | all 1,225 search files pass `gzip -t`; MD5s of the search and mapping files in `search_manifest.md5`; made read-only; `proteomes_md5` recorded in `diagnostic.json` (`5df42f6efb77659a03c6a82d6d9c7b26`) |

The largest clusters are lineage-expanded immune and sensory families: OG0000000 (5,613 sequences)
NLRC3-like NACHT/LRR proteins, most of them in cichlids (Mzeb 434, Acal 369); OG0000001 (2,696)
TRIM16-like (fish TRIM expansions); OG0000002 (2,405) olfactory receptors, most in mouse (595).

### Estimate from the real sizes

| Clusters | Number | IQ-TREE time |
|---|---|---|
| 3 to 199 sequences | 14,727 | about 5,100 h, by the pilot's fit (sizes the pilot covered; scaled to its measured total) |
| 200 to 999 | 306 | about 3,000 to 5,000 h (anchored on the pilot's 576-sequence cluster, 12.5 h; growth n to n²) |
| 1,000 to 5,613 | 19 | most of the rest; the largest alone takes about 5 days (growth n), 16 days (n^1.5) or 50 days (n²) on one core |

- **On 32 cores: about 11 to 13 days of work in total**, and no less than the largest cluster's own time.
  If the largest clusters grow like n to n^1.5 (the pilot's largest ran 2 to 3 times faster than the
  trend), **the run takes about 2 to 3 weeks**; if they grow like n², the 5,613-sequence cluster alone would
  take about 7 weeks.
- ModelFinder was about half of IQ-TREE's time in the pilot; a fixed model would roughly halve all of these.
- These replace the scaled estimate below, which assumed every cluster grows 7 times.

### Threads per tree (benchmark, 2026-10-21)

The pilot ran 32 IQ-TREE jobs at once, each on one thread (OrthoFinder's command has no `-T`). Measured on
the pilot's largest alignment (OG0000000, 598 × 1,261), its model fixed (JTT+F+I+G4), one parsimony starting
tree, 3 search iterations, seed 1; the four runs side by side (29 threads). Kept in
`pilot_runs/thread_bench/`.

| Threads | Wall time | Speed-up | Efficiency | CPU time |
|---|---|---|---|---|
| 1 | 933 s | 1.0× | 100% | 933 s |
| 2 | 461 s | 2.0× | 101% | 883 s (run alone, 2026-10-22) |
| 4 | 375 s | 2.5× | 62% | 1,307 s |
| 8 | 283 s | 3.3× | 41% | 1,745 s |
| 16 | 212 s | 4.4× | 28% | 2,239 s |

Re-extrapolated (bulk, clusters under 500 sequences, single-threaded: about 7,900 core-hours, 10 days on 32
cores; the 65 clusters of 500+ given more threads, at the measured CPU cost):

| Giant growth | Giants at 1 thread | 8 threads | 16 threads |
|---|---|---|---|
| n | about 12 days | about 14 days | about 15 days |
| n^1.5 | about 16 days | about 16 days | about 17 days |
| n² | about 50 days (the largest cluster alone) | about 20 days | about 23 days |

Threads for the giants remove the worst case (one cluster holding the run for weeks) but do not shorten the
typical case: the bulk is about 10 days of work and the machine is already full. With the agreed settings
the run is about 2 to 3 weeks; going well below 2 weeks needs a cheaper model setting for the bulk. The
speed-up comes from one alignment and side-by-side runs, and is probably conservative for larger alignments.

### Threads per tree: the run's rule (2026-10-22)

The estimates above added up the work and spread it over 32 cores. That misses the order: smallest first, the
largest clusters under 1,000 sequences (up to about 7 days on one thread) start around day 11, so the second
interim result (up to 1,000 sequences) could not come before day ~15–19 whatever the thread counts for the
giants. `pilot_runs/scripts/schedule_sim.py` simulates the tree stage in list order on 32 cores, with the
wrapper's rule (a tree waits until its threads are free) and the measured speed-ups; `threshold_search.py`
searches thread rules for the earliest second interim result, then the earliest final.

2 threads cost no extra CPU time on the benchmark (the 2-thread run ran alone; the others side by side, so its
efficiency may be a little lower on a full machine). 8 and 16 threads waste 47% and 72% of the extra cores,
which delays everything queued behind them while the machine is full. The chosen rule, now the default:
`1;600:2;1500:8;3000:16` (32 clusters on 2 threads, 9 on 8, 1 on 16). Simulated, in days from the start
(0.79–1.0× the pilot fit; worst case: giants growing as n²), with `-m TEST` for the clusters over 1,000:

| Rule | 1st result (≤300) | 2nd result (≤1,000) | Final | Worst |
|---|---|---|---|---|
| `1;500:8;2000:16` | 6.7–8.5 | 17.1–21.7 | 19.6–24.9 | 30.7 |
| `1;1000:8;2000:16` | 6.7–8.5 | 15.0–18.9 | 15.8–20.1 | 25.6 |
| `1;600:2;1500:8;3000:16` | 6.7–8.5 | 12.9–16.3 | 14.3–18.1 | 23.2 |

With the full model search for the giants too (the agreed setting) the final is 16.8–21.3 (worst 31.3). Lower
bounds with every core busy: 6.6–8.3 for the first result, 11.5–14.6 for the second. Starting 2 threads at 100
sequences gives the same dates within 0.1 day, but would apply the 598-sequence speed-up to small alignments.

### Restarting an interrupted run (tests, 2026-10-21)

Kept in `pilot_runs/resume_test/`. **A**: IQ-TREE on OG0000000 (598 sequences), killed after 150 s, rerun with
the same command. **B**: the 2,183-protein subsample (pilot 2's input), 8 threads, killed in the tree stage
with 129 trees finished and 8 running, then restarted with `-b <WorkingDirectory>`.

| | Result |
|---|---|
| IQ-TREE resumes a killed run (same prefix) | yes: "CHECKPOINT: Resuming analysis", initial tree and model parameters restored (`.ckp.gz`; ModelFinder's tested models are in `.model.gz`) |
| `-b` reuses the search results | yes |
| `-b` reuses finished trees or IQ-TREE checkpoints | **no**: it writes a new `WorkingDirectory/OrthoFinder/Results_<date>/` and rebuilt all 129 finished trees |
| Clusters on a rerun | identical (only the `# cline` comment, which holds the path, differs) |
| Alignments on a rerun | identical, 137 of 137 byte for byte |
| Trees on a rerun | **108 of 129 differ**: IQ-TREE draws a new random seed each run (257529, then 51875, for OG0000000) |

So, as OrthoFinder runs IQ-TREE, an interruption loses the whole tree stage, and the same alignments give
different trees on each run. A tree wrapper (the same one that would give the giants more threads) can fix
both without changing the method: a fixed `-seed`, recorded in `run.json`; and a cache of trees and IQ-TREE
working files keyed by the alignment's checksum, so that a restart reuses finished trees and continues
half-finished ones from IQ-TREE's checkpoint. How different the re-drawn trees are was not measured.

### Recovering an interrupted run by hand (test C, 2026-10-21)

Kept in `pilot_runs/resume_test/C/`. Input: 3 clusters of 100 sequences and 40 of 4 to 20 from the pilot
(599 proteins). Run with a multithreaded IQ-TREE entry in `~/config_orthofinder_user.json`
(`iqtree3_mt`: OrthoFinder's command plus `-T METHODTHREAD`) and the hidden option `-mt 4`. Killed with 41 of
43 trees finished, OG0000001 in its tree search and OG0000000 in ModelFinder.

1. Each unfinished tree rerun by hand with OrthoFinder's own command, in its own `Alignments_ids/` (same
   prefix): both resumed from their checkpoints. OG0000001 restored its model in 0.04 s and finished in
   4.5 min. OG0000000 chose the same model as an uninterrupted copy (VT+F+R4), but its ModelFinder took
   26.0 min against 27.2 min fresh, although about 3 min had been done before the kill: a kill during
   ModelFinder loses most of ModelFinder's progress.
2. Each finished `<OG>.treefile` moved to `Trees_ids/<OG>.txt`, as OrthoFinder does.
3. `orthofinder -fgt <Results folder>` (hidden option, "from trees"): finished the workflow in 1 s from the
   43 trees, without running IQ-TREE, into a new `Results_<date>_1/`. All 49 resolved trees match their
   final orthogroups; the hand-finished trees became final OG0000000 and OG0000002 plus OG0000019.

Notes: IQ-TREE writes `<prefix>.treefile` during the search, so "finished" must be judged from the log's
final lines, not from the tree file. `-mt 4` with `-t 16` ran 16 jobs of 4 threads (64 threads): the option
does not limit the total. Hidden options found: `-og` (stop after the orthogroups), `-fg`, `-fgt`, `-ogt`,
`-mt`. Not yet tested: whether a recovered run gives the same result as an uninterrupted one (needs a fixed
`-seed`).

### The tree wrapper (2026-10-22)

`cichlid-orthofinder-tree` (`treecache.py`), registered as OrthoFinder's tree program `cichlid_tree` in
`~/config_orthofinder_user.json`: the same IQ-TREE 3 command (`-m MFP`, the built-in ModelFinder) with a fixed
`-seed`, threads by cluster size within a shared core budget, and a cache keyed by the alignment and the
settings. `run.py` sets it up (`--seed`, `--thread-rule`, `--cores`), records every attempt in
`run_settings.json`, and restarts with `--resume` (OrthoFinder `-b` from the finished search). Kept in
`pilot_runs/wrapper_tests/`.

**Tests.** 44 offline tests pass (14 for the wrapper with a fake IQ-TREE, 2 for the budget, 9 for the progress
report). With the real IQ-TREE, on alignments of 17 sequences, all 9 cases pass: new; reused (0.0 s, IQ-TREE
not run); same seed with 1 thread and with 4 threads gives byte-identical trees (and the same tree for both);
killed in ModelFinder (162 models done) and resumed (54 s instead of about 95 to 120, identical tree); killed in
the tree search and resumed (identical tree); finished but not recorded (recorded, not rerun); two at once
(IQ-TREE runs once); IQ-TREE failure (non-zero exit, no tree).

**Inside OrthoFinder** (`integration.sh`; the 599-protein input, thread rule `1;50:4`, 16-core budget each): run U
uninterrupted; run K killed after 11 trees, restarted with `-b`, killed again during a 100-sequence tree's
search, restarted with `-b` to the end. 10 of 13 checks pass:
- no tree computed twice (44 trees, each once);
- every tree finished before a kill was only reused after it (11 after kill 1, 42 after kill 2);
- of the 16 trees running at kill 1, 15 resumed and 1 (killed before its first checkpoint) rebuilt; both trees
  running at kill 2 resumed;
- orthogroups, gene counts, alignments, hierarchical orthogroups and orthologues byte-identical to U.
The three that fail are byte-identity of the raw trees (38 of 43), and so of the resolved trees and the
duplication table: 3 are the 100-sequence trees on 4 threads (run-to-run noise of multithreaded IQ-TREE); 2 are
single-threaded trees resumed after kill 1, with the same topology and branch lengths within 0.0006.

**Conclusion.** A crash costs only the work since each running tree's last checkpoint. Single-threaded trees
that run uninterrupted are byte-identical; resumed ones keep their topology with branch lengths equal to
numerical precision; multithreaded trees are reproducible to numerical precision and can settle on an
alternative tree of near-equal likelihood.

### Reusing the saved 35-species search (2026-10-22)

`cichlid-orthofinder run --from-search <WorkingDirectory>` copies a finished search into the new run and starts
OrthoFinder from it (`-b`), only if `search_problems()` finds all 1,225 files complete and `search_mismatches()`
finds every species' `Species<i>.fa` identical, sequence by sequence and in order, to the catalog's proteome
(IDs through `SequenceIDs.txt`); otherwise the run searches again. The saved search (Oct 21) used the same
programs (OrthoFinder 3.1.5, DIAMOND 2.2.8), the same proteomes (`proteomes_md5` 5df42f6e…) and OrthoFinder's
default search settings. Checks on the real search: complete, all 35 species match (about 1 minute).

Real test: the search copied (7.5 GB, 1 minute), then `orthofinder -b <copy> -M msa -og -t 32`. Clustering
took 13 minutes, as in the saved run. `clusters_OrthoFinder_I1.2.txt` (848,665 proteins, 33,499 clusters) and
`clusters_OrthoFinder_I1.2.txt_id_pairs.txt` (126,238 lines) are identical to the saved run's except the
comment line with MCL's command (the folder path). The saved folder was not changed (no file newer than the
test). Note: `-og` did not stop OrthoFinder 3.1.5 here; it went on to the alignments and trees (FastTree, its
default with no `-T`) and was stopped by hand once the orthogroups were written.

### Resuming after a crash during the search (2026-10-22)

`--resume` restarts OrthoFinder from the run's own search (`-b`). A run killed during the search (the first
~3.5 hours) leaves some `Blast<i>_<j>.txt.gz` missing or cut short, and `-b` would cluster from partial hits.
`run.search_problems()` now checks that every file of the 35 × 35 species pairs exists and reads to the end as
gzip; if any does not, the output folder is moved aside (`output_incomplete_<time>`, kept) and the run starts
from the beginning (recorded as `search_restarted` in the attempt). Finished trees in the cache are reused
either way. Checked on the saved 35-species search: 1,225 files, complete, 47 s. Tests:
`test_run.TestSearchCheck` (complete, missing, cut short, never started, set aside). The restart branch of
`run()` itself is exercised only through these pieces, not by a killed 35-species run.

### The core budget, measured live (2026-10-22)

38 real pilot alignments (6 to 86 sequences) started at once through `cichlid-orthofinder-tree`, real IQ-TREE,
a 16-core budget and the rule `1;20:2;40:8` (the run's 1/2/8 mix, scaled to small trees). A sampler summed the
`-T` of the running IQ-TREE processes every 0.2 s (1,945 samples): peak 16 of 16, never above; mean 12.8 while
trees ran, 14 to 16 cores in use 59% of the time. All 38 trees built (24 on 1 thread, 10 on 2, 4 on 8), none
failed, 571 s. Kept in the session scratchpad (`budget_check/`).

### The model setting by size, changeable during the run (2026-10-22)

`CICHLID_TREE_MODEL` gives IQ-TREE's `-m` by size, like the thread rule; the runner's default is `hybrid`
(`MFP;1000:TEST`: ModelFinder below 1,000 sequences, `-m TEST` from 1,000), `full` is `MFP` for every tree.
The rule is kept in `<run>/tree_settings.json` (`CICHLID_TREE_SETTINGS`), read again each time a tree starts,
so `cichlid-orthofinder set-model <run> full|hybrid|<rule>` changes it while OrthoFinder runs; trees already
started keep their setting. Smallest first, the clusters of 1,000+ start around day 11, so that is the time
to decide. The model is in the cache key and in each tree's `done.json`; MFP trees keep the keys they had
before the rule existed. An unreadable file falls back to the launch setting and is recorded in
`events.jsonl` (`settings_unreadable`). Tested offline (`test_treecache.TestModelRule`, `test_run.TestModelSetting`)
and with IQ-TREE 3.1.3 on two 5-sequence pilot alignments: the first built with `-m MFP`, the file changed,
the second built with `-m TEST`, no restart.

### Interim results while the largest trees run (2026-10-22)

`cichlid-orthofinder snapshot <run> <dest>` copies the run's working files (without the search results) with
only the finished trees, and runs `orthofinder -fgt` on the copy; the running run is untouched. Clusters
without a tree are left out (their genes reported unassigned). `cichlid-orthofinder map` then pairs the
snapshot's orthogroups with the final run's, by identical genes.

**A trap: a copied OrthoFinder folder still reads the original.** `-fgt` finds its files through absolute paths
in `Log.txt` (`WorkingDirectory_Base`, `WorkingDirectory_Trees`, `FN_Orthogroups`). On a plain copy with three
trees removed, `-fgt` silently used the original's trees and gave the full output. The snapshot rewrites them.
This applies to any copy, including `native/`.

**Tests** (`pilot_runs/wrapper_tests/snapshot_test/`): (A) a copy of the uninterrupted integration run with its
three 100-sequence trees removed, as mid-run: 40 trees, 3 clusters and 300 genes left out; all 45 remaining
orthogroups have exactly the genes of a final orthogroup, and the same resolved tree (45 of 45). OrthoFinder
numbers orthogroups by size, so every ID differs (snapshot OG0000000 is final OG0000003): the mapping is
needed. (B) a snapshot of the resumed run K with all trees: orthogroups, resolved trees and duplications
byte-identical to K's own output. Offline tests cover the path rewrite and the mapping.

### The order of the tree jobs (2026-10-22)

OrthoFinder hands out its alignment and tree jobs largest first (`cmd_order = "descending"`). It has an
"ascending" order in its code, but the option is disabled, and **it crashes with a user species tree (`-s`)**:
in `tools/trees_msa.py` (line 568) the ascending branch reads `orig_commands_and_filenames`, which is defined
only when OrthoFinder infers the species tree. Set through the launcher `cichlid-orthofinder-of` (which sets
the option without changing OrthoFinder's files), the run stopped at the start of the tree stage
(`pilot_runs/wrapper_tests/order_test/A.log`). **OrthoFinder exited with code 0 after the crash**: `run.py`
judges success by "OrthoFinder run completed" in the log, not the exit code. The launcher's default is
largest first.

With the wrapper's core budget, largest first runs the clusters of 500+ sequences first, a few at a time on
8 threads, while the rest of OrthoFinder's job slots wait behind them; the small clusters finish last. An
early snapshot of the small clusters therefore needs the order changed in the wrapper, not in OrthoFinder.

### Smallest first, made safe (2026-10-22)

`cichlid-orthofinder-of` loads OrthoFinder's `tools/trees_msa.py` with the one-word fix
(`orig_commands_and_filenames` -> `commands_and_filenames` in the ascending branch; the value it computes,
`tree_tasksize`, is passed to a job runner that does not use it), **in memory**: the file on disk is not
changed, and the launcher refuses to run if its MD5 is not the checked one (`e875a68a…`). Tested
(`patch_test.sh`, 306 proteins, all trees on one thread, seed 1): plain OrthoFinder (P), the patched launcher
largest first (L) and smallest first (A) give **all 139 result files and all 33 gene trees byte-identical**;
A started its trees in order of size. Smallest first is the default (`--order ascending`).

### Interim reports, and a completeness gate (2026-10-22)

`cichlid-orthofinder run` watches the run and, once every cluster of up to each size in `--snapshot-sizes`
(default 300 and 1,000) has its tree, takes a snapshot automatically; when the run ends it maps each snapshot
to the final results. A snapshot works on a copy and never pauses OrthoFinder; a failing one is recorded and
retried, never stopping the run. Snapshots taken are kept in `run_settings.json`, so a resumed run does not
repeat them.

OrthoFinder writes `WorkingDirectory_Trees` to `Log.txt` only after the tree stage (`orthologues.py`, line
528): the snapshot reads the running run's trees from its current Results folder instead.

**OrthoFinder finishes and reports success even when tree jobs fail.** With its tree jobs killed but itself
not, it logged `ERROR: external program returned code -9`, dropped those clusters (29 of 33 trees; 154 genes
reported unassigned) and wrote "Done writing files". `run.py` now accepts an attempt only with a tree for every
cluster of 3+ sequences and no failed external program; otherwise it records the missing trees and asks for
`--resume` (the cache rebuilds only what is missing).

**End to end** (`watcher_e2e.py`; 306 proteins; smallest first; reports at 10 and 30 sequences standing in for
300 and 1,000; one session): report 1 taken mid-run; OrthoFinder killed; resumed with `-b` (30 trees reused,
4 resumed, none rebuilt, every cluster with its tree); report 2 taken; the final orthogroups, hierarchical
orthogroups, orthologues, alignments and all 33 trees byte-identical to plain OrthoFinder; both reports'
families (29 and 30) found in the final results with identical trees. 17 of 17 checks pass. (Two earlier runs
of this test overlapped with a second Claude Code process on the same conversation, between 14:16 and 14:31,
and are not counted; every other test ran before the overlap.)

### Checking the cost model on real 35-species clusters (2026-10-22)

`timing_sample.py`: 24 clusters of the 35-species clustering (8 each of about 40, 55 and 70 sequences), aligned
and trimmed with OrthoFinder's own steps (`famsa -t 1`; `trim.main(fn, fn, 10, 0.1, 500, 0.75)`), built through
the wrapper on one thread each, 30-minute cap. 18 finished; the 6 capped had long alignments (788 to 3,152
columns). Measured over predicted time (the pilot model, n^2.23): median 0.79 over all 24 (0.79 at ~40, 0.55
at ~55, at least 1.04 at ~70, where 5 of 8 were capped); single clusters vary about two-fold either way, mostly
with alignment length. ModelFinder was 46% of the time. The model is kept uncorrected (factor 1.0): it is about
right on average and not optimistic.

## Estimate for the full run (35 species, agreed settings), before the 35-species clustering

With 35 species instead of 5, a cluster holds about 7 times as many sequences (848,665 proteins against
125,879); the number of clusters changes much less. Scaled from the pilot's 275 hours:

| If IQ-TREE time grows as | Cost per cluster | IQ-TREE time | On 32 cores | Largest single cluster (one core) |
|---|---|---|---|---|
| n (lower bound) | ×7 | about 1,900 h | about 3 days | about 4 days |
| n^1.5 | ×19 | about 5,100 h | about 7 days | about 10 days |
| n^2.2 (measured) | ×77 | about 21,000 h | about 27 days | about 40 days |

- **The run is likely to take weeks**, and its end is set by the few largest clusters, each on one core,
  whatever the number of cores. In the pilot the last cluster ran for 12.5 hours; in the full run the
  largest may run for weeks on its own.
- These are extrapolations from clusters of up to 600 sequences to clusters of several thousand. The
  35-species cluster sizes (DIAMOND search and clustering on all 35, a few hours, reusable with `-b`)
  would replace the 7× assumption with the real sizes.
- **Memory:** 3.3 GB at peak in the pilot. IQ-TREE's memory grows with sequences × columns, so the largest
  35-species clusters will need several times more per job; 31 GB is probably enough, not measured.
- **Disk:** roughly 20 to 30 GB (DIAMOND results for 1,225 species pairs about 7 GB; IQ-TREE's files and the
  alignments about 7 to 10 times the pilot's). 101 GB is free on `/`.
- **Where:** this workstation can hold it, but would be occupied for weeks. That, the agreed settings, and
  whether to run the 35-species clustering first, are for the owner.

## Pilot 2: what OrthoFinder 3.1.5 writes

Finished in 26 min 45 s on 4 threads (`-t 4 -a 4`), alongside pilot 1; peak memory 165 MB; 34 MB of
output, 27 MB of it `WorkingDirectory/`. 2,183 proteins, 542 orthogroups, 336 with 4 or more genes.

### Folder structure

```
<-o>/Results_Oct20/
├── Citation.txt, Log.txt                  Log.txt: version, command line, programs, species numbering
├── Orthogroups/                           Orthogroups.txt, Orthogroups.tsv, Orthogroups.GeneCount.tsv,
│                                          Orthogroups_UnassignedGenes.tsv, Orthogroups_SingleCopyOrthologues.txt
├── Orthogroup_Sequences/OG*.fa            every orthogroup, unaligned (542)
├── MultipleSequenceAlignments/OG*.fa      final orthogroups with 2+ genes (350)
├── Resolved_Gene_Trees/Resolved_Gene_Trees.txt   ONE file, one line per final orthogroup with 4+ genes (336)
├── Species_Tree/                          SpeciesTree_rooted.txt, SpeciesTree_rooted_node_labels.txt (= the -s tree)
├── Phylogenetic_Hierarchical_Orthogroups/ N1.tsv, N2.tsv, N3.tsv (no N0.tsv here; see WorkingDirectory/N0.ids.tsv)
├── Gene_Duplication_Events/               Duplications.tsv, SpeciesTree_Gene_Duplications_0.5_Support.txt
├── Orthologues/, Putative_Xenologs/, Phylogenetically_Misplaced_Genes/, Single_Copy_Orthologue_Sequences/
├── Comparative_Genomics_Statistics/
└── WorkingDirectory/
    ├── SpeciesIDs.txt, SequenceIDs.txt    the mapping files
    ├── Species<i>.fa, diamondDBSpecies<i>.dmnd, Blast<i>_<j>.txt.gz   search inputs and results (the bulk)
    ├── clusters_OrthoFinder_I1.2.txt      the preliminary MCL clusters
    ├── Sequences_ids/, Alignments_ids/    per preliminary cluster, internal IDs; Alignments_ids also holds
    │                                      IQ-TREE's own files (.log, .iqtree, .model.gz, .ckp.gz, ...)
    ├── Trees_ids/OG*.txt                  the IQ-TREE trees as inferred, per preliminary cluster (319)
    ├── N0.ids.tsv, SpeciesTree_rooted_ids.txt, SpeciesTree_unrooted_ids.txt, pickle/, dependencies/
```

There is **no `Gene_Trees/` folder** at the top level, unlike OrthoFinder 2.

### Mapping files

- `WorkingDirectory/SpeciesIDs.txt`: `0: Acal.faa`, one line per species, numbered **in alphabetical order of
  file name**, not the catalog's order.
- `WorkingDirectory/SequenceIDs.txt`: `0_0: Acal_113014445`, one line per protein (2,183 = the input
  count). Our IDs are kept whole, underscore included.

### Labels

| File | Label |
|---|---|
| `Orthogroups/*.txt`, `*.tsv` | our ID unchanged: `Acal_113012997` |
| `Resolved_Gene_Trees.txt`, `MultipleSequenceAlignments/`, `Gene_Duplication_Events/` | **species name added**: `Acal_Acal_113012997` |
| `WorkingDirectory/Trees_ids/`, `Alignments_ids/`, `N0.ids.tsv` | internal: `0_393` |

So the sheet's guess is right: without `-X`, OrthoFinder prefixes the species name. Stripping one
`<prefix>_` where the label starts with `<prefix>_<prefix>_` gives our ID back; `-X` would stop it at the
source (not used: not an agreed setting).

### The two tree sets, and the orthogroups (the main finding)

**OrthoFinder 3 builds its IQ-TREE trees on preliminary clusters, then splits those clusters into the final
orthogroups using the rooted trees, and renumbers them.** The `OG` numbers in `WorkingDirectory/Trees_ids/`
and `Alignments_ids/` are the preliminary clusters', not the final orthogroups':

- Raw tree `OG0000000` (20 genes) became the final orthogroups `OG0000012`, `OG0000065` and `OG0000107`
  (9 + 6 + 5 genes); the final `OG0000000` (17 genes) comes from another raw tree.
- None of the 319 raw trees has exactly the gene set of the final orthogroup with the same number. 23 raw
  trees were split into 2 or more final orthogroups. Every final orthogroup lies inside a single raw tree.
- 17 final orthogroups with 4 or more genes have no IQ-TREE run of their own (they are pieces of a split
  tree); they still get a resolved tree.
- In pilot 1, the `Orthogroups/Orthogroups.txt` written during the run (22,939) is the preliminary set; the
  final count will differ.

| | `WorkingDirectory/Trees_ids/OG*.txt` (as inferred) | `Resolved_Gene_Trees/Resolved_Gene_Trees.txt` |
|---|---|---|
| Per | preliminary cluster (different numbering and gene sets) | **final orthogroup** (gene sets match `Orthogroups.txt` exactly, 336 of 336) |
| Rooted | **no** (3 children at the root, as IQ-TREE writes) | **yes** (2 children at the root; rooted with the species tree) |
| Labels | internal (`0_393`) | `Acal_Acal_113012997` |
| Branch lengths | yes | yes |
| Supports | none (IQ-TREE runs without bootstrap) | none; internal nodes are named `n0`, `n1`, ... |
| Topology | IQ-TREE's | **differs from the raw tree for 69 of 336** (21%), compared unrooted on the orthogroup's genes |

**Consequence for the open decision.** The coordinator recommends the trees as inferred. In 3.1.5 those
are not one per final orthogroup: using them means either taking the families to be the preliminary
clusters (not the orthogroups in `Orthogroups.txt` and `assign/`), or cutting each final orthogroup's subtree
out of an unrooted raw tree. The resolved trees match the final orthogroups one to one and are rooted, but
OrthoFinder has rearranged a fifth of them with the species tree. For the owner and the coordinator.

### Alignments

- `MultipleSequenceAlignments/OG*.fa` is per final orthogroup (350, for orthogroups of 2 or more genes), rows
  labelled `Acal_Acal_…`; for the orthogroups with 4 or more genes the rows match the members exactly.
- Two 3-gene orthogroups (`OG0000339`, `OG0000340`) have an alignment file with a single row named `None`.
  They are below the size that gets a tree, so no family is affected here; the conversion should check it.

### Not yet seen

- The run without `-s` (pilot 3, small): whether it also builds a species tree from a concatenated alignment.

## Differences from the sheet

- `-T iqtree` is `-T iqtree3` in 3.1.5.
- `-o` is a parent folder: the results are in `<-o>/Results_<date>/`.
- The gene trees as inferred are per preliminary cluster, not per orthogroup (above); there is no
  `Gene_Trees/` folder.
- "Keep OrthoFinder's … supports": neither tree set has support values with the agreed settings.
- The species numbering in `SpeciesIDs.txt` is alphabetical, not the catalog's order.
- `SequenceIDs.txt` and `SpeciesIDs.txt`: names as the sheet expected, in `WorkingDirectory/`.

### The status tool during the full run (2026-10-22)

At 1.1 h into the tree stage `cichlid-orthofinder status` gave "about 127.9 days": it extrapolated C * n^B
(×1.57, calibrated on the smallest trees) to the 5,613-sequence cluster (117.7 days alone on 16 threads).
Fixed to use the plan's model (above 1,000 sequences: 12.5 h at 576 × (n/576)^1.5, halved for TEST), to
calibrate on trees of 20+ sequences once 30 have finished, and to show the date of each stage, simulated in
run order on the run's cores (`status.simulate`). The runner does not import `status.py`, so the change does
not touch the running run. First real timings (3,086 trees, 3–20+ sequences): 1.2–2.5× the model's time
for the smallest (fixed overhead dominates), 1.2× at 20+; the run used about 17–25 cores' worth of IQ-TREE
work on these. Stage dates at that point: 8.6, 16.4 and 18.3 days.

## The conversion on pilot 2 (2026-10-23)

`convert.py` on pilot 2's full output (5 species, 2,183 proteins), with the defaults (tree set `raw_pruned`,
alignment source `orthogroup`, minimum family size 4), made with coordinator 0.14.1: 359 orthogroups, 336
families written, 23 orthogroups below 4 genes (no folder; their genes are in `assign/` with their
orthogroup), 0 orthogroups of 4+ genes without a tree, 183 unassigned genes (empty `best_family`).
`cichlid-data validate` passes on every `assign` file, on `families/` and on `run.json`; two conversions are
byte for byte identical.

- **Clusters, trees and orthogroups line up.** Row k of `clusters_OrthoFinder_I1.2.txt_id_pairs.txt` is
  exactly `Sequences_ids/OG<k>.fa`, and its tree's leaves are exactly its members (512 of 512). Every cluster
  of 3+ sequences has a tree. Every final orthogroup lies in exactly one cluster (542 of 542, counting the 183
  single-gene orthogroups of `Orthogroups_UnassignedGenes.tsv`, which are the unassigned genes).
- **The alignment question.** For the 293 orthogroups of 4+ genes that are a whole cluster, and 5 of the 43
  that are a piece of a split cluster, `MultipleSequenceAlignments/OGx.fa` is the cluster's alignment
  (`Alignments_ids/`) restricted to the orthogroup's rows, with the all-gap columns removed. For the other
  **38 pieces of split clusters (11% of the families) it is a new alignment**, so with the sheet's alignment
  source the pruned raw tree and `combined.aln` come from different alignments. `--alignment raw_cluster`
  uses the cluster's alignment instead. Not decided: for the owner.
- **Small orthogroups.** Every orthogroup of 4+ genes has a `MultipleSequenceAlignments` file, but 6 of the
  17 two-gene and 3 of the 6 three-gene orthogroups (all pieces of split clusters) have none. The one-row
  `None` files (`OG0000339.fa`, `OG0000340.fa`) hold a Newick tree, not an alignment. A minimum family size
  below 4 would need the cluster's alignment for these.
- **The split-cluster check** (`native/summary/split_clusters.tsv` and `split_clusters_summary.txt`): 52 of
  359 orthogroups come from 23 split clusters, 43 of them with 4+ genes; in all 43 at least one species also
  has genes in another piece of the same cluster. OrthoFinder's splits almost always separate paralogs, which
  the pruned raw trees then show in different families.
- **Models** (`native/summary/cluster_models.tsv`): all 319 cluster trees ran ModelFinder (`MFP`): OrthoFinder's
  own command gives no `-m`. In the full run the model of each cluster comes from the tree wrapper's record.
- `-m MFP` tests up to 1,232 protein models, `-m TEST` up to 224 (IQ-TREE 3.1.3, one pilot 2 alignment).
- After checking against the whole pack: `native/summary/` also holds `conversion_counts.tsv` and
  `skipped_orthogroups.tsv` (on pilot 2: the 23 orthogroups of 2 or 3 genes); a rerun into the same output
  reuses the finished families (pilot 2: 336 reused, every file unchanged); the output paths come from the
  coordinator's `ResultsLayout`.
