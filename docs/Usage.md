# Usage

### Contents: 
* [Input data](#input-data)
    * [Mandatory](#mandatory)
    * [Optional](#optional)
* [Running RCK](#running-rck)
    * [preprocessing options](#preprocessing-options)
    * [running options](#running-options)
* [Examples](#examples)

### Input data

#### Mandatory
RCK expects two mandatory peaces of the input:
* `--scnt` - clone- and allele-specific segment copy number predictions, obtained from 3rd-party tools (see [segments docs](Segments.md#converting-to-rck-format-from-clone--and-allele-specific-inference-tools) for more details).
* `--adjacencies` (unlabeled) novel adjacencies (aka SVs) obtained from 3rd-party tools (see [adjacencies docs](Adjacencies.md#converting-to-rck-format-from-sv-detection-tools) for more details).

Both inputs must be i the RCK format (refer to [segments](Segments.md#rck-segments-format) and [adjacencies](Adjacencies.md#rck-adjacency-format) docs on the formatting issues.)

#### Optional

* `--adjacency-groups` - Adjacencies groups (see [adjacencies groups docs](AdjacencyGroups.md) for more details).
* `--clone-ids` - a comma-separated list of clone ids (as present in the `scnt` input).
* **advanced** `--telomere-positions` and `--telomere-segments` -  Telomeres (either via exact locations, or via segments, for which all spanned extremities will be considered as possible additional telomeres)
* **advanced** `--fragments` - Fragments that span segments. Only works properly if no preprocessing on input `scnt` is performed.

### Running RCK

Running RCK inference algorithm is achieved through the `rck` executable (which is automatically added to your `PATH` with RCK installation).

#### preprocessing options
When running RCK, a lot of input preprocessing options, that are achieved via RCK utilities.

All preprocessing can be turned off via the `--no-pre` flag, but this is an advanced option, use with caution.

All the `--pre-scnt-xxx` flags refer to `--scnt` input clone- and allele-specific segment copy number values (in RCK format) preprocessing, that is similar to the `refine` command on the `rck-scnt-process` tool.

All the `--pre-scnb-xxx` flags refer to creating copy number boundaries for the inferred copy number values. 
By default the strategy for obtaining the copy number boundaries is the uniform min-max one, where regardless of the input clone- and allele-specific copy number, the lower bound is set to the `--pre-scnb-uniform-min` value (default 0) and the upper is set to the `--pre-scnb-uniform-max` value (default 10). 
When working with genomes with known highly amplified segments, one can think about altering the default value for the `--pre-scnb-uniform-max`.

Preprocessing of the input adjacencies concerns the reciprocality adjustments, and is similar to the `rck-adj-process reciprocal` command and the `--pre-adj-xxx` option mirror those of the `rck-adj-process reciprocal`.
Adjacency preprocessing can be turned off by specifying the `--pre-no-adj`.

#### running options

One of the main arguments in running `rck` is the `--workdir` option, specifying the working directory in which the three following directories are created:
* `raw_input` - contains exact copies of the input files
* `input` - contains fully preprocessed data
* `output` - contains inference results from the RCK algorithm  

Running options for `rck` start with the `--run-` prefix.
Running RCK without actually executing the inference algorithm can be achieved by using the `--no-run` flag.
This will prevent the actual gurobi based ilp solving and respective karyotype inference, but will preprocess (unless disabled) of all the input, and putting the preprocessed data into the `workdir/input` directory.

The `--run-g-` are the flags corresponding to setting Gurobi related options:
* `--run-g-mip-gap` - the gap between the best bound and best objective, after which the Gurobi solver will stop crunching numbers (default: 0.015, or 1.5% difference)
* `--run-g-time-limit` - the maximum time (in seconds) for gurobi to run, before stopping execution and taking the current best objective as the result (default: 28800, aka 8 hours)
* `--run-g-threads` - number of threads gurobi will use (deault: 4)
* `--run-g-allow-interrupted` - allow for gurobi run to be interrupted and still use the best obtained objective for the inference result

Other flags:
* `--run-nas-fp` - default False Positive upper bound (i.e., at most a `--run-nas-fp` fraction of input novel adjacencies can be *not* used in the inferred karyotypes). Default is 0.1
* `--run-group-m-default-fp` - default False Positive values for *molecule* adjacencies groups (unless explicitly specified in with the `fp` value in the `extra` field). Default is 0.1
* `--run-segment-length-attr` - an choice based attribute that is used to get the segments length. Default is `length_100` which means that for every segment of length `l` an `ceil(l/100)` value is used in the inference minimization. 


### Examples

The following command runs RCK inference on the clone- and allele-specific segment copy number tensor (stored in the `input.rck.scnt.tsv`)
