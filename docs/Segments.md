# Segments

### Contents: 
* [Segments overview](#segments-overview)
* [RCK segments format](#rck-segments-format)
    * [inferred clone- and haplotype-specific segment copy numbers](#inferred-clone--and-haplotype-specific-segment-copy-numbers)
    * [input clone- and allle-specific segment copy numbers](#input-allele--and-clone-specific-segment-copy-numbers)
* [Converting to RCK format from clone- and allele-specific inference tools](#converting-to-rck-format-from-clone--and-allele-specific-inference-tools)
* [Processing RCK segments](#processing-rck-segments)
* [Segments vs fragments](#segments-vs-fragments)

### Segments overview

One of the key concepts in the RCK model if the notion of *segment*.
A segment `(chr,start,end)` represents a continuous part of the reference genome's chromosome `chr` starting at `start` and ending at `end` (both inclusive).
A segment `s=(chr,start,end)` naturally has two extremities `(chr,start,-)` and `(chr,end,+)` corresponding to tail `s^t` and head `s^h` of the segment `s`.
In a diploid reference genome every segment `(chr,start,end)` (except for segment on sex chromosomes) has two haplotype-specific copies `(chr,start,end,A)` and `(chr,start,end,B)` respectively.  

Reference adjacencies correspond to pairs of adjacent extremities of consecutive segments.
For example, two consecutive segments `a=(chr,10001,20000)` and `b=(chr,20001,30000)` determine a reference adjacency `{(chr,20000,+),(chr,20001,-)}`.
Naturally for every chromosome that has two homologous copies, for every unlabeled reference adjacency `{(chr,20000,+),(chr,20001,-)}` there are two labeled reference adjacency counterparts:
`{(chr,20000,+,A),(chr,20001,-,A)}` and `{(chr,20000,+,B),(chr,20001,-,B)}`.

### RCK segments format
RCK works with segments in the following (T/C)SV (Tab/Comma Separated Values) text format (similar to that of bedpe):
````
chr	start	end	extra
````
where every entry describes a segment `(chr,start,end)`.
The `extra` field is designated to store `key=value` pairs of additional information about segments, with entries being separated by `;`.

There are several special extra fields, that RCK relies upon, when working:
* `cn` -- clone and allele/haplotype-specific copy number values of the segment (refer to the following [subsection](#inferred-clone--and-haplotype-specific-segment-copy-numbers))
* `cnb` -- clone and allele/haplotype-specific copy number boundaries of the segment (refer to the respective [subsection](#copy-number-boundaries))

#### inferred clone- and haplotype-specific segment copy numbers
The result of the main RCK algorithm (via `rck` executable) contains the `rck.scnt.tsv` file, with entries following the RCK segments format.
While the segments themselves are self-explanatory, the main important peace of information about them is the `cn` field in the `extra` column, that encode clone- and haplotype-specific copy number values.

The `cn` value is a python/JSON dict with the following structure:
```
{
  'clone_id' : {
    'A': int,
    'B'; int
  },
  ...
}
``` 
where `clone_id` corresponds to the clone, for which haplotype-specific copy numbers are provided, with the `A` and `B` entries encoding the copy number of the multiplicity of the corresponding haplotype-specific segments.

In the following example:
````
chr	start	end	extra
1	10000	20000	cn={'c1':{'A': 1, 'B': 2}, 'c2':{'A': 3 'B': 2}}
````
where for the segment `(1,10000,20000)` its `A` haplotype-specific version has `1` copy in clone `c1` and `3` copies in clone `c2`, and its `B` haplotype-specific version has `2` copies in clone `c1` and `2` copies on clone `c2`.

#### input allele- and clone-specific segment copy numbers
The input for the `RCK` method expects clone- and *allele*-specific (approximate) segment copy numbers. 
While most methods follow the notion of `major` and `minor` alleles (based on the segment copy numbers), we employ the same format of `cn` field, as was described in the previous [subsection](#inferred-clone--and-haplotype-specific-segment-copy-numbers)
The only and major difference is, that while the RCK output is haplotype-specific (i.e., `A` and `B` are matching and the same for every segment), the input is allele-specific (i.e., `A` and `B` entries do not necessarily match and can be"flipped").

### Converting to RCK format from clone- and allele-specific inference tools
RCK installation adds `rck-scnt-x2rck` segment copy number conversion executable tool to the `PATH` of your installation environment.
With the help of `rck-scnt-x2rck` one can convert clone- and allele-specific prediction from the following tools:
* **HATCHet** [[paper](https://www.biorxiv.org/content/early/2018/12/17/496174) | [code](https://github.com/raphael-group/hatchet)] (*recommended* as it has fewest limitation w.r.t. tumor heterogeneity)
* **TitanCNA** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/25060187) | [code](https://github.com/gavinha/TitanCNA)]
* **Battenberg** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/22608083) | [code](https://github.com/cancerit/cgpBattenberg)]
* **ReMixT** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/28750660) | [code](https://bitbucket.org/dranew/remixt)]

For help message of `rck-scnt-x2rck` run 
````bash
rck-scnt-x2rck --help
````

To get help in converting clone- and allele-specific predictions from a specific tool `x` run:
````bash
rck-scnt-x2rck x --help
````

### Processing RCK segments
RCK installation adds `rck-scnt-process` segment copy number processing executable tool to the `PATH` of your installation environment.
For `rck-scnt-process` the following commands are available:
* `align` -- aligning segments (and corresponding segment copy number tensors) form 1+ segment copy number tensord
* `refine` -- filling the missing spans in entries, of merging consecutive entries that have the same clone- and allele/haplotype-specific copy numbers.
This option ran by default in the main `rck` executable, unless explicitly suppressed.

Running `rck-scnt-process command --help` provides one with the help on usage of each particular command.


### Segments vs fragments
A segment and a fragment are of the same nature: a consecutive part of the reference chromosome. 
Input clone- and allele-specific copy number are usually inferred on rather large spans, which we treat as fragments. 
Such fragments are further fragmented into smaller, actual segments, based on the novel adjacencies, in such a way that novel adjacencies involve segments' extremities.
With smaller segments, we still retain information about allele-separation based on fragments, which span smaller segments. 
