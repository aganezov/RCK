# Adjacencies

### Contents: 
* [Adjacencies overview](#adjacencies-overview)
* [RCK adjacency format](#rck-adjacency-format)
    * [inferred clone- and haplotype-specific adjacency copy numbers](#inferred-clone--and-haplotype-specific-adjacency-copy-numbers)
* [Converting to RCK format from SV detection tools](#converting-to-rck-format-from-sv-detection-tools)
* [Processing RCK adjacencies](#processing-rck-adjacencies)

### Adjacencies overview
One of the key concepts in the RCK model is the notion of *adjacency*.
Adjacency is a transition between segment's extremities. 
There are two kinds of adjacencies:
* **reference** (present in the reference genome, or inherited by the derived cancer genome(s))
* **novel**  (present in derive genomes only).

Every adjacency `{(q,x,+|-),(p,y,+|-)}` describes a transition from (right|left) side of loci at coordinate `x` on chromosome `q` to (right|left) side of the loci at coordinate `y` on the chromosome `p`. 
Reference adjacencies naturally have a form of `{(chr,x,+),(chr,x+1,-)}` (i.e., on the same chromosome, neighbouring positions, and respective extremities via strands orientation).

We call two adjacencies reciprocal if some pair of their extremities resemble a reference adjacency.
For example, adjacencies `{(1,123450,+),(1,4567890,+)}` and `{(1,123451,-),(1,876534,+)}` are reciprocal, because extremity `(1,123450,+)` form first adjacency one, and extremity `(1,123451,-)` from second adjacency resemble a reference adjacency `{(1,123450,+),(1,123451,-)}`  

While sometimes novel adjacencies are classified as insertion, deletion, duplication, reversals (aka inversion), translocation, etc.
This is usually done by looking at chromosomes, coordinates, and strands of involved extremities.
For example insertion and deletion has the same *signature* (i.e., same chromosome, `+` strand on the leftmost extremity, and `-` on the rightmost extremity).
Duplication has a signature of `-` strand followed by the `-` strand on the same chromosome.
Reversal (event) usually involves two reciprocal novel adjacencies, with `+`,`+` strand signature on one, and `-`, `-` signature on another adjacency.

While, indeed, aforementioned annotation correspond to cases, where respective rearrangement events would produce such novel adjacencies, 
it can also be the case that novel adjacencies that resemble signatures described above can be produced by more complex rearrangement events, such as *chromoplexy* and *chromothripsis*.

### RCK adjacency format
RCK works with adjacencies in the following (T/C)SV (Tab/Comma Separated Values) text format:

````
aid	chr1	coord1	strand1	chr2	coord2	strand2	extra
````
where every entry thus describes an adjacency `{(chr1, coord1, strand1), (chr2, coord2, strand2)}` with an id of `aid`.
The `extra` field is designed to store `key=value` pairs of additional information about adjacencies, with entries being separated by `;`.

There are several special extra fields, that RCK relies upon, when working:
* `aid` -- copy of the `aid`. When reading adjacencies, id for the adjacency will be based on the column, not the extra field.
* `cn` -- copy number values (refer to the following [subsection](#inferred-clone--and-haplotype-specific-adjacency-copy-numbers))
* `at` -- adjacency type  (either `N` for novel (default), or `R` for reference). By default all adjacencies are considered to be noevl, unless the adjacency id starts with the lower-case `r`. 

#### inferred clone- and haplotype-specific adjacency copy numbers
The result sof the main RCK algorithm (via `rck` executable) contains the `rck.acnt.tsv` file, with entries following the RCK adjacencies format.

Both novel and reference adjacencies are output in the result, and depending on the ``--o-acnt-mix-novel-and-reference`` novel and reference adjacencies are either going to be mixed together, or separated with novel adjacencies followed by the reference ones.
While the adjacencies themselves are self-explanatory, the main important peace of information about them is the `cn` field in the `extra` column, that encodes the clone- and haplotype-specific copy number values.

The `cn` value is a python/JSON dict with the following structure:

```
{
  'clone_id': {
    'AA': int,
    'AB': int,
    'BA': int,
    'BB': int
  },
  ...
}
``` 
where `clone_id` corresponds to the clone, for which haplotype-specific copy numbers are provided, and the 
`AA`, `AB`, `BA`, `BB` entries encode the copy number of the (haplotype) labeled versions of the adjacency (where the first position is labeled with the first haplotype letter, and the second position is labeled with the second haplotype letter).

In the following example:
````
aid	chr1	coord1	strand1	chr2	coord2	strand2	extra
id1	1	123450	+	1	123760	-	cn={'c1':{'AA': 1, 'AB': 0, 'BA': 0, 'BB':0}, 'c2': {'AA': 2, 'AB': 0, 'BA': 0, 'BB':0}}
````
where for the novel adjacency `{(1,123450,+),(1,123760,-)}` with id `id1` the following labeled adjacency `{(1,123450,+,A),(1,123760,-,A)}` has a copy number 1 in clone `c1` and copy number 2 in clone `c2`.  

### Converting to RCK format from SV detection tools
RCK installation adds `rck-adj-x2rck` adjacency-conversion executable tool to the `PATH` of your installation environment.
With the help of `rck-adj-x2rck` one can convert (unlabeled) novel adjacency predictions from the following tools:

* *short-reads*
    * **Delly** [[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436805/) | [code](https://github.com/dellytools/delly)] 
    * **Manta** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/26647377) | [code](https://github.com/Illumina/manta)] 
    * **Lumpy** [[paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84) | [code](https://github.com/arq5x/lumpy-sv)] 
* *linked/barcode reads* 
    * **LongRanger** [[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4786454/) | [code](https://github.com/10XGenomics/longranger)] 
    * **GROC-SVs** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/28714986) | [code](https://github.com/grocsvs/grocsvs)]
    * **NAIBR** [[paper](https://www.ncbi.nlm.nih.gov/pubmed/29112732) | [code](https://github.com/raphael-group/NAIBR)]
* *long reads*
    * **Sniffles** [[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990442/) | [code](https://github.com/fritzsedlazeck/Sniffles)]
    * **PBSV** [paper | [code](https://github.com/PacificBiosciences/pbsv)]
    
The input is converted to the RCK adjacency format (described [above](#rck-adjacency-format)).
`rck-adj-x2rck` tries to retain as much as possible of extra information from the SV detection tools input during the conversion, and such information is stored in the `extra` column.

For help message of `rck-adj-x2rck`run:
````bash
rck-adj-x2rck --help
````

To get help in converting adjacency prediction from a specific tool `x` run:
````bash
rck-adj-x2rck x --help
```` 

The following optional command line arguments are shared for all of the input sources and can be beneficial to use:
* `-o` | `--output` -- output file (default is `stdout`)
* `--id-suffix` -- a suffix, that will be appended to every input adjacency id, when transforming to RCK format.
This can be beneficial when working with several input sources of adjacencies, and one would want to differentiate based on the source 

An example of converting adjacency prediction by `Sniffles` tool to the RCK suitable format:
````bash
rck-adj-x2rck sniffles SV_predictions.vcf --id-suffix sample-technology-sniffles -o sample-technology-sniffles.rck.adj.tsv
````
Which will convert SV prediction in `SV_predictions.vcf` produced by Sniffles in sample `sample` that was sequenced with technology `technology` into the RCK formatted adjacency calls in the file `sample-technology-sniffles.rck.adj.tsv`.

All converted adjacencies will have `id_sample-technology-sniffles`, where `id` is the VCF id provided by Sniffles.

Not that the `--id-suffix` value here is provided as an example and is not mandatory to (i) be prent at all (default value is empty string) (ii) be in the form of `sample-technology-method`, though we found it useful on several occasions. 
Depending on your needs a different suffix values may be more useful. 

###  Processing RCK adjacencies
RCK installation adds `rck-adj-process` adjacency processing executable tool to `PATH` of your installation environment.
For `rck-adj-process` the following commands are available:
* `cat` -- combining adjacencies from 1+ inputs into a single one
* `reciprocal` -- updating extremities of adjacencies in the input, so that pairs of extremities of distinct adjacencies that resemble reciprocality, but are tno exactly 1 bp apart, are brought together.
This option ran by default in the main `rck` executable, unless explicitly suppressed.

Running `rck-adj-process command --help` provides one with the help on usage of each particular command. 

