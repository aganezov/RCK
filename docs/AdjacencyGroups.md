# Adjacencies Groups

### Contents: 
* [Adjacencies Groups overview](#adjacencies-groups-overview)
* [RCK Adjacencies Groups format](#rck-adjacencies-groups-format)
* [Molecule Adjacencies Group](#molecule-adjacencies-group)


### Adjacencies Groups overview
3rd-generations sequencing experiments can produce groups of novel adjacencies for which we can infer additional, useful, information.
We assume that all (novel) adjacencies are provided in the RCK input, but then we allow additional "grouping" information.

### RCK Adjacencies Groups format
Adjacencies groups information is accepted into RCK workflow via the `--adjacency-groups` option that ust point to a file with RCK formatted adjacencies groups. 
RCK works with adjacencies groups in the following (T/C)SV (Tab/Comma Separated Values) text format:
```
gid	aids	extra
```
where every entry describes a subset of the adjacencies, that are part of the RCK input, with the group id of `gid`, and adjacencies ids in the group listed in a comma-separated fashion in the `aids` column.
Comma-separated values in the `aids` column must match entries in the `aid` column in the RCK input adjacencies file.
The extra field is designed to store `key=value` pairs of additional information about each adjacency groups, with entries being separated by `;`.

There are several special extra fields, that RCK relies upon, when working:

* `agt` -- adjacencies group type (`M` for [molecule groups](#molecule-adjacencies-group)).
* `fp` -- maximum false positive (fraction) value for the adjacencies group

### Molecule Adjacencies Group 
This type of adjacencies groups usually comes from the 3rd-generation sequencing experiments as either a group of adjacencies supported by a single long read (long read sequencing), 
or a predicted by short reads with the same barcode (10x Genomics sequencing data), or coming from a single cell experiment. 
One way or the other, all of the adjacencies in the molecule group come from a single clone: either the single (part) of the derived chromosome (long reads + barcoded cases), which is in a single cell, which represents a single clone,
or from a real single cell (single cell sequencing source), which, again, represents a single clone.

Every group in the input adjacencies groups file with the entry `agt=M;` in the `extra` column is treated as the molecule adjacency group.

So for every molecule group `U` comprised of `|U|` input adjacencies with a False Positive values of `f`, RCK forces that in at least one of the reconstructed clones, there will be at least `(1-f)*|U|` labeled  representations of novel adjacencies from `U` present.
We note, that such a constraint does not imply that only in one clone labeled realizations of adjacencies from `U` can be present. 
There may be several clones, in which different subsets of adjacencies from `U` have their labeled realizations present, but the constraints guaranties that in at least one clone, there will be `(1-f)|U|` of them.
       
