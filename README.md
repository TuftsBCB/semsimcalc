# semsimcalc.py

Requires module `networkx`. Tested with verion `1.9.1`
Requires module `pickle`
Requires modules `sys`, `time`, and `math`, which should be installed with python by default.
Requires module `numpy`

See [http://bib.oxfordjournals.org/content/13/5/569.full](http://bib.oxfordjournals.org/content/13/5/569.full) for definitions metric definitions.

## parse_go_file(go_file_name)

This function takes in the file name for a GO ontology file (obo format). The file must be of a specific format:

```
  ! comments

  [Term]
  id: GO_term
  ...
  is_a: GO_term
  is_a: GO_term

  [Term]
  id: GO_term
  ....


  [Typedef]
  ...
```

Note: The `[Typedef]` tag signals the end of GO terms, and is required (otherwise, the parser will fail to record the final GO term in the ontology file)

Please see `example_go.obo` for a full example file.


`parse_go_file` returns two objects as a tuple:
   
  1. `go_graph`:
    A networkx DiGraph object. `go_graph` represents the ontology as a DiGraph, where each `is_a` relationship is represented as an edge.

  2. `alt_ids`:
    GO ontology files provide alternate IDs for some terms (represented by `alt_id:` lines in the ontology file).
    `alt_ids` is a mapping from alternate IDs to the IDs stored in `go_graph`.
    `alt_ids` is a python dictionary, where keys are GO IDs not in `go_graph`, and values are corresponding GO IDs in `go_graph`

---

## parse_annotation_corpus(ac_file_name, alt_ids=None)

This function takes in a file name for a pre-processed annotation corpus file of a specific format:

```
  -
  protein_name
  GO_term
  GO_term
  GO_term
  -
  ...
  -
```

Note: File must both start and end with a `-`line. Please see `example_corpus.stripped` for a full example of a pre-processed annotation corpus file.

`parse_annotation_corpus` returns two objects in a tuple:

  1. `prot_to_gos`:
    This is a python dictionary mapping protein names to GO terms. Keys are protein names, values are python lists of GO terms associated with the key (from the annotation corpus)
  2. `go_to_prots`:
    This is a python dictionary mapping GO terms to protein names. Keys are GO terms, and values are python lists of protein names labeled with the key (from the annotation corpus)

If `alt_ids` is provided, then any keys in `alt_ids` that appear in the annotation corpus will be stored as their associated values in `alt_ids`.

---

# load_semsimcalc(saved_path)

This function takes in a file path to a pickled `SemSimCalculator` object.
It returns an unpickled `SemSimCalculator` object

---
---

# SemSimCalculator class

The `SemSimCalculator` class takes an ontology and an annotation corpus. It parses and uses these to calculate various semantic similarity metrics between terms, groups of terms, and proteins.

All class variables are technically public, but should be treated as private. Use the getter functions explained below to access them.
Class variables:

* `_go_graph`
  networkx DiGraph representing the GO ontology, as parsed/returned by `parse_go_file`

* `_alt_list`
  python dictionary represting alternate GO term IDs, as parsed/returned by `parse_go_file`

* `_prot_to_gos`
  python dictionary mapping protein names to their GO term labels, as parsed/returned by `parse_annotation_corpus`

* `_go_to_prots`
  python dictionary mapping GO terms to the proteins which they label, as parsed/returned by `parse_annotation_corpus`

* `_proteins`
  python list of proteins names. Contains names of all proteins that have labels

* `_num_proteins`
  integer. Size of `proteins`

* `_ic_vals`
  dictionary mapping GO term to its IC (information content) value. Initialized empty. Used for memoization

* `_go_terms`
  list of all GO terms in the graph of the ontology

* `_mica_store`
  reference to a `MicaStore` instance. Initialized as `None`, must be set manually

---

## Initialization and variable access

### Class initialization - `__init__(self, go_file_name, ac_file_name)`

Creates new instance. Call `semsimcalc.SemSimCalculator(file_name, file_name)` to use. Will return a `SemSimCalculator` object.

Takes in file names for GO ontology file (obo format) and annotation corpus file (pre-processed file of the same format that `parse_annotation_corpus` takes, as explained above).

Initializes `go_graph`, `alt_list`, `prot_to_gos`, `go_to_prots`, `proteins`, and `num_proteins`.
Creates `ic_vals` as an empty dictionary.

### `link_mica_store(self, mica_store)`

Saves a reference to a `MicaStore` instance

### `unlink_mica_store(self)`

Removes `_mica_store` reference (sets to `None`)

### `save(self, filepath)`

Pickles and saves `self` to `filepath`

### `get_go_graph(self)`

Return copy of `_go_graph`

### `get_alt_list(self)`

Return copy of `_alt_list

### `get_ptg(self)`

Return copy of `_prot_to_gos`

### `get_gtp(self)`

Return copy of `_go_to_prots`

### `get_ic_vals(self)`

Return copy of `_ic_vals`
Note: `get_ic_vals` does not inherently calculate IC values. Use `precompute_ic_vals` first if you need all IC values.

### `get_go_terms(self)`

Return copy of `_go_terms`

### `get_mica_store(self)`

Return *reference* to `_mica_store`

---

## Basic calculations

### `calc_term_probability(self, term)`

Takes in a GO term as a string.
Calculates and returns the probability of that term or any of that term's descendants (in the GO DiGraph) occuring in the annotation corpus.
That is: [number of proteins labeled with term or a descendant of term] / [number of labeled proteins in annotation corpus]

### `IC(self, term)`

Takes in a GO term as a string.
Calculates and returns the information content of that term.
Information content is defined (within this implementation) as:
```
	-ln(prob(term))
```
Where prob(term) is the same as the result of calling `calc_term_probability(term)`

Note: Once an IC is calculated, it is stored in `_ic_vals`. Subsequent calls for the IC of the same term only look up the recorded value.

### `precompute_ic_vals(self)`

Fills the `_ic_vals` dictionary used for information content memoization.
Runs `IC` on all terms in the GO ontology.

### `MICA(self, left, right)`

Takes in two GO terms as strings. Order doesn't matter.
Calculates and returns the Maximum Informative Common Ancestor. (Returns a GO term as a string)

The MICA of two terms is the common ancestor of both terms with the highest information content value.

Note: For this implementation, if left and right are the same, they are included in the list of "common ancestors."

If a `MicaStore` instance is linked (through `link_mica_store`), `MICA` first queries the `MicaStore` instance. Only if the `MicaStore` instance does not return a GO term does `MICA` calculate a result from the GO graph and annotation corpus.

---

## Comparison Metrics

### `simRes(self, left, right)`

Takes in two GO terms as strings. Order doesn't matter.
Calculates and returns the resnik score of the two terms. (Returns a float)

simRes is defined as the information content of the MICA of two terms.
See [here](http://bib.oxfordjournals.org/content/13/5/569.full) for more details.

### `simLin(self, left, right)`

Takes in two GO terms as strings. Order doesn't matter.
Calculates and returns the Lin score for the two terms. (Returns a flot)

simLin is defined as the simRes of two terms divided by the sum of the information contents for each term (left and right).
See [here](http://bib.oxfordjournals.org/content/13/5/569.full) for more details.

**Note: Currently untested**

### `simJC(self, left, right)`

Takes in two GO terms as strings. Order doesn't matter.
Calculates and returns the Jiang-Conrath score for two terms is defined as:
```
  1 - IC(left) + IC(right) - 2 * simRes(left, right)
```
See [here](http://bib.oxfordjournals.org/content/13/5/569.full) for more details.

**Note: Currently untested**

---

## Mixing methods

Example for proper function calls:

Assume `calc` is a `SemSimCalculator` instance:
```
  calc.pairwise_average_term_comp(left_term, right_term, calc.simRes)
```
Will calculate and return the average of all pairwise resnik scores for the given lists of GO terms, `left_terms` and `right_terms`.

### `pairwise_average_term_comp(self, lefts, rights, metric)`

Takes in two python lists of GO terms (`lefts`, `rights`) and a comparison metric (ex. any function from the "Comparison Metrics" section). `metric` must take in two ontology terms and return a numeric score.
Returns the average of all pairwise term comparisons, using `metric`.

### `pairwise_max_term_comp(self, lefts, rights, metric)`

Takes in two python lists of GO terms (`lefts`, `rights`) and a comparison metric (ex. any function from the "Comparison Metrics" section). `metric` must take in two ontology terms and return a numeric score.
Returns the max of all pairwise term comparisons, using `metric`.

---

## Protein comparison

Example for proper function calls:

Assume `calc` is a `SemSimCalculator` instance:
```
  calc.average_protein_comp(left_prot, right_prot, calc.simRes)
```
Will calculate and return the average of all pairwise resnik scores for the go terms associated with the two protein names, left_prot and right_prot. 

### `average_protein_comp(self, left_prot, right_prot, metric)`

Takes in two protein names as strings (left_prot, right_prot) and a reference to a comparison metric (ex. any function from the "Comparison Metrics" section). `metric` must take in two ontology terms and return a numeric score.
Returns the average of all pairwise term comparisons for the GO terms associated with the proteins `left_prot` and `right_prot` using `metric`.

### `max_protein_comp(self, left_prot, right_prot, metric)`

Takes in two protein names as strings (left_prot, right_prot) and a reference to a comparison metric (ex. any function from the "Comparison Metrics" section). `metric` must take in two ontology terms and return a numeric score.
Returns the max of all pairwise term comparisons for the GO terms associated with the proteins `left_prot` and `right_prot` using `metric`.

---
---

# MicaStore class

Wrapper class for a numpy matrix of MICA values.

Takes a numpy matrix of MICA values and an ordering of GO terms (indices in the matrix).

Class variables:

* `_micas`
  `numpy` matrix of mica values

* `_go_to_index`
  dictionary mapping GO terms to indices in the matrix (matrix must be symmetrical)

### Class initialization - `__ini__(self, matrix_filename, ordering_filename)`

Loads `numpy` matrix from `matrix_filename` into `_micas`.
Parses ordered list of GO terms from `ordering_filename` (one term per line).

### `get_micas(self)`

Returns reference to `numpy` array `_micas`

### `get_ordering(self)`

Returns copy of `_go_to_index` dictionary 

### `get_index(self, term)`

If `term` is in `_go_to_index`, return `_go_to_index[term]`, which corresponds to `term`'s index in `_micas`
If `term` is not in `_go_to_index`, return `None`

### `mica_lookup(self, left, right)`

Attempts to look up a MICA value from `_micas`.
If MICA value cannot be found (or `left` or `right` are not in `_go_to_index`), returns `None`

---
---

# strip_ac.py

Standalone script to strip down a Swiss-Prot text file (".dat").
See [http://www.uniprot.org/downloads](http://www.uniprot.org/downloads) for download location.

Only tested on Swiss-Prot currently.

To run:

```
  python strip_ac.py -i [filename] -o [filename]
```
Where:
* `-i` takes the file name for the original Swiss-Prot .dat file
* `-o` takes the file name for the output (stripped) file.

The output of this script is compatible with `semsimcal.parse_annotation_corpus`


---
---

# Example files

These files are examples of format. `example_corpus.dat` corresponds with `example_corpus.stripped`, but not with `example_go.obo`
