# Integrated analysis of bulk and single cell RNA-seq

This software performs integrated analysis of bulk DNA-seq and single cell RNA-seq data. 

## Requirements to build the software
+ gsl
+ boost
+ cmake

## Implementation notes
The main class is `TSSBState<S,P>` class, which manages the nodes of the tree in the tree structured stick breaking process. 
It works with generics and accepts the type S of the data and P of the parameters.

The nodes of the tree are represented by an abstract class `Node<S,P>`. This class is to be extended based on the problem. 

A concrete class that implement `Node<S,P>` is `CloneTreeNode` representing a node of a cancer clone tree. 

To create an instance of `TSSBState`, first create a root node.

```
CloneTreeNode *root = CloneTreeNode::create_root_node("");
root->set_nu_stick(0.0);
root->sample_node_parameters(random, model_params, 0);
```

Note that the nu-stick is set to 0 as the root represents healthy cell population and hence, we do not want/expect any data to be assigned to the root. Then, create a trivial `TSSBState`:

```
auto state = TSSBState<BulkDatum, CloneTreeNodeParam>::construct_trivial_state(root);
```

This generates a state with root and empty data. To simulate bulk data, we provide two functions in `main.cpp`:
+ `sample_tree_from_prior` and
+ `generate_bulk_data`.

`sample_tree_from_prior` takes number of data points as input and samples the tree. It generates sticks (and hence the nodes) in the process. `generate_bulk_data` accepts `TSSBState` object. It generates copy number profile along the tree for each datum to determine the variance and total copy numbers. Then, it samples the number of variant reads using these numbers to parameterize a Binomial distribution. More details to follow. The example usage is as follows:
```
    unsigned long seed = 11;
    double alpha0 = 5;
    double lambda = 0.8;
    double gamma = 2;
    double seq_err = 0.001;
    double dir_conc_mult_factor = 1;
    double birth_rate = 0.1;
    double death_rate = 0.1;
    size_t n_data = 50;
    size_t max_cn = 6;
    size_t mean_depth = 1000;

    gsl_rng *random = generate_random_object(seed);

    ModelParams params(alpha0, gamma, lambda, seq_err, dir_conc_mult_factor);
    vector<BulkDatum *> data;
    TSSBState<BulkDatum, CloneTreeNodeParam> *true_tree = sample_tree_from_prior(random, n_data, params, data);
    generate_bulk_data(random, mean_depth, max_cn, birth_rate, death_rate, data, params, true_tree);
```

Finally, we run slice sampler 

```
    size_t n_mcmc_iter = 2000;
    size_t n_mh_iter = 20;

    TSSBState<BulkDatum, CloneTreeNodeParam> sampled_tree = run_slice_sampler(random, n_mcmc_iter, n_mh_iter, params, &data);
```

`run_slice_sampler` accepts true tree as an optional argument. If it is passed in, it will compute the evaluation metrics.

## Node<S,P>
It is a representation of a node of a tree. 
+ Stores a pointer to the parent node.
+ (Optionally) Stores a pointer to `TSSBState<S,P>` that manages the node.
+ Two concrete public constructors are provided to be called by the subclasses.
    - Both constructors accept child idx, so that the name can be represented.
```
Node(size_t child_idx, Node<S,P> *parent)
Node(size_t child_idx, Node<S,P> *parent, TSSBState<S,P> *state)
```

### Data
+ nu-stick for the node
+ node name as `vector<size_t>`
    - The name matches the parent's name except for the last entry. For example: if the parent name contains `[0, 0, 1]` then the node's name is `[0, 0, 1, j]`, where `j` is the child index.
+ `unordered_map<size_t, pair<double, Node<S,P> *> > idx2child`: stores children nodes along with the psi-sticks in a pair. The child idx (of type `size_t`) serves as a key to retrieve the child.
+ `P param`: stores a reference to the parameter object.
+ `unordered_set<S *> data`: stores pointer to data assigned to the node.

## BulkDatum
Represents single nucleotide aberration. 
+ Stores total number of reads and number of variant reads.
+ It also stores observed weighted average of total copy number.
+ For testing purposes, it also accepts the weight average of variant copy number.
    - Note that this is a latent variable, and we sum over variant copies.
    - If this value is set, `is_var_cn_observed` returns `true`.

## CloneTreeNodeParam
Represents node parameters of a clone tree. There are two values, cellular prevalence and clone frequency. This class provides getters and setters.
+ `is_consistent()` function checks to ensure that cellular prevalence is at least the clone frequency.

## CloneNode
It is a concrete implementation of `Node<BulkDatum, CloneTreeNodeParam>`.

### Functions
+ `spawn_child(double psi)`: creates a new child with psi-stick and returns the child.
+ `sample_node_parameters`:
    - For the root node, initially set cellular prevalence and clone freq to 1.0
    - Other nodes will sample clone freq and cellular prevalence from uniform distribution between 0 and parent's clone frequency.
    - Parent node's clone frequency is reduced accordingly.

## Resampling stick orders
In the Python implementation of TSSB provided by Adams et. al. (2010), the resampling of stick orders draws U ~ Uniform[0, 1], then determines the interval (child) that U falls in. If U falls in the interval not covered by existing branches (i.e., not explicitly represented), then they draw new psi-sticks and represent it. In other words, they perform a lazy computation. In our implementation, we compute the size of the intervals for each branch, normalize it and performs multinomial sampling. This ensures that no new branch is created when we resample stick orders. 

**It remains to justify that our approach does not break the sampler in any way.**

# Unit testing
The basic functionalities have been tested. More tests are to be added.
