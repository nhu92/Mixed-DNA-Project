# Working in a Hierarchical Way

The take the output from the order level search into an input for the second round search. Each time, we need to generate a new reference file subset from the original huge reference or initialize a new user customized reference. This chapter is for notes on building up the looping hierarchical method, testing with multiple input data, and evaluating the pipeline reliance.

---

First, let me collect the meeting notes from last discussion of the pipeline.

1. The pipeline will be a couple of python scripts, with minor bash commands required(?)
2. The dependencies will be installed by conda (what about the python libs?)
3. 