# Working in a Hierarchical Way

The take the output from the order level search into an input for the second round search. Each time, we need to generate a new reference file subset from the original huge reference or initialize a new user customized reference. This chapter is for notes on building up the looping hierarchical method, testing with multiple input data, and evaluating the pipeline reliance.

---

First, let me collect the meeting notes from last discussion of the pipeline.

1. The pipeline will be a couple of python scripts, with minor bash commands required(?)
2. The dependencies will be installed by conda (what about the python libs?)
3. Consider NextFlow?
4. Angiosperm353 target file and reference file should be included in the final package
5. Need a well written tutorial of potential bash commands
6. Separate into stages. Reporting system. Clean up the intermediate files.
    HybPiper
    Exon extraction
    MAFFT + Fasttree
    Distance Matrix + Prediction
7. Output:
    HTML?
    Reads Mapping Info
    Merge exon trees into a single one for displaying
    Rename the NODE but with a table to track the name
    Key genes (?)
8. Consider multiprocessing using Python
9. For the **Hierarchical approach**, the first round predictions will be used as reference selection criteria
10. Naming system fix

