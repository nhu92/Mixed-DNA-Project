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
    - HybPiper
    - Exon extraction
    - MAFFT + Fasttree
    - Distance Matrix + Prediction
7. Output:
    - HTML?
    - Reads Mapping Info
    - Merge exon trees into a single one for displaying
    - Rename the NODE but with a table to track the name
    - Key genes (?)
8. Consider multiprocessing using Python
9. For the **Hierarchical approach**, the first round predictions will be used as reference selection criteria
10. Naming system fix
    - Using the Python non-conflict form? Adding trailing and tail "__" to a project name?

To build up the hierarchical prediction, I will extend the current pipeline to run again with selected species according to the initial predictions.

1. From the output, select the corresponding order based on the threshold
2. Run ref_preparation script used to make for initialize the starting reference folder
3. Run MAFFT + FastTree + Distance
4. Family level prediction

Prepare some testing data samples:

- Some random Poales(3), Asterales(3), Brassicales(3), and Lamiales(3) from the Kew
- Resample the input reads
- Make mixes: single species, 3 mixes in 1 order, 2 mixes in 2 orders, 3 mixes in 3 order

Evaluating processes:

- Evaluate the order level threshold and the family level threshold together
- The evaluation will be shown in a plot where X-Y panel is two threshold choice and color scale showed the specific evaluating parameters
- Consider the threshold choice other than Z-score from normal distribution. Maybe Kruskal-Wallis test?

First, I did a sample run for 3 Rosales species. I will use `seqkit` to thin down the reads of each species down to 1000k. Then, I selected 30 genes and run an order prediction pipeline. I will later evaluate the order prediction results based on the threshold to run different family prediction line.

---

    1. Select corresponding families (refer to github notes)
    2. Prepare the mix
    3. Make auto pipeline down to family
        3.1 Check each step and let us customize all the folder names as output and input
        3.2 Handle output from the prediction.csv
    4. Make 3 levels of Order level selection (modify the selecting logic)
        4.1 From order results - Make a selection based on the empirical PPV - Output to a list
        4.2 
    5. Run the pipeline
    6. Evaluate the results

