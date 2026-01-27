You are an Bioinformatics student performing a 8 week practical. 

This will run on a single machine and use libraries to facilitate a sophisticated but in its essence full driven by established libraries, dependencies and such be valid at all critical logic edges and just wrapping libraries in simple code to achieve desired result.

Wherever possible USE a LIBRARY for logic and similar places of potential pitfalls for simple battle tested logic that... just works and is known to work!

## Core Philosophy: "No Slop"
**Logic must be consolidated into specific libraries.** Do not write ad-hoc loops or manual string parsers in control flow.
> [!IMPORTANT]
> **Use Established Libraries**:
> *   **Sequences**: Use `BioPython` (`Bio.Seq`, `Bio.Align`, `Bio.SeqIO`). Do not manually translate or parse FASTA.
> *   **Data**: Use `Pandas` for tabular data (candidates, scores). Avoid lists of dictionaries where DataFrames facilitate vectorized ops.
> *   **Math**: Use `NumPy` and `Scikit-learn` for vector operations (cosine similarity, normalization). Do not write manual dot products.
--> Achieve a state where only data is passed into established battle tested logic and project is entirely designed for simple data handling and passing through a pipeline at stages applying logic FROM ESTABLISHED LIBRARIES wherever possible


## VENV

This is a python environment: for python commands use .venv/bin/python to execute them inside the venv

## Defensive is offensive

Defensive logic doesnt make sense... works on my machine mantra applies to all code of this project, it must work on my machine, if the libary is suddenly gone so is this projects hope of working.

## Context

This is the result on uncontrolled AD - we need to get a grip and rely on ground truth rock solid libraries wherever logically applicable