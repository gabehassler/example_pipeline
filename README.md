# example_pipeline

How to use this code:

  1. Run the `dependencies.jl` file in the `dep` folder (or just manually add the packages you need by looking for all the imported packages in `not_julia.jl`).
  2. Run the `not_julia.jl` file in the main folder in this repo.
      - If you want the pipeline to actually run, then you need to choose appropriately formated files for the __Import Data__, __Import Tree__, and __Import Plotting Labels__ buttons. To run an example, go to the `data` directory of this repository and select `yeast_continuous.csv` for __Import Data__, `yeast.txt` for __Import Tree__, and `yeast_labels.csv` for __Import Plotting Labels__.
  
That's it!

Note, to make the gui more foolproof, you would need to check that the uploaded files are properly formatted. I didn't do that for the sake of simplicity (and I was lazy), but it's probably best that the gui gets mad at the user rather than Julia raising its own errors when it can't handle improperly formatted files.
