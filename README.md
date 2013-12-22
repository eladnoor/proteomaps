proteomaps
==========

Scripts and data for generating proteome treemaps (AKA proteomaps)

------------------------------------------
File descriptions
------------------------------------------

  - ./data/KO_gene_hierarchy/KO_gene_hierarchy_changes.csv  : Table with our own gene reannotations
  - ./data/KO_gene_hierarchy/KO_gene_hierarchy_colors.csv   : Table with colors for 2nd-level categories
  - ./data/KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping/[ORGANISM]_mapping.csv :
    contains organism-specific ID mapping from KO to the relevant genomic annotations
    with the following columns [Systematic name] // [Gene name] // [KO number]
    
-------------------------------------------
Procedure for drawing a proteomap using Paver
-------------------------------------------
File -> New -> Treemap -> Simple format
    select the Hierarchy file (proteomaps/res/[ORGANISM]_hierarchy.tms)

Data -> Load Cell Sizes -> CSV Cell Size Data
    select the abundance data file (proteomaps/data/proteomic_data/[ORGANISM]_[CONDITION]_cost.csv)
    
Colors -> Load Colors
    select the color scheme table (proteomaps/data/KO_gene_hierarchy/KO_gene_hierarchy_colors.csv)
    
File -> Save Visualization
    res/[ORGANISM]_[CONDITION]_cost.xml

If the organism has more than one condition, use the first condition as a template
and load the cell sizes on top of that in order to have similarly placed cells.
