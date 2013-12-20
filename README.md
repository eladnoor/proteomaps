proteomaps
==========

Scripts and data for generating proteome treemaps (AKA proteomaps)



------------------------------------------
General input data for creating proteomaps

Please do not change! 
Original versions are still maintained by Wolf
------------------------------------------


Directory KO_gene_hierarchy: Data about the KEGG Pathways Hierarchy

  o File KO_gene_hierarchy/KO_gene_hierarchy_general.tms    Original KEGG tree in tms (proteomaps) format
  o File KO_gene_hierarchy/KO_gene_hierarchy_changes.csv    Table with our own gene reannotations
  o File KO_gene_hierarchy/KO_gene_hierarchy_colors.csv     Table with colors for 2nd-level categories


Directory KO_gene_hierarchy/KO_gene_hierarchy_organism_mapping

  o File [ORGANISM]_mapping.csv contains organism-specific ID mapping 
    [Systematic name] // [Gene name] // [KO number]
    
    
-------------------------------------------
Procedure for drawing a proteomap using Paver
-------------------------------------------
File -> New -> Treemap -> Simple format
    select the Hierarchy file (hierarchy_data/hierarchy_standardised.tms)

File -> Extend -> Path Format
    check the "Remove none extended cells" checkbox
    select KEGG to Gene ID file (hierarchy_data/*.csv)

Data -> Load Cell Sizes -> CSV Cell Size Data
    select the abundance data file (protein_data/*.csv)
    
Colors -> Load Colors
    select the color scheme table (hierarchy_data/KO_colour_table.tsv)
    
File -> Save Visualization
    *.xml
    
    
If the organism has more than one condition, use the first condition as a template
and load the cell sizes on top of that in order to have similarly placed cells.
