# Script to create a reference Word document with 0.5 inch margins
# Run this once before rendering main_figures_word.Rmd

library(officer)

base_dir <- "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions"
ref_doc_path <- file.path(base_dir, "shiny_scripts", "reference_0.5in_margins.docx")

# Create a new document with custom margins
# Set page margins to 0.5 inches on all sides
# Note: page_mar() values are in inches by default
section_props <- prop_section(
  page_size = page_size(orient = "portrait"),
  page_margins = page_mar(top = 0.5, bottom = 0.5, left = 0.5, right = 0.5)
)

ref_doc <- read_docx() %>%
  body_add_par("Reference document with 0.5 inch margins", style = "Normal")

# Set the default section properties for the document
ref_doc <- body_set_default_section(ref_doc, section_props)

print(ref_doc, target = ref_doc_path)
cat("Reference document created at:", ref_doc_path, "\n")
