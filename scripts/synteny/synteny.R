# Load required packages
library(dplyr)
library(readr)
library(geneviewer)
library(htmlwidgets)
library(chromote)

dir.create("results/synteny", recursive = TRUE, showWarnings = FALSE)

# Read and prepare the CNR1 gene-neighborhood data.
synteny_data <- read_csv("data/synteny/synteny.csv", show_col_types = FALSE) %>%
  filter(receptor == "cnr1") %>%
  rename(
    cluster = taxa,
    name = gene,
    strand = orientation
  )

# Highlighted genes
goi <- "CNR1"
goi_color <- "#e21f4c"

ortholog <- c("AKIRIN1", "SRSF12", "PNRC1")
ortholog_color <- "#f28da5"

# Create gene color map
all_genes <- unique(synteny_data$name)

gene_colors <- setNames(
  rep("#e6e7e8", length(all_genes)),
  all_genes
)

gene_colors[goi] <- goi_color
gene_colors[ortholog] <- ortholog_color

# Label styles
label_styles <- list(
  list(index = 5, fill = "black", fontWeight = "bold"),
  list(index = 3, fill = "black", fontWeight = "bold"),
  list(index = 7, fill = "black", fontWeight = "bold"),
  list(index = 1, fill = "black", fontWeight = "bold")
)

# Set dynamic plot height
dynamic_height <- paste0(
  80 * length(unique(synteny_data$cluster)),
  "px"
)

# Create chart
main_chart <- GC_chart(
  synteny_data,
  cluster = "cluster",
  group = "name",
  width = 1000,
  height = dynamic_height
) %>%
  GC_normalize(
    group = "name",
    gap = 0.01,
    preserve_gene_length = FALSE
  ) %>%
  GC_tooltip(
    formatter = "<b>Start:</b> {original_start}<br><b>End:</b> {original_end}"
  ) %>%
  GC_coordinates(
    TRUE,
    showAverage = TRUE,
    rotate = 0,
    yPositionTop = 95,
    yPositionBottom = 25,
    overlapThreshold = 20,
    ticksFormat = ".4s",
    tickStyle = list(
      strokeWidth = 0,
      lineLength = 0
    ),
    textStyle = list(
      fill = "black",
      fontSize = "10px",
      fontFamily = "Arial",
      cursor = "default",
      x = -2
    ),
    fontSize = "6px"
  ) %>%
  GC_labels(
    label = "name",
    y = 25,
    fontFamily = "Arial",
    fontStyle = "normal",
    itemStyle = label_styles
  ) %>%
  GC_genes(
    show = TRUE,
    marker = "rbox",
    marker_size = "small",
    stroke = "black",
    strokeWidth = 1
  ) %>%
  GC_legend(FALSE)

# Add links and colors
cluster_plot <- main_chart %>%
  GC_links(
    "name",
    use_group_colors = TRUE,
    linkWidth = 0.3
  ) %>%
  GC_color(
    customColors = gene_colors
  )

# Display plot
cluster_plot

# Render the widget through the browser detected by Chromote.
html_file <- tempfile("synteny-", fileext = ".html")
saveWidget(
  cluster_plot,
  html_file,
  selfcontained = TRUE
)

# Start a Chromote browser session
b = ChromoteSession$new()

# Enable the Chrome DevTools Page domain
b$Page$enable()

# Open the HTML file
b$Page$navigate(paste0("file://", normalizePath(html_file)))

# Wait until all SVG layers have been rendered
svg_ready = FALSE
for (i in 1:50) {
  
  # Check whether more than one SVG has been created
  ready = b$Runtime$evaluate(
    "document.querySelectorAll('.geneviewer svg').length > 1"
  )$result$value
  
  # Stop waiting once the figure is ready
  if (isTRUE(ready)) {
    svg_ready = TRUE
    break
  }
  
  # Pause briefly before checking again
  Sys.sleep(0.2)
}

# Stop if the SVG layers never appeared
if (!svg_ready)
  stop("Not all SVG layers appeared — increase wait time")

# JavaScript to merge all SVG layers into one SVG
merge_js = "
(function() {

  // Find the GeneViewer container
  const container = document.querySelector('.geneviewer');
  const containerRect = container.getBoundingClientRect();

  // Collect all SVG layers
  const svgs = Array.from(container.querySelectorAll('svg'));

  // Separate link layers from the other layers
  const linkSvgs = svgs.filter(s => s.querySelector('.GeneLink') !== null);
  const baseSvgs = svgs.filter(s => !linkSvgs.includes(s));

  // Preserve drawing order
  const ordered = [...linkSvgs, ...baseSvgs];

  let inner = '';

  // Copy each SVG into a translated group
  ordered.forEach(svg => {

    const r = svg.getBoundingClientRect();

    const dx = r.left - containerRect.left;
    const dy = r.top - containerRect.top;

    inner +=
      '<g transform=\"translate(' + dx + ',' + dy + ')\">' +
      svg.innerHTML +
      '</g>';

  });

  // Determine the final SVG dimensions
  const width = containerRect.width;
  const height = containerRect.height;

  // Return the merged SVG
  return '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' +
         width +
         '\" height=\"' +
         height +
         '\" viewBox=\"0 0 ' +
         width +
         ' ' +
         height +
         '\">' +
         inner +
         '</svg>';

})()
"

# Execute the JavaScript and retrieve the merged SVG
svg_merged = b$Runtime$evaluate(merge_js)$result$value

# Save the merged SVG to disk
writeLines(svg_merged, "results/synteny/synteny.svg")

# Close the browser session
b$close()
unlink(html_file)
