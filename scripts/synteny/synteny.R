# Load required packages
library(dplyr)
library(readr)
library(geneviewer)
library(htmlwidgets)
library(chromote)

# Locate synteny dataset
csv_candidates <- c(
  "data/synteny/synteny.csv",
  "synteny.csv",
  "../../data/synteny/synteny.csv"
)
csv_file <- csv_candidates[file.exists(csv_candidates)][1]
if (is.na(csv_file)) {
  stop("Could not find synteny.csv in candidate locations.")
}

synteny_data <- read_csv(csv_file, show_col_types = FALSE) %>%
  select(-starts_with("...")) %>%
  filter(receptor == "cnr1") %>%
  rename(
    cluster = taxa,
    name = gene,
    strand = orientation
  ) %>%
  mutate(
    strand = trimws(as.character(strand)),
    genomic_start = start,
    genomic_end = end,
    start_coordinate = if_else(strand == "+", start, end)
  )

stopifnot(all(synteny_data$strand %in% c("+", "-")))

goi <- "CNR1"
goi_color <- "#e21f4c"
ortholog_color <- "#f28da5"
background_color <- "#e6e7e8"

ortholog <- synteny_data %>%
  distinct(name, cluster) %>%
  count(name, name = "paralogons") %>%
  filter(paralogons > 1) %>%
  pull(name)

all_genes <- unique(synteny_data$name)

gene_colors <- setNames(
  rep(background_color, length(all_genes)),
  all_genes
)

gene_colors[ortholog] <- ortholog_color
gene_colors[goi] <- goi_color

dynamic_height <- paste0(
  80 * length(unique(synteny_data$cluster)),
  "px"
)

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
    formatter = "<b>Start:</b> {genomic_start}<br><b>End:</b> {genomic_end}"
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
    highlight = goi
  ) %>%
  GC_genes(
    show = TRUE,
    marker = "rbox",
    marker_size = "small",
    stroke = "black",
    strokeWidth = 1
  ) %>%
  GC_legend(FALSE)

cluster_plot <- main_chart %>%
  GC_links(
    "name",
    use_group_colors = TRUE,
    linkWidth = 0.3
  ) %>%
  GC_color(
    customColors = gene_colors
  )

coord_values <- unname(lapply(cluster_plot$x$series, function(s) {
  d <- as.data.frame(s$data)
  d$start_coordinate[order((d$start + d$end) / 2)]
}))

relabel_js <- r"---(
function(el, x, lookup) {

  const fmt = v => {
    const a = Math.abs(v);
    const [s, p] = a >= 1e9 ? [1e9, 'G'] : a >= 1e6 ? [1e6, 'M'] : a >= 1e3 ? [1e3, 'k'] : [1, ''];
    return (v / s).toPrecision(4) + p;
  };

  const isTick = s => /^\d+(\.\d+)?[kMG]$/.test(s);

  const rows = arr => {
    const sorted = arr.slice().sort((a, b) => a.y - b.y);
    const groups = [];
    sorted.forEach(o => {
      const g = groups[groups.length - 1];
      if (g && o.y - g[g.length - 1].y < 30) g.push(o);
      else groups.push([o]);
    });
    return groups.map(g => g.sort((a, b) => a.cx - b.cx));
  };

  const apply = () => {
    try {
      const texts = Array.from(el.querySelectorAll('text')).map(t => {
        const r = t.getBoundingClientRect();
        return { t: t, cx: r.left + r.width / 2, y: r.top, s: t.textContent.trim() };
      });

      const tickRows = rows(texts.filter(o => isTick(o.s)));
      const vals = lookup.map(v => [].concat(v));

      el.setAttribute('data-coords-status', JSON.stringify({
        clusters: vals.length,
        tickRows: tickRows.map(r => r.length),
        genes: vals.map(v => v.length)
      }));

      if (tickRows.length !== vals.length) return;
      for (let k = 0; k < vals.length; k++) {
        if (tickRows[k].length !== vals[k].length) return;
      }

      for (let k = 0; k < vals.length; k++) {
        tickRows[k].forEach((o, i) => {
          const txt = fmt(vals[k][i]);
          if (o.t.textContent !== txt) o.t.textContent = txt;
        });
      }
      el.setAttribute('data-coords-fixed', '1');
    } catch (e) {
      el.setAttribute('data-coords-status', 'ERR ' + e.message);
    }
  };

  const obs = new MutationObserver(apply);
  obs.observe(el, { childList: true, subtree: true, characterData: true });
  apply();
}
)---"

cluster_plot <- htmlwidgets::onRender(cluster_plot, relabel_js, data = coord_values)

cluster_plot

Sys.setenv(
  CHROMOTE_PORT = "9222",
  CHROMOTE_CHROME = "/usr/bin/google-chrome"
)

out_dir <- if (dir.exists("results/synteny") || dir.exists("scripts/synteny")) {
  "results/synteny"
} else {
  "."
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

temp_html <- file.path(out_dir, "temp.html")
saveWidget(
  cluster_plot,
  temp_html,
  selfcontained = TRUE
)

b = ChromoteSession$new()
b$Page$enable()
b$Page$navigate(paste0("file://", normalizePath(temp_html)))

svg_ready = FALSE
for (i in 1:50) {
  ready = b$Runtime$evaluate(
    "document.querySelectorAll('.geneviewer svg').length > 1 &&
     document.querySelector('[data-coords-fixed]') !== null"
  )$result$value
  
  if (isTRUE(ready)) {
    svg_ready = TRUE
    break
  }
  
  Sys.sleep(0.2)
}

if (!svg_ready) {
  status <- b$Runtime$evaluate(
    "Array.from(document.querySelectorAll('[data-coords-status]')).map(e => e.getAttribute('data-coords-status')).join(' | ') || 'NO STATUS: onRender did not run'"
  )$result$value
  cat("Page status:", status, "\n")
  b$close()
  unlink(temp_html)
  stop("Coordinate relabelling did not complete")
}

merge_js = "
(function() {

  const container = document.querySelector('.geneviewer');
  const containerRect = container.getBoundingClientRect();

  const svgs = Array.from(container.querySelectorAll('svg'));

  const linkSvgs = svgs.filter(s => s.querySelector('.GeneLink') !== null);
  const baseSvgs = svgs.filter(s => !linkSvgs.includes(s));

  const ordered = [...linkSvgs, ...baseSvgs];

  let inner = '';

  ordered.forEach(svg => {

    const r = svg.getBoundingClientRect();

    const dx = r.left - containerRect.left;
    const dy = r.top - containerRect.top;

    inner +=
      '<g transform=\"translate(' + dx + ',' + dy + ')\">' +
      svg.innerHTML +
      '</g>';

  });

  const width = containerRect.width;
  const height = containerRect.height;

  return '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' +
         width +
         '\" height=\"' +
         height +
         '\" viewBox=\"0 0 ' +
         width +
         ' ' +
         height +
         '\">' +
         '<rect width=\"100%\" height=\"100%\" fill=\"white\"/>' +
         inner +
         '</svg>';

})()
"

svg_merged = b$Runtime$evaluate(merge_js)$result$value

writeLines(svg_merged, file.path(out_dir, "synteny.svg"))

b$close()
unlink(temp_html)
unlink(file.path(out_dir, "temp_files"), recursive = TRUE)

