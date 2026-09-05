library(ggplot2)
library(dplyr)
library(egg)
library(tidyr)
library(stringr)
library(scales)
library(svglite)

dir.create("results/phyletic-distribution", recursive = TRUE, showWarnings = FALSE)

data <- read.csv(
  "data/phyletic-distribution/receptor_sequence_counts.csv",
  stringsAsFactors = FALSE
)

# Axis definitions
receptor_order <- c("CNR1", "CNR2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "LPAR1", "LPAR2", "LPAR3",
                    "GPR3", "GPR6", "GPR12", "GPR18", "GPR183", "GPR55", "GPR35", "GPR119", "GPR161", "CYSLTR1", "CYSLTR2", "LPAR4", "LPAR6", "LPAR5", "GPR17", "OXGR1",
                    "LTB4R1", "LTB4R2", "OXER1", "GPR31", "FFAR1", "FFAR2", "FFAR3",
                    "GPR42", "PTGDR2", "PTAFR", "GPR34", "GPR87", "GPR82", "P2RY10",
                    "GPR174", "GPR132", "PTGDR", "PTGER2", "PTGIR", "PTGER4", "PTGER1", "PTGFR",
                    "TBXA2R", "PTGER3", "GPR84", "FFAR4")

phyla_order <- c("Mammalia", "Archosauria", "Lepidosauria", "Testudines", "Amphibia",
                 "Dipnomorpha", "Coelacanthiformes", "Actinopterygii", "Chondrichthyes",
                 "Cyclostomata")

n_receptor <- length(receptor_order)

gap_after_phylum <- "Cyclostomata"
gap_size <- 1

phyla_pos <- setNames(seq_along(phyla_order), phyla_order)
if (!is.null(gap_after_phylum)) {
  gap_at <- match(gap_after_phylum, phyla_order)
  if (!is.na(gap_at) && gap_at < length(phyla_order)) {
    phyla_pos[(gap_at + 1):length(phyla_order)] <- phyla_pos[(gap_at + 1):length(phyla_order)] + gap_size
  }
}
n_phyla <- max(phyla_pos)

# 1. Separate data and identify positive singles
single_data <- data %>%
  filter(!grepl("\\|", Receptor)) %>%
  mutate(
    base_rec = sub("/LIKE$", "", Receptor),
    is_like = grepl("/LIKE$", Receptor)
  )

single_positive <- single_data %>% filter(Sequences > 0)
merged_data <- data %>% filter(grepl("\\|", Receptor) & Sequences > 0)

# 2. Process merged rectangles (layer 2)
split_tiles <- data.frame(Phyla = character(), base_rec = character(), is_split = logical())

if (nrow(merged_data) > 0) {
  user_rects <- merged_data %>%
    mutate(
      rec_start = sapply(strsplit(Receptor, "\\|"), `[`, 1),
      rec_end   = sapply(strsplit(Receptor, "\\|"), function(x) x[length(x)])
    ) %>%
    filter(rec_start %in% receptor_order & rec_end %in% receptor_order & Phyla %in% phyla_order) %>%
    rowwise() %>%
    mutate(
      covered_receptors = list(as.character(receptor_order[match(rec_start, receptor_order):match(rec_end, receptor_order)])),
      needs_split = any(covered_receptors %in% single_positive$base_rec[single_positive$Phyla == Phyla])
    ) %>%
    ungroup() %>%
    mutate(
      is_like = FALSE,
      X_idx   = unname(phyla_pos[Phyla]),
      Y_start = match(rec_start, rev(receptor_order)),
      Y_end   = match(rec_end, rev(receptor_order)),
      ymin = pmin(Y_start, Y_end) - 0.5,
      ymax = pmax(Y_start, Y_end) + 0.5,
      xmin = ifelse(needs_split, X_idx, X_idx - 0.5),
      xmax = X_idx + 0.5
    )
  
  if (nrow(user_rects) > 0) {
    split_tiles <- user_rects %>%
      unnest(covered_receptors) %>%
      filter(needs_split) %>%
      distinct(Phyla, base_rec = covered_receptors) %>%
      mutate(is_split = TRUE)
  }
} else {
  user_rects <- data.frame()
}

# 3. Process single data (layer 3)
single_processed <- single_data %>%
  filter(Sequences > 0) %>%
  filter(base_rec %in% receptor_order & Phyla %in% phyla_order) %>%
  group_by(base_rec, Phyla) %>%
  mutate(has_both = n_distinct(is_like) == 2) %>%
  ungroup() %>%
  left_join(split_tiles, by = c("Phyla", "base_rec")) %>%
  mutate(
    is_split = replace_na(is_split, FALSE),
    X_idx = unname(phyla_pos[Phyla]),
    Y_idx = match(base_rec, rev(receptor_order)),
    ymin = Y_idx - 0.5,
    ymax = Y_idx + 0.5,
    xmin = X_idx - 0.5,
    xmax = ifelse(is_split, X_idx, X_idx + 0.5)
  )

# 4. Build empty base grid (layer 1)
base_grid <- expand.grid(base_rec = receptor_order, Phyla = phyla_order, stringsAsFactors = FALSE) %>%
  mutate(
    Sequences = 0,
    is_like = FALSE,
    X_idx = unname(phyla_pos[Phyla]),
    Y_idx = match(base_rec, rev(receptor_order)),
    xmin = X_idx - 0.5,
    xmax = X_idx + 0.5,
    ymin = Y_idx - 0.5,
    ymax = Y_idx + 0.5
  )

# 5. Combine plain rectangles (painter's algorithm)
rect_data <- bind_rows(base_grid %>% select(Phyla, Sequences, is_like, xmin, xmax, ymin, ymax))

if (nrow(user_rects) > 0) {
  rect_data <- bind_rows(rect_data, user_rects %>% select(Phyla, Sequences, is_like, xmin, xmax, ymin, ymax))
}

single_rect_only <- single_processed %>% filter(!has_both & !is_like)
if (nrow(single_rect_only) > 0) {
  rect_data <- bind_rows(rect_data, single_rect_only %>% select(Phyla, Sequences, is_like, xmin, xmax, ymin, ymax))
}

rect_data <- rect_data %>%
  mutate(
    xmid = (xmin + xmax) / 2,
    ymid = (ymin + ymax) / 2,
    label_text = ifelse(Sequences > 0, paste0(Sequences, ifelse(is_like, "*", "")), "")
  )

# 5b. Diagonal-split triangles (unchanged)
triangle_data <- single_processed %>% filter(has_both | is_like)

if (nrow(triangle_data) > 0) {
  triangle_data <- triangle_data %>%
    mutate(
      group_id = paste(Phyla, base_rec, is_like, sep = "_"),
      xmid = ifelse(is_like, (xmin + xmax + xmax) / 3, (xmin + xmin + xmax) / 3),
      ymid = ifelse(is_like, (ymin + ymin + ymax) / 3, (ymin + ymax + ymax) / 3),
      label_text = paste0(Sequences, ifelse(is_like, "*", ""))
    ) %>%
    rowwise() %>%
    mutate(
      poly = list(
        if (is_like) {
          data.frame(x = c(xmin, xmax, xmax), y = c(ymin, ymin, ymax))
        } else {
          data.frame(x = c(xmin, xmin, xmax), y = c(ymin, ymax, ymax))
        }
      )
    ) %>%
    ungroup() %>%
    unnest(poly)
}

triangle_labels <- if (nrow(triangle_data) > 0) {
  triangle_data %>% distinct(group_id, xmid, ymid, label_text, Sequences)
} else {
  data.frame()
}

# 6. Color scale: monochromatic amber ramp, complements panel 1's blue
data_scope <- data %>%
  filter(Phyla %in% phyla_order) %>%
  mutate(base_rec = sub("/LIKE$", "", sapply(strsplit(Receptor, "\\|"), `[`, 1))) %>%
  filter(base_rec %in% receptor_order | grepl("\\|", Receptor))
max_val <- max(data_scope$Sequences, na.rm = TRUE)

amber_stops  <- c(0, 1, 3, 10, 30, 80, 200, 450, 600, max_val)
amber_colors <- c(
  "#ffffff",  # 0
  "#fdf3e3",  # 1
  "#fbe0b8",  # 3
  "#f6c78a",  # 10
  "#f0aa5c",  # 30
  "#d08a3e",  # 80   (was #e08a2e)
  "#aa6828",  # 200  (was #b8681a)
  "#844d1d",  # 450  (was #8f4d12)
  "#623814",  # 600  (was #6b380c)
  "#44260c"   # max  (was #4a2606)
)

# A data maximum can coincide with a predefined stop. Retain the final colour
# for duplicated stops so interpolation and the legend remain well defined.
keep_stops <- !duplicated(amber_stops, fromLast = TRUE)
amber_stops <- amber_stops[keep_stops]
amber_colors <- amber_colors[keep_stops]

# Custom transform: maps irregular amber_stops onto evenly-spaced [0,1]
# positions, so the colorbar shows equal visual gaps between labels
# regardless of how close/far the underlying values are.
even_positions <- seq(0, 1, length.out = length(amber_stops))

stop_trans <- scales::trans_new(
  name      = "even_stops",
  transform = function(x) stats::approx(x = amber_stops, y = even_positions, xout = x, rule = 2)$y,
  inverse   = function(x) stats::approx(x = even_positions, y = amber_stops, xout = x, rule = 2)$y
)

# Palette is now evenly spaced by default (no 'values' arg needed --
# the trans above already handles the uneven stop spacing)
fill_pal <- scales::gradient_n_pal(amber_colors)
get_fill <- function(x) fill_pal(stop_trans$transform(pmin(x, max_val)))

get_text_color <- function(x) {
  cols <- get_fill(x)
  rgb_mat <- grDevices::col2rgb(cols) / 255
  lum <- 0.2126 * rgb_mat["red", ] + 0.7152 * rgb_mat["green", ] + 0.0722 * rgb_mat["blue", ]
  ifelse(lum > 0.5, "grey15", "white")
}

rect_data <- rect_data %>% mutate(text_color = get_text_color(Sequences))
if (nrow(triangle_labels) > 0) {
  triangle_labels <- triangle_labels %>% mutate(text_color = get_text_color(Sequences))
}

# 7. Build heatmap
main_plot <- ggplot() +
  geom_rect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Sequences),
    color = "grey80", linewidth = 0.2
  )

if (nrow(triangle_data) > 0) {
  main_plot <- main_plot +
    geom_polygon(
      data = triangle_data,
      aes(x = x, y = y, group = group_id, fill = Sequences),
      color = "grey80", linewidth = 0.2
    )
}

main_plot <- main_plot +
  geom_text(
    data = rect_data,
    aes(x = xmid, y = ymid, label = label_text, color = text_color),
    size = 4, family = "sans"
  )

if (nrow(triangle_labels) > 0) {
  main_plot <- main_plot +
    geom_text(
      data = triangle_labels,
      aes(x = xmid, y = ymid, label = label_text, color = text_color),
      size = 3.2, family = "sans"
    )
}

main_plot <- main_plot +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = amber_colors,
    trans   = stop_trans,
    limits  = c(0, max_val),
    oob     = scales::squish,
    na.value = "white",
    name    = "Sequences",
    breaks  = amber_stops,
    labels  = as.character(round(amber_stops)),
    guide   = guide_colorbar(
      direction   = "horizontal",
      barwidth    = unit(6, "cm"),
      barheight   = unit(0.35, "cm"),
      frame.colour = "grey40",
      frame.linewidth = 0.3,
      ticks.colour = "grey40",
      title.position = "left",
      title.vjust = 0.9,
      label.position = "bottom"
    )
  ) +
  scale_x_continuous(
    breaks = unname(phyla_pos),
    labels = names(phyla_pos),
    expand = c(0, 0),
    limits = c(0.5, n_phyla + 0.5),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = 1:n_receptor,
    labels = rev(receptor_order),
    expand = c(0, 0),
    limits = c(0.5, n_receptor + 0.5)
  ) +
  coord_fixed(ratio = 0.8) +
  theme_classic(base_family = "sans", base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 12, colour = "grey20"),
    axis.text.y = element_text(size = 12, colour = "grey20", hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold", colour = "grey20"),
    legend.text = element_text(size = 8, colour = "grey20"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 10, 5, 5)
  )

# 8. Save
ggsave(
  "results/phyletic-distribution/sequence_counts_vertebrates.svg",
  plot = main_plot,
  device = svglite,
  width = 11,
  height = 20,
  units = "in"
)
