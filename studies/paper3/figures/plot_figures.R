suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(jsonlite)
  library(digest)
  library(grid)
})
options(warn = 2)
# Grobs are built on a null PDF device (same metrics as the default device), so no Rplots.pdf is written to the working directory.
grDevices::pdf(NULL)

# Paper 3 publication figures from one figure-source directory (contract publication_figure_sources_v1, written by the
# figure_sources stage). Main figures: figure1 (framework and evidence states), figure2 (detection and breadth),
# figure3 (positional context), figure4 (a priori testability and record cases); supplementary figures: figureS1
# (query properties and alignment quality), figureS2 (Monte Carlo precision), figureS3 (positional robustness).
# The record-case figure (figureS4) is drawn by plot_record_cases.R. No panel contains a table: lookup-style content is
# exported as CSV (tables/ inside the output directory) and every panel is drawn from the exported CSV.
# House style (paper3_figure_style.R): one width per profile (170 mm by default), at most 200 mm tall, Arial 7 pt, panel tags
# A, B, C without panel titles. Method text, filters, transformations and symbol definitions belong in the caption: they are
# recorded per figure as legend_notes in the artwork manifest. Hue marks platforms only. Evidence states use grey value,
# fill and outline style. Proportions carry 95% Clopper-Pearson intervals; control means carry the relative-entropy Hoeffding
# interval of the Methods. UK spelling, units in parentheses, thousands separators, en dashes in ranges, true minus signs.
# R object names stay figure1-figure4; the manuscript numbers them Figures 2-5 (the graphical abstract is Figure 1).
# Exported CSV names carry the manuscript numbers (e.g. fig2B_* for figure1 panel B, figS4B_* for figureS2, tableS6_* for Table S6).
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% c(3, 4))
input <- normalizePath(arguments[1])
style_file <- normalizePath(arguments[2])
output <- arguments[3]
selection <- if (length(arguments) == 4) arguments[4] else "all"
manifest <- fromJSON(file.path(input, "manifest.json"))
stopifnot(manifest$contract == "publication_figure_sources_v1")
for (filename in names(manifest$files)) {
  stopifnot(!grepl("(^/|(^|/)\\.\\.(/|$))", filename),
            digest(file.path(input, filename), algo = "sha256", file = TRUE) == manifest$files[[filename]])
}
source(style_file)
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
tables_directory <- file.path(output, "tables")
if (!dir.exists(tables_directory)) dir.create(tables_directory, recursive = TRUE)
figure_width_mm <- FIG$width
figure_dpi <- 300
groups <- c("ONT", "Illumina", "BGI/DNBSEQ", "PacBio")
colours <- pal_platform[groups]
draws_per_parent <- manifest$counts$draws_per_parent
draw_prefixes <- as.numeric(unlist(manifest$draw_prefixes))
final_prefix <- max(draw_prefixes)
target_limit <- manifest$target_limit
parent_label <- function(id) gsub("_", " ", sub("_forward$", "", id))
base_pt <- FIG$base
dense_pt <- FIG$dense
load_table <- function(name) read_tsv(file.path(input, paste0(name, ".tsv")),
                                     show_col_types = FALSE, na = c("", "not_estimated"))
group_order <- function(data) mutate(data, group = factor(group, levels = groups))
serial <- function(values) {
  values <- as.character(values)
  if (length(values) <= 2) paste(values, collapse = " and ") else paste0(paste(values[-length(values)], collapse = ", "), ", and ", values[length(values)])
}
comma <- function(x) scales::comma(x, accuracy = 1)
# Text column outside a plot area: no clipping and room in the right margin.
no_clip <- coord_cartesian(clip = "off")

# Write a table to tables/ (never overwriting different content) and return the values read back,
# so that the panels are built from exactly what is exported.
export_csv <- function(data, name) {
  path <- file.path(tables_directory, name)
  temporary <- tempfile(fileext = ".csv")
  write_csv(data, temporary, na = "NA")
  if (file.exists(path)) {
    stopifnot(digest(temporary, algo = "sha256", file = TRUE) == digest(path, algo = "sha256", file = TRUE))
  } else {
    stopifnot(file.copy(temporary, path))
  }
  unlink(temporary)
  read_csv(path, show_col_types = FALSE, na = "NA")
}

# Schematic panels drawn in millimetre coordinates (x 0 to the figure width, y downwards from the top of the panel).
canvas_plot <- function(height, layers) {
  ggplot() + layers +
    scale_x_continuous(limits = c(0, figure_width_mm), expand = c(0, 0)) +
    scale_y_reverse(limits = c(height, 0), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_void(base_family = FONT_SANS) +
    theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_rect(fill = "white", colour = NA),
          plot.tag = element_text(size = FIG$tag, face = "bold", hjust = 0, vjust = 1, family = FONT_SANS),
          plot.tag.position = c(0, 1), plot.tag.location = "plot")
}
mm_text <- function(x, y, label, size = base_pt, face = "plain", hjust = 0, colour = ink, family = FONT_SANS, vjust = 0.5, parse = FALSE) {
  annotate("text", x = x, y = y, label = label, size = fs(size), fontface = face, hjust = hjust,
           vjust = vjust, colour = colour, family = family, lineheight = 0.95, parse = parse)
}
# Species names (Latin binomials, flag from the record-case verification) are set in italics with plotmath. A leading comma of the upright remainder is
# set with the italic name, so the plotmath italic correction falls inside the name.
mm_text_species <- function(x, y, organism, italic, rest, size = base_pt) {
  if (!isTRUE(italic)) return(mm_text(x, y, paste0(organism, rest), size))
  quote <- function(text) paste0('"', gsub('"', '\\\\"', text), '"')
  comma <- startsWith(rest, ",")
  label <- paste0("italic(", quote(paste0(organism, if (comma) "," else "")), ")*", quote(if (comma) substring(rest, 2) else rest))
  annotate("text", x = x, y = y, label = label, parse = TRUE, size = fs(size), hjust = 0, vjust = 0.5, colour = ink, family = FONT_SANS)
}
mm_rect <- function(xmin, xmax, ymin, ymax, fill = "white", colour = NA, linewidth = lw(0.5), linetype = "solid", alpha = 1) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = colour,
           linewidth = linewidth, linetype = linetype, alpha = alpha)
}
mm_segment <- function(x, y, xend, yend, colour = grey_dark, linewidth = lw(0.5), head = FALSE, linetype = "solid", lineend = "butt") {
  arguments <- list("segment", x = x, y = y, xend = xend, yend = yend, colour = colour, linewidth = linewidth,
                    linetype = linetype, lineend = lineend)
  if (head) arguments$arrow <- arrow(length = unit(1.2, "mm"), type = "closed", angle = 22)
  do.call(annotate, arguments)
}

saved <- list()
# Manuscript numbering: the graphical abstract is Figure 1, so the R figures figure1-figure4 are manuscript Figures 2-5.
# Supplementary figures are numbered in the order of first citation: R figureS2 is Figure S4 and the record cases
# (R figureS4, plot_record_cases.R) are Figure S2. File names keep the R names.
manuscript_numbers <- c(figure1 = "Figure 2", figure2 = "Figure 3", figure3 = "Figure 4", figure4 = "Figure 5",
                        figureS1 = "Figure S1", figureS2 = "Figure S4", figureS3 = "Figure S3")
save_figure <- function(name, plot, height, panels, source_tables, note, csv_tables = character(), legend_notes = character()) {
  if (selection != "all" && selection != name) return(invisible(NULL))
  if (height > FIG$max_h + 1e-6) stop(name, ": ", height, " mm exceeds the ", FIG$max_h, "-mm limit of profile ", FIG_PROFILE)
  stopifnot(name %in% names(manuscript_numbers))
  extensions <- c("pdf", "png", if (FIG$tiff) "tiff")
  for (extension in extensions) {
    filename <- file.path(output, paste0("paper3_", name, ".", extension))
    stopifnot(!file.exists(filename))
    device <- switch(extension, pdf = grDevices::cairo_pdf, png = ragg::agg_png, tiff = ragg::agg_tiff)
    arguments <- list(filename, plot, width = figure_width_mm, height = height, units = "mm", bg = "white",
                      dpi = if (extension == "tiff") 600 else figure_dpi, device = device)
    if (extension == "tiff") arguments$compression <- "lzw"
    do.call(ggplot2::ggsave, arguments)
  }
  saved[[name]] <<- list(manuscript_figure = unname(manuscript_numbers[name]), panels = panels, source_tables = source_tables,
                         exported_csv = csv_tables, note = note, legend_notes = legend_notes, profile = FIG_PROFILE,
                         width_mm = figure_width_mm, height_mm = height, png_dpi = figure_dpi)
  message("Rendered ", name)
}

# ---------------------------------------------------------------------------------------------
# Figure 1 (manuscript Figure 2). Framework and evidence states
# ---------------------------------------------------------------------------------------------
counts <- manifest$counts
testability <- load_table("length_testability")
production_evalue <- unique(testability$production_evalue)
stopifnot(length(production_evalue) == 1)
informative_limit <- unique(testability$informativeness_limit)
testable_limit <- unique(testability$testable_chance_limit)
stopifnot(length(informative_limit) == 1, length(testable_limit) == 1)
# Step numbers: a black disc with a white numeral (Arial has no circled digits).
step_marker <- function(x, y, number) list(
  annotate("point", x = x, y = y, shape = 21, size = 3.4, fill = ink, colour = ink, stroke = 0),
  mm_text(x, y, as.character(number), dense_pt, "bold", hjust = 0.5, colour = "white"))

# A: four steps in one row; exclusions as dashed boxes with black text.
step_width <- (figure_width_mm - 3 * 6) / 4
step_left <- (seq_len(4) - 1) * (step_width + 6)
step_top <- 4.5
step_bottom <- 14.5
steps <- data.frame(
  title = c("Query panel", "Matched search", "Single-HSP support", "Record context"),
  detail = c(paste0(comma(counts$physical_references), " queries, ", comma(counts$canonical_parents), " parents\n",
                    draws_per_parent, " draws per parent"),
             paste0("E = ", production_evalue, ", target limit ", comma(target_limit), "\ntestability by query length"),
             "strict or tolerant\none HSP per match", "exact version and length\nexplicit eligibility"))
exclusions <- data.frame(step = c(2, 4), label = c("Untestable, chance-limited or\nincomplete: excluded from counts",
                                                  "Excluded or unresolved records:\noutside the linear-genomic set"))
workflow_height <- 28
workflow_panel <- canvas_plot(workflow_height, c(
  list(mm_rect(step_left, step_left + step_width, step_top, step_bottom, colour = grey_dark),
       mm_segment(step_left[-4] + step_width + 0.8, (step_top + step_bottom) / 2, step_left[-1] - 0.8, (step_top + step_bottom) / 2,
                  colour = ink, head = TRUE),
       mm_text(step_left + 6, step_top + 2.6, steps$title, base_pt, "bold"),
       mm_text(step_left + 2.2, step_top + 6.6, steps$detail, dense_pt),
       mm_segment(step_left[exclusions$step] + step_width / 2, step_bottom, step_left[exclusions$step] + step_width / 2, 17.2,
                  colour = ink, head = TRUE),
       mm_rect(step_left[exclusions$step], step_left[exclusions$step] + step_width, 17.4, 24.6, colour = grey_dark, linetype = "dashed"),
       mm_text(step_left[exclusions$step] + step_width / 2, 21, exclusions$label, dense_pt, hjust = 0.5)),
  unlist(lapply(1:4, function(index) step_marker(step_left[index] + 3, step_top + 2.6, index)), recursive = FALSE)))
invisible(export_csv(data.frame(step = 1:4, title = steps$title, detail = gsub("\n", " ", steps$detail)), "fig2A_workflow_steps.csv"))

# B: matching rules as a drawn schematic (one 24-nt query bar above six subject alignments), achromatic.
examples <- load_table("matching_examples")
stopifnot(nrow(examples) == 6, length(unique(examples$reference)) == 1, all(nchar(examples$reference) == 24))
matching_rules <- export_csv(data.frame(
  example = examples$example, query = examples$reference, aligned_subject = examples$subject,
  strict = examples$strict, tolerant = examples$tolerant), "fig2B_matching_rules.csv")
query_sequence <- matching_rules$query[1]
query_length <- nchar(query_sequence)
parse_alignment <- function(reference, subject) {
  if (grepl(" / ", subject, fixed = TRUE)) {
    parts <- strsplit(subject, " / ", fixed = TRUE)[[1]]
    ends <- cumsum(nchar(parts))
    stopifnot(length(parts) == 2, ends[2] == nchar(reference))
    return(data.frame(kind = "aligned", start = c(0, ends[1]), end = ends, lane = 1:2))
  }
  subject_chars <- strsplit(subject, "")[[1]]
  reference_chars <- strsplit(reference, "")[[1]]
  stopifnot(length(subject_chars) == length(reference_chars))
  kind <- ifelse(subject_chars == ".", "unaligned", ifelse(subject_chars == "-", "gap", "aligned"))
  runs <- rle(kind)
  bounds <- cumsum(runs$lengths)
  parts <- data.frame(kind = runs$values, start = bounds - runs$lengths, end = bounds, lane = 0)
  mismatch <- which(kind == "aligned" & subject_chars != reference_chars)
  if (length(mismatch)) parts <- rbind(parts, data.frame(kind = "mismatch", start = mismatch - 1, end = mismatch, lane = 0))
  parts
}
rule_x0 <- 50
rule_scale <- 3
rule_rows_y <- 23 + 6.2 * (seq_len(nrow(matching_rules)) - 1)
alignment_parts <- bind_rows(lapply(seq_len(nrow(matching_rules)), function(row) {
  mutate(parse_alignment(matching_rules$query[row], matching_rules$aligned_subject[row]), row = row)
})) %>% mutate(
  ycentre = rule_rows_y[row] + ifelse(lane == 0, 0, (lane - 1.5) * 2.7),
  half = ifelse(lane == 0, 1.5, 1.05),
  xmin = rule_x0 + start * rule_scale, xmax = rule_x0 + end * rule_scale)
bars <- filter(alignment_parts, kind %in% c("aligned", "unaligned"))
bars <- bars[order(bars$kind != "unaligned"), ]
mismatches <- filter(alignment_parts, kind == "mismatch")
gaps <- filter(alignment_parts, kind == "gap")
symbol_columns <- c(strict = 134, tolerant = 156)
draw_outcome <- function(x, y, accepted) {
  if (accepted) {
    mm_segment(c(x - 1.5, x - 0.4), c(y + 0.1, y + 1.3), c(x - 0.4, x + 1.7), c(y + 1.3, y - 1.5), colour = ink, linewidth = lw(1), lineend = "round")
  } else {
    mm_segment(c(x - 1.2, x - 1.2), c(y - 1.2, y + 1.2), c(x + 1.2, x + 1.2), c(y + 1.2, y - 1.2), colour = grey_mid, linewidth = lw(1), lineend = "round")
  }
}
rules_height <- 69
legend_y <- 63
rules_panel <- canvas_plot(rules_height, c(
  list(mm_text(0, 9, "Example", base_pt, "bold"), mm_text(rule_x0, 9, "Alignment to the query", base_pt, "bold"),
    mm_text(symbol_columns[["strict"]], 6.2, "Strict", base_pt, "bold", hjust = 0.5),
    mm_text(symbol_columns[["strict"]], 10.4, "100% identity\nfull length", dense_pt, hjust = 0.5),
    mm_text(symbol_columns[["tolerant"]], 6.2, "Tolerant", base_pt, "bold", hjust = 0.5),
    mm_text(symbol_columns[["tolerant"]], 10.4, "\u226595% identity\n\u226590% span, ungapped", dense_pt, hjust = 0.5),
    mm_text(0, 16.5, paste0("Query (", query_length, " nt)"), base_pt),
    mm_rect(rule_x0, rule_x0 + query_length * rule_scale, 16.5 - 1.5, 16.5 + 1.5, fill = ink),
    mm_text(0, rule_rows_y, matching_rules$example, base_pt),
    mm_rect(bars$xmin, bars$xmax, bars$ycentre - bars$half, bars$ycentre + bars$half,
      fill = ifelse(bars$kind == "aligned", grey_fill, "white"), colour = grey_dark,
      linetype = ifelse(bars$kind == "aligned", "solid", "dashed")),
    mm_rect(mismatches$xmin, mismatches$xmax, mismatches$ycentre - mismatches$half, mismatches$ycentre + mismatches$half, fill = ink),
    mm_segment(rule_x0 + gaps$start * rule_scale, gaps$ycentre, rule_x0 + gaps$end * rule_scale, gaps$ycentre, colour = grey_dark, linewidth = lw(0.75))),
  lapply(seq_len(nrow(matching_rules)), function(row) list(
    draw_outcome(symbol_columns[["strict"]], rule_rows_y[row], matching_rules$strict[row]),
    draw_outcome(symbol_columns[["tolerant"]], rule_rows_y[row], matching_rules$tolerant[row]))),
  list(
    mm_rect(0, 3.5, legend_y - 1.3, legend_y + 1.3, fill = grey_fill, colour = grey_dark), mm_text(5, legend_y, "Aligned"),
    mm_rect(22, 25.5, legend_y - 1.3, legend_y + 1.3, colour = grey_dark, linetype = "dashed"), mm_text(27, legend_y, "Unaligned query bases"),
    mm_rect(60, 62.6, legend_y - 1.3, legend_y + 1.3, fill = ink), mm_text(64.5, legend_y, "Mismatch"),
    mm_segment(80, legend_y, 84.5, legend_y, colour = grey_dark, linewidth = lw(0.75)), mm_text(86.5, legend_y, "Gap within one HSP"),
    draw_outcome(118, legend_y, TRUE), mm_text(121, legend_y, "Accepted"),
    draw_outcome(138, legend_y, FALSE), mm_text(141, legend_y, "Rejected"))))

# C: evidence states as a CONSORT-style decision tree with the canonical parents (strict / tolerant) in every end state.
# Every parent has exactly one state; states that the tree does not draw must be empty.
state_counts <- load_table("evidence_states")
drawn_states <- c("untestable", "chance_limited", "failed_or_incomplete", "raw_zero", "endpoint_non_positive", "cap_limited_non_positive",
                  "positive_unresolved_breadth", "identified_positive")
stopifnot(all(state_counts$parents[!state_counts$state %in% drawn_states] == 0),
          all(tapply(state_counts$parents, state_counts$rule, sum) == counts$canonical_parents))
state_n <- function(state, rule) state_counts$parents[state_counts$state == state & state_counts$rule == rule]
testable_states <- c("raw_zero", "endpoint_non_positive", "cap_limited_non_positive", "positive_unresolved_breadth", "identified_positive")
testable_n <- c(strict = sum(state_counts$parents[state_counts$rule == "strict" & state_counts$state %in% testable_states]),
                tolerant = sum(state_counts$parents[state_counts$rule == "tolerant" & state_counts$state %in% testable_states]))
stopifnot(testable_n[["strict"]] + state_n("failed_or_incomplete", "strict") == counts$strict_testable_parents,
          testable_n[["tolerant"]] + state_n("failed_or_incomplete", "tolerant") == counts$tolerant_testable_parents)
pair_n <- function(state) paste0(comma(state_n(state, "strict")), " / ", comma(state_n(state, "tolerant")))
evidence_states <- export_csv(data.frame(
  evidence_state = c("Untestable", "Chance-limited", "Failed or incomplete", "Raw zero", "Endpoint non-positive", "Cap-limited non-positive",
                     "Positive, unresolved breadth", "Identified positive"),
  state = drawn_states,
  counted_as = c("excluded", "excluded", "excluded", "0", "0", "0", "1", "1"),
  strict_parents = vapply(drawn_states, state_n, numeric(1), rule = "strict"),
  tolerant_parents = vapply(drawn_states, state_n, numeric(1), rule = "tolerant"),
  search_complete = c("n/a", "n/a", "No", "Yes", "Yes", "Yes", "Yes", "Yes"),
  raw_hsps = c("n/a", "n/a", "Unknown", "0", ">0", ">0", ">0", ">0"),
  accepted_support = c("n/a", "n/a", "Unknown", "No", "No", "No", "Yes, additional breadth unresolved", "Yes"),
  target_cap_reached = c("n/a", "n/a", "Unknown", "No", "No", "Yes", "Either", "Either")), "fig2C_evidence_states.csv")
state_label <- function(name, state) paste0(name, " (n = ", pair_n(state), ")")
question_x <- 33
question_w <- 44
cap_x <- 77.5
cap_w <- 28
state_x <- 132
state_w <- 68
tree_nodes <- data.frame(
  id = c("S0", "T1", "T2", "Q1", "Q2", "Q3", "Q4", "Q5", "U1", "U2", "L1", "L2", "L3", "L4", "L5", "L6"),
  label = c(paste0("Canonical parents (n = ", comma(counts$canonical_parents), ")"),
            paste0("Score retainable at E = ", production_evalue, "\nand chance bound \u2264 ", informative_limit, "?"),
            paste0("Chance bound \u2264 ", testable_limit, "?"),
            "Search complete?", "Raw HSPs > 0?", "Accepted support?", "Reached target cap?", "Additional breadth resolved?",
            state_label("Untestable, excluded", "untestable"), state_label("Chance-limited, excluded", "chance_limited"),
            state_label("Failed or incomplete, excluded", "failed_or_incomplete"), state_label("Raw zero", "raw_zero"),
            state_label("Endpoint non-positive", "endpoint_non_positive"), state_label("Cap-limited non-positive", "cap_limited_non_positive"),
            state_label("Positive, unresolved breadth", "positive_unresolved_breadth"), state_label("Identified positive", "identified_positive")),
  kind = c("start", rep("question", 7), rep("excluded", 3), rep("zero", 3), rep("positive", 2)),
  x = c(question_x, question_x, question_x, question_x, question_x, question_x, cap_x, question_x,
        state_x, state_x, state_x, state_x, state_x, state_x, state_x, question_x),
  y = c(7, 17.5, 28, 37, 46, 57, 57, 69, 17.5, 28, 37, 46, 52.5, 61.5, 69, 79.5),
  w = c(question_w, question_w, question_w, question_w, question_w, question_w, cap_w, question_w, rep(state_w, 7), question_w),
  h = c(6, 8, 6, 6, 6, 6, 6, 6, rep(6, 7), 6))
stopifnot(!anyNA(tree_nodes$label))
tree_edges <- export_csv(data.frame(
  from = c("S0", "T1", "T1", "T2", "T2", "Q1", "Q1", "Q2", "Q2", "Q3", "Q3", "Q4", "Q4", "Q5", "Q5"),
  to = c("T1", "U1", "T2", "U2", "Q1", "L1", "Q2", "L2", "Q3", "Q4", "Q5", "L3", "L4", "L6", "L5"),
  answer = c("", "no", "yes", "no", "yes", "no", "yes", "no", "yes", "no", "yes", "no", "yes", "yes", "no")), "fig2C_decision_tree_edges.csv")
tree_nodes$fill <- ifelse(tree_nodes$kind == "positive", grey_fill, "white")
tree_nodes$border <- ifelse(tree_nodes$kind == "question", grey_mid, grey_dark)
tree_nodes$linetype <- ifelse(tree_nodes$kind == "excluded", "dashed", "solid")
node_at <- function(id) tree_nodes[match(id, tree_nodes$id), ]
edge_label <- function(x, y, label, hjust = 0) mm_text(x, y, label, dense_pt, hjust = hjust, colour = grey_text)
tree_edge_layers <- lapply(seq_len(nrow(tree_edges)), function(index) {
  edge <- tree_edges[index, ]
  from <- node_at(edge$from)
  to <- node_at(edge$to)
  if (edge$from == "Q4") {
    x_branch <- from$x + from$w / 2 + 3
    return(list(
      mm_segment(from$x + from$w / 2, from$y, x_branch, from$y),
      mm_segment(x_branch, from$y, x_branch, to$y),
      mm_segment(x_branch, to$y, to$x - to$w / 2, to$y, head = TRUE),
      edge_label(x_branch + 1, to$y + ifelse(edge$answer == "no", -1.8, 1.8), edge$answer)))
  }
  if (from$x == to$x) {
    label <- if (edge$from == "T2") paste0("yes, testable (n = ", comma(testable_n[["strict"]]), " / ", comma(testable_n[["tolerant"]]), ")") else edge$answer
    list(mm_segment(from$x, from$y + from$h / 2, to$x, to$y - to$h / 2, head = TRUE),
         edge_label(from$x + 1.2, (from$y + from$h / 2 + to$y - to$h / 2) / 2, label))
  } else {
    x_start <- from$x + from$w / 2
    x_end <- to$x - to$w / 2
    list(mm_segment(x_start, from$y, x_end, to$y, head = TRUE),
         edge_label((x_start + x_end) / 2, from$y - 1.6, edge$answer, hjust = 0.5))
  }
})
# Step brackets on the left edge link the tree to the steps of panel A.
bracket <- function(top, bottom, number) list(
  mm_segment(c(4.5, 4.5, 4.5), c(top, top, bottom), c(6.5, 4.5, 6.5), c(top, bottom, bottom), colour = grey_dark),
  step_marker(2, (top + bottom) / 2, number))
tree_height <- 84
tree_panel <- canvas_plot(tree_height, c(
  list(mm_text(state_x - state_w / 2, 7, "End states, canonical parents (strict / tolerant)", dense_pt, colour = grey_text),
       mm_rect(tree_nodes$x - tree_nodes$w / 2, tree_nodes$x + tree_nodes$w / 2, tree_nodes$y - tree_nodes$h / 2,
               tree_nodes$y + tree_nodes$h / 2, fill = tree_nodes$fill, colour = tree_nodes$border, linetype = tree_nodes$linetype),
       mm_text(tree_nodes$x, tree_nodes$y, tree_nodes$label, base_pt, hjust = 0.5)),
  bracket(13.5, 40, 2), bracket(43, 82.5, 3), tree_edge_layers))
figure1_legend_notes <- c(
  "A, The four steps. Untestable, chance-limited and incomplete searches leave every count. Record context defines the linear-genomic set.",
  "B, Conceptual 24-nt examples. HSP: high-scoring segment pair. Two separate bars are two HSPs, each assessed alone. Strict support is nested within tolerant support.",
  paste0("C, Evidence states in the order of assessment, with the canonical parents (strict / tolerant) per end state. Dashed outlines mark excluded states, ",
         "white boxes counted non-positive states and grey boxes counted positive states. E: expected chance count of the search (E-value). ",
         "Step numbers refer to panel A."))
figure1_heights <- c(workflow_height, rules_height, tree_height)
figure1 <- assemble(workflow_panel, rules_panel, tree_panel, heights = figure1_heights) & theme(plot.margin = margin(0, 0, 0, 0))
save_figure("figure1", figure1, sum(figure1_heights), c("Workflow steps", "Matching rules schematic", "Evidence-state flow with parent counts"),
            c("matching_examples.tsv", "manifest.json", "length_testability.tsv", "evidence_states.tsv"),
            "Conceptual method illustration with the canonical-parent count of every evidence state",
            c("fig2A_workflow_steps.csv", "fig2B_matching_rules.csv", "fig2C_evidence_states.csv", "fig2C_decision_tree_edges.csv"),
            legend_notes = figure1_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure 2 (manuscript Figure 3). Detection and breadth
# ---------------------------------------------------------------------------------------------
# Native testable cohorts per endpoint (testable, complete and resolved parents). Excluded testability states stay out of every count.
# Platform and global unions cover every condition.
detection <- load_table("detection") %>% filter(population == "all_records")
endpoint_values <- c("eligible_parents", "reference_positive", "reference_rate", "control_draws", "control_positive", "control_rate",
                     "censored", "chance_limited", "partially_censored", "unresolved", "incomplete", "other_inapplicable")
native <- detection %>% filter(cohort == "native")
stopifnot(all(native %>% group_by(order) %>% summarise(n = n_distinct(planned_parents), .groups = "drop") %>% pull(n) == 1))
native <- native %>% select(scope, order, block, group, label, min_length, max_length, planned_parents, endpoint, all_of(endpoint_values)) %>%
  pivot_wider(names_from = endpoint, values_from = all_of(endpoint_values))
common <- detection %>% filter(cohort == "common_strict_tolerant") %>%
  select(order, endpoint, eligible_parents, reference_positive, control_positive, control_draws) %>%
  pivot_wider(names_from = endpoint, values_from = c(eligible_parents, reference_positive, control_positive, control_draws), names_prefix = "common_")
stopifnot(identical(common$eligible_parents_common_strict, common$eligible_parents_common_tolerant))
block_levels <- c("Reference classes", "Platform unions")
classes <- native %>% left_join(common, by = "order") %>%
  mutate(block = factor(block, levels = block_levels)) %>%
  arrange(block, match(group, c(groups, "Global union")), order)
stopifnot(!anyNA(classes$block))
stopifnot(all(classes$eligible_parents_strict + classes$censored_strict + classes$chance_limited_strict + classes$partially_censored_strict +
                classes$unresolved_strict + classes$incomplete_strict + classes$other_inapplicable_strict == classes$planned_parents))
pct <- function(positive, n) ifelse(n > 0, 100 * positive / n, NA_real_)
table1 <- export_csv(data.frame(
  block = as.character(classes$block), platform = classes$group, reference_class = classes$label,
  length_nt = ifelse(classes$min_length == classes$max_length, as.character(classes$min_length), paste0(classes$min_length, "-", classes$max_length)),
  planned_parents_n = classes$planned_parents,
  strict_testable_n = classes$eligible_parents_strict, strict_censored_n = classes$censored_strict, strict_chance_limited_n = classes$chance_limited_strict,
  strict_unresolved_n = classes$unresolved_strict, strict_incomplete_n = classes$incomplete_strict,
  ref_pos_strict = classes$reference_positive_strict, ref_pct_strict = pct(classes$reference_positive_strict, classes$eligible_parents_strict),
  ctrl_draws_strict = classes$control_draws_strict, ctrl_pos_strict = classes$control_positive_strict,
  ctrl_pct_strict = pct(classes$control_positive_strict, classes$control_draws_strict),
  tolerant_testable_n = classes$eligible_parents_tolerant, tolerant_censored_n = classes$censored_tolerant,
  tolerant_chance_limited_n = classes$chance_limited_tolerant, tolerant_partially_censored_n = classes$partially_censored_tolerant,
  tolerant_unresolved_n = classes$unresolved_tolerant, tolerant_incomplete_n = classes$incomplete_tolerant,
  ref_pos_tolerant = classes$reference_positive_tolerant, ref_pct_tolerant = pct(classes$reference_positive_tolerant, classes$eligible_parents_tolerant),
  ctrl_draws_tolerant = classes$control_draws_tolerant, ctrl_pos_tolerant = classes$control_positive_tolerant,
  ctrl_pct_tolerant = pct(classes$control_positive_tolerant, classes$control_draws_tolerant),
  common_parents_n = classes$eligible_parents_common_strict,
  common_ref_pos_strict = classes$reference_positive_common_strict, common_ref_pos_tolerant = classes$reference_positive_common_tolerant,
  common_ctrl_draws = classes$control_draws_common_strict,
  common_ctrl_pos_strict = classes$control_positive_common_strict, common_ctrl_pos_tolerant = classes$control_positive_common_tolerant),
  "table1_reference_classes.csv")
stopifnot(isTRUE(all.equal(table1$ref_pct_strict / 100, ifelse(table1$strict_testable_n > 0, classes$reference_rate_strict, NA_real_))),
          isTRUE(all.equal(table1$ctrl_pct_tolerant / 100, ifelse(table1$tolerant_testable_n > 0, classes$control_rate_tolerant, NA_real_))))

# Panel A: every display class and union, grouped by platform. Strict and tolerant sit in vertically offset sub-rows.
# References: filled markers with 95% Clopper-Pearson intervals. Control draws: open markers on a grey band, the
# relative-entropy Hoeffding interval of the control mean (figure_sources). A rule without testable parents is written out.
control_limits <- load_table("detection") %>% filter(population == "all_records", cohort == "native") %>%
  select(order, endpoint, control_lower, control_upper)
detection_rows <- table1 %>% mutate(order = classes$order,
  facet = ifelse(platform == "Global union", "All platforms", platform),
  row_name = ifelse(block == "Platform unions", "All classes", reference_class)) %>%
  mutate(facet = factor(facet, levels = platform_levels)) %>%
  arrange(facet, block == "Platform unions", order) %>% mutate(y = rev(seq_len(n())))
detection_points <- bind_rows(lapply(c("strict", "tolerant"), function(rule) {
  data <- detection_rows %>% transmute(order, facet, platform, row_name, y, rule = rule,
    n = .data[[paste0(rule, "_testable_n")]], k = .data[[paste0("ref_pos_", rule)]],
    control_draws = .data[[paste0("ctrl_draws_", rule)]], control_positive = .data[[paste0("ctrl_pos_", rule)]])
  data %>% left_join(filter(control_limits, endpoint == rule) %>% select(-endpoint), by = "order")
})) %>% mutate(reference_pct = ifelse(n > 0, 100 * k / n, NA_real_),
  control_pct = ifelse(control_draws > 0, 100 * control_positive / control_draws, NA_real_))
detection_points <- bind_cols(detection_points, ci_binom(detection_points$k, detection_points$n)) %>%
  mutate(reference_lower_pct = 100 * lo, reference_upper_pct = 100 * hi, control_lower_pct = 100 * control_lower,
         control_upper_pct = 100 * control_upper) %>% select(-lo, -hi, -control_lower, -control_upper)
detection_points <- export_csv(detection_points %>% arrange(desc(y), rule) %>%
  transmute(facet = as.character(facet), row_name, rule, testable_parents = n, reference_positive = k, reference_pct, reference_lower_pct,
            reference_upper_pct, control_draws, control_positive, control_pct, control_lower_pct, control_upper_pct, y), "fig3A_detection_points.csv")
stopifnot(all(is.na(detection_points$control_pct) == (detection_points$testable_parents == 0)),
          all(detection_points$control_lower_pct <= detection_points$control_pct + 1e-9, na.rm = TRUE),
          all(detection_points$control_upper_pct >= detection_points$control_pct - 1e-9, na.rm = TRUE))
detection_points <- detection_points %>% mutate(facet = factor(facet, levels = platform_levels),
  ypos = y + ifelse(rule == "strict", 0.2, -0.2), shape = ifelse(rule == "strict", 21, 24),
  colour = unname(pal_platform[as.character(facet)]))
row_info <- detection_rows %>% transmute(facet, y, row_name, strict_n = strict_testable_n, tolerant_n = tolerant_testable_n,
  n_label = paste0(comma(strict_testable_n), " / ", comma(tolerant_testable_n), " (", comma(planned_parents_n), ")"))
untestable_rows <- filter(row_info, strict_n == 0 & tolerant_n == 0)
untestable_rules <- detection_points %>% filter(testable_parents == 0, !(y %in% untestable_rows$y))
label_lookup <- setNames(row_info$row_name, row_info$y)
top_y <- max(row_info$y)
n_column_x <- 104
legend_keys <- c("Strict, reference (95% CI)", "Tolerant, reference (95% CI)", "Strict, control draws", "Tolerant, control draws",
                 "Control interval")
legend_dummy <- data.frame(facet = factor(platform_levels[1], levels = platform_levels), x = 0, y = top_y,
                           key = factor(legend_keys, levels = legend_keys))
detection_plot <- ggplot(filter(detection_points, testable_parents > 0)) +
  geom_rect(aes(xmin = control_lower_pct, xmax = pmax(control_upper_pct, control_lower_pct + 0.25), ymin = ypos - 0.13, ymax = ypos + 0.13),
            fill = grey_band, colour = NA) +
  geom_point(aes(control_pct, ypos, shape = I(shape), colour = I(colour)), fill = "white", size = 1.2, stroke = lw(0.6)) +
  geom_linerange(aes(xmin = reference_lower_pct, xmax = reference_upper_pct, y = ypos, colour = I(colour)), linewidth = lw(0.75)) +
  geom_point(aes(reference_pct, ypos, shape = I(shape), colour = I(colour), fill = I(colour)), size = 1.5, stroke = lw(0.5)) +
  geom_point(data = legend_dummy, aes(x, y, shape = key), alpha = 0, inherit.aes = FALSE) +
  geom_text(data = untestable_rules, aes(x = 3, y = ypos, label = "not testable"), hjust = 0, size = fs(dense_pt), fontface = "italic",
            colour = grey_text) +
  geom_text(data = untestable_rows, aes(x = 1, y = y, label = "not testable at this database size"), hjust = 0, size = fs(dense_pt),
            fontface = "italic", colour = grey_text) +
  geom_text(data = row_info, aes(x = n_column_x, y = y, label = n_label), hjust = 0, size = fs(dense_pt)) +
  geom_text(data = data.frame(facet = factor(platform_levels[1], levels = platform_levels), y = top_y + 1.15),
            aes(x = n_column_x, y = y, label = "Testable parents\nstrict / tolerant (planned)"), hjust = 0, vjust = 0.5, lineheight = 0.95,
            size = fs(dense_pt), fontface = "bold") +
  geom_blank(data = data.frame(facet = factor(platform_levels[1], levels = platform_levels), y = top_y + 1.7), aes(x = 0, y = y)) +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_shape_manual(values = c(21, 24, 1, 2, 15), limits = legend_keys, name = NULL,
    guide = guide_legend(nrow = 2, byrow = FALSE, override.aes = list(alpha = 1, size = c(1.5, 1.5, 1.2, 1.2, 2.6),
      colour = c(ink, ink, ink, ink, grey_band), fill = c(ink, ink, "white", "white", grey_band)))) +
  scale_x_continuous("Parents with an accepted match (%)", breaks = seq(0, 100, 25), expand = expansion(add = c(1, 1))) +
  scale_y_continuous(NULL, breaks = function(limits) row_info$y[row_info$y >= limits[1] & row_info$y <= limits[2]],
                     labels = function(breaks) unname(label_lookup[as.character(breaks)]), expand = expansion(add = 0.55)) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold"), strip.placement = "outside",
        panel.spacing.y = unit(1.2, "mm"), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "top", legend.justification = "left", legend.key.width = unit(2.5, "mm"),
        panel.grid.major.x = element_line(colour = "#EEEEEE", linewidth = lw(0.25)),
        plot.margin = margin(3, 30, 1.5, 1.5, "mm"))

# Panel B: retained accessions per detected parent (strict breadth cohort, reference count above 0), one point per parent.
breadth <- load_table("paired_breadth") %>% filter(group != "Global union") %>% group_order()
breadth_summary <- breadth %>% group_by(group, .drop = FALSE) %>%
  summarise(parents = n(), detected = sum(reference > 0), zero_pairs = sum(reference == 0 & control_mean == 0),
            highest_control_mean = if (n()) max(control_mean) else NA_real_, .groups = "drop")
breadth_points <- export_csv(breadth %>% filter(reference > 0) %>% arrange(group, reference, parent_id) %>%
  transmute(group = as.character(group), parent_id, reference_accessions = reference, control_mean), "fig3B_breadth_per_detected_parent.csv")
breadth_points <- breadth_points %>% mutate(group = factor(group, levels = groups))
breadth_axis <- setNames(paste0(breadth_summary$group, "\nn detected = ", comma(breadth_summary$detected),
                                "\nboth zero = ", comma(breadth_summary$zero_pairs)), as.character(breadth_summary$group))
breadth_breaks <- c(1, 10, 100, 1000)
breadth_breaks <- breadth_breaks[breadth_breaks <= max(10, max(breadth_points$reference_accessions, 1) * 1.5)]
breadth_plot <- ggplot(breadth_points, aes(group, reference_accessions, colour = group)) +
  ggbeeswarm::geom_quasirandom(width = 0.32, size = 0.9, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_discrete(NULL, labels = breadth_axis, drop = FALSE) +
  scale_y_continuous("Retained accessions\nper detected parent", transform = pseudo_log, breaks = breadth_breaks, labels = comma,
                     limits = c(0.8, NA)) +
  theme(axis.text.x = element_text(size = dense_pt, lineheight = 0.95), plot.margin = margin(3, 30, 1.5, 1.5, "mm"))
highest_control <- max(breadth$control_mean)
figure2_legend_notes <- c(
  paste0("A, Detection per display class and per platform union. Circles: strict. Triangles: tolerant. Filled markers: references with 95% ",
         "Clopper-Pearson intervals. Open markers: the mean of the ", draws_per_parent, " control draws per parent, on a grey band that ",
         "shows the relative-entropy Hoeffding interval of the control mean (Methods). The n column gives the testable, complete and resolved ",
         "parents per rule and, in parentheses, the planned parents of the row. Classes share parents, so class rows need not sum to the union row."),
  paste0("B, Retained accessions per detected parent in the strict breadth cohort, one point per parent, on a pseudo-logarithmic axis. ",
         "n detected: parents with a reference count above 0. Both zero: parents with zero reference and control breadth. The highest ",
         "control mean is ", format(signif(highest_control, 2), nsmall = 0), " accessions."))
figure2_heights <- c(4.6 * max(row_info$y) + 30, 46)
figure2 <- assemble(detection_plot, breadth_plot, heights = figure2_heights)
save_figure("figure2", figure2, sum(figure2_heights), c("Detection by class and union with intervals", "Breadth per detected parent"),
  c("detection.tsv", "paired_breadth.tsv"),
  "Descriptive construct-level support in rule-specific testable cohorts. Untestable and chance-limited parents are excluded",
  c("table1_reference_classes.csv", "fig3A_detection_points.csv", "fig3B_breadth_per_detected_parent.csv"),
  legend_notes = figure2_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure 3 (manuscript Figure 4). Positional context (strict, a priori testable events of every condition)
# ---------------------------------------------------------------------------------------------
# A: eligibility of the accessions with strict matches of testable queries. Reasons are grouped by state and, for exclusions,
# listed in the order of application. A single-reason step without records (targeting, RNA, circular) is shown at 0.
reason_levels <- c(
  linear_genomic_dna_known_exclusions_removed = "Linear-genomic DNA",
  record_type_transcript = "Transcript record", record_type_amplicon = "Amplicon record",
  record_type_synthetic = "Synthetic record", record_type_vector = "Vector record",
  origin_synthetic = "Synthetic origin", origin_vector = "Vector origin",
  origin_artificial_other = "Other artificial origin", origin_engineered = "Engineered origin",
  targeted_record = "Targeted record", rna_record = "RNA", circular_record = "Circular",
  molecule_not_established = "Molecule not established", record_type_not_established = "Record type not established",
  topology_not_established = "Topology not established", molecule_conflicting = "Conflicting molecule",
  record_type_conflicting = "Conflicting record type", topology_conflicting = "Conflicting topology",
  origin_conflicting = "Conflicting origin", targeted_conflicting = "Conflicting targeting",
  version_mismatch = "Version mismatch", length_conflict = "Length conflict",
  conflicting_record_evidence = "Conflicting evidence", unresolved_id = "Unresolved identifier",
  identifier_not_normalised = "Identifier not normalised", unresolved_identifier = "Unresolved identifier")
state_labels <- c(included = "Included", excluded = "Excluded\n(in order of\napplication)", unresolved = "Unresolved")
flow <- load_table("eligibility_flow")
stopifnot(all(flow$reason %in% names(reason_levels)), all(flow$state %in% names(state_labels)))
stopifnot(sum(flow$accessions[!(flow$state == "unresolved" & flow$reason == "unresolved_identifier")]) == counts$testable_strict_accessions)
single_steps <- c("targeted_record", "rna_record", "circular_record")
flow <- bind_rows(flow, data.frame(state = "excluded", reason = setdiff(single_steps, flow$reason), accessions = 0))
flow <- export_csv(flow %>% arrange(match(state, names(state_labels)), match(reason, names(reason_levels))) %>%
  transmute(state, reason, reason_label = unname(reason_levels[reason]), accessions), "fig4A_eligibility_flow.csv")
flow_totals <- flow %>% group_by(state) %>% summarise(accessions = sum(accessions), .groups = "drop")
flow_count <- function(name) comma(sum(flow_totals$accessions[flow_totals$state == name]))
flow <- flow %>% mutate(state_label = factor(unname(state_labels[state]), levels = unname(state_labels)),
  reason_label = factor(reason_label, levels = rev(unique(reason_label))))
flow_plot <- ggplot(flow, aes(accessions, reason_label)) +
  geom_segment(aes(x = 0, xend = accessions, yend = reason_label, colour = state, linetype = state), linewidth = lw(1)) +
  geom_point(aes(colour = state, fill = state), shape = 21, size = 1.6, stroke = lw(0.75)) +
  geom_text(aes(label = comma(accessions)), hjust = -0.35, size = fs(dense_pt)) +
  facet_grid(state_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c(included = grey_dark, excluded = grey_mid, unresolved = grey_dark), guide = "none") +
  scale_fill_manual(values = c(included = grey_dark, excluded = grey_mid, unresolved = "white"), guide = "none") +
  scale_linetype_manual(values = c(included = "solid", excluded = "solid", unresolved = "dotted"), guide = "none") +
  scale_x_continuous("Distinct accession.version records", labels = comma, expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(NULL) +
  coord_cartesian(clip = "off") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold", size = dense_pt), axis.ticks.y = element_blank(),
        panel.spacing.y = unit(1.2, "mm"), axis.text.y = element_text(size = dense_pt))

# B: distance to the nearest stored record end, eligible linear-genomic records of at least the primary minimum length.
# Groups with more than five events: empirical cumulative distribution. Sparse groups: every event as a point in a strip.
position_policy <- manifest$position_policy
main_window <- position_policy$primary_window
main_minimum <- position_policy$primary_minimum_length
position <- load_table("position_events") %>% group_order()
main_position <- filter(position, record_length >= main_minimum)
position_totals <- main_position %>% group_by(group) %>%
  summarise(events = n(), accessions = n_distinct(accession_version), .groups = "drop")
dense_groups <- as.character(position_totals$group[position_totals$events > 5])
sparse_groups <- as.character(position_totals$group[position_totals$events <= 5])
invisible(export_csv(main_position %>% arrange(group, distance, accession_version) %>%
  transmute(group = as.character(group), accession_version, record_length, distance,
            display = ifelse(as.character(group) %in% dense_groups, "ecdf", "points")), "fig4B_end_distances.csv"))
sparse_y <- setNames(-0.17 * seq_along(sparse_groups), sparse_groups)
sparse_points <- main_position %>% filter(as.character(group) %in% sparse_groups) %>%
  mutate(part = "sparse", y = unname(sparse_y[as.character(group)]))
dense_points <- main_position %>% filter(as.character(group) %in% dense_groups) %>% mutate(part = "ecdf")
event_label <- function(group) {
  row <- position_totals[as.character(position_totals$group) == group, ]
  paste0(group, " (n = ", comma(row$events), ")")
}
axis_lookup <- setNames(vapply(sparse_groups, event_label, character(1)), as.character(sparse_y))
distance_max <- max(c(10, main_position$distance))
distance_breaks <- c(0, 10^(1:ceiling(log10(distance_max))))
# Direct labels: the earliest-rising curve is labelled right of the W line, below its plateau; every other curve to the right of
# the point where it reaches 40%, below the curve.
curve_labels <- bind_rows(lapply(dense_groups, function(name) {
  values <- dense_points$distance[as.character(dense_points$group) == name]
  data.frame(group = name, rise = unname(quantile(values, 0.4, type = 1)), q85 = unname(quantile(values, 0.85, type = 1)),
             q40 = unname(quantile(values, 0.4, type = 1)), label = event_label(name))
}))
if (nrow(curve_labels)) {
  curve_labels <- curve_labels %>% arrange(rise) %>%
    mutate(first = row_number() == 1, x = ifelse(first, main_window * 1.4, pmax(q40, 1) * 1.8), y = ifelse(first, 0.8, 0.26),
           hjust = 0, vjust = 0.5, part = "ecdf")
}
ecdf_breaks <- seq(0, 1, 0.25)
distance_plot <- ggplot() +
  geom_vline(xintercept = main_window, linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  stat_ecdf(data = dense_points, aes(distance, colour = group), geom = "step", pad = FALSE, linewidth = lw(0.75)) +
  geom_point(data = sparse_points, aes(distance, y, colour = group), size = 1.3, alpha = 0.9) +
  geom_text(data = data.frame(part = "ecdf"), aes(x = main_window * 1.4, y = 0.55, label = paste0("W = ", comma(main_window), " nt")),
            hjust = 0, size = fs(dense_pt)) +
  (if (nrow(curve_labels)) geom_text(data = curve_labels, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
                                     size = fs(dense_pt))) +
  facet_grid(factor(part, levels = c("ecdf", "sparse")) ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous("Distance to nearest record end (nt)", transform = pseudo_log, breaks = distance_breaks, labels = label_power_expr,
                     limits = c(0, distance_max * 1.05)) +
  scale_y_continuous("Cumulative events (%)", breaks = function(limits) if (limits[2] > 0.5) ecdf_breaks else sort(unname(sparse_y)),
                     labels = function(breaks) ifelse(breaks >= 0, as.character(round(100 * breaks)), unname(axis_lookup[as.character(breaks)])),
                     expand = expansion(add = 0.07)) +
  coord_cartesian(clip = "off") +
  theme(strip.text = element_blank(), panel.spacing.y = unit(1.5, "mm"))

# C: share of events within W of a record end, observed with 95% Clopper-Pearson interval, against the share expected under
# uniform placement (vertical tick). One event per unit is chosen by a position-independent hash.
sensitivity <- load_table("position_sensitivity") %>% filter(group != "Global union") %>% group_order()
variant_labels <- c(all_events = "All events", one_per_accession = "One per accession", one_per_bioproject = "One per BioProject",
                    excluding_targeted_insdc = "Without INSDC-targeted records")
variant_labels <- variant_labels[names(variant_labels) %in% unique(load_table("position_sensitivity")$variant)]
variant_shapes <- c("All events" = 16, "One per accession" = 1, "One per BioProject" = 2, "Without INSDC-targeted records" = 5)[unname(variant_labels)]
stopifnot(all(sensitivity$variant %in% names(variant_labels)))
terminal <- export_csv(sensitivity %>% filter(window == main_window, minimum_length == main_minimum) %>%
  arrange(match(as.character(group), groups), match(variant, names(variant_labels))) %>%
  transmute(group = as.character(group), variant, unit_rule, events, terminal, expected, ratio, units_without_project,
            observed_fraction = ifelse(events > 0, terminal / events, NA_real_), expected_fraction = ifelse(events > 0, expected / events, NA_real_),
            observed_lower = ci_binom(terminal, events)$lo, observed_upper = ci_binom(terminal, events)$hi,
            ratio_upper_limit = ifelse(expected > 0, events / expected, NA_real_)), "fig4C_terminal_ratios.csv")
terminal <- terminal %>% mutate(group = factor(group, levels = groups),
  variant = factor(unname(variant_labels[variant]), levels = rev(unname(variant_labels))),
  kn_label = paste0(comma(terminal), "/", comma(events)),
  oe_label = ifelse(is.na(ratio), "n.e.", ifelse(expected < 1, "n.e.", paste0(formatC(ratio, format = "f", digits = 1), "\u00d7"))),
  oe_colour = ifelse(is.na(ratio) | expected < 1, grey_text, ink))
header_group <- factor(groups[groups %in% terminal$group][1], levels = groups)
terminal_header <- data.frame(group = header_group, variant = factor(rev(levels(terminal$variant))[1], levels = levels(terminal$variant)))
terminal_plot <- ggplot(terminal, aes(y = variant)) +
  geom_linerange(aes(xmin = 100 * observed_lower, xmax = 100 * observed_upper, colour = group), linewidth = lw(0.75), na.rm = TRUE) +
  geom_point(aes(x = 100 * observed_fraction, shape = variant, colour = group), fill = "white", size = 1.4, na.rm = TRUE) +
  geom_point(aes(x = 100 * expected_fraction), shape = 124, size = 2.2, colour = ink, na.rm = TRUE) +
  geom_text(aes(x = 106, label = kn_label), hjust = 0, size = fs(dense_pt)) +
  geom_text(aes(x = 138, label = oe_label, colour = I(oe_colour)), hjust = 0, size = fs(dense_pt)) +
  geom_text(data = terminal_header, aes(x = 106, y = 4.85, label = "k/n"), hjust = 0, size = fs(dense_pt), fontface = "bold") +
  geom_text(data = terminal_header, aes(x = 138, y = 4.85, label = "O/E"), hjust = 0, size = fs(dense_pt), fontface = "bold") +
  facet_grid(group ~ ., switch = "y") +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(values = variant_shapes, breaks = unname(variant_labels), name = NULL,
                     guide = guide_legend(ncol = 2, override.aes = list(colour = ink))) +
  scale_x_continuous(paste0("Events within ", comma(main_window), " nt of a record end (%)"), breaks = seq(0, 100, 25),
                     expand = expansion(add = 2)) +
  scale_y_discrete(NULL, expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold", size = dense_pt), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), panel.spacing.y = unit(1.2, "mm"), legend.position = "bottom",
        legend.text = element_text(size = dense_pt), legend.key.height = unit(2.6, "mm"), plot.margin = margin(3, 26, 1.5, 1.5, "mm"))
figure3_legend_notes <- c(
  paste0("A, Eligibility of the ", flow_count("included"), " + ", flow_count("excluded"), " + ", flow_count("unresolved"), " = ",
         comma(sum(flow_totals$accessions)), " accession.version records with strict matches of testable queries. ",
         "Exclusion reasons are applied in a fixed order (record type, origin, targeting, RNA, circular), and each record counts under its ",
         "first reason. Excluded records keep their matches in the all-record population."),
  paste0("B, Distance from strict, orientation-deduplicated events of testable queries to the nearest stored end of eligible linear-genomic ",
         "records of at least ", comma(main_minimum / 1000), " kb, on a pseudo-logarithmic axis. Groups with more than five events are shown as ",
         "empirical cumulative distributions. Sparse groups show every event as a point. Dashed line: W = ", comma(main_window), " nt."),
  paste0("C, Share of events within W = ", comma(main_window), " nt of a record end in the same population. Points: observed share k/n with ",
         "95% Clopper-Pearson interval. Vertical tick: the share expected under uniform placement, E/n, with the terminal probability ",
         "min(1, 2(W + 1)/(L \u2212 a + 1)) per event. O/E: observed / expected terminal events. n.e.: not estimated because E < 1. ",
         "One event per accession or BioProject is chosen by a position-independent hash. A stored record end is not necessarily a ",
         "chromosome end."))
figure3_heights <- c(42, 70)
figure3 <- assemble(flow_plot, distance_plot, terminal_plot, design = "AA\nBC", heights = figure3_heights)
save_figure("figure3", figure3, sum(figure3_heights) + 3, c("Accession eligibility", "End-distance distribution", "Terminal share with intervals"),
  c("eligibility_flow.tsv", "position_events.tsv", "position_sensitivity.tsv"),
  "Strict, a priori testable events of every condition. Ordered eligibility reasons. Descriptive geometry with unit-level sensitivities",
  c("fig4A_eligibility_flow.csv", "fig4B_end_distances.csv", "fig4C_terminal_ratios.csv"), legend_notes = figure3_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure 4 (manuscript Figure 5). Short-query calibration and record context
# ---------------------------------------------------------------------------------------------
# A: a-priori testability from the production calibration (length_testability) and benchmark planted exact-spike recovery at E = 10
# for three database sizes. E-value at the minimum qualifying score decides score retainability; the iid chance bound decides
# informativeness. Both are expected chance counts on one axis; the limits are read from the table.
class_labels <- c(testable = "Testable", chance_limited = "Chance-limited", censored_untestable = "Untestable")
curves <- export_csv(testability %>% transmute(rule = endpoint, length_nt = query_length, testability, testability_reason,
    score_retainable, partially_score_censored, minimum_qualifying_score, word_size, effective_search_space,
    evalue_at_minimum_score, expected_chance_bound, production_evalue, informativeness_limit, testable_chance_limit) %>%
  arrange(rule, length_nt), "fig5A_testability_curves.csv")
stopifnot(all(curves$testability %in% names(class_labels)))
spikes <- load_table("short_query_detectability")
database_labels <- c(core_nt_div10 = "core_nt / 10", core_nt = "core_nt", core_nt_x10 = "core_nt \u00d7 10")
stopifnot(all(spikes$database %in% names(database_labels)), all(spikes$permissive_diagnostic_recovered == 1))
spike_summary <- export_csv(spikes %>% group_by(database, length_nt = length) %>%
  summarise(parents = n_distinct(parent_id), planned_spikes = n(), calibrated_recovered = sum(production_calibrated_recovered),
    diagnostic_recovered = sum(permissive_diagnostic_recovered), exact_evalue = unique(calibrated_exact_evalue),
    effective_search_space = unique(effective_search_space), .groups = "drop") %>%
  arrange(match(database, names(database_labels)), length_nt), "fig5A_spike_recovery.csv")
length_breaks <- c(8, 10, 16, 24, 32, 40, 61, 90)
strict_label <- "Strict (exact full length)"
tolerant_label <- "Tolerant (minimum qualifying score)"
# The benchmark is named "benchmark" in visible text, with a multiplication sign for "x 10".
benchmark_series <- "Exact E-value at core_nt / 10 and \u00d7 10 (benchmark)"
missed_marker <- "Benchmark exact spikes of testable lengths missed at E = 10"
testable_label <- "Testable panel length"
y_floor <- 1e-12
y_ceiling <- 1e9
curve_data <- curves %>% mutate(rule = factor(recode(rule, strict = strict_label, tolerant = tolerant_label), levels = c(strict_label, tolerant_label)),
  class = factor(unname(class_labels[testability]), levels = unname(class_labels)),
  below = evalue_at_minimum_score < y_floor, y_plot = pmax(evalue_at_minimum_score, y_floor))
# Points mark the lengths that are testable under the rule of the facet. The curves cover every panel length.
testable_points <- curve_data %>% filter(testability == "testable") %>% mutate(point_label = factor(testable_label))
stopifnot(all(c(strict_label, tolerant_label) %in% testable_points$rule))
curve_lines <- bind_rows(
  curve_data %>% transmute(rule, length_nt, series = "E-value at minimum score (core_nt)", expectation = evalue_at_minimum_score),
  curve_data %>% transmute(rule, length_nt, series = "Chance bound (core_nt)", expectation = expected_chance_bound),
  spike_summary %>% filter(database != "core_nt") %>%
    transmute(rule = factor(strict_label, levels = levels(curve_data$rule)), length_nt, database,
              series = benchmark_series, expectation = exact_evalue)) %>%
  mutate(database = coalesce(database, "core_nt"))
core_spikes <- spike_summary %>% filter(database == "core_nt") %>% inner_join(curves %>% filter(rule == "strict"), by = "length_nt")
stopifnot(nrow(core_spikes) == n_distinct(spike_summary$length_nt), isTRUE(all.equal(core_spikes$exact_evalue, core_spikes$evalue_at_minimum_score)))
strict_testable_lengths <- curves$length_nt[curves$rule == "strict" & curves$testability == "testable"]
all_missed <- spike_summary %>% filter(calibrated_recovered < planned_spikes)
missed <- all_missed %>% filter(length_nt %in% strict_testable_lengths) %>%
  mutate(rule = factor(strict_label, levels = levels(curve_data$rule)), y_plot = pmin(pmax(exact_evalue, y_floor), y_ceiling),
         marker = missed_marker)
production_e <- unique(curves$production_evalue)
stopifnot(all(all_missed$exact_evalue > production_e),
          all(with(spike_summary, calibrated_recovered[exact_evalue <= production_e] == planned_spikes[exact_evalue <= production_e])))
spike_text <- if (nrow(missed)) {
  paste0("Crosses: benchmark exact spikes of testable lengths missed at E = 10 (", serial(sort(unique(missed$length_nt))), " nt).")
} else "All benchmark exact spikes of testable lengths were recovered at E = 10."
size_labels <- spike_summary %>% filter(database != "core_nt") %>% group_by(database) %>%
  filter(exact_evalue <= y_ceiling) %>% slice_min(length_nt, n = 1, with_ties = TRUE) %>% ungroup() %>%
  mutate(rule = factor(strict_label, levels = levels(curve_data$rule)), label = recode(database, core_nt_div10 = "/10", core_nt_x10 = "\u00d710"))
limit_lines <- data.frame(y = c(production_e, unique(curves$informativeness_limit), unique(curves$testable_chance_limit)),
  label = c(paste0("E = ", production_e), paste0("chance ", unique(curves$informativeness_limit)),
            paste0("chance ", unique(curves$testable_chance_limit))),
  colour = c(accent_alert, grey_dark, grey_dark), linetype = c("solid", "dotted", "dashed"), vjust = c(-0.35, -0.35, 1.35))
# Every panel length is a point, coded by its a priori state under the rule of the facet (same coding as Figure 2C).
state_names <- c(testable = "testable", chance_limited = "chance-limited", censored_untestable = "untestable")
state_keys <- c(testable = "Testable length", `chance-limited` = "Chance-limited length", untestable = "Untestable length")
length_points <- curve_data %>% mutate(state = factor(unname(state_names[testability]), levels = names(state_point_fill)))
first_testable_length <- curve_data %>% filter(testability == "testable") %>% group_by(rule) %>%
  summarise(length_nt = min(length_nt), .groups = "drop") %>%
  mutate(label = paste0(ifelse(rule == strict_label, "strict", "tolerant"), " \u2265 ", length_nt, " nt"))
stopifnot(nrow(first_testable_length) == 2)
floor_label <- function(breaks) lapply(breaks, function(value) {
  exponent <- round(log10(value))
  power <- if (exponent < 0) bquote(10^.(paste0("\u2212", -exponent))) else bquote(10^.(as.character(exponent)))
  if (abs(value - y_floor) < 1e-24) bquote("<" * .(power)) else power
})
curve_types <- c("E-value at minimum score" = "solid", "Chance bound" = "longdash")
curve_lines <- curve_lines %>% mutate(series = recode(series, "E-value at minimum score (core_nt)" = "E-value at minimum score",
                                                      "Chance bound (core_nt)" = "Chance bound"))
evalue_plot <- ggplot(curve_data, aes(length_nt)) +
  geom_hline(data = limit_lines, aes(yintercept = y, colour = I(colour), linetype = I(linetype)), linewidth = lw(0.5)) +
  geom_vline(data = first_testable_length, aes(xintercept = length_nt), linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  geom_text(data = first_testable_length, aes(x = length_nt * 1.04, y = 3e8, label = label), hjust = 0, size = fs(dense_pt)) +
  geom_line(data = filter(curve_lines, series != benchmark_series), aes(y = expectation, linetype = series, group = interaction(series, database)),
            colour = grey_mid, linewidth = lw(0.75)) +
  geom_line(data = filter(curve_lines, series == benchmark_series), aes(y = expectation, group = interaction(series, database)),
            colour = grey_mid, linewidth = lw(0.6), linetype = "dotted") +
  geom_point(data = length_points, aes(y = y_plot, fill = state), shape = 21, colour = grey_dark, size = 1.5, stroke = lw(0.5)) +
  (if (nrow(missed)) geom_point(data = missed, aes(y = y_plot), shape = 4, colour = accent_alert, size = 1.8, stroke = lw(1))) +
  geom_text(data = size_labels, aes(x = length_nt * 0.93, y = exact_evalue, label = label), hjust = 1, size = fs(dense_pt), colour = grey_text) +
  geom_text(data = mutate(limit_lines, rule = factor(tolerant_label, levels = levels(curve_data$rule))),
    aes(x = 125, y = y, label = label, colour = I(colour), vjust = vjust), hjust = 1, size = fs(dense_pt)) +
  facet_wrap(~rule, nrow = 1) +
  scale_fill_manual(values = state_point_fill, labels = state_keys, name = NULL, drop = FALSE,
                    guide = guide_legend(order = 1, override.aes = list(size = 1.6))) +
  scale_linetype_manual(values = curve_types, name = NULL, guide = guide_legend(order = 2, override.aes = list(colour = grey_mid))) +
  scale_x_log10("Query length (nt)", breaks = length_breaks, limits = c(6, 125), expand = expansion(mult = 0)) +
  scale_y_log10("Expected chance count", breaks = 10^seq(-12, 9, 3), labels = function(breaks) as.expression(floor_label(breaks))) +
  coord_cartesian(ylim = c(y_floor, y_ceiling), expand = FALSE, clip = "on") +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.key.width = unit(5, "mm"), legend.spacing.x = unit(1, "mm"),
        strip.text = element_text(hjust = 0.5), panel.spacing.x = unit(4, "mm"))

case_headings <- c(terminal = "Terminal match", synthetic_or_vector = "Synthetic or vector record", internal_or_ambiguous = "Internal match")
# Query names follow the terminology of the text; the identifier of the source record stays in the exported CSV and the caption.
query_names <- c(ILLUMINA_TRUSEQ_R2_adapter = "TruSeq R2 adapter", RAPID_adapter_forward = "Rapid adapter")
display_query <- function(id) ifelse(id %in% names(query_names), unname(query_names[id]), parent_label(id))
case_flank_shown <- 25
# Record cases come from the record_cases stage (deterministic selection, 0-3 cases); every case is verified at sequence level.
cases <- load_table("record_cases")
verification <- load_table("record_case_verification")
cases <- cases[!is.na(cases$accession_version), ]
case_columns <- c("panel", "category", "accession_version", "organism", "organism_italic", "topology", "molecule", "record_type", "record_context",
  "record_length_nt", "match_start", "match_end", "match_length_nt", "distance_to_nearest_end_nt", "group", "query_id", "query_label",
  "display_orientation", "query_sequence", "subject_sequence", "display_left_flank", "display_right_flank", "segment_roles")
if (nrow(cases)) {
  stopifnot(nrow(cases) <= 3, identical(verification$accession_version, cases$accession_version), all(cases$category %in% names(case_headings)),
    !anyDuplicated(cases$category), all(verification$query_parent_matches), all(verification$searched_sequence_found),
    all(verification$hsp_full_length_exact), all(verification$record_window_equals_hsp), all(verification$displayed_window_equals_record),
    all(verification$packet_identical | is.na(verification$packet_identical)), all(cases$qseq == cases$sseq),
    all(cases$upper - cases$lower + 1 == nchar(cases$qseq)),
    all(lengths(strsplit(cases$display_segment_roles, ";", fixed = TRUE)) == nchar(cases$qseq)))
  cases <- cases[order(match(cases$category, names(case_headings))), ]
  verification <- verification[match(cases$accession_version, verification$accession_version), ]
  cases$display_left_flank[is.na(cases$display_left_flank)] <- ""
  cases$display_right_flank[is.na(cases$display_right_flank)] <- ""
  case_data <- export_csv(data.frame(
    panel = LETTERS[1 + seq_len(nrow(cases))], category = cases$category, accession_version = cases$accession_version,
    organism = cases$organism, organism_italic = verification$organism_italic, topology = cases$topology, molecule = cases$molecule,
    record_type = cases$record_type, record_context = cases$context, record_length_nt = cases$record_length, match_start = cases$lower,
    match_end = cases$upper, match_length_nt = nchar(cases$qseq),
    distance_to_nearest_end_nt = ifelse(cases$topology == "circular", NA, pmin(cases$lower - 1, cases$record_length - cases$upper)),
    group = cases$group, query_id = cases$query_id, query_label = parent_label(cases$query_id),
    display_orientation = cases$display_orientation, query_sequence = cases$qseq, subject_sequence = cases$sseq,
    display_left_flank = cases$display_left_flank, display_right_flank = cases$display_right_flank,
    segment_roles = cases$display_segment_roles), "fig5BD_record_cases.csv")
  stopifnot(all(case_data$distance_to_nearest_end_nt[case_data$topology == "linear"] == cases$distance_to_end[cases$topology == "linear"]))
  case_data$display_left_flank[is.na(case_data$display_left_flank)] <- ""
  case_data$display_right_flank[is.na(case_data$display_right_flank)] <- ""
} else {
  case_data <- export_csv(setNames(data.frame(matrix(character(), nrow = 0, ncol = length(case_columns))), case_columns), "fig5BD_record_cases.csv")
}
invisible(export_csv(data.frame(
  panel = case_data$panel, accession_version = case_data$accession_version, organism = case_data$organism, topology = case_data$topology,
  record_length_nt = case_data$record_length_nt, match_start = case_data$match_start, match_end = case_data$match_end,
  distance_to_nearest_end_nt = case_data$distance_to_nearest_end_nt, query = case_data$query_label,
  match_length_nt = case_data$match_length_nt, case_category = unname(case_headings[case_data$category])), "tableS8_record_case_summary.csv"))

# B-D: one record diagram per case. Linear records: the record bar with a scale bar and, when the match lies too close to an end
# to be resolved at panel scale, a zoom of the record end. Circular records: a ring with position 1 at the top.
# Sequences: the record strand that carries the query in its searched orientation (minus strand for minus-strand matches),
# a fixed window of the same width for every case, coordinates at both ends, the match bold on a light tint of the platform
# colour of the query, barcode or index segments on a darker tint.
track_x0 <- 20
track_x1 <- if (any(case_data$topology == "linear")) 112 else figure_width_mm - 3
inset_x0 <- 126
inset_x1 <- figure_width_mm - 3
minimum_marker_mm <- 0.9
sequence_x0 <- 36
# Courier New advances 0.6 em, so 6-pt bases need 1.27 mm and fit the 1.29-mm pitch of a 91-nt window.
sequence_pt <- 6
case_window <- if (nrow(case_data)) max(case_data$match_length_nt) + 40 else 0
sequence_dx <- if (nrow(case_data)) min(1.33, (figure_width_mm - sequence_x0 - 17) / case_window) else 1.33
nice_length <- function(target) {
  candidates <- as.vector(outer(c(1, 2, 5), 10^(0:9)))
  max(candidates[candidates <= target])
}
length_text <- function(nt) {
  if (nt >= 1e6) paste0(nt / 1e6, " Mb") else if (nt >= 1e3) paste0(nt / 1e3, " kb") else paste0(nt, " nt")
}
case_header <- function(row, y) {
  quote <- function(text) paste0('"', gsub('"', '\\\\"', text), '"')
  type_parts <- c(row$topology, if (row$record_type == "genomic") "genomic" else if (row$record_context == "vector") "vector", row$molecule)
  strand <- ifelse(row$display_orientation == "-", "minus", "plus")
  rest <- paste0(", ", paste(type_parts, collapse = " "), ", ", comma(row$record_length_nt), " nt")
  organism <- if (isTRUE(row$organism_italic)) paste0("italic(", quote(row$organism), ")") else quote(row$organism)
  label <- paste0("bold(", quote(paste0(case_headings[[row$category]], ", ", row$accession_version, ".  ")), ")*", organism, "*", quote(rest))
  list(mm_text(5, y, label, base_pt, parse = TRUE),
       mm_text(5, y + 3.4, paste0(display_query(row$query_id), ", ", row$match_length_nt, "-nt exact match on the ", strand, " strand"), base_pt))
}
record_case_panel <- function(row) {
  circular <- row$topology == "circular"
  marker_colour <- unname(pal_platform[row$group])
  stopifnot(!is.na(marker_colour))
  length_nt <- row$record_length_nt
  layers <- case_header(row, 3.2)
  if (circular) {
    centre_x <- track_x0 + 8; centre_y <- 16; radius <- 6.2
    angle <- seq(0, 2 * pi, length.out = 361)
    match_angle <- 2 * pi * ((row$match_start - 1 + row$match_end) / 2) / length_nt
    marker <- c(sin(match_angle), -cos(match_angle))
    layers <- c(layers, list(
      annotate("path", x = centre_x + radius * sin(angle), y = centre_y - radius * cos(angle), colour = grey_light, linewidth = lw(3)),
      mm_segment(centre_x, centre_y - radius + 1.1, centre_x, centre_y - radius - 1.1, colour = ink, linewidth = lw(0.75)),
      mm_segment(centre_x + (radius - 2) * marker[1], centre_y + (radius - 2) * marker[2],
                 centre_x + (radius + 2) * marker[1], centre_y + (radius + 2) * marker[2], colour = marker_colour, linewidth = lw(3)),
      mm_text(centre_x, centre_y, paste0(comma(length_nt), "\nnt"), dense_pt, hjust = 0.5),
      mm_text(centre_x + radius + 4, centre_y - 2.2, paste0("Match ", comma(row$match_start), "\u2013", comma(row$match_end)), dense_pt),
      mm_text(centre_x + radius + 4, centre_y + 2.2, "Position 1 at the top tick", dense_pt)))
  } else {
    position_x <- function(position, x0 = track_x0, x1 = track_x1, from = 0, to = length_nt) x0 + (x1 - x0) * (position - from) / (to - from)
    marker_centre <- position_x((row$match_start - 1 + row$match_end) / 2)
    marker_width <- max(position_x(row$match_end) - position_x(row$match_start - 1), minimum_marker_mm)
    nearest_left <- (row$match_start - 1) <= (length_nt - row$match_end)
    distance <- row$distance_to_nearest_end_nt
    scale_nt <- nice_length(length_nt / 4)
    layers <- c(layers, list(
      mm_rect(track_x0, track_x1, 12.8, 15.2, fill = "#EEEEEE", colour = grey_mid),
      mm_rect(marker_centre - marker_width / 2, marker_centre + marker_width / 2, 11.6, 16.4, fill = marker_colour),
      mm_text(track_x0, 18.2, "1", dense_pt, hjust = 0.5), mm_text(track_x1, 18.2, comma(length_nt), dense_pt, hjust = 1),
      mm_segment(track_x0, 21.6, position_x(scale_nt), 21.6, colour = ink, linewidth = lw(1)),
      mm_text(position_x(scale_nt) + 1.5, 21.6, length_text(scale_nt), dense_pt),
      mm_text(marker_centre + ifelse(marker_centre > (track_x0 + track_x1) / 2, marker_width / 2, -marker_width / 2), 9.6,
              paste0("Match ", comma(row$match_start), "\u2013", comma(row$match_end)), dense_pt,
              hjust = ifelse(marker_centre > (track_x0 + track_x1) / 2, 1, 0))))
    # Zoom of the nearest record end: window of 1, 2 or 5 x 10^k nt, at least 100 nt and at least 1.3 times distance plus match.
    zoom_nt <- max(100, min(as.vector(outer(c(1, 2, 5), 10^(2:9)))[as.vector(outer(c(1, 2, 5), 10^(2:9))) >= 1.3 * (distance + row$match_length_nt)]))
    zoom_from <- if (nearest_left) 0 else length_nt - zoom_nt
    zoom_to <- if (nearest_left) zoom_nt else length_nt
    zoom_x <- function(position) position_x(position, inset_x0, inset_x1, zoom_from, zoom_to)
    zoom_start <- zoom_x(row$match_start - 1)
    zoom_end <- zoom_x(row$match_end)
    end_x <- if (nearest_left) inset_x0 else inset_x1
    edge_x <- if (nearest_left) zoom_start else zoom_end
    layers <- c(layers, list(
      mm_segment(position_x(zoom_from), 15.2, inset_x0, 15.2, colour = grey_mid, linewidth = lw(0.25), linetype = "dotted"),
      mm_segment(position_x(zoom_to), 12.8, inset_x1, 12.8, colour = grey_mid, linewidth = lw(0.25), linetype = "dotted"),
      mm_rect(inset_x0, inset_x1, 12.8, 15.2, fill = "#EEEEEE", colour = grey_mid),
      mm_rect(zoom_start, max(zoom_end, zoom_start + minimum_marker_mm), 11.6, 16.4, fill = marker_colour),
      mm_text((inset_x0 + inset_x1) / 2, 9.6, paste0(if (nearest_left) "First " else "Last ", length_text(zoom_nt), " of the record"),
              dense_pt, hjust = 0.5),
      mm_text(inset_x0, 18.2, comma(zoom_from + 1), dense_pt, hjust = 0), mm_text(inset_x1, 18.2, comma(zoom_to), dense_pt, hjust = 1),
      annotate("segment", x = min(end_x, edge_x), xend = max(end_x, edge_x), y = 21.6, yend = 21.6, colour = ink, linewidth = lw(0.5),
               arrow = arrow(length = unit(1, "mm"), ends = "both", type = "closed", angle = 22)),
      mm_text(if (nearest_left) inset_x0 else inset_x1, 24, paste0(comma(distance), " nt to the record end"), dense_pt,
              hjust = if (nearest_left) 0 else 1)))
  }
  # Sequence window: case_window bases with the match centred where the record allows it.
  left_flank <- row$display_left_flank
  right_flank <- row$display_right_flank
  flanks <- case_window - row$match_length_nt
  left_n <- min(nchar(left_flank), max(ceiling(flanks / 2), flanks - nchar(right_flank)))
  right_n <- min(nchar(right_flank), flanks - left_n)
  left_shown <- if (left_n) substring(left_flank, nchar(left_flank) - left_n + 1) else ""
  right_shown <- if (right_n) substring(right_flank, 1, right_n) else ""
  subject <- strsplit(paste0(left_shown, row$subject_sequence, right_shown), "")[[1]]
  query <- strsplit(row$query_sequence, "")[[1]]
  roles <- strsplit(row$segment_roles, ";", fixed = TRUE)[[1]]
  stopifnot(length(roles) == length(query))
  minus <- row$display_orientation == "-"
  first_position <- if (minus) row$match_end + left_n else row$match_start - left_n
  last_position <- if (minus) row$match_start - right_n else row$match_end + right_n
  x <- sequence_x0 + (seq_along(subject) - 0.5) * sequence_dx
  match_index <- left_n + seq_along(query)
  emphasis <- roles %in% c("barcode", "index")
  record_y <- 29.6
  query_y <- 33.6
  record_name <- if (minus) "Record (\u2212)" else "Record (+)"
  coordinate_x0 <- sequence_x0 - 4.2
  coordinate_x1 <- x[length(x)] + sequence_dx / 2 + 4.2
  tint_runs <- rle(emphasis)
  tint_end <- cumsum(tint_runs$lengths)
  tint_start <- tint_end - tint_runs$lengths + 1
  layers <- c(layers, list(
    mm_rect(sequence_x0 + (match_index[tint_start] - 1) * sequence_dx, sequence_x0 + match_index[tint_end] * sequence_dx,
            record_y - 1.9, query_y + 1.9, fill = marker_colour, alpha = ifelse(tint_runs$values, 0.45, 0.2)),
    mm_text(0, record_y, record_name, dense_pt), mm_text(0, query_y, "Query", dense_pt),
    mm_text(coordinate_x0 - 1.8, record_y, comma(first_position), dense_pt, hjust = 1),
    mm_text(coordinate_x0, record_y, "5\u2032", dense_pt, hjust = 0.5), mm_text(coordinate_x1, record_y, "3\u2032", dense_pt, hjust = 0.5),
    mm_text(coordinate_x1 + 1.8, record_y, comma(last_position), dense_pt, hjust = 0),
    mm_text(x[-match_index], record_y, subject[-match_index], sequence_pt, hjust = 0.5, family = FONT_MONO, colour = grey_text),
    mm_text(x[match_index], record_y, subject[match_index], sequence_pt, "bold", hjust = 0.5, family = FONT_MONO),
    mm_text(x[match_index], query_y, query, sequence_pt, "bold", hjust = 0.5, family = FONT_MONO),
    mm_text(x[min(match_index)] - sequence_dx / 2 - 2.4, query_y, "5\u2032", dense_pt, hjust = 0.5),
    mm_text(x[max(match_index)] + sequence_dx / 2 + 2.4, query_y, "3\u2032", dense_pt, hjust = 0.5)))
  canvas_plot(case_height, layers)
}
case_height <- 37
case_panels <- lapply(seq_len(nrow(case_data)), function(index) record_case_panel(case_data[index, ]))
case_heights <- rep(case_height, length(case_panels))
if (!length(case_panels)) {
  case_panels <- list(canvas_plot(12, list(mm_text(5, 5, "No record case was selected. The record-case selection tables list the category status."))))
  case_heights <- 12
}
case_ids <- if (nrow(case_data)) paste0(case_data$panel, ", ", case_data$query_id, collapse = "; ") else ""
figure4_legend_notes <- c(
  paste0("A, A priori testability against query length for the accepted core_nt dimensions at E = ", production_e, ". Points: every panel length ",
         "at the E-value of its minimum qualifying score, coded by its testability state. Solid curve: E-value at the minimum qualifying score. ",
         "Long-dashed curve: chance bound. Dotted curves (\u00d710, /10): exact E-values at ten times and one tenth of the database size ",
         "(benchmark). Horizontal lines: E = ", production_e, " and the chance limits ", unique(curves$informativeness_limit), " and ",
         unique(curves$testable_chance_limit), ". Vertical dashed lines: the shortest testable length per rule. ", spike_text),
  "B\u2013D, The record cases, selected by a fixed rule from strict matches of testable queries and verified against the local database.",
  paste0("Linear records are drawn to scale with a scale bar; the zoom shows the record end nearest to the match. Position 1 of a circular ",
         "record is its stored origin. Match markers and sequence tints use the platform colour of the query; the darker tint marks the ",
         "barcode or index segment. Sequences show ", case_window, " nt of the record strand that carries the query in its searched orientation, ",
         "with record coordinates at both ends. Query identifiers: ", case_ids, "."))
figure4_heights <- c(64, case_heights)
# The record diagrams are millimetre canvases over the full width, so they are not aligned with the plot area of panel A.
figure4 <- assemble(c(list(evalue_plot), lapply(case_panels, function(panel) wrap_elements(full = panel))), heights = figure4_heights)
save_figure("figure4", figure4, sum(figure4_heights),
  c("A priori testability and exact-spike recovery against query length", paste("Record diagram:", unname(case_headings[case_data$category]))),
  c("length_testability.tsv", "short_query_detectability.tsv", "record_cases.tsv", "record_case_verification.tsv"),
  "Testability classes from the production calibration. Benchmark exact spikes at three database sizes. Record cases from the record_cases stage",
  c("fig5A_testability_curves.csv", "fig5A_spike_recovery.csv", "fig5BD_record_cases.csv", "tableS8_record_case_summary.csv"),
  legend_notes = figure4_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure S1. Query properties and matching diagnostics
# ---------------------------------------------------------------------------------------------
condition_names <- c(
  ont_barcode = "ONT: Native barcode", ont_native_barcode_context_k40 = "ONT: K40 context",
  ont_native_barcode_context_k46 = "ONT: K46 context", ont_rapid_barcode_context_90 = "ONT: Rapid barcode context",
  ont_ligation_adapter = "ONT: Ligation adapter", ont_rapid_adapter = "ONT: Rapid adapter", ont_pcr_adapter = "ONT: PCR adapter",
  ont_cdna_adapter = "ONT: cDNA adapter",
  illumina_i7_index = "Illumina: bare i7", illumina_i5_index = "Illumina: bare i5",
  illumina_i7_index_adapter = "Illumina: embedded i7", illumina_i5_index_adapter = "Illumina: embedded i5",
  illumina_truseq = "Illumina: TruSeq", illumina_nextera = "Illumina: Nextera", illumina_smallrna = "Illumina: Small RNA",
  illumina_bgi = "BGI/DNBSEQ: adapters", pacbio_barcode = "PacBio: barcode", pacbio_barcode_adapter = "PacBio: barcode-hairpin",
  pacbio_smrtbell_adapter = "PacBio: SMRTbell / C2", pacbio_c2_primer = "PacBio: C2 primer")
# Panel A shows the strict-testable parents only.
properties <- load_table("query_properties") %>% filter(strict_testable) %>% group_order()
stopifnot(nrow(properties) > 0)
# Conditions without a declared label are shown under their condition name, after the declared ones.
extra_conditions <- sort(setdiff(unique(properties$condition), names(condition_names)))
condition_names <- c(condition_names, setNames(extra_conditions, extra_conditions))
properties <- properties %>%
  mutate(condition_label = factor(unname(condition_names[condition]), levels = rev(unname(condition_names[names(condition_names) %in% condition]))))
length_counts <- properties %>% count(group, condition_label, length, name = "parents")
composition_counts <- properties %>% count(group, condition_label, gc_fraction, name = "parents")
property_count_limit <- max(c(length_counts$parents, composition_counts$parents))
property_size <- scale_size_continuous(limits = c(1, property_count_limit), range = c(1, 3),
  breaks = c(1, 100, 200), name = "Parents at value")
first_testable <- function(rule) min(testability$query_length[testability$endpoint == rule & testability$testability == "testable"])
strict_first <- first_testable("strict")
tolerant_first <- first_testable("tolerant")
stopifnot(all(testability$testability[testability$endpoint == "strict" & testability$query_length < strict_first] != "testable"),
          all(testability$testability[testability$endpoint == "tolerant" & testability$query_length < tolerant_first] != "testable"))
stopifnot(min(properties$length) >= strict_first)
# The dotted line marks the first tolerant-testable panel length when strict-testable parents lie below it.
tolerant_line <- any(properties$length < tolerant_first)
length_plot <- ggplot(length_counts, aes(length, condition_label, colour = group)) +
  (if (tolerant_line) geom_vline(xintercept = tolerant_first - 0.5, colour = grey_dark, linetype = "dotted", linewidth = lw(0.5))) +
  geom_point(aes(size = parents), alpha = 0.6, stroke = 0) + property_size +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous("Length (nt)", limits = c(0, 120), breaks = sort(unique(c(0, strict_first, if (tolerant_line) tolerant_first, 60, 90, 120)))) +
  labs(y = NULL) + theme(axis.text.y = element_text(size = dense_pt))
composition_plot <- ggplot(composition_counts, aes(gc_fraction * 100, condition_label, colour = group)) +
  geom_point(aes(size = parents), alpha = 0.6, stroke = 0) + property_size +
  scale_colour_manual(values = colours, guide = "none") + scale_x_continuous("GC (%)", limits = c(0, 100)) +
  labs(y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
properties_panel <- wrap_elements(full = ((length_plot | composition_plot) + plot_layout(widths = c(1.5, 0.75), guides = "collect")) &
  theme(legend.position = "bottom", plot.margin = margin(4.5, 1.5, 0, 1.5, "mm")))

# Independent benchmark. The by-length pair recall is reported in the text and exported for Table S9 (no panel).
benchmark_all <- load_table("benchmark_metrics")
benchmark <- benchmark_all %>% filter(database == "core_nt", population == "all_lengths") %>%
  mutate(arm_label = recode(search_space, local_fixture = "Local fixture", production_calibrated = "Calibrated E=10",
    permissive_diagnostic = "Diagnostic E=1e12")) %>%
  arrange(match(search_space, c("local_fixture", "production_calibrated", "permissive_diagnostic")), endpoint)
stopifnot(nrow(benchmark) == 6, all(benchmark_all$false_positive == 0),
  all(benchmark$positive_pairs + benchmark$negative_pairs == manifest$benchmarks$core_nt$pair_universe))
negative_counts <- benchmark %>% distinct(endpoint, negative_pairs)
stopifnot(nrow(negative_counts) == 2)
negative_note <- paste0("No truth-negative pairs accepted: ", paste0(negative_counts$endpoint, " 0/",
  scales::comma(negative_counts$negative_pairs, accuracy = 1), collapse = ", "), " (each arm).")
benchmark_length <- benchmark_all %>%
  filter(population == "length", search_space %in% c("production_calibrated", "permissive_diagnostic")) %>%
  mutate(arm_label = recode(search_space, production_calibrated = "E = 10", permissive_diagnostic = "E = 1e12"))
stopifnot(!anyNA(benchmark_length$recall), all(benchmark_length$false_positive == 0))
pair_recall <- export_csv(benchmark_length %>%
  arrange(match(database, c("core_nt", "core_nt_div10", "core_nt_x10")), match(search_space, c("production_calibrated", "permissive_diagnostic")), endpoint, length) %>%
  transmute(database, arm = search_space, arm_label, rule = endpoint, length_nt = length, truth_positive_pairs = positive_pairs,
    recovered_pairs = true_positive, missed_pairs = false_negative, pair_recall = recall), "tableS9_pair_recall_by_length.csv")
# At production scale an exact spike is recovered exactly when its exact E-value passes the production threshold.
strict_calibrated <- filter(pair_recall, database == "core_nt", arm == "production_calibrated", rule == "strict") %>%
  inner_join(testability %>% filter(endpoint == "strict") %>% select(length_nt = query_length, evalue_at_minimum_score), by = "length_nt")
stopifnot(all(strict_calibrated$pair_recall[strict_calibrated$evalue_at_minimum_score > production_evalue] == 0),
  all(strict_calibrated$pair_recall[strict_calibrated$evalue_at_minimum_score <= production_evalue] == 1),
  all(pair_recall$pair_recall[pair_recall$arm == "permissive_diagnostic" & pair_recall$rule == "strict"] == 1))
# Returned HSPs of the strict-testable reference sequences (figure_sources.alignment_quality).
quality <- load_table("quality_marginal") %>% group_order()
quality_hsps <- manifest$counts$strict_testable_reference_hsps
stopifnot(sum(load_table("quality_acceptance")$hsp_rows) >= quality_hsps)
quality_guides <- data.frame(measure = c("Identity", "Query span"), threshold = c(95, 90))
quality_plot <- ggplot(quality, aes(value, cumulative, colour = group)) +
  geom_step(linewidth = lw(0.75)) + geom_point(data = quality %>% group_by(group, measure) %>% filter(n() == 1) %>% ungroup(), size = 1.3) +
  geom_vline(data = quality_guides, aes(xintercept = threshold), inherit.aes = FALSE, linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  facet_wrap(~measure, nrow = 1, scales = "free_x") + scale_colour_manual(values = colours, name = NULL) +
  scale_x_continuous("Identity or query span (%)", limits = function(limits) c(if (min(limits) >= 80) 80 else 0, 100),
    breaks = function(limits) if (min(limits) >= 75) seq(80, 100, 5) else seq(0, 100, 25)) +
  scale_y_continuous("HSPs (cumulative %)", limits = c(0, 1), labels = function(x) round(100 * x)) +
  theme(legend.position = "right", strip.text = element_text(hjust = 0.5))
joint <- load_table("quality_joint") %>% group_order() %>%
  group_by(group, identity_bin, span_bin) %>% summarise(hsp_rows = sum(hsp_rows), .groups = "drop")
joint_identity_start <- min(80, floor(min(joint$identity_bin) / 10) * 10)
joint_plot <- ggplot(joint, aes(identity_bin + 0.5, span_bin + 0.5, fill = group, alpha = hsp_rows)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = 95, linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  facet_wrap(~group, nrow = 1) + scale_fill_manual(values = colours, guide = "none") +
  scale_alpha_continuous(transform = "log10", range = c(0.4, 1), name = "HSPs per bin",
    breaks = c(1, 100, 10000, 100000), labels = label_power_expr,
    guide = guide_legend(override.aes = list(fill = ink), nrow = 1)) +
  scale_x_continuous("Identity bin (%)", breaks = seq(joint_identity_start, 100, 5) + 0.5, labels = seq(joint_identity_start, 100, 5)) +
  scale_y_continuous("Span bin (%)", breaks = seq(0, 100, 25) + 0.5, labels = seq(0, 100, 25)) +
  coord_cartesian(xlim = c(joint_identity_start, 101), ylim = c(0, 101), expand = FALSE) +
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0.5))
figureS1_legend_notes <- c(
  paste0("A, Length and GC content of the strict-testable parents per condition. Strict matching is testable from ", strict_first, " nt.",
         if (tolerant_line) paste0(" Dotted line: tolerant matching is testable from ", tolerant_first, " nt (first testable panel length).") else ""),
  paste0("B, Cumulative distributions of identity and query span over all ", comma(quality_hsps), " returned HSPs of strict-testable reference ",
         "queries, rejected HSPs included. Dashed lines: 95% identity and 90% query span."),
  "C, The same HSPs in 1%-wide bins labelled by their lower bound (100 denotes exactly 100%). Colour intensity: HSPs per bin on a logarithmic scale.")
figureS1_heights <- c(74, 46, 64)
figure_s1 <- assemble(properties_panel, quality_plot, joint_plot, heights = figureS1_heights)
save_figure("figureS1", figure_s1, sum(figureS1_heights) + 4, c("Query properties", "Marginal quality", "Joint quality"),
  c("query_properties.tsv", "benchmark_metrics.tsv", "quality_marginal.tsv", "quality_joint.tsv", "quality_acceptance.tsv"),
  "Independent current-version fixture QA",
  "tableS9_pair_recall_by_length.csv", legend_notes = figureS1_legend_notes)

# ---------------------------------------------------------------------------------------------
# Target-limit comparison (500 versus 2,000 targets and production) and cap-reaching searches. No panel: the values are
# exported for Table S10.
# ---------------------------------------------------------------------------------------------
parent_names <- load_table("logical_queries") %>% filter(role == "real") %>% distinct(parent_id, parent_logical_id = logical_id)
parent_display <- parent_names %>% group_by(parent_id) %>%
  summarise(reference_parent = paste(parent_label(sort(parent_logical_id)), collapse = " / "), .groups = "drop")
target_summary <- load_table("target_limit_summary")
if (nrow(target_summary)) target_summary <- target_summary %>% arrange(match(group, c(groups, "Global union")), endpoint, population, gate)
invisible(export_csv(target_summary, "tableS10_target_limit_summary.csv"))
cap_membership <- load_table("cap_memberships")
caps <- load_table("cap_sequences") %>%
  mutate(display_label = ifelse(is.na(aliases), paste0("Draw ", substr(sequence_id, 5, 12)), parent_label(aliases)),
    role_label = recode(roles, real = "Reference", control = "Control", .default = "Reference and control")) %>%
  left_join(select(cap_membership, sequence_id, membership_parent_id = parent_id, draw_index, membership_role = role), by = "sequence_id") %>%
  left_join(parent_display, by = c("membership_parent_id" = "parent_id")) %>%
  mutate(row_label = paste0(groups, " ", reference_parent, ", ", ifelse(membership_role == "real", "reference", paste0("draw ", draw_index)),
                            ", ", length, " nt"))
stopifnot(!anyNA(caps$reference_parent), all(caps$targets >= target_limit))
invisible(export_csv(data.frame(
  display_label = caps$display_label, group = caps$groups, role = caps$role_label, length_nt = caps$length,
  raw_hsp_rows = caps$hsp_rows, accepted_strict_hsps = caps$strict_retained, accepted_tolerant_hsps = caps$tolerant_retained,
  search_sequence_id = caps$sequence_id, full_aliases = caps$aliases, logical_memberships = caps$logical_memberships,
  reference_parent = caps$reference_parent, control_draw_index = ifelse(caps$membership_role == "control", caps$draw_index, NA),
  search_label = caps$row_label, target_cap = rep(target_limit, nrow(caps))), "tableS10_cap_reaching_searches.csv"))
target_parents <- load_table("target_limit_parents")
if (nrow(target_parents)) {
  target_parents <- target_parents %>% left_join(parent_display, by = "parent_id") %>% relocate(reference_parent) %>%
    arrange(match(sub(";.*$", "", groups), groups), endpoint, population, reference_parent)
  stopifnot(!anyNA(target_parents$reference_parent))
}
invisible(export_csv(target_parents, "tableS10_target_limit_parents.csv"))

# ---------------------------------------------------------------------------------------------
# Figure S2. Control quality (Table S6) and Monte Carlo precision
# ---------------------------------------------------------------------------------------------
qc <- load_table("control_qc_summary")
qc_parents <- load_table("control_parent_qc") %>% group_order() %>%
  mutate(testable_identical = parent_identical_canonical * strict_testable, testable_repeated = repeated_canonical * strict_testable)
draw_qc <- load_table("control_draw_qc") %>% select(parent_id, invariants_passed) %>%
  left_join(select(qc_parents, parent_id, group), by = "parent_id")
stopifnot(!anyNA(draw_qc$group))
quality_by_group <- bind_rows(lapply(groups[groups %in% qc_parents$group], function(name) {
  parents <- filter(qc_parents, group == name)
  draw_rows <- filter(draw_qc, group == name)
  cbind(group = name, summarise(parents, parents = n(), distinct_sequences_median = median(distinct_canonical),
    distinct_sequences_min = min(distinct_canonical), distinct_sequences_max = max(distinct_canonical),
    parent_identical_canonical = sum(parent_identical_canonical), repeated_canonical = sum(repeated_canonical),
    parent_identical_exact = sum(parent_identical_exact), repeated_exact = sum(repeated_exact),
    draws_passing_invariants = sum(draw_rows$invariants_passed), draws_total = nrow(draw_rows),
    testable_parents = sum(strict_testable), testable_parent_identical_canonical = sum(testable_identical),
    testable_repeated_canonical = sum(testable_repeated)))
}))
quality_total <- cbind(group = "All groups", summarise(qc_parents, parents = n(), distinct_sequences_median = median(distinct_canonical),
  distinct_sequences_min = min(distinct_canonical), distinct_sequences_max = max(distinct_canonical),
  parent_identical_canonical = sum(parent_identical_canonical), repeated_canonical = sum(repeated_canonical),
  parent_identical_exact = sum(parent_identical_exact), repeated_exact = sum(repeated_exact),
  draws_passing_invariants = sum(draw_qc$invariants_passed), draws_total = nrow(draw_qc),
  testable_parents = sum(strict_testable), testable_parent_identical_canonical = sum(testable_identical),
  testable_repeated_canonical = sum(testable_repeated)))
control_quality <- export_csv(rbind(quality_by_group, quality_total), "tableS6_control_quality.csv")
total_row <- nrow(control_quality)
stopifnot(control_quality$group[total_row] == "All groups",
  control_quality$parent_identical_canonical[total_row] == qc$parent_identical_canonical, control_quality$repeated_canonical[total_row] == qc$repeated_canonical,
  control_quality$parent_identical_exact[total_row] == qc$parent_identical_exact, control_quality$repeated_exact[total_row] == qc$repeated_exact,
  control_quality$draws_passing_invariants[total_row] == qc$selected_draws, control_quality$draws_total[total_row] == qc$selected_draws,
  control_quality$testable_parents[total_row] == qc$testable_parents,
  control_quality$testable_parent_identical_canonical[total_row] == qc$testable_parent_identical_canonical,
  control_quality$testable_repeated_canonical[total_row] == qc$testable_repeated_canonical)

# Control draws of the strict-testable parents (legend note; all parents: Table S6).
testable_control_note <- paste0(scales::comma(qc$testable_parents), " strict-testable parents have ", scales::comma(qc$testable_draws),
  " control draws. ", if (qc$testable_parent_identical_canonical == 0 && qc$testable_repeated_canonical == 0)
    "None of these draws repeats within its parent or equals its parent." else
    paste0(scales::comma(qc$testable_parent_identical_canonical), " of these draws equal their parent and ",
           scales::comma(qc$testable_repeated_canonical), " repeat within their parent."),
  " Reverse complements count as identical.")

precision <- load_table("precision")
# A: fixed-prefix estimates with Monte Carlo intervals, detection (linear, %) and breadth (pseudo-logarithmic) in separate rows.
precision_selected <- precision %>% filter(condition == "all_conditions", population == "all_records", group != "Global union") %>%
  group_order() %>% mutate(endpoint = factor(endpoint, levels = c("strict", "tolerant")),
    measure = factor(measure, levels = c("detection", "breadth")),
    scale = ifelse(measure == "detection", 100, 1),
    display_estimate = estimate * scale, display_lower = lower * scale, display_upper = upper * scale,
    marker_fill = ifelse(status == "conditional_bound", unname(colours[as.character(group)]), "white"),
    rule_label = factor(recode(as.character(endpoint), strict = "Strict", tolerant = "Tolerant"), levels = c("Strict", "Tolerant")))
fixed_precision <- load_table("precision_fixed_cohort") %>% filter(population == "all_records", group != "Global union") %>%
  mutate(group = factor(group, levels = groups), measure = factor(measure, levels = c("detection", "breadth")),
    display_estimate = ifelse(measure == "detection", estimate * 100, estimate),
    rule_label = factor(recode(endpoint, strict = "Strict", tolerant = "Tolerant"), levels = c("Strict", "Tolerant")))
precision_row <- function(measure_name, y_title, transform = "identity") {
  ggplot(filter(precision_selected, measure == measure_name),
         aes(draws_per_parent, display_estimate, colour = group, shape = rule_label, group = rule_label)) +
    geom_linerange(aes(ymin = display_lower, ymax = display_upper), linewidth = lw(0.75), position = position_dodge(width = 1.4)) +
    geom_point(aes(fill = I(marker_fill)), position = position_dodge(width = 1.4), size = 1.5) +
    geom_line(data = filter(fixed_precision, measure == measure_name), aes(draws_per_parent, display_estimate, colour = group,
      linetype = rule_label, group = rule_label), inherit.aes = FALSE, linewidth = lw(0.5)) +
    facet_wrap(~group, nrow = 1, scales = "free_y") +
    scale_x_continuous("Draws per parent", breaks = draw_prefixes) +
    scale_y_continuous(y_title, transform = transform, breaks = scales::breaks_pretty(n = 3)) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_shape_manual(values = c(Strict = 21, Tolerant = 24), name = NULL) +
    scale_linetype_manual(values = c(Strict = "dotted", Tolerant = "dashed"), name = NULL) +
    guides(shape = guide_legend(override.aes = list(colour = ink, fill = ink)), linetype = guide_legend()) +
    theme(strip.text = element_text(hjust = 0.5))
}
precision_detection <- precision_row("detection", "Control detection (%)") + theme(axis.title.x = element_blank())
precision_breadth <- precision_row("breadth", "Control breadth\n(accessions)", pseudo_log) + theme(strip.text = element_blank())
precision_plot <- wrap_elements(full = (precision_detection / precision_breadth) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.margin = margin(4.5, 1.5, 0.5, 1.5, "mm")))
precision_final <- precision %>% filter(draws_per_parent == final_prefix, condition == "all_conditions") %>%
  mutate(group = factor(group, levels = rev(c(groups, "Global union"))),
    endpoint = factor(endpoint, levels = c("strict", "tolerant")),
    population = recode(population, all_records = "All records", linear_genomic = "Linear-genomic"),
    cell = factor(paste(endpoint, measure, sep = " / "), levels = c("strict / detection", "strict / breadth", "tolerant / detection", "tolerant / breadth")),
    status_label = case_when(precision_met %in% "yes" ~ "Met", precision_met %in% "no" ~ "Missed",
      status == "diagnostic_resolved_subset_no_coverage_claim" ~ "Subset", TRUE ~ "No cohort"),
    ratio = half_width / target_half_width)
precision_counts <- precision %>% filter(draws_per_parent == final_prefix)
precision_count_label <- function(data) {
  counts <- c(met = sum(data$precision_met %in% "yes"), missed = sum(data$precision_met %in% "no"),
    subset = sum(is.na(data$precision_met) & data$status == "diagnostic_resolved_subset_no_coverage_claim"),
    empty = sum(is.na(data$precision_met) & data$status == "no_resolved_matched_population"))
  stopifnot(sum(counts) == nrow(data))
  parts <- c(paste(counts[["met"]], "targets met"), paste(counts[["missed"]], "missed"), paste(counts[["subset"]], "resolved-subset diagnostics"),
    if (counts[["empty"]]) paste(counts[["empty"]], "without a resolved cohort"))
  paste0(nrow(data), " rows: ", serial(parts), ".")
}
precision_tiles <- export_csv(precision_final %>% transmute(population, group = as.character(group), rule = as.character(endpoint), measure,
  half_width, target_half_width, half_width_over_target = ratio, status = status_label) %>%
  arrange(population, match(group, c(groups, "Global union")), rule, measure), "figS4B_precision_tiles.csv")
grey_ramp <- c("#F7F7F7", "#D9D9D9", "#BDBDBD", "#969696", "#636363", "#252525")
precision_final <- precision_final %>% mutate(text_colour = ifelse(log10(ratio) > 1.2, "white", "black"),
  tile_label = paste0(status_label, "\n", vapply(signif(ratio, 3), format, character(1), scientific = FALSE, trim = TRUE)))
precision_status <- ggplot(precision_final, aes(cell, group, fill = ratio)) +
  geom_tile(colour = "white", linewidth = lw(1)) + geom_text(aes(label = tile_label, colour = I(text_colour)), size = fs(dense_pt), lineheight = 0.95) +
  facet_wrap(~population, nrow = 1) +
  scale_fill_gradientn(colours = grey_ramp, transform = "log10", name = "Half-width / target",
    breaks = c(1, 3, 10, 30, 100, 300), labels = c("1", "3", "10", "30", "100", "300"),
    guide = guide_colourbar(barwidth = unit(45, "mm"), barheight = unit(2.5, "mm"))) +
  scale_x_discrete(labels = c("strict / breadth" = "Strict\nbreadth", "strict / detection" = "Strict\ndetection",
    "tolerant / breadth" = "Tolerant\nbreadth", "tolerant / detection" = "Tolerant\ndetection")) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom", strip.text = element_text(hjust = 0.5))
figureS2_legend_notes <- c(testable_control_note,
  paste0("A, Control means at each declared draw prefix with their Monte Carlo intervals. Filled markers: conditional bound. Open markers: ",
         "resolved-subset diagnostic. Lines: the final-parent cohort followed back through the prefixes (dotted: strict, dashed: tolerant). ",
         "Breadth uses a pseudo-logarithmic axis. Y scales differ between panels."),
  paste0("B, Final precision targets. ", precision_count_label(filter(precision_counts, condition == "all_conditions")),
         " Met: half-width \u2264 target. Tile text: status and half-width / target. Grey value: the same ratio on a logarithmic scale. Values: Table S6."))
figureS2_heights <- c(92, 82)
figure_s3 <- assemble(precision_plot, precision_status, heights = figureS2_heights)
save_figure("figureS2", figure_s3, sum(figureS2_heights) + 4, c("Prefix precision", "Precision targets"),
  c("control_draw_qc.tsv", "control_parent_qc.tsv", "control_qc_summary.tsv", "precision.tsv", "precision_fixed_cohort.tsv"),
  paste0("All final-stage ", precision_count_label(precision_counts), " Intervals describe conditional Monte Carlo error."),
  c("tableS6_control_quality.csv", "figS4B_precision_tiles.csv"), legend_notes = figureS2_legend_notes)

# ---------------------------------------------------------------------------------------------
# Figure S3. Positional robustness
# ---------------------------------------------------------------------------------------------
dense_position <- filter(position, group %in% c("ONT", "Illumina"))
sparse_position <- filter(position, group %in% c("BGI/DNBSEQ", "PacBio"))
all_position_totals <- position %>% group_by(group) %>% summarise(events = n(), accessions = n_distinct(accession_version), .groups = "drop")
all_position_labels <- setNames(paste0(all_position_totals$group, "\n", comma(all_position_totals$events),
  ifelse(all_position_totals$events == 1, " event\n", " events\n"), comma(all_position_totals$accessions),
  ifelse(all_position_totals$accessions == 1, " accession", " accessions")), as.character(all_position_totals$group))
position_distance_breaks <- c(0, 10^seq(2, ceiling(log10(max(100, position$distance))), 2))
position_scatter <- ggplot(position, aes(distance, record_length)) +
  geom_hex(data = dense_position, aes(fill = group, alpha = after_stat(count)), bins = 22) +
  geom_point(data = sparse_position, aes(colour = group), size = 1.3, alpha = 0.8) +
  geom_vline(xintercept = main_window, linetype = "dashed", linewidth = lw(0.5), colour = grey_dark) +
  geom_hline(yintercept = main_minimum, linetype = "dashed", linewidth = lw(0.5), colour = grey_dark) +
  facet_wrap(~group, nrow = 1, labeller = labeller(group = all_position_labels)) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_alpha_continuous(transform = "log10", range = c(0.15, 1), breaks = c(1, 10, 100, 1000), labels = label_power_expr,
    name = "Events per bin", guide = guide_legend(override.aes = list(fill = ink))) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous("Distance to nearest record end (nt)", transform = pseudo_log, breaks = position_distance_breaks, labels = label_power_expr) +
  scale_y_log10("Record length (nt)", breaks = 10^(2:7), labels = label_power_expr) +
  theme(strip.text = element_text(hjust = 0.5), legend.position = "right")
density <- load_table("density") %>% group_order()
density_summary <- density %>% group_by(group) %>% summarise(accessions = n(), mean = mean(strict_orientation_deduplicated_events),
  median = median(strict_orientation_deduplicated_events), .groups = "drop")
single_event_share <- density %>% group_by(group) %>%
  summarise(share = mean(strict_orientation_deduplicated_events == 1), .groups = "drop")
# Nested queries (Illumina embedded adapters; ONT 24-nt barcode inside K40/K46) count once per query (nesting report), so the
# single-event share can fall below 50%; the y-axis floor follows the data.
ecdf_floor <- floor(10 * min(single_event_share$share)) / 10
density_plot <- ggplot(density, aes(strict_orientation_deduplicated_events, colour = group)) +
  stat_ecdf(geom = "step", pad = FALSE, linewidth = lw(0.75)) +
  scale_x_log10("Match events per accession", labels = comma) +
  scale_y_continuous("Accessions\n(cumulative %)", breaks = seq(ecdf_floor, 1, 0.1), labels = function(x) round(100 * x)) +
  coord_cartesian(ylim = c(ecdf_floor, 1)) +
  scale_colour_manual(values = colours, labels = setNames(paste0(as.character(density_summary$group), " (n = ",
    comma(density_summary$accessions), ")"), as.character(density_summary$group)), name = NULL) +
  theme(legend.position = "right")
# C: every window x minimum-length cell for every unit variant (no cell selection); the 25-nt window is diagnostic.
grid_export <- export_csv(load_table("position_sensitivity") %>%
  arrange(match(variant, names(variant_labels)), match(group, c(groups, "Global union")), minimum_length, window) %>%
  transmute(variant, unit_rule, group, minimum_length_nt = minimum_length, window_nt = window, events, terminal, expected, ratio,
            diagnostic_window = window == 25, units_without_project), "figS3C_positional_sensitivity.csv")
grid <- grid_export %>% filter(group != "Global union") %>%
  mutate(group = factor(group, levels = groups), variant = factor(unname(variant_labels[variant]), levels = unname(variant_labels)),
    cell = paste0(minimum_length_nt / 1000, " kb / ", window_nt, " nt", ifelse(diagnostic_window, " (diagnostic)", "")),
    x = ifelse(is.na(ratio), 0, ratio), state = ifelse(is.na(ratio), "No eligible events", "Ratio"))
cell_levels <- grid %>% distinct(minimum_length_nt, window_nt, cell) %>% arrange(desc(minimum_length_nt), desc(window_nt)) %>% pull(cell)
grid <- grid %>% mutate(cell = factor(cell, levels = cell_levels))
ratio_breaks <- c(0, 10, 100, 1000)
sensitivity_plot <- ggplot(grid, aes(x, cell, colour = group, shape = variant)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = grey_dark, linewidth = lw(0.5)) +
  geom_point(data = filter(grid, state == "Ratio"), size = 1.2, position = position_dodge(width = 0.7)) +
  geom_point(data = filter(grid, state != "Ratio"), shape = 4, colour = grey_mid, size = 1.2, position = position_dodge(width = 0.7)) +
  facet_wrap(~group, nrow = 1) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_shape_manual(values = variant_shapes, name = NULL, guide = guide_legend(override.aes = list(colour = ink))) +
  scale_x_continuous("Observed / expected terminal events", transform = pseudo_log, breaks = ratio_breaks,
                     labels = function(x) comma(x), limits = c(0, max(1500, max(grid$x, na.rm = TRUE) * 1.2))) +
  labs(y = "Minimum record length / window") +
  theme(axis.text.y = element_text(size = dense_pt), legend.position = "bottom", strip.text = element_text(hjust = 0.5))
figureS3_legend_notes <- c(
  paste0("A, Distance to the nearest record end against record length for strict matches of testable queries in eligible linear-genomic ",
         "records of all lengths. Both axes are logarithmic (distance pseudo-logarithmic). Hexagons: events per bin. Points: groups with few ",
         "events. Dashed lines: records of ", comma(main_minimum / 1000), " kb and the distance W = ", comma(main_window), " nt."),
  paste0("B, Match events per hit-positive accession (strict matches of testable queries, all records), orientation-deduplicated per query. ",
         "The y-axis starts at ", round(100 * ecdf_floor), "%. Nested queries (for example K40 and K46 contexts) add events at one site."),
  paste0("C, Observed / expected terminal events for every minimum-length and window cell and every unit variant, on a pseudo-logarithmic axis. ",
         "Dashed line: ratio 1. Grey crosses: no eligible events. The 25-nt window is diagnostic. A record without a BioProject is its own unit."))
figureS3_heights <- c(54, 40, 92)
figure_s4 <- assemble(position_scatter, density_plot, sensitivity_plot, heights = figureS3_heights)
save_figure("figureS3", figure_s4, sum(figureS3_heights) + 4, c("Distance versus length", "Conditional HSP density", "Positional sensitivity by unit"),
  c("position_events.tsv", "density.tsv", "position_sensitivity.tsv"),
  "Strict, a priori testable events of every condition. Sparse events are drawn as exact points",
  "figS3C_positional_sensitivity.csv", legend_notes = figureS3_legend_notes)

# ---------------------------------------------------------------------------------------------
# Supplementary Tables re-exported from the same source (Word tables are regenerated from these)
# ---------------------------------------------------------------------------------------------
precision_export <- function(measure_name) {
  precision %>% filter(draws_per_parent == final_prefix, condition == "all_conditions", measure == measure_name) %>%
    mutate(group = factor(group, levels = c(groups, "Global union"))) %>% arrange(group, endpoint, population) %>%
    transmute(group = as.character(group), rule = endpoint, population, planned_parents, matched_parents,
      mean = estimate, range_lower = lower, range_upper = upper, half_width, target_half_width,
      control_draws, unique_seeds, target_limited_parents, precision_met, status)
}
invisible(export_csv(precision_export("detection"), "tableS6_detection_precision.csv"))
invisible(export_csv(precision_export("breadth"), "tableS6_breadth_precision.csv"))
invisible(export_csv(data.frame(arm = benchmark$search_space, arm_label = benchmark$arm_label, rule = benchmark$endpoint,
  truth_positive = benchmark$positive_pairs, recovered = benchmark$true_positive, missed = benchmark$false_negative,
  truth_negative = benchmark$negative_pairs, accepted_negative = benchmark$false_positive, recall = benchmark$recall),
  "tableS9_independent_benchmark.csv"))
invisible(export_csv(spike_summary %>% transmute(database, length_nt, canonical_parents = parents, exact_spikes = planned_spikes,
    calibrated_e10_recovered = calibrated_recovered, diagnostic_e1e12_recovered = diagnostic_recovered,
    exact_evalue_at_database_scale = exact_evalue, effective_search_space), "tableS9_exact_spike_detectability.csv"))

if (length(saved)) {
  suffix <- if (selection == "all") "" else paste0("_", selection)
  # PNG files are byte-deterministic; cairo writes the current time into every PDF, so PDFs are listed by name and their
  # content digest (creation time excluded) is added by build_figures.sh.
  filenames <- paste0("paper3_", names(saved), ".png")
  csv_files <- sort(list.files(tables_directory, pattern = "\\.csv$"))
  result <- list(contract = "publication_figure_artwork_v1",
    source_manifest_sha256 = digest(file.path(input, "manifest.json"), algo = "sha256", file = TRUE),
    plot_script_sha256 = digest(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1]), algo = "sha256", file = TRUE),
    style_sha256 = digest(style_file, algo = "sha256", file = TRUE), profile = FIG_PROFILE, fonts = c(sans = FONT_SANS, mono = FONT_MONO),
    figures = saved,
    R_version = R.version.string,
    package_versions = setNames(lapply(c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales", "jsonlite", "digest", "ragg", "hexbin", "ggbeeswarm", "systemfonts"),
      function(package) as.character(packageVersion(package))), c("ggplot2", "dplyr", "tidyr", "readr", "patchwork", "scales", "jsonlite", "digest", "ragg", "hexbin", "ggbeeswarm", "systemfonts")),
    files = setNames(lapply(filenames, function(name) digest(file.path(output, name), algo = "sha256", file = TRUE)), filenames),
    pdf_files = paste0("paper3_", names(saved), ".pdf"),
    tables = setNames(lapply(csv_files, function(name) digest(file.path(tables_directory, name), algo = "sha256", file = TRUE)), csv_files))
  write_json(result, file.path(output, paste0("paper3_figure_manifest", suffix, ".json")), auto_unbox = TRUE, pretty = TRUE)
}
