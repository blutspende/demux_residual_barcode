suppressPackageStartupMessages({
  library(grid)
  library(readr)
  library(jsonlite)
  library(digest)
})
options(warn = 2)
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3)
input <- normalizePath(arguments[1])
style_file <- normalizePath(arguments[2])
output <- arguments[3]
manifest <- fromJSON(file.path(input, "manifest.json"))
# Record cases (figureS4, manuscript Figure S2) from one figure-source directory (contract publication_figure_sources_v1): the
# cases selected by the record_cases stage with their sequence-level verification. House style of paper3_figure_style.R: Arial
# 7 pt, Courier New for sequences, panel tags A-C, no title or footer (the caption text is recorded as legend_notes), the match
# on a light tint of the platform colour of the query and barcode or index segments on a darker tint, UK spelling.
stopifnot(manifest$contract == "publication_figure_sources_v1")
for (name in names(manifest$files)) {
  stopifnot(!grepl("(^/|(^|/)\\.\\.(/|$))", name),
    digest(file.path(input, name), algo = "sha256", file = TRUE) == manifest$files[[name]])
}
source(style_file)
cases <- read_tsv(file.path(input, "record_cases.tsv"), show_col_types = FALSE, na = c("", "not_estimated"))
cases <- cases[!is.na(cases$accession_version), ]
# Sequence-level verification of every case; species names (Latin binomials) are set in italics.
verification <- read_tsv(file.path(input, "record_case_verification.tsv"), show_col_types = FALSE, na = c("", "not_estimated"))
stopifnot(nrow(cases) <= 3, !anyDuplicated(cases$accession_version))
if (nrow(cases)) {
  stopifnot(identical(verification$accession_version, cases$accession_version), all(verification$displayed_window_equals_record),
    all(verification$record_window_equals_hsp), all(verification$hsp_full_length_exact),
    all(verification$packet_identical | is.na(verification$packet_identical)))
  cases$organism_italic <- verification$organism_italic
  cases$packet_identical <- verification$packet_identical
}
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
width_mm <- FIG$width
case_pitch <- 61
case_start <- 1
height_mm <- max(20, case_start + case_pitch * max(1, nrow(cases)))
stopifnot(height_mm <= FIG$max_h)
query_names <- c(ILLUMINA_TRUSEQ_R2_adapter = "TruSeq R2 adapter", RAPID_adapter_forward = "Rapid adapter")
display_query <- function(id) if (id %in% names(query_names)) unname(query_names[id]) else gsub("_", " ", sub("_forward$", "", id))
platform_colour <- function(group) {
  colour <- pal_platform[group]
  stopifnot(!anyNA(colour))
  unname(colour)
}
tint <- function(colour, alpha) grDevices::adjustcolor(colour, alpha.f = alpha)

write_text <- function(text, xpos, top, size = FIG$base, bold = FALSE, colour = ink, family = FONT_SANS, just = "left") {
  grid.text(text, x = unit(xpos, "mm"), y = unit(height_mm - top, "mm"), just = c(just, "centre"),
    gp = gpar(fontsize = size, fontface = if (bold) "bold" else "plain", col = colour, fontfamily = family))
}
text_width <- function(text, size, face = "plain", family = FONT_SANS) {
  convertWidth(grobWidth(textGrob(text, gp = gpar(fontsize = size, fontface = face, fontfamily = family))), "mm", valueOnly = TRUE)
}

draw_case <- function(row, index, top) {
  category_label <- switch(row$category, terminal = "Terminal match", synthetic_or_vector = "Synthetic or vector record",
    internal_or_ambiguous = "Internal match", row$category)
  write_text(LETTERS[index], 0, top + 2, FIG$tag, TRUE)
  write_text(paste0(category_label, ", ", row$accession_version), 5, top + 2, FIG$base, TRUE)
  roles <- strsplit(row$display_segment_roles, ";", fixed = TRUE)[[1]]
  stopifnot(length(roles) == nchar(row$qseq), row$qseq == row$sseq)
  emphasis <- roles %in% c("barcode", "index")
  strand <- if (row$display_orientation == "-") "minus" else "plus"
  write_text(paste0("Query: ", display_query(row$query_id), " (", row$query_id, "), ", nchar(row$qseq), "-nt strict exact match on the ",
                    strand, " strand"), 5, top + 6.2)
  title_lines <- strwrap(row$title, width = 110)
  for (line_index in seq_along(title_lines)) {
    line <- title_lines[line_index]
    line_top <- top + 9.8 + 3.2 * (line_index - 1)
    if (line_index == 1 && isTRUE(row$organism_italic) && startsWith(line, row$organism)) {
      gp <- gpar(fontsize = FIG$base, fontface = "italic", fontfamily = FONT_SANS)
      grid.text(row$organism, x = unit(5, "mm"), y = unit(height_mm - line_top, "mm"), just = c("left", "centre"), gp = gp)
      write_text(substring(line, nchar(row$organism) + 1), 5 + text_width(row$organism, FIG$base, "italic"), line_top)
    } else {
      write_text(line, 5, line_top)
    }
  }
  metadata <- paste0(row$molecule, ", ", row$topology, ", ", row$record_type, ". Origin: ", row$origin, ". Linear-genomic set: ",
                     row$linear_eligibility, ". Taxid ", row$taxids, ".")
  write_text(metadata, 5, top + 16.6, FIG$dense)
  coordinate_label <- paste0("Record positions ", format(row$lower, big.mark = ",", scientific = FALSE), "\u2013",
    format(row$upper, big.mark = ",", scientific = FALSE), " of ", format(row$record_length, big.mark = ",", scientific = FALSE), " nt, ",
    if (row$topology == "circular") "circular record (position 1 is the stored origin)" else
      paste(format(row$distance_to_end, big.mark = ",", scientific = FALSE), "nt from the nearest record end"))
  write_text(coordinate_label, 5, top + 20.2, FIG$dense)
  track_top <- top + 25
  grid.lines(x = unit(c(22, width_mm - 2), "mm"), y = unit(rep(height_mm - track_top, 2), "mm"), gp = gpar(col = grey_mid, lwd = lw(1) * .pt))
  match_x <- 22 + (width_mm - 24) * ((row$lower + row$upper) / 2 - 1) / max(1, row$record_length - 1)
  grid.lines(x = unit(rep(match_x, 2), "mm"), y = unit(height_mm - track_top - c(-1.8, 1.8), "mm"),
             gp = gpar(col = platform_colour(row$group), lwd = lw(2.5) * .pt))
  write_text("1", 22, track_top + 3, FIG$dense, just = "centre")
  write_text(format(row$record_length, big.mark = ",", scientific = FALSE), width_mm - 2, track_top + 3, FIG$dense, just = "right")
  left <- if (is.na(row$display_left_flank)) "" else row$display_left_flank
  right <- if (is.na(row$display_right_flank)) "" else row$display_right_flank
  # Two alignment lines of 60 bases: longer queries show correspondingly shorter flanks.
  flank <- max(0, floor((120 - nchar(row$qseq)) / 2))
  left <- if (flank) substring(left, max(1, nchar(left) - flank + 1)) else ""
  right <- if (flank) substring(right, 1, flank) else ""
  subject <- strsplit(paste0(left, row$sseq, right), "", fixed = TRUE)[[1]]
  query <- c(rep(".", nchar(left)), strsplit(row$qseq, "", fixed = TRUE)[[1]], rep(".", nchar(right)))
  in_match <- c(rep(FALSE, nchar(left)), rep(TRUE, nchar(row$qseq)), rep(FALSE, nchar(right)))
  strong <- c(rep(FALSE, nchar(left)), emphasis, rep(FALSE, nchar(right)))
  stopifnot(length(subject) == length(query), length(subject) <= 120)
  marker <- platform_colour(row$group)
  pitch <- 2.4
  record_label <- if (row$display_orientation == "-") "Record (\u2212)" else "Record (+)"
  for (offset in seq(1, length(subject), by = 60)) {
    indices <- offset:min(length(subject), offset + 59)
    line_index <- (offset - 1) / 60
    subject_top <- top + 32 + 10 * line_index
    write_text(record_label, 0, subject_top, FIG$dense)
    write_text("Query", 0, subject_top + 4, FIG$dense)
    for (position in seq_along(indices)) {
      base <- indices[position]
      xpos <- 22 + pitch * (position - 0.5)
      if (in_match[base]) {
        grid.rect(x = unit(xpos, "mm"), y = unit(height_mm - subject_top - 2, "mm"), width = unit(pitch, "mm"), height = unit(8, "mm"),
          gp = gpar(fill = tint(marker, if (strong[base]) 0.45 else 0.2), col = NA))
      }
      write_text(subject[base], xpos, subject_top, FIG$dense, bold = in_match[base], colour = if (in_match[base]) ink else grey_text,
                 family = FONT_MONO, just = "centre")
      write_text(query[base], xpos, subject_top + 4, FIG$dense, bold = in_match[base], colour = if (in_match[base]) ink else grey_light,
                 family = FONT_MONO, just = "centre")
    }
  }
  # Verification result of this case (every check is asserted above from record_case_verification.tsv).
  lines <- paste("Verified: exact full-length HSP of the searched query. The record sequence at these coordinates equals the HSP.",
    if (isTRUE(row$packet_identical)) "The archived record sequence is identical." else "")
  write_text(trimws(lines), 0, top + 55, FIG$dense)
}

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  if (nrow(cases)) {
    for (index in seq_len(nrow(cases))) draw_case(cases[index, ], index, case_start + case_pitch * (index - 1))
  } else {
    write_text("No record case was selected. The record-case selection tables list the category status.", 0, 6)
  }
}

for (extension in c("pdf", "png", if (FIG$tiff) "tiff")) {
  path <- file.path(output, paste0("paper3_figureS4.", extension))
  stopifnot(!file.exists(path))
  if (extension == "pdf") grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4)
  else if (extension == "png") ragg::agg_png(path, width = width_mm, height = height_mm, units = "mm", res = 300)
  else ragg::agg_tiff(path, width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw")
  draw_figure()
  grDevices::dev.off()
}
legend_notes <- c(
  "Deterministically selected strict matches of testable queries with exact-version local sequence verification.",
  paste0("Record (\u2212) or (+): the record strand that carries the query in its searched orientation, with up to 50 nt of flank, clipped at ",
         "stored record boundaries. Query: the searched query; dots mark positions outside the match. The match is set on a light tint of the ",
         "platform colour of the query, and the barcode or index segment on a darker tint. Track: the record drawn to scale with the match ",
         "marker. Selection rules and complete case evidence: Table S8."))
files <- "paper3_figureS4.png"
report <- list(contract = "publication_record_case_artwork_v1",
  source_manifest_sha256 = digest(file.path(input, "manifest.json"), algo = "sha256", file = TRUE),
  figures = list(figureS4 = list(width_mm = width_mm, height_mm = height_mm, selected_cases = nrow(cases), profile = FIG_PROFILE,
    source_tables = c("record_cases.tsv", "record_case_verification.tsv"), legend_notes = legend_notes)),
  code_sha256 = digest(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1]), algo = "sha256", file = TRUE),
  style_sha256 = digest(style_file, algo = "sha256", file = TRUE),
  files = setNames(lapply(files, function(name) digest(file.path(output, name), algo = "sha256", file = TRUE)), files),
  pdf_files = "paper3_figureS4.pdf")
write_json(report, file.path(output, "paper3_case_figure_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
