# paper3_figure_style.R: the house style of every Paper 3 figure (main, supplementary, record cases).
#
# One profile switch sets width, maximum height, font sizes and panel tags (environment variable PAPER3_FIG_PROFILE):
#   default   170 mm wide, at most 200 mm tall, 7-pt text, 8-pt bold tags A, B, C (preprint, Cell Press, OUP, BMC width)
#   nature    183 mm wide, at most 170 mm tall, 7-pt text, 8-pt bold tags a, b, c
#   plos_bmc  170 mm wide, at most 200 mm tall, 8-pt text, 10-pt bold tags A, B, C (also writes 600-dpi LZW TIFF files)
# Text is Arial (sequences Courier New), black, axis titles and tick labels 20% above their previous sizes, units in parentheses, UK spelling.
# Hue is reserved for platforms. Evidence states, rules and thresholds use grey value, fill, outline style and shape.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

FIG_PROFILES <- list(
  default = list(width = 170, max_h = 200, base = 7, dense = 6.5, tag = 8, tags = "A", tiff = FALSE),
  nature = list(width = 183, max_h = 170, base = 7, dense = 6.5, tag = 8, tags = "a", tiff = FALSE),
  plos_bmc = list(width = 170, max_h = 200, base = 8, dense = 8, tag = 10, tags = "A", tiff = TRUE))
# Axis titles are set 20% larger than the base text of the profile (8.4 pt at base 7), so they stand out from tick labels.
AXIS_TITLE_SCALE <- 1.2
AXIS_TEXT_SCALE <- 1.2
FIG_PROFILE <- Sys.getenv("PAPER3_FIG_PROFILE", "default")
if (!FIG_PROFILE %in% names(FIG_PROFILES)) stop("Unknown PAPER3_FIG_PROFILE: ", FIG_PROFILE)
FIG <- FIG_PROFILES[[FIG_PROFILE]]

# Fonts must resolve to the real files: a silent fallback would change metrics and the embedded font.
FONT_SANS <- "Arial"
FONT_MONO <- "Courier New"
for (family in c(FONT_SANS, FONT_MONO)) {
  path <- systemfonts::match_fonts(family)$path
  if (!grepl(tolower(gsub(" ", "", if (family == FONT_SANS) "arial" else "cour")), tolower(basename(path)), fixed = TRUE)) {
    stop("Font '", family, "' is not installed (matched ", path, "). Install the Microsoft core font before rendering.")
  }
}

# Units: lw() converts points to ggplot linewidth (1 lwd = 1/96 inch = 0.75 pt), fs() points to text size (mm).
lw <- function(pt) pt / (.pt * 0.75)
fs <- function(pt) pt / .pt

# Platforms (Okabe-Ito). The order ONT, Illumina, BGI/DNBSEQ, PacBio holds in every figure.
pal_platform <- c(ONT = "#0072B2", Illumina = "#D55E00", `BGI/DNBSEQ` = "#CC79A7", PacBio = "#009E73", `All platforms` = "#1A1A1A")
platform_levels <- names(pal_platform)
scale_colour_platform <- function(...) scale_colour_manual(values = pal_platform, breaks = platform_levels, ...)
scale_fill_platform <- function(...) scale_fill_manual(values = pal_platform, breaks = platform_levels, ...)

# Evidence states, identical in every figure.
ink <- "#000000"
grey_dark <- "#404040"     # testable / included points and bars, outlines of counted states
grey_mid <- "#9E9E9E"      # excluded bars, rejected marks
grey_light <- "#BDBDBD"    # chance-limited points
grey_fill <- "#E0E0E0"     # counted positive boxes, aligned segments
grey_band <- "#D9D9D9"     # control intervals
grey_text <- "#595959"     # secondary text (flanks, notes inside panels); contrast 7:1 on white
accent_alert <- "#8E2043"  # thresholds and missed spikes only, never in a panel with BGI/DNBSEQ
state_point_fill <- c(testable = grey_dark, `chance-limited` = grey_light, untestable = "white")

theme_paper3 <- function(base = FIG$base, family = FONT_SANS) {
  theme_classic(base_size = base, base_family = family, base_line_size = lw(0.5), base_rect_size = lw(0.5)) +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(size = base * AXIS_TEXT_SCALE, colour = ink),
      axis.title = element_text(size = base * AXIS_TITLE_SCALE, face = "plain"),
      axis.line = element_line(colour = ink, linewidth = lw(0.5)),
      axis.ticks = element_line(colour = ink, linewidth = lw(0.5)),
      axis.ticks.length = unit(1, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base, face = "bold", hjust = 0, margin = margin(0.5, 0.5, 0.5, 0.5, "mm")),
      strip.placement = "outside",
      legend.title = element_text(size = base),
      legend.text = element_text(size = base),
      legend.key.size = unit(3, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.5, "mm"),
      legend.background = element_blank(),
      panel.spacing = unit(2.5, "mm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(size = FIG$tag, face = "bold", hjust = 0, vjust = 1, family = family),
      plot.tag.position = c(0, 1),
      plot.tag.location = "plot",
      plot.margin = margin(3, 1.5, 1.5, 1.5, "mm"))
}
theme_set(theme_paper3())
update_geom_defaults("point", list(size = 1.5, stroke = lw(0.5)))
update_geom_defaults("text", list(family = FONT_SANS, size = fs(FIG$base), colour = ink))
update_geom_defaults("line", list(linewidth = lw(0.75)))
update_geom_defaults("step", list(linewidth = lw(0.75)))
update_geom_defaults("segment", list(linewidth = lw(0.5)))
update_geom_defaults("linerange", list(linewidth = lw(0.75)))

# Panel tags come from patchwork, one level only (A, B, C or a, b, c).
assemble <- function(..., heights = NULL, widths = NULL, design = NULL) {
  wrap_plots(..., heights = heights, widths = widths, design = design) + plot_annotation(tag_levels = FIG$tags)
}

# Clopper-Pearson interval of k successes in n trials (vectorised).
ci_binom <- function(k, n, conf = 0.95) {
  a <- (1 - conf) / 2
  data.frame(lo = ifelse(n == 0, NA_real_, ifelse(k == 0, 0, qbeta(a, k, n - k + 1))),
             hi = ifelse(n == 0, NA_real_, ifelse(k == n, 1, qbeta(1 - a, k + 1, n - k))))
}

# Pseudo-logarithmic axes for data with zeros (asinh-type transform, sigma 1, base 10).
pseudo_log <- transform_pseudo_log(sigma = 1, base = 10)
# Tick labels 0, 1, 10, 10^2, ...: plotmath expressions with a true superscript.
label_power <- function(x) {
  lapply(x, function(value) {
    if (is.na(value)) return(NA)
    if (value == 0) return(quote(0))
    exponent <- round(log10(value))
    if (abs(value - 10^exponent) > 1e-9 * value) return(format(value, big.mark = ",", trim = TRUE, scientific = FALSE))
    if (exponent == 0) quote(1) else if (exponent == 1) quote(10) else bquote(10^.(exponent))
  })
}
label_power_expr <- function(x) as.expression(label_power(x))
count_label <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
