# ---- paths: this repo contains no absolute paths; see scripts/_paths.R ----
# Set MIGRATE_DATA_DIR to your data root. Run scripts from the repository root.
if (!exists(".migrate_paths_loaded")) source(
  if (file.exists("scripts/_paths.R")) "scripts/_paths.R"
  else if (file.exists("upstream/_paths.R")) "upstream/_paths.R" else "_paths.R")

# ============================================================
# CONCEPTUAL SIGNATURE FIGURE
# ------------------------------------------------------------
# The paper's argument as a flow diagram: regulation cleans the air in
# place (left, equalizing) while residential mobility and slower
# minority-area cleanup (right, inequitable) work against convergence,
# so ABSOLUTE gaps narrow but RELATIVE gaps persist.
# Pure ggplot2, no external rendering engine. Writes a PNG.
#
# Output: pm25_avoidance_lag3_olag1/figures/fig_conceptual_signature.png
# ============================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(stringr) })
base_dir <- migrate_data_dir()
fig_dir  <- file.path(base_dir, "pm25_avoidance_lag3_olag1", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# palette
NAVY <- "#1F3864"; GREEN <- "#1b7837"; GREENF <- "#d9f0d3"; RED <- "#b2182b"; REDF <- "#fddbc7"; GREYF <- "#eeeeee"
wrap <- function(s, w=26) str_wrap(s, width=w)

# nodes: id, x, y (center), w, h, fill, border, textcol, label
nodes <- rbindlist(list(
 data.table(id="top", x=5.0, y=9.1, w=6.6, h=0.95, fill=NAVY,  border=NAVY,  tcol="white",
            label="Air pollution fell nationally, 2011-2019 (Clean Air Act regulation)"),
 # LEFT / equalizing track
 data.table(id="L1",  x=2.4, y=6.9, w=3.7, h=1.05, fill=GREYF, border="grey55", tcol="black",
            label="Neighborhoods clean up IN PLACE (nearly all of the decline)"),
 data.table(id="L2",  x=2.4, y=4.5, w=3.7, h=1.20, fill=GREENF, border=GREEN, tcol="black",
            label="Minority areas' ABSOLUTE cleanup was larger (they started dirtier)"),
 data.table(id="Lout",x=2.4, y=1.7, w=3.7, h=1.10, fill=GREEN, border=GREEN, tcol="white",
            label="ABSOLUTE Black-white gap NARROWS (-19%)"),
 # RIGHT / inequitable track
 data.table(id="R1",  x=7.6, y=6.9, w=3.9, h=1.05, fill=REDF,  border=RED, tcol="black",
            label="Residential mobility redistributes millions across the pollution surface"),
 data.table(id="R2",  x=7.7, y=4.6, w=3.9, h=1.80, fill=RED,   border=RED, tcol="white",
            label="Conditional on baseline, minority-origin movers reach dirtier destinations (avoidance gap), and dirtier places pull minority in-migration (sorting)"),
 data.table(id="Rout",x=7.6, y=1.7, w=3.9, h=1.10, fill=RED,   border=RED, tcol="white",
            label="RELATIVE disparities PERSIST"),
 # CENTER / siting (pollution dynamics, inequitable)
 data.table(id="SIT", x=4.5, y=3.05, w=3.2, h=1.30, fill=REDF,  border=RED, tcol="black",
            label="But net of starting pollution, cleanup was SLOWER where minorities live (unequal cleanup)")
), fill=TRUE)

# edges: from id -> to id
edges <- rbindlist(lapply(list(c("top","L1"),c("top","R1"),c("top","SIT"),
                               c("L1","L2"),c("L2","Lout"),
                               c("R1","R2"),c("R2","Rout"),
                               c("SIT","Rout")), function(e){
  a <- nodes[id==e[1]]; b <- nodes[id==e[2]]
  data.table(x=a$x, y=a$y-a$h/2, xend=b$x, yend=b$y+b$h/2,
             col=ifelse(e[2] %in% c("L2","Lout"), GREEN, ifelse(e[1]=="top", "grey45", RED)))
}))

p <- ggplot() +
  geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend, color=I(col)),
               arrow=arrow(length=unit(.14,"in"), type="closed"), linewidth=.8) +
  geom_rect(data=nodes, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2),
            fill=nodes$fill, color=nodes$border, linewidth=.7) +
  geom_text(data=nodes, aes(x=x, y=y, label=wrap(label)), color=nodes$tcol,
            size=3.0, lineheight=.95, fontface=ifelse(nodes$id %in% c("Lout","Rout","top"),"bold","plain")) +
  annotate("text", x=2.4, y=8.05, label="What equalized", fontface="bold", color=GREEN, size=3.4) +
  annotate("text", x=7.6, y=8.05, label="What did not", fontface="bold", color=RED, size=3.4) +
  labs(title="Residential mobility reproduces environmental inequality even as pollution falls",
       subtitle="Regulation cleaned the air in place (absolute gaps narrow); unequal mobility and slower minority-area cleanup keep relative gaps open") +
  coord_cartesian(xlim=c(0,10), ylim=c(0.9,9.8)) +
  theme_void(base_size=12) +
  theme(plot.title=element_text(face="bold", hjust=.5),
        plot.subtitle=element_text(hjust=.5, size=9.5, color="grey30"),
        plot.margin=margin(10,10,10,10))

ggsave(file.path(fig_dir, "fig_conceptual_signature.png"), p, width=10, height=7.2, dpi=200)
message("DONE. Conceptual figure: ", file.path(fig_dir, "fig_conceptual_signature.png"))

# ============================================================
# GRAPHICAL ABSTRACT (separate file; not a manuscript figure)
# ------------------------------------------------------------
# IOP displays graphical abstracts at 500 px wide (120 px on the issue
# page) and asks for >= 12 pt type, so this is a deliberately reduced
# version of the conceptual figure: three beats, large text, wide format.
# Upload separately as a Source File named with the "GA" prefix; do NOT
# place it in the manuscript PDF or the SI.
#
# Output: figures/GA_conceptual_signature.png  (2400 x 1110 px)
# ============================================================
ga <- rbindlist(list(
  data.table(x=1.85, w=2.7, fill=GREENF, border=GREEN, tcol="black",
             label="Regulation cleaned\nthe air in place"),
  data.table(x=5.00, w=2.7, fill=REDF,   border=RED,   tcol="black",
             label="Mobility reallocated\nit unequally"),
  data.table(x=8.15, w=2.7, fill=NAVY,   border=NAVY,  tcol="white",
             label="Absolute gaps closed\nRelative gaps did not")
))
ga[, `:=`(y=2.45, h=1.75)]
ga_arrows <- data.table(x=c(3.25, 6.40), xend=c(3.60, 6.75), y=2.45, yend=2.45)

pga <- ggplot() +
  # three panels
  geom_rect(data=ga, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2),
            fill=ga$fill, color=ga$border, linewidth=1.1) +
  geom_text(data=ga, aes(x=x, y=y, label=label), color=ga$tcol,
            size=5.4, lineheight=1.05, fontface="bold") +
  # connectors
  geom_segment(data=ga_arrows, aes(x=x, xend=xend, y=y, yend=yend),
               arrow=arrow(length=unit(.16,"in"), type="closed"),
               linewidth=1.0, color="grey35") +
  # headline and the single key number
  annotate("text", x=5, y=4.15, label="Residential mobility reproduces environmental-opportunity inequality",
           size=5.4, fontface="bold", color=NAVY) +
  annotate("text", x=5, y=3.55,
           label="National block-group migration flows \u00d7 eight pollutants, 2011\u20132019",
           size=4.6, color="grey25") +
  annotate("text", x=5, y=0.80,
           label="Movers from majority-Black neighborhoods reached destinations 0.118 \u00b5g m\u207b\u00b3\ndirtier than comparable white-origin movers leaving equally polluted places",
           size=4.6, fontface="bold", color=RED, lineheight=1.15) +
  coord_cartesian(xlim=c(0.3, 9.7), ylim=c(0.1, 4.6), expand=FALSE) +
  theme_void()

ggsave(file.path(fig_dir, "GA_conceptual_signature.png"), pga,
       width=8, height=3.7, dpi=300)
cat("Graphical abstract written: ", file.path(fig_dir, "GA_conceptual_signature.png"), "\n", sep="")
