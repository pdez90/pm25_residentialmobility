# ============================================================
# DAG: how housing markets translate baseline inequality into
# unequal exposure through residential mobility (reviewer request)
# ------------------------------------------------------------
# Pure ggplot2 (no extra packages). Solid arrows = pathway studied;
# dashed = confounding/selection paths the design must reckon with.
# Output: figures/figS_dag.png
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out_dir <- "/Users/priyanka/MIGRATE_dir/pm25_avoidance_lag3_olag1"
fig_dir <- file.path(out_dir,"figures"); dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

nodes <- data.table(
  id=c("hist","orig","mkt","move","dest","gap","ses","pref"),
  x  =c(0.6,  2.2,  4.4,  6.6,  8.4,  8.4,  4.4,  6.6),
  y  =c(2.6,  2.6,  2.6,  2.6,  2.6,  1.05, 4.15, 4.15),
  lab=c("Historical siting,\nredlining, zoning",
        "Baseline neighborhood\npollution + composition",
        "Housing market\n(prices, discrimination,\nsteering, segregation)",
        "Move\ndestination",
        "Destination\nexposure",
        "Observed avoidance gap\n(conditional on origin pollution)",
        "Origin SES\n(income, poverty,\neducation)  [controlled]",
        "Unobserved preferences,\nselection into moving\n[bounded by robustness]"),
  main=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE))
edges <- data.table(
  from=c("hist","orig","mkt","move","dest","ses","ses","pref","orig"),
  to  =c("orig","mkt","move","dest","gap","mkt","move","move","gap"),
  style=c("solid","solid","solid","solid","solid","dashed","dashed","dashed","dashed"))
E <- merge(merge(edges, nodes[, .(from=id, x0=x, y0=y)], by="from"),
           nodes[, .(to=id, x1=x, y1=y)], by="to")
# shrink arrows so they start/end at box edges
E[, `:=`(dx=x1-x0, dy=y1-y0)]; E[, L:=sqrt(dx^2+dy^2)]
E[, `:=`(xs=x0+dx/L*0.78, ys=y0+dy/L*0.42, xe=x1-dx/L*0.78, ye=y1-dy/L*0.42)]

p <- ggplot() +
  geom_segment(data=E, aes(x=xs,y=ys,xend=xe,yend=ye,linetype=style),
               arrow=arrow(length=unit(0.16,"cm"), type="closed"),
               linewidth=.55, color="grey25") +
  geom_label(data=nodes, aes(x,y,label=lab, fill=main),
             size=2.9, lineheight=.95, label.padding=unit(0.18,"cm"),
             label.r=unit(0.12,"cm"), color="black") +
  scale_fill_manual(values=c(`TRUE`="#e8f0e8",`FALSE`="#f5f5f5"), guide="none") +
  scale_linetype_identity() +
  annotate("text", x=0.6, y=0.9, hjust=0, size=2.7, color="grey35",
           label="Solid: pathway studied.  Dashed: confounding/selection paths (controlled, or bounded by robustness checks).") +
  coord_cartesian(xlim=c(-0.1,9.6), ylim=c(0.6,4.9), expand=FALSE) +
  theme_void()
ggsave(file.path(fig_dir,"figS_dag.png"), p, width=9.5, height=4.2, dpi=300)
cat("DONE. DAG written to figures/figS_dag.png\n")
