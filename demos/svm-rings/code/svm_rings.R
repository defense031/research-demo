# Support vector machines and the kernel trick on two concentric rings.
# No line separates the rings in 2D. Lifting each point with
# z = x1^2 + x2^2 lets a flat plane separate them. This script fits both
# SVMs and writes every figure, table, and number that main.tex reads.
library(e1071)
library(ggplot2)

set.seed(376)
n_per <- 150
ring <- function(n, r, label) {
  theta <- runif(n, 0, 2 * pi)
  rad   <- r + rnorm(n, 0, 0.18)
  data.frame(x1 = rad * cos(theta), x2 = rad * sin(theta), ring = label)
}
dat <- rbind(ring(n_per, 1.0, "inner"), ring(n_per, 2.1, "outer"))
dat$ring <- factor(dat$ring)
dat$z <- dat$x1^2 + dat$x2^2                      # the lift

test_idx <- unlist(lapply(split(seq_len(nrow(dat)), dat$ring),
                          function(i) sample(i, round(0.3 * length(i)))))
train <- dat[-test_idx, ]
test  <- dat[test_idx, ]

# Linear SVM in the original 2D space
fit2d <- svm(ring ~ x1 + x2, data = train, kernel = "linear", cost = 1, scale = FALSE)
acc2d <- mean(predict(fit2d, test) == test$ring)

# Linear SVM in the lifted 3D space: its boundary is a flat plane
fit3d <- svm(ring ~ x1 + x2 + z, data = train, kernel = "linear", cost = 1, scale = FALSE)
pred3d <- predict(fit3d, test)
acc3d <- mean(pred3d == test$ring)
w   <- drop(t(fit3d$coefs) %*% fit3d$SV)          # w1, w2, w3 (on x1, x2, z)
rho <- fit3d$rho                                  # plane: w . (x1,x2,z) = rho
plane_z <- function(x1, x2) (rho - w[1] * x1 - w[2] * x2) / w[3]

cols <- c(inner = "#B5985A", outer = "#2B2B28")
lim  <- c(-2.9, 2.9)
base_theme <- theme_minimal(base_size = 28) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        legend.text = element_text(size = 28))

# ---- 2D plots ---------------------------------------------------------------
p2d <- ggplot(dat, aes(x1, x2, colour = ring)) +
  geom_point(size = 2.8, alpha = 0.8) +
  scale_colour_manual(values = cols, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression(x[1]), y = expression(x[2])) + base_theme
ggsave("../figures/rings_2d.pdf", p2d, width = 6, height = 6)

# best linear boundary in 2D: w2d . x = rho2d
w2d <- drop(t(fit2d$coefs) %*% fit2d$SV)
grid2 <- expand.grid(x1 = seq(lim[1], lim[2], length.out = 300),
                     x2 = seq(lim[1], lim[2], length.out = 300))
grid2$region <- predict(fit2d, grid2)
dat$pred2d <- predict(fit2d, dat)
wrong2d <- dat[dat$pred2d != dat$ring, ]
p2d_line <- ggplot() +
  geom_raster(data = grid2, aes(x1, x2, fill = region), alpha = 0.18) +
  geom_point(data = dat, aes(x1, x2, colour = ring), size = 2.2, alpha = 0.8) +
  geom_point(data = wrong2d, aes(x1, x2), shape = 21, size = 4.2, stroke = 1.2,
             colour = "#B22222", fill = NA) +
  geom_abline(intercept = fit2d$rho / w2d[2], slope = -w2d[1] / w2d[2],
              linewidth = 1.4, colour = "#B22222", linetype = "dashed") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression(x[1]), y = expression(x[2])) + base_theme
ggsave("../figures/rings_2d_line.pdf", p2d_line, width = 6, height = 6)

# ---- 3D flipbook: tilt from a top view into 3D ------------------------------
proj <- function(x, y, z, tilt, az, zscale = 0.55) {
  X <- x * cos(az) - y * sin(az)
  Y <- x * sin(az) + y * cos(az)
  data.frame(sx = X, sy = Y * cos(tilt) + zscale * z * sin(tilt),
             depth = Y * sin(tilt) - zscale * z * cos(tilt))
}
n_tilt <- 10
tilts <- seq(0, 68, length.out = n_tilt) * pi / 180
azs   <- seq(0, 25, length.out = n_tilt) * pi / 180
frames <- c(seq_len(n_tilt), n_tilt + 1)          # the last frame adds the plane
sq <- expand.grid(x1 = c(-2.9, 2.9), x2 = c(-2.9, 2.9))[c(1, 2, 4, 3), ]

for (k in frames) {
  i <- min(k, n_tilt); show_plane <- k > n_tilt
  pts <- cbind(dat, proj(dat$x1, dat$x2, dat$z, tilts[i], azs[i]))
  pts$above <- dat$z > plane_z(dat$x1, dat$x2)
  floor <- proj(sq$x1, sq$x2, 0, tilts[i], azs[i])
  pl    <- proj(sq$x1, sq$x2, plane_z(sq$x1, sq$x2), tilts[i], azs[i])
  g <- ggplot() +
    geom_polygon(data = floor, aes(sx, sy), fill = NA, colour = "grey80")
  if (show_plane) {
    g <- g + geom_point(data = pts[!pts$above, ], aes(sx, sy, colour = ring), size = 2.8, alpha = 0.8) +
      geom_polygon(data = pl, aes(sx, sy), fill = "#B22222", alpha = 0.22, colour = "#B22222") +
      geom_point(data = pts[pts$above, ], aes(sx, sy, colour = ring), size = 2.8, alpha = 0.8)
  } else {
    pts <- pts[order(-pts$depth), ]
    g <- g + geom_point(data = pts, aes(sx, sy, colour = ring), size = 2.8, alpha = 0.8)
  }
  g <- g + scale_colour_manual(values = cols, name = NULL) +
    coord_equal(xlim = c(-3.3, 3.3), ylim = c(-3.2, 5.0)) +
    theme_void(base_size = 28) + theme(legend.position = "top",
                                        legend.text = element_text(size = 28))
  ggsave(sprintf("../figures/lift_%02d.pdf", k), g, width = 6, height = 6.6)
}
writeLines(sprintf("\\includegraphics<%d>[height=0.70\\textheight]{figures/lift_%02d.pdf}%%", frames, frames),
           "../tables/lift_frames.tex")

# ---- back in 2D: the plane meets the bowl in a circle -----------------------
grid <- expand.grid(x1 = seq(lim[1], lim[2], length.out = 300),
                    x2 = seq(lim[1], lim[2], length.out = 300))
grid$z <- grid$x1^2 + grid$x2^2
grid$region <- predict(fit3d, grid)
grid$dv <- attr(predict(fit3d, grid, decision.values = TRUE), "decision.values")[, 1]
pback <- ggplot() +
  geom_raster(data = grid, aes(x1, x2, fill = region), alpha = 0.18) +
  geom_contour(data = grid, aes(x1, x2, z = dv), breaks = 0, colour = "#B22222", linewidth = 1.4) +
  geom_point(data = train, aes(x1, x2, colour = ring), size = 1.8, alpha = 0.45) +
  geom_point(data = test, aes(x1, x2, colour = ring), size = 2.8, shape = 21, stroke = 1.3, fill = "white") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_colour_manual(values = cols, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression(x[1]), y = expression(x[2])) + base_theme
ggsave("../figures/rings_back.pdf", pback, width = 6, height = 6)

# ---- tables and numbers -----------------------------------------------------
pct <- function(a) sprintf("%.0f\\%%", 100 * a)
writeLines(c(
  "\\begin{tabular}{lc}", "\\toprule",
  "\\textbf{Linear SVM on} & \\textbf{Held-out correct} \\\\", "\\midrule",
  sprintf("$(x_1, x_2)$ & %s \\\\", pct(acc2d)),
  sprintf("$(x_1, x_2, x_1^2 + x_2^2)$ & %s \\\\", pct(acc3d)),
  "\\bottomrule", "\\end{tabular}"), "../tables/accuracy.tex")
writeLines(c(
  sprintf("\\newcommand{\\npts}{%d}", nrow(dat)),
  sprintf("\\newcommand{\\ntrain}{%d}", nrow(train)),
  sprintf("\\newcommand{\\ntest}{%d}", nrow(test)),
  sprintf("\\newcommand{\\acctwo}{%s}", pct(acc2d)),
  sprintf("\\newcommand{\\accthree}{%s}", pct(acc3d)),
  sprintf("\\newcommand{\\nsv}{%d}", fit3d$tot.nSV)), "../tables/numbers.tex")

cat("2D linear acc:", acc2d, " lifted linear acc:", acc3d, " SVs:", fit3d$tot.nSV, "\n")
cat("plane w:", round(w, 3), " rho:", round(rho, 3), "\n")
print(table(Prediction = pred3d, Truth = test$ring))

# ---- the animated HTML deck (index.html) ----
# e1071: a positive decision value means the first class level
pos_class   <- fit3d$levels[1]
pos_class2d <- fit2d$levels[1]
dat$set <- "train"; dat$set[test_idx] <- "test"
web <- list(
  points = data.frame(x1 = round(dat$x1, 4), x2 = round(dat$x2, 4),
                      ring = as.character(dat$ring), set = dat$set),
  plane = list(w = round(w, 6), rho = round(rho, 6), pos_class = pos_class),
  line2d = list(w = round(w2d, 6), rho = round(fit2d$rho, 6), pos_class = pos_class2d),
  n = nrow(dat), ntrain = nrow(train), ntest = nrow(test),
  acc2d = round(100 * acc2d), acc3d = round(100 * acc3d), nsv = fit3d$tot.nSV)
# One self-contained page: data and logos go inline, so the file plays on its
# own after a download and on GitHub Pages.
data_js <- paste0("window.SVM = ", jsonlite::toJSON(web, dataframe = "columns",
                  auto_unbox = TRUE, digits = NA), ";")
img64 <- function(f) paste0("data:image/png;base64,",
                            jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size)))
page <- paste(readLines("deck_template.html", warn = FALSE), collapse = "\n")
page <- sub("<!-- Template.*?-->\n", "<!-- Built by code/svm_rings.R from code/deck_template.html. Edit the template, not this file. -->\n", page)
page <- sub("/*SVM_DATA*/", data_js, page, fixed = TRUE)
page <- sub("{{CREST}}", img64("../figures/usma_crest.png"), page, fixed = TRUE)
page <- sub("{{MULE}}", img64("../figures/math_mule.png"), page, fixed = TRUE)
writeLines(page, "../index.html")
