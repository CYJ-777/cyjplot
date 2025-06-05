
# cyjplot

🧬 `cyjplot` is a custom R package that provides easy-to-use ggplot2-based plotting functions.

---

## 📦 Installation

You can install the package directly from GitHub:

```r
install.packages("devtools")  # if not installed yet
devtools::install_github("CYJ-777/cyjplot")

library(cyjplot)
```

---

## 📊 Available Functions

### 1. `boxplot()`
Custom boxplot that:
- Automatically chooses between t-test or Wilcoxon test
- Ensures control group is blue, altered group is orange
- Saves both PNG and PDF

```r
boxplot(df, x = "Group", y = "Response", filename = "boxplot.png")
```

---

### 2. `correlation()`
Generates a scatterplot with Pearson correlation coefficient, p-value, and n.

```r
correlation(df, xvar = "GeneA", yvar = "DrugResponse", outfile_prefix = "GeneA_vs_Response")
```

---

### 3. `survival()`
Creates Kaplan-Meier survival curves with log-rank test and risk table.

```r
fit <- survfit(Surv(time, status) ~ group, data = df)
survival(fit, data = df, filename = "KM_curve")
```

---

### 4. `longboxplot()`
Faceted boxplot across multiple genes/conditions with optional strip background and annotation.

```r
longboxplot(expr_merged, pvals, title_text = "Expression Comparison", filename = "longplot")
```

---

## 🧠 Author

Developed by [CYJ-777](https://github.com/CYJ-777) (島田研)

---

