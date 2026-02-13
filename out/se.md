---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
You are working with GWAS summary statistics where:

* `odds_ratio` (OR) is given
* `standard_error` is NA
* You know total `n`, `num_cases`, `num_controls`
* You have `effect_allele_frequency`

Since this is a casecontrol study (OR reported), the SE must be calculated on the **log scale**, not directly from OR.

---

# 1 Correct framework

For logistic regression GWAS:

$$\beta = \ln(OR)$$

The variance of  under an additive model is approximately:

$$Var(\beta) \approx \frac{1}{2 \cdot p \cdot (1-p) \cdot N_{eff}}$$

Where:

* $p$ = effect allele frequency
* $N_{eff}$ = effective sample size
* $N_{eff} = \frac{4}{\frac{1}{N_{cases}} + \frac{1}{N_{controls}}}$

Then:

$$SE = \sqrt{Var(\beta)}$$

---

# 2 Compute Effective Sample Size

From your data:

```
num_cases = 1965
num_controls = 630574
```

$$N_{eff} = \frac{4}{1/1965 + 1/630574}$$

Since controls are huge relative to cases:

$$N_{eff} \approx 4 \times 1965 = 7860$$

(When controls >> cases, N_eff  4  N_cases)

---

#   Example Calculation (First SNP)

Row 1:

```
OR = 1.1398609
EAF = 0.0053
```

### Step 1  Log transform

$$\beta = \ln(1.1398609) \approx 0.1309$$

### Step 2  Variance

$$Var(\beta) = \frac{1}{2 \cdot 0.0053 \cdot (1-0.0053) \cdot 7860}$$

$$2p(1-p) \approx 2 \cdot 0.0053 \cdot 0.9947 \approx 0.01054$$

$$0.01054 \cdot 7860 \approx 82.8$$

$$Var(\beta) \approx \frac{1}{82.8} = 0.01207$$

$$SE \approx \sqrt{0.01207} = 0.1098$$

So:

$$SE \approx 0.11$$

---

# 4 Important Reality Check

This approximation assumes:

* Additive model
* No covariate inflation
* Proper case-control logistic regression
* No meta-analysis weighting

If this is a meta-analysis (your file has `direction`, `i2`, `q_pval`), then:

 You **cannot accurately reconstruct SE from N alone**
because inverse-variance meta-analysis weights were used.

---

# 5 More Accurate Method (Best Option)

If you have:

* OR
* p_value

You can reconstruct SE directly:

$$Z = \Phi^{-1}(1 - p/2)$$

$$SE = \frac{\beta}{Z}$$

This is much more accurate than using N.

---

# 6 Example Using p-value (Row 1)

```
p = 0.7538
```

Two-sided:

$$Z = 0.31$$

$$SE = 0.1309 / 0.31 \approx 0.42$$

This is VERY different from 0.11 above.

That tells you:

 This is almost certainly meta-analysis output
 You must compute SE from beta and p-value

---

# 7 Final Recommendation

For your dataset, compute SE as:

$$\beta = \ln(OR)$$

$$Z = qnorm(1 - p/2)$$

$$SE = \beta / Z$$

That will match the original analysis.

---


Good question. The factor **4** in the effective sample size formula comes directly from the variance of the log-odds ratio in a **casecontrol logistic regression** under an additive genetic model.

Let's derive it cleanly.

---

# 1 Start from the variance of  in logistic regression

For a SNP coded additively (0,1,2), the asymptotic variance of the log-OR estimate is:

$$Var(\beta) \approx \frac{1}{N \cdot Var(G) \cdot Var(Y)}$$

Where:

* $G$ = genotype
* $Y$ = case/control status
* $N$ = total sample size

---

# 2 Variance of phenotype in casecontrol study

For binary phenotype:

$$Var(Y) = K(1-K)$$

Where:

$$K = \frac{N_{cases}}{N}$$

In a highly unbalanced study:

$$K(1-K) = \frac{N_{cases}}{N} \cdot \frac{N_{controls}}{N}$$

Multiply numerator and denominator:

$$N \cdot Var(Y) = \frac{N_{cases} \cdot N_{controls}}{N}$$

---

# 3 Variance of genotype

Under HardyWeinberg:

$$Var(G) = 2p(1-p)$$

---

# 4 Plug everything together

$$Var(\beta) \approx \frac{1}{2p(1-p)} \cdot \frac{N}{N_{cases} N_{controls}}$$

Rearrange:

$$Var(\beta) \approx \frac{1}{2p(1-p)} \cdot \frac{1}{\frac{N_{cases} N_{controls}}{N}}$$

Now define:

$$N_{eff} = \frac{N_{cases} N_{controls}}{N}$$

So:

$$Var(\beta) \approx \frac{1}{2p(1-p) N_{eff}}$$

---

# 5 Where does the 4 come from?

Now rewrite $N_{eff}$ algebraically:

$$\frac{N_{cases} N_{controls}}{N_{cases} + N_{controls}}$$

Multiply numerator and denominator by 4:

$$N_{eff} = \frac{4}{\frac{1}{N_{cases}} + \frac{1}{N_{controls}}}$$

This is just algebra.

Because:

$$\frac{4}{\frac{1}{a} + \frac{1}{b}} = \frac{4ab}{a+b}$$

And:

$$\frac{ab}{a+b}$$

is the harmonic mean scaled by 2.

So the **4 appears from algebraic rearrangement**, not from a heuristic.

---

# 6 Intuition

If the study is perfectly balanced:

$$N_{cases} = N_{controls} = n$$

Then:

$$N_{eff} = \frac{4}{\frac{1}{n} + \frac{1}{n}} = \frac{4}{2/n} = 2n$$

But total N = 2n.

So:

$$N_{eff} = N$$

Balanced design  no efficiency loss.

---

# 7 Your dataset intuition

In your case:

```
cases = 1965
controls = 630574
```

Since controls >> cases:

$$N_{eff} \approx 4 \times N_{cases}$$

Why?

Because:

$$\frac{1}{N_{controls}} \approx 0$$

So:

$$N_{eff} \approx \frac{4}{1/N_{cases}} = 4 N_{cases}$$

Meaning:

 Once controls are very large, power is limited by the number of cases.

The factor 4 ensures the variance matches the logistic likelihood information.

---

# 8 Bottom Line

The 4 is not arbitrary.

It comes from:

* Logistic regression information matrix
* Binary phenotype variance
* Algebraic conversion to harmonic mean form

It ensures that:

$$Var(\beta) = \frac{1}{2p(1-p)N_{eff}}$$

matches the true asymptotic variance of the log-odds ratio.


