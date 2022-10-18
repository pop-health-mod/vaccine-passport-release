# One-covariate models -----
# create formula as a string
predictors <- c(
  # base ITS model with level and slope change
  "ns(week, knots = knot_pos)", "pass_anno", "pass_anno:week_anno",
  # intercept for prior vaccine coverage
  "vaxcov_grp",
  ## intercept and slope for age
  # pre-intervention trend and impact of intervention
  "age", "age:ns(week, knots = knot_pos)",
  "age:pass_anno", "age:pass_anno:week_anno"
)

fmla_age <- paste("rate_1dose ~ ", paste(predictors, collapse = " + "), sep = "")

# create formula by replacing age with income
fmla_income <- gsub("age", "quin_income", fmla_age)

# create formula by replacing age with vismin
fmla_vismin <- gsub("age", "quin_vismin", fmla_age)

# Interaction models ----
predictors_age.income <- c(
  # base ITS model with level and slope change
  "ns(week, knots = knot_pos)", "pass_anno", "pass_anno:week_anno",
  
  # intercept for prior vaccine coverage
  "vaxcov_grp",
  
  ## intercept and slope for age
  # pre-intervention trend and impact of intervention
  "age", "age:ns(week, knots = knot_pos)", 
  "age:pass_anno", "age:pass_anno:week_anno",
  
  ## intercept and slope for income quintile
  # pre-intervention trend and impact of intervention
  "quin_income", "quin_income:ns(week, knots = knot_pos)", 
  "quin_income:pass_anno", "quin_income:pass_anno:week_anno",
  
  # age x income interactions
  "age:quin_income", "age:quin_income:ns(week, knots = knot_pos)",
  "age:quin_income:pass_anno", "age:quin_income:pass_anno:week_anno"
)

fmla_age.income <- paste("rate_1dose ~ ", paste(predictors_age.income, collapse = " + "), sep = "")

# create formula by replacing income with vismin
fmla_age.vismin <- gsub("quin_income", "quin_vismin", fmla_age.income)
