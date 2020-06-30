---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Use the `saveRDS` function to save the R object causing the bug.
2. Upload the RDS file.
For example:
```r
paths <- lineagePath(tree)
fixations <- fixationSites(paths) # Error occurs

saveRDS("paths.rds", paths) # Upload "paths.rds"
```

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
